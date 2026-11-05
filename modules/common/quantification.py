##########
#IMPORTED#
##########

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, peak_widths
from scipy.ndimage import uniform_filter1d

import modules.common.transients as tr
import modules.common.preprocess as pp

def extract_dff_summary(fiberbehav_df, mouse, batch, group,
                         zone_cols, behav_cols,
                         baseline_col=None,
                         merge_into=None,
                         dff_col='dFF',
                         fps=20,
                         use_zscore=False,
                         shock_col=None,
                         shock_exclude_cols=None,
                         shock_exclude_seconds=2,
                         stim_sniff_pairs=None):
    """
    Extract mean dFF and AUC dFF during each zone/behaviour for one animal.

    fiberbehav_df : pd.DataFrame
            Combined fiber + behavior dataframe (fiberbehavnotderived.csv).
    mouse, batch, group : str
            Animal identifiers.
    zone_cols : list of str
            Binary zone columns. Time spent in each → mean dFF + AUC dFF.
            E.g. ['Closed arm', 'Open arm', 'Center'].
    behav_cols : list of str
            Binary behaviour columns. Same metrics computed independently.
            E.g. ['Head dipping'].
    baseline_col : str, optional
            Binary column in fiberbehav_df used as the z-score baseline mask
            (1 = baseline period) when use_zscore=True. E.g. 'Closed arm'.
            If None, the baseline is instead defined as the frames where none
            of the behav_cols are 1 (i.e. animal not engaged in any listed
            behaviour).
    merge_into : dict, optional
            {behaviour: zone} pairs where frames of a behaviour should ALSO be
            included in the zone's signal. 
            E.g. {'Head dipping': 'Open arm'} ensures head dipping frames
            are counted in Open arm dFF in addition to their own column.
    dff_col : str
            Column name for raw dFF signal.
    fps : float
            Frame rate, used to convert frame count to seconds for AUC 
            (AUC = sum(dFF) / fps, i.e. dFF integrated over time in seconds).
    use_zscore : bool
            If True, z-score the dFF trace using the median/std computed only
            from the baseline mask (see baseline_col), then apply that baseline
            to the entire trace.
    shock_col : str, optional
        Binary column marking shock delivery (1 = shock frame). If provided
        together with `shock_exclude_cols`, any bout of a column listed in
        `shock_exclude_cols` that contains a shock will have its last
        `shock_exclude_seconds` seconds excluded from the dFF metrics, so
        shock-evoked dFF doesn't artificially inflate that column's signal
        (e.g. CS+ trials that end in a shock vs CS- trials that never do).
    shock_exclude_cols : list of str, optional
        Columns (usually a subset of zone_cols/behav_cols, e.g. ['CS+'])
        for which shock-containing bouts get trimmed as described above.
        Bouts of these columns that do NOT contain a shock are left intact.
    shock_exclude_seconds : float
        Duration (in seconds) to trim from the end of shock-containing
        bouts. Default 2s.
    stim_sniff_pairs : dict, optional
        {stim_col: sniff_col} pairs, e.g. {'Stim Novel': 'Sniff Novel'}.
        For each pair, computes mean/AUC dFF over frames where the stim
        is active AND the corresponding sniff is NOT active (i.e. the
        'passive' stim-exposure period, excluding active sniffing bouts).
        Results are added to the main record as
        '{stim_col} no sniff mean dFF' / '{stim_col} no sniff AUC dFF'.
    """

    merge_into = merge_into or {}
    shock_exclude_cols = shock_exclude_cols or []
    stim_sniff_pairs = stim_sniff_pairs or {}

    # ── Signal preparation ────────────────────────────────────────────────────
    if dff_col not in fiberbehav_df.columns:
        raise ValueError(f"Column '{dff_col}' not found in dataframe.")

    raw_signal = fiberbehav_df[dff_col].values.copy()
    signal = raw_signal.copy()

    # ── Baseline mask (used for z-scoring if requested, and for baseline_record) ─
    if baseline_col is not None:
        if baseline_col not in fiberbehav_df.columns:
            raise ValueError(f"Baseline column '{baseline_col}' not found in dataframe.")
        baseline_mask = fiberbehav_df[baseline_col].values.astype(bool)
    else:
        baseline_mask = np.ones(len(fiberbehav_df), dtype=bool)
        for col in behav_cols:
            if col in fiberbehav_df.columns:
                baseline_mask &= ~fiberbehav_df[col].values.astype(bool)
            else:
                print(f"  Warning: behaviour column '{col}' not found, "
                      f"skipping it when building baseline mask.")

    baseline_signal_raw = raw_signal[baseline_mask]

    if use_zscore:
        if baseline_signal_raw.size < 2:
            raise ValueError(
                "Not enough baseline samples to compute z-score baseline "
                "(need at least 2)."
            )
        signal = (signal - np.median(baseline_signal_raw)) / np.std(baseline_signal_raw, ddof=1)

    record = {
        'Mouse' : mouse,
        'Batch' : batch,
        'Group' : group,
    }

    # ── Separate baseline record (always on raw dFF, regardless of use_zscore) ──
    baseline_record = {
        'Mouse' : mouse,
        'Batch' : batch,
        'Group' : group,
    }

    if baseline_signal_raw.size >= 2:
        baseline_record['Baseline std dFF']  = round(float(np.std(baseline_signal_raw, ddof=1)), 6)
        baseline_record['Baseline var dFF']  = round(float(np.var(baseline_signal_raw, ddof=1)), 6)
        baseline_record['Baseline mean dFF'] = round(float(np.mean(baseline_signal_raw)), 6)
        baseline_record['Baseline AUC dFF']  = round(float(np.sum(baseline_signal_raw) / fps), 6)
    else:
        print("  Warning: not enough baseline frames, baseline metrics set to NaN.")
        baseline_record['Baseline std dFF']  = float('nan')
        baseline_record['Baseline var dFF']  = float('nan')
        baseline_record['Baseline mean dFF'] = float('nan')
        baseline_record['Baseline AUC dFF']  = float('nan')

    # ── Whole-trace variability (raw dFF, includes stims and sniffs) ──────────
    valid = np.isfinite(raw_signal)
    if 'dFF ExclusionMask' in fiberbehav_df.columns:       # drop corrupted regions
        valid &= ~fiberbehav_df['dFF ExclusionMask'].values.astype(bool)

    whole = raw_signal[valid]
    if whole.size >= 2:
        baseline_record['Whole trace std dFF'] = round(float(np.std(whole, ddof=1)), 6)
        baseline_record['Whole trace var dFF'] = round(float(np.var(whole, ddof=1)), 6)
    else:
        baseline_record['Whole trace std dFF'] = float('nan')
        baseline_record['Whole trace var dFF'] = float('nan')

    # ── Helper: find contiguous True runs (bouts) in a boolean mask ───────────
    def _get_bouts(mask):
        bouts = []
        in_bout = False
        start = None
        for i, v in enumerate(mask):
            if v and not in_bout:
                start = i
                in_bout = True
            elif not v and in_bout:
                bouts.append((start, i - 1))
                in_bout = False
        if in_bout:
            bouts.append((start, len(mask) - 1))
        return bouts

    # ── Helper: trim last N seconds off shock-containing bouts ────────────────
    def _exclude_shock_tail(mask, col_label):
        if shock_col is None or col_label not in shock_exclude_cols:
            return mask
        if shock_col not in fiberbehav_df.columns:
            print(f"  Warning: shock column '{shock_col}' not found, "
                  f"skipping shock exclusion for '{col_label}'.")
            return mask

        shock_mask = fiberbehav_df[shock_col].values.astype(bool)
        n_trim = int(round(shock_exclude_seconds * fps))
        mask = mask.copy()

        for start, end in _get_bouts(mask):
            if shock_mask[start:end + 1].any():
                trim_start = max(start, end + 1 - n_trim)
                mask[trim_start:end + 1] = False

        return mask

    # ── Helper: compute metrics for a boolean mask ────────────────────────────
    def _metrics(mask, label):
        if mask.sum() == 0:
            print(f"  Warning: no frames found for '{label}', filling with NaN.")
            record[f'{label} mean dFF']  = float('nan')
            record[f'{label} AUC dFF']   = float('nan')
            return
        sig_masked             = signal[mask]
        record[f'{label} mean dFF'] = round(float(np.mean(sig_masked)), 6)
        record[f'{label} AUC dFF']  = round(float(np.sum(sig_masked) / fps), 6)

    # ── Zone metrics ──────────────────────────────────────────────────────────
    for col in zone_cols:
        if col not in fiberbehav_df.columns:
            print(f"  Warning: zone column '{col}' not found, filling with NaN.")
            record[f'{col} mean dFF'] = float('nan')
            record[f'{col} AUC dFF']  = float('nan')
            continue

        mask = fiberbehav_df[col].values.astype(bool)

        for behav, target_zone in merge_into.items():
            if target_zone == col and behav in fiberbehav_df.columns:
                behav_mask = fiberbehav_df[behav].values.astype(bool)
                mask = mask | behav_mask

        mask = _exclude_shock_tail(mask, col)

        _metrics(mask, col)

    # ── Behaviour metrics (independent, no merging) ───────────────────────────
    for col in behav_cols:
        if col not in fiberbehav_df.columns:
            print(f"  Warning: behaviour column '{col}' not found, filling with NaN.")
            record[f'{col} mean dFF'] = float('nan')
            record[f'{col} AUC dFF']  = float('nan')
            continue

        mask = fiberbehav_df[col].values.astype(bool)
        mask = _exclude_shock_tail(mask, col)
        _metrics(mask, col)

    # ── Stim-without-sniff metrics (plethysmograph experiments) ──────────────
    for stim_col, sniff_col in stim_sniff_pairs.items():
        label = f'{stim_col} no sniff'
        if stim_col not in fiberbehav_df.columns or sniff_col not in fiberbehav_df.columns:
            print(f"  Warning: '{stim_col}' or '{sniff_col}' not found, "
                  f"filling '{label}' with NaN.")
            record[f'{label} mean dFF'] = float('nan')
            record[f'{label} AUC dFF']  = float('nan')
            continue

        stim_mask  = fiberbehav_df[stim_col].values.astype(bool)
        sniff_mask = fiberbehav_df[sniff_col].values.astype(bool)
        mask = stim_mask & ~sniff_mask
        _metrics(mask, label)

    return record, baseline_record

def plot_transients(time, sig, valid, peaks, amplitudes, width_heights,
                    left_t, right_t, median_level, prominence, title,
                    zoom_start_s=0.0, zoom_len_s=120.0, save_path=None):
    """Whole trace + zoomed panel with threshold, detected peaks, amplitudes and widths."""
    fig, axes = plt.subplots(2, 1, figsize=(14, 7))
    t0 = time[0]
    windows = [(time[0], time[-1]), (t0 + zoom_start_s, t0 + zoom_start_s + zoom_len_s)]

    # contiguous excluded regions (for shading)
    inv = (~valid).astype(int)
    edges = np.diff(np.concatenate([[0], inv, [0]]))
    starts, ends = np.where(edges == 1)[0], np.where(edges == -1)[0] - 1

    for ax, (xmin, xmax) in zip(axes, windows):
        ax.plot(time, sig, color='k', lw=0.7, label='dFF')
        for s, e in zip(starts, ends):
            ax.axvspan(time[s], time[e], color='grey', alpha=0.3, lw=0)
        ax.axhline(median_level, color='tab:blue', ls=':', lw=1, label='median')
        ax.axhline(median_level + prominence, color='tab:red', ls='--', lw=1,
                   label=f'median + prominence thr ({prominence:.3f})')
        if len(peaks):
            ax.plot(time[peaks], sig[peaks], 'v', color='tab:red', ms=6, label='transients')
            ax.vlines(time[peaks], sig[peaks] - amplitudes, sig[peaks],
                      color='tab:red', lw=1.2)                           # amplitude
            ax.hlines(width_heights, left_t, right_t, color='tab:green', lw=1.5)  # width
        ax.set_xlim(xmin, xmax)
        # y-limits adapted to the visible window
        m = (time >= xmin) & (time <= xmax) & np.isfinite(sig)
        if m.any():
            pad = 0.05 * (sig[m].max() - sig[m].min() + 1e-12)
            ax.set_ylim(sig[m].min() - pad, sig[m].max() + pad)
        ax.set_ylabel('dFF')

    axes[0].set_title(f'{title} — {len(peaks)} transients')
    axes[0].legend(loc='upper right', fontsize=8, ncol=4)
    axes[1].set_xlabel('Time (s)')
    plt.tight_layout()
    if save_path is not None:
        fig.savefig(save_path, dpi=150)
        plt.close(fig)       # free memory when looping over mice
    return fig

def extract_transients(fiberbehav_df, mouse, batch, group,
                       zone_cols=None, behav_cols=None,
                       baseline_col=None,
                       baseline_exclude_cols=None,
                       merge_into=None,
                       stim_sniff_pairs=None,
                       dff_col='dFF',
                       fps=10,
                       threshold_mode='mad',
                       mad_k=3.0,
                       abs_prominence=None,
                       min_distance_s=0.5,
                       min_width_s=0.0,
                       smooth_s=0.0,
                       rel_height=0.5,
                       exclusion_col='dFF ExclusionMask',
                       plot=False, plot_path=None,
                       zoom_start_s=0.0, zoom_len_s=120.0):
    """
    Detect transients on the whole trace, then split metrics by period.

    Periods
    -------
    'Whole trace'  : all valid frames
    'Baseline'     : valid frames where none of `baseline_exclude_cols` is active
                     (default: zone_cols + behav_cols), or `baseline_col` if given
    each zone_col / behav_col : frames where that column is active
                     (zone masks also include merged behaviours, as in extract_dff_summary)
    '{stim} no sniff' : stim active AND matching sniff inactive (stim_sniff_pairs, {stim: sniff})

    A transient belongs to a period if its PEAK falls inside it. Periods can
    overlap (a sniff inside a stim counts in both), which is why the events table
    has one boolean column per period instead of a single label.

    Returns
    -------
    record    : dict, one row per mouse, columns '{period} N transients',
                '{period} duration (min)', '{period} frequency (per min)',
                '{period} mean amplitude', '{period} mean width (s)'
    events_df : one row per transient, with a boolean column per period
    """
    zone_cols = zone_cols or []
    behav_cols = behav_cols or []
    merge_into = merge_into or {}
    stim_sniff_pairs = stim_sniff_pairs or {}

    if dff_col not in fiberbehav_df.columns:
        raise ValueError(f"Column '{dff_col}' not found in dataframe.")

    n = len(fiberbehav_df)
    sig = fiberbehav_df[dff_col].values.astype(float)
    valid = np.isfinite(sig)
    if exclusion_col in fiberbehav_df.columns:
        valid &= ~fiberbehav_df[exclusion_col].values.astype(bool)

    time = (fiberbehav_df['Time(s)'].values.astype(float)
            if 'Time(s)' in fiberbehav_df.columns else np.arange(n) / fps)

    def _col(c):
        return fiberbehav_df[c].values.astype(bool)

    # ── Period masks ──────────────────────────────────────────────────────────
    periods = {'Whole trace': np.ones(n, dtype=bool)}

    if baseline_col is not None:
        periods['Baseline'] = _col(baseline_col)
    else:
        excl = baseline_exclude_cols if baseline_exclude_cols is not None else zone_cols + behav_cols
        base = np.ones(n, dtype=bool)
        for c in excl:
            if c in fiberbehav_df.columns:
                base &= ~_col(c)
            else:
                print(f"  Warning: '{c}' not found, ignored for baseline mask.")
        periods['Baseline'] = base

    for c in zone_cols:
        if c not in fiberbehav_df.columns:
            print(f"  Warning: zone column '{c}' not found."); continue
        m = _col(c)
        for behav, target in merge_into.items():
            if target == c and behav in fiberbehav_df.columns:
                m = m | _col(behav)
        periods[c] = m

    for c in behav_cols:
        if c not in fiberbehav_df.columns:
            print(f"  Warning: behaviour column '{c}' not found."); continue
        periods[c] = _col(c)

    for stim, sniff in stim_sniff_pairs.items():
        if stim in fiberbehav_df.columns and sniff in fiberbehav_df.columns:
            periods[f'{stim} no sniff'] = _col(stim) & ~_col(sniff)
        else:
            print(f"  Warning: '{stim}' or '{sniff}' not found.")

    record = {'Mouse': mouse, 'Batch': batch, 'Group': group}
    event_cols = ['Mouse', 'Batch', 'Group', 'Time(s)', 'Amplitude', 'Width (s)', 'Peak height'] + list(periods)
    empty = pd.DataFrame(columns=event_cols)

    def _fill_nan():
        for p in periods:
            record[f'{p} N transients'] = np.nan
            record[f'{p} duration (min)'] = np.nan
            record[f'{p} frequency (per min)'] = np.nan
            record[f'{p} mean amplitude'] = np.nan
            record[f'{p} mean width (s)'] = np.nan

    if valid.sum() < 2:
        _fill_nan()
        return record, empty

    # ── Detection on the whole trace ──────────────────────────────────────────
    median_level = np.median(sig[valid])
    sig_f = np.where(valid, sig, median_level)
    n_smooth = int(round(smooth_s * fps))
    if n_smooth > 1:
        sig_f = uniform_filter1d(sig_f, size=n_smooth)

    sigma = 1.4826 * np.median(np.abs(sig_f[valid] - median_level))

    if threshold_mode == 'mad':
        prominence = mad_k * sigma
    elif threshold_mode == 'absolute':
        if abs_prominence is None:
            raise ValueError("abs_prominence must be set when threshold_mode='absolute'.")
        prominence = abs_prominence
    else:
        raise ValueError("threshold_mode must be 'mad' or 'absolute'.")

    peaks, props = find_peaks(sig_f, prominence=prominence,
                              distance=max(1, int(round(min_distance_s * fps))))

    amplitudes = np.array([]); widths_s = np.array([])
    width_heights = np.array([]); left_t = np.array([]); right_t = np.array([])

    if len(peaks) > 0:
        widths, w_heights, left_ips, right_ips = peak_widths(
            sig_f, peaks, rel_height=rel_height,
            prominence_data=(props['prominences'], props['left_bases'], props['right_bases']))
        widths_all_s = widths / fps

        keep = np.zeros(len(peaks), dtype=bool)
        for k, p in enumerate(peaks):
            lb, rb = props['left_bases'][k], props['right_bases'][k]
            keep[k] = valid[lb:rb + 1].all() and widths_all_s[k] >= min_width_s

        peaks         = peaks[keep]
        amplitudes    = props['prominences'][keep]
        widths_s      = widths_all_s[keep]
        width_heights = w_heights[keep]
        idx = np.arange(n)
        left_t  = np.interp(left_ips[keep],  idx, time)
        right_t = np.interp(right_ips[keep], idx, time)

    # ── QC plot ───────────────────────────────────────────────────────────────
    if plot:
        plot_transients(time, sig_f, valid, peaks, amplitudes, width_heights,
                        left_t, right_t, median_level, prominence,
                        title=f'{batch} {mouse} ({group}) — {dff_col} — thr={prominence:.3f}',
                        zoom_start_s=zoom_start_s, zoom_len_s=zoom_len_s,
                        save_path=plot_path)

    # ── Metrics per period ────────────────────────────────────────────────────
    for p, m in periods.items():
        dur_min = (valid & m).sum() / fps / 60
        in_p = m[peaks] if len(peaks) else np.array([], dtype=bool)
        n_p = int(in_p.sum())
        record[f'{p} N transients']        = n_p
        record[f'{p} duration (min)']      = round(float(dur_min), 6)
        record[f'{p} frequency (per min)'] = round(n_p / dur_min, 6) if dur_min > 0 else np.nan
        record[f'{p} mean amplitude']      = round(float(np.mean(amplitudes[in_p])), 6) if n_p else np.nan
        record[f'{p} mean width (s)']      = round(float(np.mean(widths_s[in_p])), 6) if n_p else np.nan

    record['Noise sigma (MAD)']    = round(float(sigma), 6)
    record['Prominence threshold'] = round(float(prominence), 6)

    if len(peaks) == 0:
        return record, empty

    events_df = pd.DataFrame({
        'Mouse': mouse, 'Batch': batch, 'Group': group,
        'Time(s)': time[peaks],
        'Amplitude': amplitudes,
        'Width (s)': widths_s,
        'Peak height': sig[peaks],
    })
    for p, m in periods.items():
        events_df[p] = m[peaks]
    return record, events_df

def compute_PETH_bybout_metrics(PETH_list, subject_list, group_list,
                                 timewindow, BOI, step=1, max_bouts=None,
                                 dff_column='465'):
    """
    Compute mean and max dFF before and after the event for each bout group,
    per mouse. Returns a long-format DataFrame (one row per mouse × bout group).

    Parameters
    ----------
    PETH_list    : list of np.ndarray, each shape (n_bouts_i, timepoints)
    subject_list : list of str
    group_list   : list of str
    timewindow   : [PRE_TIME, POST_TIME]
    BOI          : str
    step         : int  — bouts per group (must match the value used in PETH_by_bout)
    max_bouts    : int, optional
    dff_column   : str  — channel label used in column names ('465' or '560')
    """
    PRE_TIME, POST_TIME = float(timewindow[0]), float(timewindow[1])

    if max_bouts is None:
        max_bouts = max(len(p) for p in PETH_list)

    n_positions = max_bouts // step
    timepoints  = PETH_list[0].shape[1]

    # Index of time-0 in the linspace(-PRE_TIME, POST_TIME, timepoints) vector
    pre_idx = round(PRE_TIME * (timepoints - 1) / (PRE_TIME + POST_TIME))

    records = []

    for mouse, group, peth in zip(subject_list, group_list, PETH_list):
        for pos_idx in range(n_positions):
            bout_indices = range(pos_idx * step, (pos_idx + 1) * step)

            # Collect bouts this mouse actually has within this group
            traces = [peth[b] for b in bout_indices if len(peth) > b]
            if not traces:
                continue

            # Average across bouts in the group for this mouse → (timepoints,)
            mouse_trace = np.mean(traces, axis=0)

            pre  = mouse_trace[:pre_idx]
            post = mouse_trace[pre_idx:]

            if step == 1:
                bout_label = f'Bout #{pos_idx + 1}'
            else:
                start = pos_idx * step + 1
                end   = (pos_idx + 1) * step
                bout_label = f'Bouts {start}–{end}'

            records.append({
                'Subject'   : mouse,
                'Group'     : group,
                'Bout_group': bout_label,
                'Bout_group_idx': pos_idx + 1,
                f'{dff_column} Mean dFF before {BOI}': float(np.mean(pre)),
                f'{dff_column} Mean dFF after {BOI}' : float(np.mean(post)),
                f'{dff_column} Max dFF before {BOI}' : float(np.max(pre)),
                f'{dff_column} Max dFF after {BOI}'  : float(np.max(post)),
            })

    return pd.DataFrame(records)