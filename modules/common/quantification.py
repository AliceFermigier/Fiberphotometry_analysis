##########
#IMPORTED#
##########

import pandas as pd
import numpy as np

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
                         shock_exclude_seconds=2):
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
    """
    merge_into = merge_into or {}
    shock_exclude_cols = shock_exclude_cols or []

    # ── Signal preparation ────────────────────────────────────────────────────
    if dff_col not in fiberbehav_df.columns:
        raise ValueError(f"Column '{dff_col}' not found in dataframe.")

    signal = fiberbehav_df[dff_col].values.copy()

    if use_zscore:
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

        baseline_signal = signal[baseline_mask]
        if baseline_signal.size < 2:
            raise ValueError(
                "Not enough baseline samples to compute z-score baseline "
                "(need at least 2)."
            )
        signal = (signal - np.median(baseline_signal)) / np.std(baseline_signal, ddof=1)

    record = {
        'Mouse' : mouse,
        'Batch' : batch,
        'Group' : group,
    }

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

    return record

def variance_transients(fiberbehav_df, list_BOI, mouse, group, exp, batch, threshold):
    """
    Calculates variance, transient frequency, and amplitude during whole trace and during behaviours
    
    Output:
    - A dataframe with variance, transient frequency and amplitude.
    """
    
    # Calculate variance during whole trace, baseline and post-baseline periods
    variance = np.var(fiberbehav_df['dFF'])
    
    # Calculate transients for whole trace, baseline and post-baseline periods
    peaks_df, peak_frequency, peak_amplitude, transients_fig = tr.transients(fiberbehav_df, threshold)

    # Store the results in a dataframe
    results_df = pd.DataFrame({
        'Batch': batch,
        'Subject': mouse,
        'Group': group,
        'Variance': variance,
        'Transients Frequency': peak_frequency,
        'Transients Amplitude': peak_amplitude,
    }, index=[0])

    if list_BOI != []:
        for behavior in list_BOI:
            if behavior not in fiberbehav_df.columns:
                print(f"Warning: Behavior '{behavior}' not found in DataFrame.")
                continue

            # Find rows where behavior is active
            behavior_mask = fiberbehav_df[behavior] == 1

            # Peaks during that behavior
            behavior_peaks = peaks_df[(peaks_df['Peaks'] == 1) & behavior_mask]

            # Duration of behavior period in seconds
            behavior_time = fiberbehav_df.loc[behavior_mask, 'Time(s)']
            if behavior_time.empty:
                duration = np.nan
            else:
                duration = behavior_time.iloc[-1] - behavior_time.iloc[0]

            # Calculate frequency and mean amplitude during behavior
            if duration and len(behavior_peaks) > 0:
                freq = len(behavior_peaks) / duration
                amp = behavior_peaks['Filtered dFF'].mean()
            else:
                freq = np.nan
                amp = np.nan

            # Add to results
            results_df[f'{behavior} Transients Frequency'] = freq
            results_df[f'{behavior} Transients Amplitude'] = amp
    
    return results_df, transients_fig

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