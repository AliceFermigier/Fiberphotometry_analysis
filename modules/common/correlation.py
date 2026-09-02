from scipy.signal import correlate, correlation_lags, fftconvolve, savgol_filter
import numpy as np
import matplotlib.pyplot as plt
from statsmodels.tsa.stattools import grangercausalitytests
from scipy.stats import combine_pvalues
import pandas as pd
import importlib

#import functions
import modules.common.preprocess as pp
importlib.reload(pp)
import modules.common.behavplot as bp
importlib.reload(bp)

def compute_peth_crosscorr(peth_465_list, peth_560_list, sr,
                            max_lag_s=5, min_bouts=1):
    """
    Compute normalized cross-correlation between 465nm and 560nm PETH traces,
    pooled across bouts and mice.

    Parameters
    ----------
    peth_465_list, peth_560_list : list of np.ndarray, shape (n_bouts, timepoints)
        Matched per-mouse PETH arrays for the two channels.
    sr : float
        Sampling rate (Hz).
    max_lag_s : float
        Maximum lag to return (seconds).
    min_bouts : int
        Minimum number of bouts a mouse must have to be included.

    Returns
    -------
    lags_s      : np.ndarray  — lag axis in seconds
    mean_xcorr  : np.ndarray  — mean cross-correlation
    sem_xcorr   : np.ndarray  — SEM across bouts
    peak_lag_s  : float       — lag at peak correlation
    xcorr_matrix: np.ndarray  — all individual bout cross-correlations (n_bouts_total, lags)
    valid_mouse_idx : list of int — indices into the input lists that were kept
    """

    n_tp      = peth_465_list[0].shape[1]
    lags_s    = correlation_lags(n_tp, n_tp, mode='full') / sr
    mask      = np.abs(lags_s) <= max_lag_s
    lags_trim = lags_s[mask]

    per_mouse = []
    valid_mouse_idx = []
    for idx, (peth_465, peth_560) in enumerate(zip(peth_465_list, peth_560_list)):
        n_bouts = min(len(peth_465), len(peth_560))
        if n_bouts < min_bouts:
            continue
        bout_xcorrs = []
        for s1, s2 in zip(peth_465, peth_560):
            s1 = s1 - s1.mean();  s2 = s2 - s2.mean()
            norm = np.sqrt(np.dot(s1, s1) * np.dot(s2, s2))
            if norm < 1e-10:
                continue
            bout_xcorrs.append(correlate(s2, s1, mode='full')[mask] / norm)
        if bout_xcorrs:
            per_mouse.append(np.mean(bout_xcorrs, axis=0))
            valid_mouse_idx.append(idx)

    if not per_mouse:
        raise ValueError("No valid mice found.")

    per_mouse_arr = np.array(per_mouse)
    mean_xcorr    = per_mouse_arr.mean(axis=0)
    sem_xcorr     = per_mouse_arr.std(axis=0) / np.sqrt(len(per_mouse_arr))
    peak_lag_s    = float(lags_trim[np.argmax(mean_xcorr)])

    return lags_trim, mean_xcorr, sem_xcorr, peak_lag_s, per_mouse_arr, valid_mouse_idx


def compute_crosscorr_significance(per_mouse_arr, lags_s,
                                    peth_465_list, peth_560_list, sr,
                                    method='phase_randomization',
                                    n_shuffles=1000, ci=95):

    mean_xcorr = per_mouse_arr.mean(axis=0)
    sem_xcorr  = per_mouse_arr.std(axis=0) / np.sqrt(len(per_mouse_arr))
    peak_lag_s = float(lags_s[np.argmax(mean_xcorr)])

    # Derive the slice indices once from lags_s length — avoids all float comparison
    n_tp   = peth_465_list[0].shape[1]
    center = n_tp - 1              
    half_n = len(lags_s) // 2     

    def _mouse_mean_xcorr(p465, p560):
        xcorrs = []
        for s1, s2 in zip(p465, p560):
            s1 = s1 - s1.mean();  s2 = s2 - s2.mean()
            norm = np.sqrt(np.dot(s1, s1) * np.dot(s2, s2))
            if norm < 1e-10:
                continue
            xcorr_full = correlate(s2, s1, mode='full') / norm
            # Slice exactly len(lags_s) elements centred at lag=0
            xcorrs.append(xcorr_full[center - half_n : center + half_n + 1])
        return np.mean(xcorrs, axis=0) if xcorrs else None

    shuffle_group_means = []
    for _ in range(n_shuffles):
        shuf_mouse_means = []
        for p465, p560 in zip(peth_465_list, peth_560_list):
            if method == 'phase_randomization':
                p560_shuf = np.array([_phase_randomize(row) for row in p560])
            elif method == 'trial_permutation':
                p560_shuf = p560[np.random.permutation(len(p560))]
            else:
                raise ValueError(f"Unknown method '{method}'.")
            m = _mouse_mean_xcorr(p465, p560_shuf)
            if m is not None:
                shuf_mouse_means.append(m)
        if shuf_mouse_means:
            shuffle_group_means.append(np.mean(shuf_mouse_means, axis=0))

    shuffle_arr = np.array(shuffle_group_means)
    tail_pct    = (100 - ci) / 2
    ci_low      = np.percentile(shuffle_arr, tail_pct,       axis=0)
    ci_high     = np.percentile(shuffle_arr, 100 - tail_pct, axis=0)
    is_sig      = (mean_xcorr > ci_high) | (mean_xcorr < ci_low)

    return mean_xcorr, sem_xcorr, peak_lag_s, ci_low, ci_high, is_sig

def _phase_randomize(signal):
    """
    Destroy cross-signal temporal structure while preserving the
    signal's own power spectrum (autocorrelation).
    """
    f   = np.fft.rfft(signal)
    phi = np.random.uniform(0, 2 * np.pi, len(f))
    phi[0] = 0                              # keep DC real
    if len(signal) % 2 == 0:
        phi[-1] = 0                         # keep Nyquist real for even-length signals
    return np.fft.irfft(np.abs(f) * np.exp(1j * phi), n=len(signal))

def compute_baseline_crosscorr(fiberbehav_df, behaviours_excluded_list,
                                sr, pad_s=5, max_lag_s=1,
                                min_segment_s=5,
                                exclusion_col='dFF ExclusionMask'):
    """
    Compute cross-correlation between 465nm and 560nm during baseline periods,
    excluding all behavioral bouts (with padding) and corrupted regions.

    Parameters
    ----------
    behaviours_excluded_list : list of str
        Derived behavior columns (1=onset, -1=offset) to exclude.
    pad_s : float
        Padding in seconds around each excluded bout.
    min_segment_s : float
        Minimum length of a contiguous baseline segment to use.
    """
    pad_samples     = int(pad_s * sr)
    min_seg_samples = int(min_segment_s * sr)
    n               = len(fiberbehav_df)
    baseline_mask   = np.ones(n, dtype=bool)

    # ── Exclude each behavior bout + padding ──────────────────────────────────
    for behav in behaviours_excluded_list:
        if behav not in fiberbehav_df.columns:
            print(f"  [!] '{behav}' not found, skipping.")
            continue
        onsets  = fiberbehav_df.index[fiberbehav_df[behav] == 1].tolist()
        offsets = fiberbehav_df.index[fiberbehav_df[behav] == -1].tolist()
        for onset, offset in zip(onsets, offsets):
            start = max(0, onset - pad_samples)
            end   = min(n - 1, offset + pad_samples)
            baseline_mask[start : end + 1] = False

    # ── Exclude corrupted regions ─────────────────────────────────────────────
    if exclusion_col in fiberbehav_df.columns:
        corrupted = fiberbehav_df[exclusion_col].fillna(False).astype(bool).values
        baseline_mask &= ~corrupted

    print(f"  Baseline: {baseline_mask.sum()} / {n} frames retained "
          f"({100 * baseline_mask.sum() / n:.1f}%)")

    # ── Find contiguous baseline segments ─────────────────────────────────────
    segments = []
    in_seg, start = False, 0
    for i, val in enumerate(baseline_mask):
        if val and not in_seg:
            start, in_seg = i, True
        elif not val and in_seg:
            if i - start >= min_seg_samples:
                segments.append((start, i))
            in_seg = False
    if in_seg and n - start >= min_seg_samples:
        segments.append((start, n))

    print(f"  Found {len(segments)} baseline segments ≥ {min_segment_s}s")
    if not segments:
        print("  [!] No usable baseline segments.")
        return None, None

    # ── Cross-correlate within each segment, then average ────────────────────
    xcorrs = []
    for start, end in segments:
        seg = fiberbehav_df.iloc[start:end][['dFF', '560 dFF']].dropna()
        if len(seg) < min_seg_samples:
            continue
        s1 = seg['dFF'].values;      s1 -= s1.mean()
        s2 = seg['560 dFF'].values;  s2 -= s2.mean()
        norm = np.sqrt(np.dot(s1, s1) * np.dot(s2, s2))
        if norm < 1e-10:
            continue
        xc   = correlate(s2, s1, mode='full') / norm
        lags = correlation_lags(len(s1), len(s1), mode='full') / sr
        mask = np.abs(lags) <= max_lag_s
        xcorrs.append(xc[mask])
        lags_trimmed = lags[mask]

    if not xcorrs:
        return None, None

    mean_xcorr = np.mean(xcorrs, axis=0)
    sem_xcorr  = np.std(xcorrs,  axis=0) / np.sqrt(len(xcorrs))
    return lags_trimmed, mean_xcorr, sem_xcorr

def plot_crosscorr_with_significance(lags_s, mean_xcorr, sem_xcorr,
                                      peak_lag_s, ci_low, ci_high, is_sig,
                                      BOI, exp, group, n_bouts,
                                      color='cornflowerblue', sig_style='bar',
                                      baseline_xcorr=None,   
                                      baseline_sem=None):    
    
    fig, ax = plt.subplots(figsize=(8, 4))

    ax.fill_between(lags_s, ci_low, ci_high,
                    color='grey', alpha=0.25, label='Shuffle null (95% CI)')
    ax.fill_between(lags_s, mean_xcorr - sem_xcorr, mean_xcorr + sem_xcorr,
                    color=color, alpha=0.3)
    ax.plot(lags_s, mean_xcorr, color=color, linewidth=2,
            label=f'Event xcorr  (n={n_bouts} mice)')

    # ── Baseline overlay ──────────────────────────────────────────────────────
    if baseline_xcorr is not None:
        if baseline_sem is not None:
            ax.fill_between(lags_s,
                            baseline_xcorr - baseline_sem,
                            baseline_xcorr + baseline_sem,
                            color='orange', alpha=0.2)
        ax.plot(lags_s, baseline_xcorr, color='orange', linewidth=1.5,
                linestyle='--', label='Baseline xcorr')

    if sig_style == 'bar':
        ax.fill_between(lags_s, 0.97, 1.0, where=is_sig,
                        transform=ax.get_xaxis_transform(),
                        color='firebrick', alpha=0.85, linewidth=0, label='p < 0.05')
    elif sig_style == 'overlay':
        sig_xcorr = np.where(is_sig, mean_xcorr, np.nan)
        ax.plot(lags_s, sig_xcorr, color='firebrick', linewidth=4,
                solid_capstyle='round', alpha=0.75, zorder=4, label='p < 0.05')

    ax.axvline(x=0, color='black', linewidth=1, linestyle='--')
    ax.axvline(x=peak_lag_s, color='firebrick', linewidth=1.5, linestyle=':',
               label=f'Peak lag = {peak_lag_s:.3f} s')
    ax.axhline(y=0, color='grey', linewidth=0.8, linestyle=':')
    ax.set_xlabel('Lag (s)', fontsize=12)
    ax.set_ylabel('Normalized cross-correlation', fontsize=12)
    ax.set_title(f'465→560 cross-correlation — {BOI}, {exp}, {group}', fontsize=12)
    ax.legend(fontsize=9)
    ax.margins(0, 0.05)
    plt.tight_layout()
    return fig

def deconvolve_rgeco(signal, sr, tau_rise=0.5, tau_decay=2.0):
    """
    Remove R-GECO indicator kinetics from a 560nm dFF trace by
    Wiener deconvolution with a double-exponential impulse response.

    Parameters
    ----------
    signal    : np.ndarray  — raw 560nm dFF trace
    sr        : float       — sampling rate (Hz)
    tau_rise  : float       — R-GECO rise time constant (s), default 0.5
    tau_decay : float       — R-GECO decay time constant (s), default 2.0

    Returns
    -------
    deconv : np.ndarray  — deconvolved trace (same length)
    """
    from numpy.fft import fft, ifft

    t = np.arange(len(signal)) / sr
    # Double-exponential kernel
    kernel = (np.exp(-t / tau_decay) - np.exp(-t / tau_rise))
    kernel = np.where(kernel < 0, 0, kernel)
    kernel /= kernel.sum()

    # Wiener deconvolution (SNR regularisation avoids noise explosion)
    SNR = 10.0
    H   = fft(kernel, n=len(signal))
    S   = fft(signal)
    deconv = np.real(ifft(S * np.conj(H) / (np.abs(H)**2 + 1 / SNR)))
    return deconv

def test_granger_causality(fiberbehav_df,behavior_col,max_lag_s=1,sr=None,
    alpha=0.05,exclude_corrupted=True,exclusion_col='dFF ExclusionMask',exclusion_padding_s=1,
    min_bout_s=2,min_std=0.01,return_raw=False):
    """
    Test whether dFF Granger-causes 560 dFF within behavior bouts.
    Robust version:
    ----------------
    - rejects bouts overlapping excluded regions
    - rejects bouts too close to excluded regions
    - rejects short bouts
    - rejects low-variance bouts
    - combines p-values using Fisher's method
    Parameters
    ----------
    fiberbehav_df : pd.DataFrame
    behavior_col : str
        Behavioral column containing:
            1  = onset
           -1  = offset
    max_lag_s : float 
        Maximum lag tested in seconds.
    exclusion_padding_s : float
        Reject bouts within ±padding seconds of excluded regions.
    min_bout_s : float
        Minimum bout duration required.
    min_std : float
        Minimum signal standard deviation required.
    Returns
    -------
    pd.DataFrame
        lag_s
        mean_F
        combined_p
        significant
        n_bouts
    """
    if sr is None:
        sr = round(pp.samplerate(fiberbehav_df))

    max_lag_samples = int(max_lag_s * sr)
    exclusion_padding = int(exclusion_padding_s * sr)
    min_bout_samples = int(min_bout_s * sr)

    # FIND BOUTS
    onsets = fiberbehav_df.index[fiberbehav_df[behavior_col] == 1].tolist()
    offsets = fiberbehav_df.index[fiberbehav_df[behavior_col] == -1].tolist()

    lag_results = {}
    for lag in range(1, max_lag_samples + 1):
        lag_results[lag] = {
            'F': [],
            'p': []
        }
    accepted_bouts = 0
    rejected_bouts = 0

    for onset, offset in zip(onsets, offsets):

        window_pre_s = 2
        window_post_s = 2
        start_idx = max(0, onset - int(window_pre_s * sr))
        end_idx   = min(len(fiberbehav_df) - 1,offset + int(window_post_s * sr))
        bout_len = end_idx - start_idx

        if bout_len < min_bout_samples:
            rejected_bouts += 1
            continue

        # EXCLUSION MASK CHECK
        if exclude_corrupted and exclusion_col in fiberbehav_df.columns:
            padded_start = max(0, start_idx - exclusion_padding)
            padded_end = min(len(fiberbehav_df) - 1,end_idx + exclusion_padding)
            exclusion_window = fiberbehav_df.loc[padded_start:padded_end,exclusion_col]
            if exclusion_window.any():
                rejected_bouts += 1
                continue

        # EXTRACT SEGMENT
        seg = fiberbehav_df.loc[start_idx:end_idx,['dFF', '560 dFF']].dropna()

        # enough data for Granger
        if len(seg) < max_lag_samples * 3:
            rejected_bouts += 1
            continue

        # LOW VARIANCE REJECTION
        if seg['dFF'].std() < min_std:
            rejected_bouts += 1
            continue
        if seg['560 dFF'].std() < min_std:
            rejected_bouts += 1
            continue

        # GRANGER TEST
        try:
            gc = grangercausalitytests(
                seg.values,
                maxlag=max_lag_samples,
                verbose=False)
            for lag, res in gc.items():
                f_stat = res[0]['ssr_ftest'][0]
                p_val = res[0]['ssr_ftest'][1]

                lag_results[lag]['F'].append(f_stat)
                lag_results[lag]['p'].append(p_val)

            accepted_bouts += 1
        except Exception:
            rejected_bouts += 1
            continue

    # COMBINE RESULTS
    output = []
    for lag in lag_results:
        F_vals = lag_results[lag]['F']
        p_vals = lag_results[lag]['p']

        if len(F_vals) == 0:
            continue
        # Fisher method for p-values
        combined_p = combine_pvalues(p_vals,method='fisher')[1]

        output.append({

            'lag_samples': lag,
            'lag_s': lag / sr,

            'mean_F': np.mean(F_vals),
            'std_F': np.std(F_vals),

            'combined_p': combined_p,
            'significant': combined_p < alpha,
 

            'accepted_bouts': accepted_bouts,
            'rejected_bouts': rejected_bouts
        })
        
    if return_raw:
        return pd.DataFrame(output), lag_results, accepted_bouts, rejected_bouts
    return pd.DataFrame(output)

def compute_joint_psth(peth_465, peth_560):
    n_bouts, n_tp = peth_465.shape

    jpsth_raw = sum(
        np.outer(peth_465[i], peth_560[i]) for i in range(n_bouts)
    ) / n_bouts

    mean_465  = peth_465.mean(axis=0)
    mean_560  = peth_560.mean(axis=0)
    predictor = np.outer(mean_465, mean_560)

    # Normalisation: correlation-coefficient scale
    predictor_std = np.outer(
        peth_465.std(axis=0),
        peth_560.std(axis=0)
    )
    predictor_std = np.where(predictor_std < 1e-10, np.nan, predictor_std)

    jpsth_corrected  = jpsth_raw - predictor
    jpsth_normalized = jpsth_corrected/predictor_std

    coincidence = np.diag(jpsth_normalized)
    coincidence_raw = np.diag(jpsth_raw)
    predictor_diag = np.diag(predictor)

    return jpsth_raw, predictor, jpsth_normalized, coincidence, coincidence_raw, predictor_diag

def compute_baseline_jpsth(fiberbehav_df, behaviours_excluded_list,
                             sr, timewindow, pad_s=2,
                             min_segment_s=None,
                             exclusion_col='dFF ExclusionMask'):
    """
    Compute a JPSTH from randomly sampled baseline windows,
    matched in length to the PETH time window.
    Used as a noise-floor baseline to subtract from event-aligned JPS THs.
    """
    from itertools import islice

    PRE_TIME, POST_TIME = float(timewindow[0]), float(timewindow[1])
    win_samples = int((PRE_TIME + POST_TIME) * sr) + 1
    pad_samples = int(pad_s * sr)
    if min_segment_s is None:
        min_segment_s = PRE_TIME + POST_TIME + 2 * pad_s
    min_seg_samples = int(min_segment_s * sr)
    n = len(fiberbehav_df)

    # Build exclusion mask
    baseline_mask = np.ones(n, dtype=bool)
    for behav in behaviours_excluded_list:
        if behav not in fiberbehav_df.columns:
            continue
        onsets  = fiberbehav_df.index[fiberbehav_df[behav] == 1].tolist()
        offsets = fiberbehav_df.index[fiberbehav_df[behav] == -1].tolist()
        for onset, offset in zip(onsets, offsets):
            baseline_mask[max(0, onset - pad_samples):
                          min(n, offset + pad_samples + 1)] = False
    if exclusion_col in fiberbehav_df.columns:
        baseline_mask &= ~fiberbehav_df[exclusion_col].fillna(False).astype(bool).values

    # Find contiguous segments
    segments, in_seg, start = [], False, 0
    for i, val in enumerate(baseline_mask):
        if val and not in_seg:
            start, in_seg = i, True
        elif not val and in_seg:
            if i - start >= min_seg_samples:
                segments.append((start, i))
            in_seg = False
    if in_seg and n - start >= min_seg_samples:
        segments.append((start, n))

    # Extract windows from baseline segments
    sig_465 = fiberbehav_df['dFF'].values
    sig_560 = fiberbehav_df['560 dFF'].values

    windows_465, windows_560 = [], []
    for start, end in segments:
        for t in range(start, end - win_samples, win_samples // 2):
            w465 = sig_465[t : t + win_samples]
            w560 = sig_560[t : t + win_samples]
            if np.any(np.isnan(w465)) or np.any(np.isnan(w560)):
                continue
            # Z-score each window
            std465 = w465.std()
            std560 = w560.std()
            if std465 < 1e-10 or std560 < 1e-10:
                continue
            windows_465.append((w465 - w465.mean()) / std465)
            windows_560.append((w560 - w560.mean()) / std560)

    if not windows_465:
        print("  [!] No usable baseline windows found.")
        return None, None

    peth_465_bl = np.stack(windows_465)
    peth_560_bl = np.stack(windows_560)
    _, _, jpsth_bl, coinc_bl, _, _ = compute_joint_psth(peth_465_bl, peth_560_bl)

    print(f"  Baseline JPSTH from {len(windows_465)} windows.")
    return jpsth_bl, coinc_bl

def stationary_bootstrap(signal, average_block_size=120):
    """
    Resample a time series using stationary bootstrap (Politis & Romano 1994).
    Blocks of geometrically distributed length are resampled with replacement,
    preserving autocorrelation at lags < average_block_size.
    NaNs are temporarily filled with linear interpolation to avoid propagation.
    """
    signal  = np.asarray(signal, dtype=float)
    nan_mask = np.isnan(signal)

    # Temporarily fill NaNs so they don't contaminate blocks
    if nan_mask.any():
        filled = pd.Series(signal).interpolate(
            method='linear', limit_direction='both'
        ).values
    else:
        filled = signal.copy()

    l = len(filled)
    p = 1 / average_block_size
    resampled = np.empty_like(filled)
    j = np.random.randint(l)
    resampled[0] = filled[j]
    for i in range(1, l):
        if np.random.random() < 1 - p:
            j = (j + 1) % l
        else:
            j = np.random.randint(l)
        resampled[i] = filled[j]

    # Re-apply the original NaN mask
    resampled[nan_mask] = np.nan
    return resampled

def compute_shuffle_jpsth(dfiberbehav_df, peth_465_real,
                           BOI, event, timewindow,
                           behav_cols,
                           event_time_threshold, sr,
                           n_shuffles=100, baseline=False,
                           sig_to_shuffle = '560 dFF',
                           maxboutsnumber=None,
                           exclusion_col='dFF ExclusionMask'):
    """
    Build a shuffle null distribution for the JPSTH by bootstrapping
    the 560nm signal before PETH extraction, then compute the z-score JPSTH.

    Parameters
    ----------
    dfiberbehav_df   : pd.DataFrame  — per-mouse aligned dataframe
    peth_465_real    : np.ndarray (n_bouts, timepoints) — real 465 PETH
    n_shuffles       : int           — number of bootstrappings
    exclusion_col    : str           — column marking corrupted frames

    Returns
    -------
    jpsth_shuffle_mean : np.ndarray (n_tp, n_tp)
    jpsth_shuffle_std  : np.ndarray (n_tp, n_tp)  — floored at 1e-10 → nan
    """
    n            = len(dfiberbehav_df)

    sig_560_orig = dfiberbehav_df[sig_to_shuffle].values.copy()

    if exclusion_col in dfiberbehav_df.columns:
        excl_mask = dfiberbehav_df[exclusion_col].fillna(False).astype(bool).values
    else:
        excl_mask = np.zeros(n, dtype=bool)

    n_bouts_real = len(peth_465_real)
    jpsth_shuffles = []

    for _ in range(n_shuffles):
        sig_shifted = stationary_bootstrap(sig_560_orig,average_block_size=3*sr) #bootstrapping in 3s blocks
        sig_shifted[excl_mask] = np.nan          

        df_shuf = dfiberbehav_df.copy()
        df_shuf[sig_to_shuffle] = sig_shifted

        peth_560_shuf = bp.PETH(
            df_shuf, BOI, event, timewindow,
            behav_cols,
            baselinewindow=baseline, maxboutsnumber=maxboutsnumber,
            dFF_column=sig_to_shuffle,
            baseline_start_stop_s=[timewindow[0],1.0],
            baseline_method='median'
        )

        # Keep only as many bouts as the real PETH has
        n_bouts = min(n_bouts_real, len(peth_560_shuf))
        if n_bouts == 0:
            continue

        p465 = peth_465_real[:n_bouts]
        p560 = peth_560_shuf[:n_bouts]

        jpsth_shuf = sum(
            np.outer(p465[i], p560[i]) for i in range(n_bouts)
        ) / n_bouts
        jpsth_shuffles.append(jpsth_shuf)

    if len(jpsth_shuffles) < n_shuffles // 2:
        print(f"  [!] Only {len(jpsth_shuffles)}/{n_shuffles} valid shuffles.")

    stack = np.stack(jpsth_shuffles)                     # (n_valid, n_tp, n_tp)
    mean  = stack.mean(axis=0)
    std   = stack.std(axis=0)
    std   = np.where(std < 1e-10, np.nan, std)

    return mean, std

def _smooth_coincidence(coincidence, sr_eff, smooth_window_s=0.5):
    """Savitzky-Golay smoothing, returns smoothed array."""
    win = max(5, int(smooth_window_s * sr_eff))
    if win % 2 == 0:
        win += 1                          # must be odd
    win = min(win, len(coincidence) - 1)  # must be < signal length
    if win < 5:
        return coincidence.copy()
    return savgol_filter(coincidence, window_length=win, polyorder=3)

def extract_coincidence_dynamics(coincidence, timewindow,
                                  smooth_window_s=0.5,
                                  onset_fraction=0.2):
    PRE_TIME, POST_TIME = float(timewindow[0]), float(timewindow[1])
    n_tp      = len(coincidence)
    sr_eff    = (n_tp - 1) / (PRE_TIME + POST_TIME)
    peri_time = np.linspace(-PRE_TIME, POST_TIME, n_tp)
    pre_idx   = round(PRE_TIME * sr_eff)
    pre_idx_max = round((PRE_TIME+0.1) * sr_eff)

    smooth = _smooth_coincidence(coincidence, sr_eff, smooth_window_s)

    # ── 1. Maximum: search only from t=0 to end ───────────────────────────────
    post_smooth = smooth[pre_idx_max:]
    idx_max     = pre_idx_max + int(np.nanargmax(post_smooth))
    t_max       = float(peri_time[idx_max])
    val_max     = float(smooth[idx_max])

    # ── 2. Trough: search only from t=0 to peak ──────────────────────────────
    idx_trough = pre_idx + int(np.nanargmin(smooth[pre_idx : idx_max + 1]))

    # ── 3. Inflection point: max first derivative between trough and peak ─────
    deriv      = np.gradient(smooth, peri_time)
    rise_deriv = deriv[idx_trough : idx_max + 1]

    if len(rise_deriv) > 1:
        idx_inflection = idx_trough + int(np.argmax(rise_deriv))
        t_inflection   = float(peri_time[idx_inflection])
    else:
        t_inflection = np.nan

    # ── 4. Onset: threshold crossing on the rising flank ─────────────────────
    baseline  = float(np.nanmean(smooth[:pre_idx]))
    threshold = baseline + onset_fraction * (val_max - baseline)

    t_onset = np.nan
    for idx in range(idx_trough, idx_max + 1):
        if smooth[idx] >= threshold:
            t_onset = float(peri_time[idx])
            break

    return {
        't_max_coincidence' : t_max,
        'max_coincidence'   : val_max,
        't_inflection'      : t_inflection,
        't_onset'           : t_onset,
    }

def extract_coincidence_metrics(coincidence, timewindow, step_s=0.2, smooth_window_s=0.5, onset_fraction=0.2):
    """
    Extract mean and max co-activation from the coincidence diagonal.

    Returns global metrics (before / after / total) plus binned metrics
    of width step_s radiating outward from t=0 in both directions.

    Parameters
    ----------
    coincidence : np.ndarray (n_timepoints,)
    timewindow  : [PRE_TIME, POST_TIME]
    step_s      : float — bin width in seconds (default 0.2)

    Returns
    -------
    dict of scalar metrics, one row per mouse in the output DataFrame.
    """
    PRE_TIME, POST_TIME = float(timewindow[0]), float(timewindow[1])
    n_tp    = len(coincidence)
    sr_eff  = (n_tp - 1) / (PRE_TIME + POST_TIME)   # samples per second
    pre_idx = round(PRE_TIME * sr_eff)

    def t_to_idx(t):
        return int(np.clip(round(pre_idx + t * sr_eff), 0, n_tp - 1))

    metrics = {}

    # ── Global metrics (unchanged) ────────────────────────────────────────────
    pre  = coincidence[:pre_idx]
    post = coincidence[pre_idx:]
    metrics['mean_coincidence_before'] = float(np.nanmean(pre))
    metrics['max_coincidence_before']  = float(np.nanmax(pre))
    metrics['mean_coincidence_after']  = float(np.nanmean(post))
    metrics['max_coincidence_after']   = float(np.nanmax(post))
    metrics['mean_coincidence_total']  = float(np.nanmean(coincidence))
    metrics['max_coincidence_total']   = float(np.nanmax(coincidence))

    # ── Bins after event: [0, step_s], [step_s, 2*step_s], … ─────────────────
    bin_edges_after = np.arange(0, POST_TIME, step_s)
    for t0 in bin_edges_after:
        t1  = t0 + step_s
        i0  = t_to_idx(t0)
        i1  = n_tp if t1 >= POST_TIME else t_to_idx(t1)   # last bin reaches end
        if i1 <= i0:
            continue
        window = coincidence[i0:i1]
        label  = f'{t0:+.2f}_to_{t1:+.2f}s'
        metrics[f'mean_coinc_{label}'] = float(np.nanmean(window))
        metrics[f'max_coinc_{label}']  = float(np.nanmax(window))

    # ── Bins before event: [0, -step_s], [-step_s, -2*step_s], … ────────────
    bin_edges_before = np.arange(0, -PRE_TIME, -step_s)
    for t1 in bin_edges_before:
        t0  = t1 - step_s
        i0  = 0 if t0 <= -PRE_TIME else t_to_idx(t0)      # last bin reaches start
        i1  = t_to_idx(t1)
        if i1 <= i0:
            continue
        window = coincidence[i0:i1]
        label  = f'{t0:+.2f}_to_{t1:+.2f}s'
        metrics[f'mean_coinc_{label}'] = float(np.nanmean(window))
        metrics[f'max_coinc_{label}']  = float(np.nanmax(window))

    # ── Dynamics metrics ──────────────────────────────────────────────────────
    dynamics = extract_coincidence_dynamics(
        coincidence, timewindow,
        smooth_window_s=smooth_window_s,
        onset_fraction=onset_fraction
    )
    metrics.update(dynamics)

    return metrics

def plot_joint_psth(jpsth_corrected, coincidence, timewindow,
                    BOI, event, exp, group,
                    n_bouts, cmap='RdBu_r', dff_column='465/560',
                    mouse=None,
                    vmin=None, vmax=None,           
                    coincidence_sem=None):         
    PRE_TIME, POST_TIME = float(timewindow[0]), float(timewindow[1])
    n_tp      = jpsth_corrected.shape[0]
    peri_time = np.linspace(-PRE_TIME, POST_TIME, n_tp)

    fig = plt.figure(figsize=(10, 9))
    gs  = fig.add_gridspec(1, 2, width_ratios=[10, 2], wspace=0.08)
    ax_jpsth = fig.add_subplot(gs[0, 0])
    ax_coinc = fig.add_subplot(gs[0, 1], sharey=ax_jpsth)

    # ── Colour scale ──────────────────────────────────────────────────────────
    if vmin is None and vmax is None:
        abs_max = np.nanmax(np.abs(jpsth_corrected))
        vmin, vmax = -abs_max, abs_max
    elif vmin is None:
        vmin = -vmax
    elif vmax is None:
        vmax = -vmin

    extent = [-PRE_TIME, POST_TIME, POST_TIME, -PRE_TIME]
    im = ax_jpsth.imshow(
        jpsth_corrected,
        cmap=cmap, aspect='auto', interpolation='bilinear',
        extent=extent, vmin=vmin, vmax=vmax
    )
    ax_jpsth.plot(peri_time, peri_time,
                  color='black', linewidth=1, linestyle='--', alpha=0.6,
                  label='Zero lag (t₁ = t₂)')
    ax_jpsth.axvline(x=0, color='white', linewidth=2, linestyle='--', alpha=0.8)
    ax_jpsth.axhline(y=0, color='white', linewidth=2, linestyle='--', alpha=0.8)
    ax_jpsth.set_xlabel('465nm time re. event (s)', fontsize=13)
    ax_jpsth.set_ylabel('560nm time re. event (s)', fontsize=13)
    ax_jpsth.set_title(
        f'Joint PSTH — {BOI} {event.capitalize()}\n{exp}, {mouse} {group}  (n={n_bouts} bouts)',
        fontsize=13
    )
    ax_jpsth.legend(loc='upper left', fontsize=9, framealpha=0.6)
    cbar = fig.colorbar(im, ax=ax_jpsth, fraction=0.046, pad=0.04)
    cbar.set_label(f'Normalised co-activation\n({dff_column} ΔF/F)', fontsize=10)

    # ── Coincidence panel ─────────────────────────────────────────────────────
    if coincidence_sem is not None:
        ax_coinc.fill_betweenx(
            peri_time,
            coincidence - coincidence_sem,
            coincidence + coincidence_sem,
            color='slategray', alpha=0.25
        )
    ax_coinc.plot(coincidence, peri_time, color='slategray', linewidth=1.5)
    ax_coinc.fill_betweenx(peri_time, 0, coincidence,
                            where=(coincidence > 0),
                            color='firebrick', alpha=0.4, label='Positive')
    ax_coinc.fill_betweenx(peri_time, 0, coincidence,
                            where=(coincidence < 0),
                            color='steelblue', alpha=0.4, label='Negative')
    ax_coinc.axhline(y=0, color='white', linewidth=1, linestyle=':')
    ax_coinc.axvline(x=0, color='black', linewidth=0.8, linestyle='--')
    ax_coinc.set_xlabel('Co-act.', fontsize=10)
    ax_coinc.set_title('Coinc.', fontsize=10)
    plt.setp(ax_coinc.get_yticklabels(), visible=False)

    plt.tight_layout()
    return fig

def plot_coincidence(coincidence, timewindow, BOI, event, exp, group,
                     n_bouts, ylim=None, coincidence_sem=None, zscore=False):
    """
    Standalone coincidence histogram with optional SEM shading.
    """
    PRE_TIME, POST_TIME = float(timewindow[0]), float(timewindow[1])
    peri_time = np.linspace(-PRE_TIME, POST_TIME, len(coincidence))

    fig, ax = plt.subplots(figsize=(7, 4))

    # Shuffle null reference bands (only meaningful for z-scored data)
    if zscore:
        try:
            ax.axhspan(-1, 1, color='grey',     alpha=0.12, label='Shuffle ±1 SD')
            ax.axhspan(-2, 2, color='grey',     alpha=0.07, label='Shuffle ±2 SD')
            ax.axhline(y= 2, color='grey', linewidth=0.8, linestyle='--')
            ax.axhline(y=-2, color='grey', linewidth=0.8, linestyle='--')
        except :
            print('Range too small to draw significance z-score lines')

    if coincidence_sem is not None:
        ax.fill_between(peri_time,
                        coincidence - coincidence_sem,
                        coincidence + coincidence_sem,
                        color='slategray', alpha=0.25)
    ax.plot(peri_time, coincidence, color='slategray', linewidth=2)
    ax.fill_between(peri_time, 0, coincidence,
                    where=(coincidence > 0),
                    color='firebrick', alpha=0.35, label='Positive')
    ax.fill_between(peri_time, 0, coincidence,
                    where=(coincidence < 0),
                    color='steelblue', alpha=0.35, label='Negative')
    ax.axvline(x=0, color='black', linewidth=1.5, linestyle='--',
               label=f'{event.capitalize()} {BOI}')
    ax.axhline(y=0, color='grey', linewidth=0.8, linestyle=':')

    ax.set_xlabel('Time re. event (s)', fontsize=12)
    ax.set_ylabel('Normalised co-activation', fontsize=12)
    ax.set_title(
        f'Coincidence — {BOI} {event.capitalize()}\n{exp}, {group}  (n={n_bouts} bouts)',
        fontsize=12
    )
    ax.legend(fontsize=10)
    ax.margins(0, 0.1)

    if ylim is not None:
        ax.set_ylim(ylim)
        
    plt.tight_layout()
    return fig

def plot_granger_results(gc_df, BOI, exp, group, alpha=0.05):
    """
    Two-panel publication figure for Granger causality results.

    Top    : mean F-statistic ± std across bouts, with individual bout
             F-values shown as a swarm to convey bout-to-bout variability.
    Bottom : −log₁₀(combined p-value) with significance threshold line.
    """
    fig, (ax_f, ax_p) = plt.subplots(2, 1, figsize=(7, 6), sharex=True,
                                       gridspec_kw={'hspace': 0.08})

    lags    = gc_df['lag_s'].values
    mean_f  = gc_df['mean_F'].values
    std_f   = gc_df['std_F'].values
    log_p   = -np.log10(gc_df['combined_p'].values)
    sig     = gc_df['significant'].values
    n_bouts = int(gc_df['accepted_bouts'].iloc[0])

    # ── F-statistic panel ─────────────────────────────────────────────────────
    ax_f.fill_between(lags, mean_f - std_f, mean_f + std_f,
                      color='cornflowerblue', alpha=0.25)
    ax_f.plot(lags, mean_f,
              color='cornflowerblue', linewidth=2,
              label=f'Mean ± SD  (n = {n_bouts} bouts)')
    ax_f.axhline(y=1, color='grey', linewidth=0.8, linestyle=':',
                 label='F = 1 (no effect)')

    # Highlight significant lags on the trace
    sig_f = np.where(sig, mean_f, np.nan)
    ax_f.plot(lags, sig_f,
              color='firebrick', linewidth=4,
              solid_capstyle='round', alpha=0.75,
              label=f'p < {alpha}  (Fisher)')

    ax_f.set_ylabel('F-statistic', fontsize=12)
    ax_f.legend(fontsize=9, loc='upper right')
    ax_f.margins(0, 0.15)
    plt.setp(ax_f.get_xticklabels(), visible=False)

    # ── −log₁₀(p) panel ──────────────────────────────────────────────────────
    ax_p.plot(lags, log_p,
              color='slategrey', linewidth=2)
    ax_p.fill_between(lags, 0, log_p,
                      where=sig,
                      color='firebrick', alpha=0.3)
    ax_p.axhline(y=-np.log10(alpha),
                 color='firebrick', linewidth=1.2, linestyle='--',
                 label=f'p = {alpha}')

    # Annotate the significant bar using the same trick as cross-correlation
    ax_p.fill_between(lags, 0.95, 1.0,
                      where=sig,
                      transform=ax_p.get_xaxis_transform(),
                      color='firebrick', alpha=0.85,
                      linewidth=0, label=f'p < {alpha}')

    ax_p.set_xlabel('Lag (s)', fontsize=12)
    ax_p.set_ylabel(r'$-\log_{10}(p)$', fontsize=12)
    ax_p.legend(fontsize=9, loc='upper right')
    ax_p.margins(0, 0.15)

    plt.suptitle(
        f'Granger causality: 465 nm → 560 nm\n{BOI} — {exp}, {group}',
        fontsize=12
    )
    plt.tight_layout()
    return fig