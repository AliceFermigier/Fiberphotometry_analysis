from scipy.signal import correlate, correlation_lags, fftconvolve
import numpy as np
import matplotlib.pyplot as plt
from statsmodels.tsa.stattools import grangercausalitytests
from scipy.stats import combine_pvalues
import pandas as pd
import importlib

#import functions
import modules.common.preprocess as pp
importlib.reload(pp)

def compute_peth_crosscorr(peth_465_list, peth_560_list, sr,
                            max_lag_s=5, min_bouts =1):
    """
    Compute normalized cross-correlation between 465nm and 560nm PETH traces,
    pooled across bouts and mice.

    Convention: positive lag means 560 follows 465 (consistent with ACh
    driving calcium with a delay).

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
    """
    all_xcorrs = []

    for peth_465, peth_560 in zip(peth_465_list, peth_560_list):
        n_bouts = min(len(peth_465), len(peth_560))
        if n_bouts < min_bouts:
            continue

        for i in range(n_bouts):
            sig1 = peth_465[i] - peth_465[i].mean()   # zero-mean
            sig2 = peth_560[i] - peth_560[i].mean()

            norm = np.sqrt(np.dot(sig1, sig1) * np.dot(sig2, sig2))
            if norm < 1e-10:
                continue

            # correlate(sig2, sig1): positive lag → sig2 (560) follows sig1 (465)
            xcorr = correlate(sig2, sig1, mode='full') / norm
            all_xcorrs.append(xcorr)

    if not all_xcorrs:
        raise ValueError("No valid bouts found.")

    n_tp  = peth_465_list[0].shape[1]
    lags  = correlation_lags(n_tp, n_tp, mode='full')
    lags_s = lags / sr

    mask          = np.abs(lags_s) <= max_lag_s
    lags_s        = lags_s[mask]
    xcorr_matrix  = np.array(all_xcorrs)[:, mask]

    mean_xcorr = xcorr_matrix.mean(axis=0)
    sem_xcorr  = xcorr_matrix.std(axis=0) / np.sqrt(len(xcorr_matrix))
    peak_lag_s = float(lags_s[np.argmax(mean_xcorr)])

    return lags_s, mean_xcorr, sem_xcorr, peak_lag_s, xcorr_matrix

def plot_peth_crosscorr(lags_s, mean_xcorr, sem_xcorr, peak_lag_s,
                         BOI, exp, group, n_bouts,
                         color='cornflowerblue', fill_alpha=0.25):
    """
    Plot mean ± SEM cross-correlogram with peak lag annotation.
    """
    fig, ax = plt.subplots(figsize=(8, 4))

    ax.fill_between(lags_s, mean_xcorr - sem_xcorr, mean_xcorr + sem_xcorr,
                    color=color, alpha=fill_alpha)
    ax.plot(lags_s, mean_xcorr, color=color, linewidth=2,
            label=f'Mean xcorr  (n={n_bouts} bouts)')

    ax.axvline(x=0, color='black', linewidth=1, linestyle='--', label='Zero lag')
    ax.axvline(x=peak_lag_s, color='firebrick', linewidth=1.5, linestyle=':',
               label=f'Peak lag = {peak_lag_s:.2f} s')
    ax.axhline(y=0, color='grey', linewidth=0.8, linestyle=':')

    ax.set_xlabel('Lag (s)  [positive = 560 follows 465]', fontsize=13)
    ax.set_ylabel('Normalized cross-correlation', fontsize=13)
    ax.set_title(f'465 → 560 cross-correlation\n{BOI} — {exp}, {group}', fontsize=13)
    ax.legend(fontsize=10)
    ax.margins(0, 0.05)
    plt.tight_layout()
    return fig

def compute_crosscorr_significance(peth_465_list, peth_560_list, sr,
                                    max_lag_s=5, n_shuffles=1000, ci=95):
    """
    Assess cross-correlation significance via trial-shuffle permutation test.

    Null hypothesis: the co-modulation between 465nm and 560nm traces is no
    stronger than expected from randomly paired (unpaired) bouts.

    Parameters
    ----------
    peth_465_list, peth_560_list : list of np.ndarray (n_bouts, timepoints)
    sr          : float   — sampling rate (Hz)
    n_shuffles  : int     — number of permutations
    ci          : float   — confidence interval width (%)

    Returns
    -------
    lags_s      : np.ndarray
    mean_xcorr  : np.ndarray   — real cross-correlation
    sem_xcorr   : np.ndarray
    peak_lag_s  : float
    ci_low      : np.ndarray   — lower bound of shuffle null
    ci_high     : np.ndarray   — upper bound of shuffle null
    is_sig      : np.ndarray   — boolean mask, True where real > null CI
    """
    def _xcorr_matrix(traces_465, traces_560, n_tp):
        xcorrs = []
        for s1, s2 in zip(traces_465, traces_560):
            s1 = s1 - s1.mean();  s2 = s2 - s2.mean()
            norm = np.sqrt(np.dot(s1, s1) * np.dot(s2, s2))
            if norm < 1e-10:
                continue
            xcorrs.append(correlate(s2, s1, mode='full') / norm)
        return np.array(xcorrs) if xcorrs else None

    # Flatten all bouts across mice
    all_465 = np.concatenate(peth_465_list, axis=0)
    all_560 = np.concatenate(peth_560_list, axis=0)
    n_tp    = all_465.shape[1]

    lags   = correlation_lags(n_tp, n_tp, mode='full')
    lags_s = lags / sr
    mask   = np.abs(lags_s) <= max_lag_s

    # ── Real cross-correlation ────────────────────────────────────────────────
    real_mat   = _xcorr_matrix(all_465, all_560, n_tp)
    mean_xcorr = real_mat[:, mask].mean(axis=0)
    sem_xcorr  = real_mat[:, mask].std(axis=0) / np.sqrt(len(real_mat))
    peak_lag_s = float(lags_s[mask][np.argmax(mean_xcorr)])

    # ── Shuffle null distribution ─────────────────────────────────────────────
    shuffle_means = []
    for _ in range(n_shuffles):
        perm         = np.random.permutation(len(all_465))
        shuf_mat     = _xcorr_matrix(all_465, all_560[perm], n_tp)
        if shuf_mat is not None:
            shuffle_means.append(shuf_mat[:, mask].mean(axis=0))

    shuffle_arr = np.array(shuffle_means)          # (n_shuffles, n_lags)
    alpha       = (100 - ci) / 2
    ci_low      = np.percentile(shuffle_arr, alpha,       axis=0)
    ci_high     = np.percentile(shuffle_arr, 100 - alpha, axis=0)
    is_sig      = (mean_xcorr > ci_high) | (mean_xcorr < ci_low)

    return lags_s[mask], mean_xcorr, sem_xcorr, peak_lag_s, ci_low, ci_high, is_sig

def plot_crosscorr_with_significance(lags_s, mean_xcorr, sem_xcorr,
                                      peak_lag_s, ci_low, ci_high, is_sig,
                                      BOI, exp, group, n_bouts,
                                      color='cornflowerblue',
                                      sig_style='bar'):
    """
    Parameters
    ----------
    sig_style : str
        'bar'     — thin coloured bar at the top of the axes where significant.
        'overlay' — significant portion of the trace replotted thicker on top.
    """
    fig, ax = plt.subplots(figsize=(8, 4))

    # Shuffle null envelope
    ax.fill_between(lags_s, ci_low, ci_high,
                    color='grey', alpha=0.25, label='Shuffle null (95% CI)')

    # SEM ribbon + mean trace
    ax.fill_between(lags_s,
                    mean_xcorr - sem_xcorr,
                    mean_xcorr + sem_xcorr,
                    color=color, alpha=0.3)
    ax.plot(lags_s, mean_xcorr,
            color=color, linewidth=2,
            label=f'Cross-correlation  (n={n_bouts} bouts)')

    # ── Significance indicator ────────────────────────────────────────────────
    if sig_style == 'bar':
        # Thin coloured bar pinned to the top of the axes in axes coordinates
        # (x = data coords, y = axes fraction → independent of y-axis scale)
        ax.fill_between(lags_s, 0.97, 1.0,
                        where=is_sig,
                        transform=ax.get_xaxis_transform(),
                        color='firebrick', alpha=0.85,
                        linewidth=0, label='p < 0.05')

    elif sig_style == 'overlay':
        # Replot only the significant samples as a thicker line on top
        sig_xcorr = np.where(is_sig, mean_xcorr, np.nan)
        ax.plot(lags_s, sig_xcorr,
                color='firebrick', linewidth=4,
                solid_capstyle='round',
                alpha=0.75, zorder=4,
                label='p < 0.05')

    # Reference lines
    ax.axvline(x=0,          color='black',     linewidth=1,   linestyle='--')
    ax.axvline(x=peak_lag_s, color='firebrick', linewidth=1.5, linestyle=':',
               label=f'Peak lag = {peak_lag_s:.3f} s')
    ax.axhline(y=0,          color='grey',      linewidth=0.8, linestyle=':')

    ax.set_xlabel('Lag (s)  [positive = 560 follows 465]', fontsize=12)
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
    alpha=0.05,exclude_corrupted=True,exclusion_col='dFF ExclusionMask',exclusion_padding_s=2,
    min_bout_s=2,min_std=0.01):
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
        bout_len = offset - onset

        if bout_len < min_bout_samples:
            rejected_bouts += 1
            continue

        # EXCLUSION MASK CHECK
        if exclude_corrupted and exclusion_col in fiberbehav_df.columns:
            padded_start = max(0, onset - exclusion_padding)
            padded_end = min(
                len(fiberbehav_df) - 1,
                offset + exclusion_padding)
            exclusion_window = fiberbehav_df.loc[padded_start:padded_end,exclusion_col]
            if exclusion_window.any():
                rejected_bouts += 1
                continue

        # EXTRACT SEGMENT
        seg = fiberbehav_df.loc[onset:offset,['dFF', '560 dFF']].dropna()

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

    return pd.DataFrame(output)

def compute_joint_psth(peth_465, peth_560):
    """
    Compute the joint PSTH between two fiber photometry channels,
    event-aligned.

    Parameters
    ----------
    peth_465, peth_560 : np.ndarray, shape (n_bouts, timepoints)
        Z-scored PETH arrays (from the existing PETH function).

    Returns
    -------
    jpsth_raw       : np.ndarray (timepoints, timepoints)
        Mean across bouts of the outer product z465 ⊗ z560.
    predictor       : np.ndarray (timepoints, timepoints)
        Outer product of the individual mean PSTHs — expected correlation
        from signal-driven modulation alone.
    jpsth_corrected : np.ndarray (timepoints, timepoints)
        Shuffle-corrected JPSTH: residual correlation beyond the predictor.
    coincidence     : np.ndarray (timepoints,)
        Main diagonal of jpsth_corrected — instantaneous co-activation.
    """
    n_bouts, n_tp = peth_465.shape

    # Raw JPSTH: mean outer product across bouts
    jpsth_raw = sum(
        np.outer(peth_465[i], peth_560[i]) for i in range(n_bouts)
    ) / n_bouts

    # Predictor: outer product of individual mean PSTHs
    predictor = np.outer(peth_465.mean(axis=0), peth_560.mean(axis=0))

    # Corrected JPSTH and coincidence histogram
    jpsth_corrected = jpsth_raw - predictor
    coincidence     = np.diag(jpsth_corrected)

    return jpsth_raw, predictor, jpsth_corrected, coincidence

def plot_joint_psth(jpsth_corrected, coincidence, timewindow,
                    BOI, event, exp, group,
                    n_bouts, cmap='RdBu_r', dff_column='465/560'):
    """
    Plot the shuffle-corrected JPSTH with marginal PSTHs and coincidence
    histogram.

    Layout
    ------
    ┌──────────────┬───┐
    │  JPSTH 2D    │   │  ← right strip: coincidence histogram
    │  (heatmap)   │   │
    └──────────────┴───┘
    The diagonal dashed line marks zero lag (t1 = t2).
    A ridge below the diagonal means 465 leads 560.

    Parameters
    ----------
    jpsth_corrected : np.ndarray (timepoints, timepoints)
    coincidence     : np.ndarray (timepoints,)  — diagonal of jpsth_corrected
    timewindow      : [PRE_TIME, POST_TIME]
    n_bouts         : int  — number of bouts used
    """
    PRE_TIME, POST_TIME = float(timewindow[0]), float(timewindow[1])
    n_tp      = jpsth_corrected.shape[0]
    peri_time = np.linspace(-PRE_TIME, POST_TIME, n_tp)

    # ── Layout: main heatmap + narrow right strip for coincidence ─────────────
    fig = plt.figure(figsize=(10, 9))
    gs  = fig.add_gridspec(
        1, 2,
        width_ratios = [10, 2],
        wspace       = 0.08,
    )
    ax_jpsth = fig.add_subplot(gs[0, 0])
    ax_coinc = fig.add_subplot(gs[0, 1], sharey=ax_jpsth)

    # ── JPSTH heatmap ─────────────────────────────────────────────────────────
    abs_max = np.nanmax(np.abs(jpsth_corrected))
    extent  = [-PRE_TIME, POST_TIME, POST_TIME, -PRE_TIME]   # (left, right, bottom, top)

    im = ax_jpsth.imshow(
        jpsth_corrected,
        cmap        = cmap,
        aspect      = 'auto',
        interpolation = 'bilinear',
        extent      = extent,
        vmin        = -abs_max,
        vmax        =  abs_max,
        origin      = 'upper',
    )

    # Diagonal: zero-lag reference line
    ax_jpsth.plot(peri_time, peri_time,
                  color='black', linewidth=1, linestyle='--',
                  label='Zero lag (t₁ = t₂)', alpha=0.6)

    # Event lines
    ax_jpsth.axvline(x=0, color='white', linewidth=2, linestyle='--',
                     alpha=0.8, label=f'{event.capitalize()} (465 axis)')
    ax_jpsth.axhline(y=0, color='white', linewidth=2, linestyle='--',
                     alpha=0.8, label=f'{event.capitalize()} (560 axis)')

    ax_jpsth.set_xlabel('465nm time re. event (s)', fontsize=13)
    ax_jpsth.set_ylabel('560nm time re. event (s)', fontsize=13)
    ax_jpsth.set_title(
        f'Joint PSTH — {BOI} {event.capitalize()}\n{exp}, {group}  (n={n_bouts} bouts)',
        fontsize=13
    )
    ax_jpsth.legend(loc='upper left', fontsize=9, framealpha=0.6)

    cbar = fig.colorbar(im, ax=ax_jpsth, fraction=0.046, pad=0.04)
    cbar.set_label(f'Corrected co-activation\n(z-scored {dff_column} ΔF/F)', fontsize=10)

    # ── Coincidence histogram (main diagonal) ─────────────────────────────────
    ax_coinc.plot(coincidence, peri_time,
                  color='slategray', linewidth=1.5)
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