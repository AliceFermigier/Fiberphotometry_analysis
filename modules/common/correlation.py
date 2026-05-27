from scipy.signal import correlate, correlation_lags, fftconvolve
import numpy as np
import matplotlib.pyplot as plt
from statsmodels.tsa.stattools import grangercausalitytests
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

def test_granger_causality(fiberbehav_df, behavior_col, max_lag_s=1,
                            sr=None, alpha=0.05):
    """
    Test whether 465nm Granger-causes 560nm within behavior bouts.

    Returns a DataFrame with F-statistic and p-value for each lag tested,
    averaged across bouts.
    """
    if sr is None:
        sr = round(pp.samplerate(fiberbehav_df))
    max_lag_samples = int(max_lag_s * sr)

    onsets  = fiberbehav_df.index[fiberbehav_df[behavior_col] == 1].tolist()
    offsets = fiberbehav_df.index[fiberbehav_df[behavior_col] == -1].tolist()

    results = []
    for onset, offset in zip(onsets, offsets):
        seg = fiberbehav_df.loc[onset:offset, ['dFF', '560 dFF']].dropna()
        if len(seg) < max_lag_samples * 3:
            continue
        try:
            gc = grangercausalitytests(seg.values, maxlag=max_lag_samples)
            for lag, res in gc.items():
                f_stat = res[0]['ssr_ftest'][0]
                p_val  = res[0]['ssr_ftest'][1]
                results.append({'lag_samples': lag,
                                 'lag_s': lag / sr,
                                 'F': f_stat, 'p': p_val})
        except Exception:
            continue

    return pd.DataFrame(results).groupby('lag_s').mean().reset_index()

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