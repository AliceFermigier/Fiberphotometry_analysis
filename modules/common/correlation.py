from scipy.signal import correlate, correlation_lags
import numpy as np
import matplotlib.pyplot as plt

def compute_peth_crosscorr(peth_465_list, peth_560_list, sr,
                            max_lag_s=5, min_bouts=1):
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