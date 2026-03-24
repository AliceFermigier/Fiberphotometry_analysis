import matplotlib.pyplot as plt
from scipy.signal import medfilt
import numpy as np

import modules.common.preprocess as pp

def smooth(x: np.ndarray, window: int) -> np.ndarray:
    """Moving-average smoothing"""
    if window < 1:
        return x.copy()
    return np.convolve(x, np.ones(window) / window, mode="same")

def iterative_median_filter(fiberpho_df, column_name, n = 1):
    """
    fiberpho_df = fiberphototry data with column names
    n = Median-filter granularity: n=1 -> 1 s steps, n=2 -> 0.5 s steps, etc.
        (1/n must produce an integer number of samples)
    """
    n_windows = 20*n
    time = fiberpho_df['Time(s)'].to_numpy()

    sr = pp.samplerate(fiberpho_df)

    raw_465 = fiberpho_df['465 Deinterleaved'].to_numpy()

    col = plt.cm.jet(np.linspace(0, 1, n_windows))
    
    medfit      = np.zeros((len(time), n_windows))
    neg         = np.zeros(n_windows)
    pos         = np.zeros(n_windows)
    rationegpos = np.zeros(n_windows)
    
    hist_bins = np.arange(-0.0025, 0.0025 + 0.0001, 0.0001)
    n_bins    = len(hist_bins) - 1
    DeltaH    = np.zeros((n_bins, n_windows))
    
    start_idx = sr * 20
    
    fig = plt.figure(figsize=(18, 10))
    ax_top   = fig.add_subplot(2, 4, (1, 3))
    ax_res   = fig.add_subplot(2, 4, (5, 7))
    ax_hist  = fig.add_subplot(2, 4, 4)
    ax_ratio = fig.add_subplot(2, 4, 8)
    
    ax_top.set_title("Raw + filtered baselines")
    ax_res.set_title("Residual signals")
    ax_hist.set_title("Residual distributions")
    ax_ratio.set_title("Neg / Pos ratio")
    
    ax_top.plot(time[start_idx:], raw_465[start_idx:], "m", label="raw")
    ax_top.text(time[-1], raw_465[start_idx:].max(), "raw",
                color="m", fontweight="bold", fontsize=6)
    
    p_min  = raw_465[start_idx:].min()
    p_max  = raw_465[start_idx:].max()
    stepz  = (p_max - p_min) / n_windows
    minmax = np.arange(p_min, p_max + stepz, stepz)
    
    for i in range(1, n_windows + 1):
        idx0 = i - 1   # 0-based index
    
        # scipy.signal.medfilt requires an odd kernel length
        kernel = int(i * sr / n)
        if kernel % 2 == 0:
            kernel += 1
        base = medfilt(raw_465, kernel)
        medfit[:, idx0] = base
    
        d = raw_465 - base

        Delta = np.zeros((len(time), n_windows)) 
        Delta[:, idx0] = d
    
        neg[idx0] = d[d < 0].sum()
        pos[idx0] = d[d >= 0].sum()
    
        frac_neg = (smooth(d, int(2 * sr)) < 0).sum() / len(d)
        if frac_neg > 0 and neg[idx0] != 0:
            rationegpos[idx0] = (pos[idx0] / abs(neg[idx0])) / frac_neg
    
        counts, _ = np.histogram(d, bins=hist_bins)
        DeltaH[:, idx0] = counts
    
        # Filtered baseline
        ax_top.plot(time[start_idx:], base[start_idx:], color=col[idx0])
        if idx0 < len(minmax):
            ax_top.text(time[-1], minmax[idx0],
                        f"{i/n:.0f}s", color=col[idx0], fontsize=6)
    
        # Stacked residual
        offset = i * 0.001
        ax_res.axhline(offset, color=[0.5, 0.5, 0.5], linestyle=":")
        sc = ax_res.scatter(time[start_idx:], offset + d[start_idx:],
                            s=4, c=d[start_idx:], cmap="jet", vmin=-0.005, vmax=0.005)
        ax_res.text(time[-1], offset, f"{i*n:.0f}s", color=col[idx0], fontsize=6)
    
        ax_hist.plot(counts, color=col[idx0])
    
        ax_ratio.plot(i, pos[idx0], "o", color=col[idx0])
        ax_ratio.plot(i, neg[idx0], "o", color=col[idx0])
        ax_ratio.plot(i, rationegpos[idx0], "*", color=col[idx0])
    
    ax_top.axis("tight")
    ax_res.axis("tight")
 
    best_idx   = int(np.argmax(rationegpos))
    best_val   = rationegpos[best_idx]
    best_win_s = (best_idx + 1) / n

    ax_ratio.plot(best_idx + 1, best_val, "ko")

    ax_top.plot(time[start_idx:], medfit[start_idx:, best_idx],
                "k", linewidth=2)
    if best_idx < len(minmax):
        ax_top.text(time[-1], minmax[best_idx],
                    f"{best_win_s:.0f}s (best)", color="k",
                    fontweight="bold", fontsize=6)

    ax_res.plot(time[start_idx:],
                (best_idx + 1) * 0.001 + Delta[start_idx:, best_idx], "k")
    ax_res.text(time[-1], (best_idx + 1) * 0.001,
                f"{best_win_s:.0f}s (best)", color="k", fontweight="bold", fontsize=6)
 
def hybrid_median_filter(fiberpho_df):
# Use best-window baseline except where the long (20-s) filter is lower,
# which avoids over-subtracting genuine negative transients.
    filt_diff = np.where(medfit[:, -1] - medfit[:, best_idx] < 0)[0]
    
    physiofit          = medfit[:, best_idx].copy()
    physiofit[filt_diff] = medfit[filt_diff, -1]
    physiofit          = hampel_filter(physiofit, 5)
    
def median_filter(fiberpho_df, column_name = '405 Deinterleaved'):
    iso_kernel = int((best_idx + 1) * sr / n)
    if iso_kernel % 2 == 0:
        iso_kernel += 1
    isofit = medfilt(iso, iso_kernel), 5

    return filtered_df

def plot_median_filter_results(fiberpho_df, column_name = '405 Deinterleaved'):
    ax_top.plot(time[start_idx:], physiofit[start_idx:], "g",
                linewidth=2, label="hybrid fit")
    ax_top.plot(time[start_idx:], iso[start_idx:] / 2.05, "k",
                label="iso / 2.05")
    ax_top.plot(time[start_idx:], isofit[start_idx:] / 2.05, "g",
                linewidth=1, label="isofit / 2.05")
    
    ax_top.legend(fontsize=6)
    plt.tight_layout()
    plt.show()