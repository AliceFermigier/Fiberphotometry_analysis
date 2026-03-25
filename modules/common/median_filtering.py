import matplotlib.pyplot as plt
from scipy.signal import medfilt
import numpy as np
import pandas as pd

import modules.common.preprocess as pp
import modules.common.clean_signal as cs

def _force_odd_kernel(kernel: int) -> int:
    """scipy.signal.medfilt requires an odd kernel length."""
    return kernel if kernel % 2 != 0 else kernel + 1

def _hampel_filter(signal: np.ndarray, half_win: int, n_sigma: float = 3.0) -> np.ndarray:
    """
    Hampel identifier: replace outliers with the local median.
 
    Parameters
    ----------
    signal   : 1-D array
    half_win : number of samples on each side of the centre point
    n_sigma  : detection threshold (default 3 σ, same as MATLAB's hampel)
    """
    signal = signal.copy().astype(float)
    n = len(signal)
    k = 1.4826  # MAD → σ consistency constant for Gaussian data
 
    for i in range(n):
        lo = max(0, i - half_win)
        hi = min(n, i + half_win + 1)
        neighbourhood = signal[lo:hi]
        local_median = np.median(neighbourhood)
        local_mad = k * np.median(np.abs(neighbourhood - local_median))
        if np.abs(signal[i] - local_median) > n_sigma * local_mad:
            signal[i] = local_median
    return signal

def median_filter(signal, sample_rate, window_s) -> np.ndarray:
    """
    Apply a median filter to a 1-D signal.
 
    Parameters
    ----------
    signal      : 1-D NumPy array or pandas Series
    sample_rate : sampling frequency in Hz
    window_s    : filter window duration in seconds
 
    Returns
    -------
    filtered_signal : np.ndarray of the same length as `signal`
    """
    if isinstance(signal, pd.Series):
        signal = signal.to_numpy()
    
    kernel = _force_odd_kernel(round(window_s * sample_rate))
    # Pad by half kernel to avoid edge zero-padding artifacts
    half = kernel // 2
    padded = np.pad(signal, half, mode='edge')
    filtered = medfilt(padded.astype(float), kernel)

    return filtered[half:-half]  

def median_filter_dff(fiberpho_df, column_name, window_s, sample_rate = None) -> pd.DataFrame:
    """
    Apply a median filter to one channel of a fiber-photometry DataFrame.
 
    Parameters
    ----------
    fiberpho_df : DataFrame with a 'Time(s)' column and signal columns
    column_name : name of the column to filter
                  (e.g. '465 Deinterleaved')
    window_s    : filter window in seconds
    sample_rate : sampling frequency in Hz; inferred from 'Time(s)' if None
 
    Returns
    -------
    medianfilter_fiberpho_df : DataFrame with columns ['Time(s)', column_name]
                               where column_name contains the filtered signal
    """
    if sample_rate is None:
        sample_rate = pp.samplerate(fiberpho_df)
 
    raw_signal = fiberpho_df[column_name].to_numpy()
    filtered_signal = median_filter(raw_signal, sample_rate, window_s)
 
    return pd.DataFrame({
        "Time(s)": fiberpho_df["Time(s)"].to_numpy(),
        f"{column_name}" : filtered_signal,
    })

def iterative_median_filter(fiberpho_df, column_name, step_size = 1.0) -> pd.DataFrame:
    """
    Sweep median-filter window sizes to find the optimal baseline, then
    return a hybrid-fit baseline as a DataFrame.
 
    The sweep tests windows from `step_size` seconds up to 20 seconds in
    steps of `step_size`.  For each window the residual (signal - baseline)
    is computed and a neg/pos ratio is derived:
 
        ratio = (sum_positive / |sum_negative|) / fraction_of_time_below_zero
 
    The window that maximises this ratio best isolates fast positive
    transients while tracking slow baseline drift.  A hybrid fit is then
    built: the best-window baseline is used for positive segments (usually outside of high transients); 
    the longest (20-s) baseline is substituted wherever it lies below the
    best-window baseline (usually during high physiological transients).
 
    Parameters
    ----------
    fiberpho_df : DataFrame with 'Time(s)' and signal columns
    column_name : column to process (e.g. '465 Deinterleaved')
    step_size   : window step in seconds (default 1.0)
                  Smaller values → finer sweep but longer runtime.
                  Must divide evenly into the sample period.
 
    Returns
    -------
    result_df : DataFrame with columns:
                  'Time(s)'          — original time vector
                  'raw_signal'       — input signal unchanged
                  'best_baseline'    — median-filtered baseline (best window)
                  'hybrid_baseline'  — hybrid fit (best + long-filter blend)
                  'residual'         — raw_signal − hybrid_baseline
    """
    # ── Setup ────────────────────────────────────────────────────────────
    sample_rate = pp.samplerate(fiberpho_df)
    time        = fiberpho_df["Time(s)"].to_numpy()
    raw_signal  = fiberpho_df[column_name].to_numpy().astype(float)
 
    # Number of windows tested: e.g. step_size=1 → 20 windows (1 s … 20 s)
    #                                step_size=0.5 → 40 windows (0.5 s … 20 s)
    n_windows   = round(20.0 / step_size)

    # Skip the first 10 s (recording artefacts / LED stabilisation)
    warmup_idx  = round(sample_rate * 10)
 
    # Pre-allocate result matrices
    all_baselines = np.zeros((len(time), n_windows))   # one column per window
    all_residuals = np.zeros((len(time), n_windows))
    neg_sum       = np.zeros(n_windows)
    pos_sum       = np.zeros(n_windows)
    negpos_ratio  = np.zeros(n_windows)
 
    hist_edges = np.arange(-0.0025, 0.0025 + 0.0001, 0.0001)
    n_bins     = len(hist_edges) - 1
    residual_histograms = np.zeros((n_bins, n_windows))
 
    # ── Plotting setup ────────────────────────────────────────────────────
    colour_map   = plt.cm.jet(np.linspace(0, 1, n_windows))
    signal_min   = raw_signal[warmup_idx:].min()
    signal_max   = raw_signal[warmup_idx:].max()
    label_step   = (signal_max - signal_min) / n_windows
    label_levels = np.arange(signal_min, signal_max + label_step, label_step)
 
    fig = plt.figure(figsize=(18, 10))
    ax_signals    = fig.add_subplot(2, 4, (1, 3))   # raw + all baselines
    ax_residuals  = fig.add_subplot(2, 4, (5, 7))   # stacked residuals
    ax_histograms = fig.add_subplot(2, 4, 4)         # residual distributions
    ax_ratio      = fig.add_subplot(2, 4, 8)         # neg/pos ratio per window
 
    ax_signals.set_title(f"Raw signal + filtered baselines  [{column_name}]")
    ax_residuals.set_title("Residuals (signal − baseline), stacked by window")
    ax_histograms.set_title("Residual distributions")
    ax_ratio.set_title("Neg / Pos ratio  (★ = best window)")
 
    # Plot raw signal
    ax_signals.plot(time[warmup_idx:], raw_signal[warmup_idx:],
                    color="m", linewidth=0.8, label="raw")
    ax_signals.text(time[-1], raw_signal[warmup_idx:].max(), "raw",
                    color="m", fontweight="bold", fontsize=6)
 
    # ── Window sweep ──────────────────────────────────────────────────────
    for win_idx in range(n_windows):           # 0-based
        window_s = (win_idx + 1) * step_size   # actual window in seconds
 
        # --- Median filter baseline ---
        kernel = _force_odd_kernel(round(window_s * sample_rate))
        baseline = medfilt(raw_signal, kernel)
        all_baselines[:, win_idx] = baseline
 
        # --- Residual ---
        residual = raw_signal - baseline
        all_residuals[:, win_idx] = residual
 
        # --- Neg / Pos statistics ---
        neg_sum[win_idx] = residual[residual < 0].sum()
        pos_sum[win_idx] = residual[residual >= 0].sum()
 
        smoothed_residual = cs.smoothing_moving_average(residual, int(2 * sample_rate))
        fraction_below_zero = (smoothed_residual < 0).sum() / len(residual)
 
        if fraction_below_zero > 0 and neg_sum[win_idx] != 0:
            negpos_ratio[win_idx] = (
                pos_sum[win_idx] / abs(neg_sum[win_idx])
            ) / fraction_below_zero
 
        # --- Residual histogram ---
        counts, _ = np.histogram(residual, bins=hist_edges)
        residual_histograms[:, win_idx] = counts
 
        # --- Plot filtered baseline ---
        colour = colour_map[win_idx]
        ax_signals.plot(time[warmup_idx:], baseline[warmup_idx:],
                        color=colour, linewidth=0.6)
        if win_idx < len(label_levels):
            ax_signals.text(time[-1], label_levels[win_idx],
                            f"{window_s:.1f}s", color=colour, fontsize=6)
 
        # --- Stacked residual plot ---
        stack_offset = (win_idx + 1) * 0.01
        ax_residuals.axhline(stack_offset, color=[0.5, 0.5, 0.5], linestyle=":")
        ax_residuals.scatter(
            time[warmup_idx:],
            stack_offset + residual[warmup_idx:],
            s=4, c=residual[warmup_idx:], cmap="jet", vmin=-0.005, vmax=0.005,
        )
        ax_residuals.text(time[-1], stack_offset,
                          f"{window_s:.1f}s", color=colour, fontsize=6)
 
        ax_histograms.plot(counts, color=colour)
 
        ax_ratio.plot(win_idx + 1, pos_sum[win_idx],  "o", color=colour)
        ax_ratio.plot(win_idx + 1, neg_sum[win_idx],  "o", color=colour)
        ax_ratio.plot(win_idx + 1, negpos_ratio[win_idx], "*", color=colour)
 
    ax_signals.axis("tight")
    ax_residuals.axis("tight")
 
    # ── Best window selection ─────────────────────────────────────────────
    best_idx   = int(np.argmax(negpos_ratio))
    best_ratio = negpos_ratio[best_idx]
    best_win_s = (best_idx + 1) * step_size
 
    ax_ratio.plot(best_idx + 1, best_ratio, "ko", markersize=10,
                  label=f"best: {best_win_s:.1f}s")
    ax_ratio.legend(fontsize=7)
 
    # Overlay best baseline
    ax_signals.plot(time[warmup_idx:], all_baselines[warmup_idx:, best_idx],
                    "k", linewidth=2.0)
    if best_idx < len(label_levels):
        ax_signals.text(time[-1], label_levels[best_idx],
                        f"{best_win_s:.1f}s (best)", color="k",
                        fontweight="bold", fontsize=6)
 
    # Overlay best residual
    best_stack_offset = (best_idx + 1) * 0.001
    ax_residuals.plot(
        time[warmup_idx:],
        best_stack_offset + all_residuals[warmup_idx:, best_idx],
        "k", linewidth=1.2,
    )
    ax_residuals.text(time[-1], best_stack_offset,
                      f"{best_win_s:.1f}s (best)", color="k",
                      fontweight="bold", fontsize=6)
 
    # ── Hybrid fit ────────────────────────────────────────────────────────
    # Rule: use the best-window baseline everywhere EXCEPT where the long
    # (20-s) baseline dips below it — those dips correspond to real negative
    # transients that we do not want to subtract away.
    long_baseline = all_baselines[:, -1]        # last column = longest window
    best_baseline = all_baselines[:, best_idx]
 
    hybrid_baseline = best_baseline.copy()
    substitute_mask = long_baseline < best_baseline   # where long is lower
    hybrid_baseline[substitute_mask] = long_baseline[substitute_mask]
 
    # Light Hampel pass to smooth any discontinuities at the join
    hybrid_baseline = _hampel_filter(hybrid_baseline, half_win=5)
 
    # Final residual against the hybrid baseline
    final_residual = raw_signal - hybrid_baseline
 
    # Overlay hybrid baseline on signal plot
    ax_signals.plot(time[warmup_idx:], hybrid_baseline[warmup_idx:],
                    color="g", linewidth=2.0, label="hybrid baseline")
    ax_signals.legend(fontsize=6)
 
    plt.tight_layout()
    plt.show()
 
    # ── Return result DataFrame ───────────────────────────────────────────
    result_df = pd.DataFrame({
        "Time(s)":          time,
        f"{column_name}": hybrid_baseline,
    })
 
    print(
        f"[iterative_median_filter] column='{column_name}'  "
        f"step_size={step_size}s  "
        f"best_window={best_win_s:.1f}s  "
        f"best_ratio={best_ratio:.4f}"
    )
 
    return result_df, best_win_s, fig