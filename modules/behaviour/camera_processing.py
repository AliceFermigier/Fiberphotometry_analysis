# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 15:51:56 2023

Functions for preprocessing fiberphotometry data

@author: alice fermigier
"""

#%%
##########
#IMPORTED#
##########

import pandas as pd
import numpy as np
import h5py
import matplotlib.pyplot as plt
import warnings
import scipy
from pathlib import Path

import modules.common.genplot as gp

#%%
###################
#DEFINED FUNCTIONS#
###################

## For BONSAI setup

def extract_sync_channel(raw_path, sync_channel = "DIO04"):
    '''
    Extracts the onset times of TTLs (0→1 transitions) of fiber/miniscope sync to Bonsai.
    '''
    with h5py.File(raw_path, 'r') as f:
        base = "DataAcquisition/FPConsole/Signals/Series0001/"

        time_sync = f[base + "DigitalIO/Time"][:]
        ttl_sync = f[base + f"DigitalIO/{sync_channel}"][:]

    # Find rising edges
    rising_edges = np.where(np.diff(ttl_sync) == 1)[0] + 1

    ttl_sync_df = pd.DataFrame({'Time(s)': time_sync[rising_edges]})

    return ttl_sync_df

def get_timestamps_from_bonsai_csv(file_path):
    df = pd.read_csv(file_path)
    # Ensure correct column names
    df.columns = ['Time(s)', 'Event']
    # Filter rows where Event == True
    output_df = df[df['Event'] == True][['Time(s)']].reset_index(drop=True)
    return output_df

def get_and_clean_sync_timestamps_from_bonsai_csv(file_path, ttl_sync_df):
    df = pd.read_csv(file_path)
    df.columns = ['Time(s)', 'Event']
    true_times = df[df['Event'] == True]['Time(s)'].values
    
    expected_ipi = np.median(np.diff(ttl_sync_df['Time(s)'].values))
    threshold = expected_ipi * 0.5  # half-period: nothing real can be shorter

    iteration = 0
    while True:
        ipi = np.diff(true_times)
        short = ipi < threshold

        if short.sum() == 0:
            break

        keep = np.ones(len(true_times), dtype=bool)
        short_idx = np.where(short)[0]

        for i in short_idx:
            # Skip if already marked for removal
            if not keep[i] or not keep[i + 1]:
                continue
            # A spurious insertion creates short IPIs on BOTH sides:
            # check if the pulse AFTER the short gap also has a short gap before it
            # i.e. both ipi[i] and ipi[i-1] are short → pulse i is the spurious one
            # OR both ipi[i] and ipi[i+1] are short → pulse i+1 is the spurious one
            prev_short = (i > 0) and (ipi[i - 1] < threshold)
            next_short = (i < len(ipi) - 1) and (ipi[i + 1] < threshold)

            if prev_short:
                # pulse i is sandwiched: remove it
                keep[i] = False
            elif next_short:
                # pulse i+1 is sandwiched: remove it
                keep[i + 1] = False
            else:
                # isolated short IPI: remove the second (more conservative)
                keep[i + 1] = False

        true_times = true_times[keep]
        iteration += 1

    if iteration > 0:
        print(f"  Removed spurious Bonsai pulses in {iteration} pass(es), "
              f"{len(true_times)} pulses remaining.")

    return pd.DataFrame({'Time(s)': true_times})

def get_start_stop_timestamps_from_bonsai_csv(file_path):
    df = pd.read_csv(file_path)
    # Ensure correct column names
    df.columns = ['Time(s)', 'Event']
    # Filter rows where Event == True
    start_df = df[df['Event'] == True][['Time(s)']].reset_index(drop=True)
    stop_df = df[df['Event'] == False][['Time(s)']].reset_index(drop=True)
    start_time = start_df.values[0][0]
    stop_time = stop_df.values[0][0]
    output_df = pd.DataFrame({'Time(s)':[start_time,stop_time]})
    return output_df

def time_gap(deinterleaved_df, led_df):
    time_led = led_df['Time(s)']
    time_fiber = deinterleaved_df['Time(s)']

    print(f"Bonsai : {time_led.iloc[0]}-{time_led.iloc[-1]}s ; Start Doric : {time_fiber.iloc[0]}-{time_fiber.iloc[-1]}s")

    time_gap = time_led[0]-time_fiber[0]
    slope = 1.0 #sets slope to 1 by default
    return slope, time_gap

def time_mapping(ttl_sync_df, led_df):
    ttl_times_doric  = ttl_sync_df['Time(s)'].values
    ttl_times_bonsai = led_df['Time(s)'].values

    # Trim to the same number of pulses in case of mismatches
    n = min(len(ttl_times_doric), len(ttl_times_bonsai))
    ttl_times_doric  = ttl_times_doric[:n]
    ttl_times_bonsai = ttl_times_bonsai[:n]

    # Work in relative time to avoid large offset absorbing the slope
    t0_doric  = ttl_times_doric[0]
    t0_bonsai = ttl_times_bonsai[0]
    doric_rel  = ttl_times_doric  - t0_doric
    bonsai_rel = ttl_times_bonsai - t0_bonsai

    # Linear regression on relative times: drift only
    slope, intercept_rel, r_value, _, _ = scipy.stats.linregress(doric_rel, bonsai_rel)

    # Reproject intercept back to absolute Bonsai time
    # bonsai = slope * (doric - t0_doric) + t0_bonsai + intercept_rel
    #        = slope * doric + (t0_bonsai - slope * t0_doric + intercept_rel)
    intercept = t0_bonsai - slope * t0_doric + intercept_rel

    session_duration = ttl_times_doric[-1] - ttl_times_doric[0]
    drift_ms = (slope - 1.0) * session_duration * 1000

    print(f"Pulses used       : {n}")
    print(f"R²                : {r_value**2:.8f}")
    print(f"Slope             : {slope:.6f}")
    print(f"Intercept         : {intercept:.4f} s")
    print(f"Accumulated drift : {drift_ms:.1f} ms")
    if abs(drift_ms) > 50:
        print(f"Large drift detected (>50ms).")

    return slope, intercept

def diagnose_sync(ttl_sync_df, led_df):
    d = ttl_sync_df['Time(s)'].values
    b = led_df['Time(s)'].values

    print(f"Doric  : {len(d)} pulses, span={d[-1]-d[0]:.1f}s")
    print(f"Bonsai : {len(b)} pulses, span={b[-1]-b[0]:.1f}s")
    print(f"Pulse count difference: {len(d) - len(b)}")

    # Compare the span covered by each — should match if alignment is correct
    ipi_d = np.diff(d)
    ipi_b = np.diff(b)
    print(f"\nDoric  IPI: mean={ipi_d.mean():.4f}s  std={ipi_d.std():.5f}s")
    print(f"Bonsai IPI: mean={ipi_b.mean():.4f}s  std={ipi_b.std():.5f}s")

    # The key check: do the spans match?
    span_diff = (d[-1] - d[0]) - (b[-1] - b[0])
    print(f"\nSpan difference (Doric - Bonsai): {span_diff:.3f}s")
    if abs(span_diff) > 0.5:
        print("WARNING: spans differ significantly — likely a pulse count mismatch at start or end")

def diagnose_bonsai_sync(file_path, ttl_sync_df):
    """
    Diagnose spurious pulses in Bonsai sync CSV.
    expected_ipi: expected inter-pulse interval in seconds (0.05s = 20Hz)
    """
    expected_ipi = np.median(np.diff(ttl_sync_df['Time(s)'].values))

    df = pd.read_csv(file_path)
    df.columns = ['Time(s)', 'Event']
    
    true_times = df[df['Event'] == True]['Time(s)'].values
    ipi = np.diff(true_times)
    
    print(f"Total True events : {len(true_times)}")
    print(f"Expected IPI      : {expected_ipi*1000:.3f} ms")
    print(f"IPI mean          : {ipi.mean()*1000:.3f} ms")
    print(f"IPI std           : {ipi.std()*1000:.3f} ms")
    print(f"IPI min           : {ipi.min()*1000:.3f} ms")
    print(f"IPI max           : {ipi.max()*1000:.3f} ms")
    
    # Flag suspiciously short IPIs (duplicates = same pulse logged twice)
    short = ipi < expected_ipi * 0.5
    # Flag suspiciously long IPIs (missed pulse)
    long  = ipi > expected_ipi * 1.5
    
    print(f"\nSuspiciously short IPIs (<{expected_ipi*0.5*1000:.0f}ms, likely duplicates): {short.sum()}")
    print(f"Suspiciously long  IPIs (>{expected_ipi*1.5*1000:.0f}ms, likely missed):     {long.sum()}")
    
    if short.sum() > 0:
        print("\nDuplicate candidate times:")
        dup_indices = np.where(short)[0]
        for i in dup_indices[:20]:  # show first 20
            print(f"  idx {i+1}: {true_times[i]:.6f}s → {true_times[i+1]:.6f}s  "
                  f"(gap={ipi[i]*1000:.2f}ms)")
    
    if long.sum() > 0:
        print("\nMissed pulse candidate times:")
        long_indices = np.where(long)[0]
        for i in long_indices[:20]:
            print(f"  idx {i+1}: {true_times[i]:.6f}s → {true_times[i+1]:.6f}s  "
                  f"(gap={ipi[i]*1000:.2f}ms)")
    
    return true_times, ipi, short.sum(), long.sum()

def correct_behav_timestamps(behaviour_timestamps_df, slope, intercept, time_col='Time(s)', verbose=False):
    behaviour_timestamps_df = behaviour_timestamps_df.copy()
    behaviour_timestamps_df[time_col] = (behaviour_timestamps_df[time_col] - intercept) / slope
    time = behaviour_timestamps_df[time_col].values
    if verbose==True and len(time)>=1:
        print(f'Behavioural timestamps between {time[0]} and {time[-1]}')
    return behaviour_timestamps_df

def align_behav_timestamps(fiberpho_df, behaviour_timestamps_df, behavior_col, time_col='Time(s)'):    
    """
    Align event timestamps to fiber photometry time points.
    Adds a new binary column (1 if event time matches closest fp time, else 0).
    
    Parameters:
        fp_df (pd.DataFrame): Fiber photometry data with a 'Time(s)' column.
        event_df (pd.DataFrame): Event timestamps with a 'Time(s)' column.
        behavior_col (str): Name for the new behavioral column to add.
        time_col (str): Name of the time column (default 'Time(s)').
    
    Returns:
        pd.DataFrame: Original fp_df with an added binary column.
    """
    # Initialize behavioral column with zeros
    fiberpho_df[behavior_col] = 0

    fp_times = fiberpho_df[time_col].values
    event_times = behaviour_timestamps_df[time_col].values
    
    # For each event timestamp, find the closest fiberphotometry time
    indices = np.searchsorted(fp_times, event_times)

        # Adjust indices if necessary (to ensure closest match)
    for i, idx in enumerate(indices):
        if idx == 0:
            closest_idx = 0
        elif idx >= len(fp_times):
            closest_idx = len(fp_times) - 1
        else:
            # Compare distances to find closest
            left = fp_times[idx - 1]
            right = fp_times[idx]
            closest_idx = idx - 1 if abs(event_times[i] - left) < abs(event_times[i] - right) else idx
        fiberpho_df.at[closest_idx, behavior_col] = 1

    return fiberpho_df

def align_camera_flashes(coordinates_df, frame_times_df, mouse=None, batch=None, save_dir_QC=None, expected_fps=20, method='tail'):
    """
    Adds a real timestamp to each DLC / Boris frame.
    
    coordinates_df: DLC / Boris coordinates after filtering
    frame_times_df: df of times for each camera frame (len = n_frames)
    """

    frame_times = frame_times_df.values
    n_dlc = len(coordinates_df)
    n_cam = len(frame_times)
    diff = n_cam-n_dlc

    if diff>1:
        if method == 'even':
            print(f"[!] DLC has {n_dlc}, camera flashes {n_cam} "
                f"({diff} extra flashes). Resampling evenly across recording ")
            idx = np.linspace(0, n_cam - 1, n_dlc).round().astype(int)
            frame_times = frame_times[idx]
        elif method == 'tail':
            print(f"[!] DLC has {n_dlc}, camera flashes {n_cam} ({diff} extra). "
            f"Trimming from the end — video recording likely stopped before flash logger.")
            frame_times = frame_times[:n_dlc]

    elif diff!=0:
        print(f"[!] Truncating tail: DLC has {n_dlc}, camera flashes {n_cam}")
        min_len = min(n_dlc, n_cam)
        # truncate both so they match
        coordinates_df = coordinates_df.iloc[:min_len].copy()
        frame_times = frame_times[:min_len]

    camera_alignment_QC(frame_times_df, coordinates_df=coordinates_df, mouse=mouse, batch=batch,
                    expected_fps=expected_fps, figsize=(14, 8), save_dir=save_dir_QC, show=True)

    coordinates_df["Time(s)"] = frame_times

    return coordinates_df

## For Julien's setup

def get_camera_flashes_from_csv(file_path):
    time_df = pd.read_csv(file_path, header=None)
    timestamps = time_df.loc[time_df[1]==0][0].values
    return pd.DataFrame({'Time(s)': timestamps})

def get_led_flashes_from_csv(file_path):
    timestamps_df = pd.read_csv(file_path, header=None)
    timestamps = timestamps_df[0].values
    return pd.DataFrame({'Time(s)': timestamps})

def align_fiber_with_led_flashes(deinterleaved_df, led_df):
    time_led = led_df['Time(s)'][::2]
    time_fiber = deinterleaved_df['Time(s)']

    min_len = min(len(time_led), len(time_fiber))

    if len(time_led) != len(time_fiber):
        warnings.warn(f'Mismatched lengths: LED = {len(time_led)}, Fiber = {len(time_fiber)}. Truncating to {min_len} samples.')

    aligned_df = pd.DataFrame({
        'Time(s)': time_led.iloc[:min_len].values,
        '405 Deinterleaved': deinterleaved_df['405 Deinterleaved'].iloc[:min_len].values,
        '470 Deinterleaved': deinterleaved_df['465 Deinterleaved'].iloc[:min_len].values
    })

    return aligned_df

## For original setup (deprecated)

def timestamp_camera(rawdata_df) : #deprecated
    """
    Function to extract the timestamps where the camera starts and stops
    --> Parameters
        camera : pd dataframe, camera I/O with sample rate = 12kSps
    --> Returns
        (camera_start, camera_stop) = timestamp when camera starts and stops in seconds (truncated to 0,1s) #camera_stop à enlever si pas besoin
    """
    ind_list = np.where(rawdata_df['DI/O-3'] == 1)[0].tolist()
    (ind_start, ind_stop) = (ind_list[0],ind_list[len(ind_list)-1])
    return (gp.truncate(rawdata_df.at[ind_start, 'Time(s)'], 1),
            gp.truncate(rawdata_df.at[ind_stop, 'Time(s)'], 1))

def load_camera_df_doric(file_path, plot=False):
    with h5py.File(file_path, 'r') as f:
        base = "DataAcquisition/FPConsole/Signals/Series0001/"
        dio_path = base + "DigitalIO/DIO03"
        time_path = base + "DigitalIO/Time"

        if dio_path not in f or time_path not in f:
            print(f"Missing DIO03 or Time path in {file_path}")
            return pd.DataFrame(columns=['Time(s)', 'Camera flashes'])

        camera = f[dio_path][:]
        time = f[time_path][:]

    if len(camera) == 0 or len(time) == 0:
        print(f"Empty camera or time array in {file_path}")
        return pd.DataFrame(columns=['Time(s)', 'Camera flashes'])

    camera_df = pd.DataFrame({
        'Time(s)': time,
        'Camera flashes': camera,
    })

    if plot:
        plt.figure(figsize=(12, 4))
        plt.plot(camera_df['Time(s)'], camera_df['Camera flashes'], drawstyle='steps-post')
        plt.title(f'Camera Flashes - {file_path.split("/")[-1]}')
        plt.xlabel('Time (s)')
        plt.ylabel('Flash Signal')
        plt.tight_layout()
        plt.grid(True)
        plt.show()

    return camera_df

def get_camera_flashes(file_path):
    camera_df = load_camera_df_doric(file_path)
    if camera_df.empty:
        print(f"Camera dataframe is empty for file: {file_path}")
        return pd.DataFrame(columns=['Time(s)'])

    camera_diff = camera_df['Camera flashes'].diff()
    starts = np.where(camera_diff == 1)[0].tolist()
    stops = np.where(camera_diff == -1)[0].tolist()

    if len(starts) == 0 or len(stops) == 0:
        print(f"No flash events detected in {file_path}")
        return pd.DataFrame(columns=['Time(s)'])

    timestamps = []
    for start, stop in zip(starts, stops):
        flash_index = round((start + stop) / 2)
        timestamps.append(camera_df.loc[flash_index]['Time(s)'])

    return pd.DataFrame({'Time(s)': timestamps})

def camera_alignment_QC(frame_times_df, coordinates_df=None, mouse=None, batch=None,
                         expected_fps=20, figsize=(14, 8), save_dir=None, show=True):
    """
    Diagnoses camera TTL frame drops / gaps / extra pulses, and (optionally)
    compares against DLC frame count to flag systematic camera-vs-DLC mismatches.

    Parameters
    ----------
    frame_times_df : pd.DataFrame or Series
        Timestamps of camera TTL flashes.
    coordinates_df : pd.DataFrame, optional
        DLC/Boris coordinates df, used to report n_dlc vs n_cam mismatch.
    mouse, batch : str, optional
        Used for figure titles/filenames if saving.
    expected_fps : float
        Expected camera frame rate, used to flag deviation from nominal IFI.
    figsize : tuple
        Figure size.
    save_dir : str or Path, optional
        If provided, saves the figure as PNG to this directory.
    show : bool
        If False, closes the figure after saving instead of displaying it
        (useful for batch QC over many mice without popping up dozens of plots).

    Returns
    -------
    dict
        Summary stats, useful for logging/aggregating QC across mice.
    """
    flash_times = np.asarray(frame_times_df).flatten()
    diffs = np.diff(flash_times)
    expected_ifi = 1 / expected_fps

    n_cam = len(flash_times)
    n_dlc = len(coordinates_df) if coordinates_df is not None else None
    diff_n = (n_cam - n_dlc) if n_dlc is not None else None

    drop_mask = diffs > 2 * expected_ifi
    n_drops = drop_mask.sum()

    # ── Figure: 2 panels ─────────────────────────────────────────────────
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, sharex=True)

    ax1.plot(diffs, color='black', linewidth=0.8)
    ax1.axhline(expected_ifi, color='grey', ls='--', linewidth=1, label=f'Expected IFI ({expected_ifi*1000:.1f} ms)')
    if n_drops > 0:
        drop_idx = np.where(drop_mask)[0]
        ax1.scatter(drop_idx, diffs[drop_idx], color='red', s=15, zorder=5, label=f'Drops (n={n_drops})')
    ax1.set_ylabel('Inter-frame interval (s)')
    ax1.legend(loc='upper right', fontsize=9)
    title = 'Camera alignment QC'
    if mouse is not None:
        title += f' — {batch} {mouse}' if batch is not None else f' — {mouse}'
    ax1.set_title(title)

    # Cumulative dropped-frame count over the recording — shows WHERE drops
    # concentrate (flat = uniform drops, steepening = drops increasing over time)
    cumulative_drops = np.cumsum(drop_mask)
    ax2.plot(cumulative_drops, color='darkred', linewidth=1.2)
    ax2.set_xlabel('Frame index')
    ax2.set_ylabel('Cumulative dropped frames')

    plt.tight_layout()

    # ── Summary stats ────────────────────────────────────────────────────
    stats = {
        'mouse': mouse,
        'batch': batch,
        'n_flashes': n_cam,
        'n_dlc': n_dlc,
        'diff_n_cam_dlc': diff_n,
        'mean_ifi_ms': diffs.mean() * 1000,
        'std_ifi_ms': diffs.std() * 1000,
        'max_ifi_ms': diffs.max() * 1000,
        'min_ifi_ms': diffs.min() * 1000,
        'n_drops': int(n_drops),
        'pct_drops_first_half': float(drop_mask[:len(drop_mask)//2].mean() * 100),
        'pct_drops_second_half': float(drop_mask[len(drop_mask)//2:].mean() * 100),
    }

    print(f"Total flashes   : {stats['n_flashes']}")
    if n_dlc is not None:
        print(f"DLC frames      : {stats['n_dlc']}  (diff = {stats['diff_n_cam_dlc']})")
    print(f"Mean IFI        : {stats['mean_ifi_ms']:.2f} ms  ({1000/stats['mean_ifi_ms']:.2f} Hz)")
    print(f"Std IFI         : {stats['std_ifi_ms']:.2f} ms")
    print(f"Max IFI         : {stats['max_ifi_ms']:.2f} ms  ← potential gap")
    print(f"Min IFI         : {stats['min_ifi_ms']:.2f} ms  ← potential burst")
    print(f"Frames > 2x IFI : {stats['n_drops']}  ← dropped frames")
    print(f"Drops 1st half  : {stats['pct_drops_first_half']:.2f}%   "
          f"Drops 2nd half : {stats['pct_drops_second_half']:.2f}%  "
          f"{'⚠️  drops concentrated later' if stats['pct_drops_second_half'] > 2*stats['pct_drops_first_half'] else ''}")

    # ── Saving ───────────────────────────────────────────────────────────
    if save_dir is not None:
        save_dir = Path(save_dir)
        save_dir.mkdir(parents=True, exist_ok=True)
        fname = f"camera_QC_{batch}_{mouse}.png" if mouse is not None else "camera_QC.png"
        fig_path = save_dir / fname
        fig.savefig(fig_path, dpi=150, bbox_inches='tight')
        print(f"Saved QC figure to {fig_path}")

    if show:
        plt.show()
    else:
        plt.close(fig)

    return stats