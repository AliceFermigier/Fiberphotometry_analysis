import pandas as pd
import numpy as np
from scipy.ndimage import uniform_filter1d
import json
import importlib
import matplotlib.pyplot as plt
from dash import Dash, dcc, html, Input, Output, State
import plotly.graph_objects as go

import modules.behaviour.mouse_position as mp
importlib.reload(mp)


def _sanitize_time_column(df, time_col='T1'):
    """
    Convert time column (ms) to numeric seconds.
    Returns new df with df[time_col] in seconds (float), and raises a helpful error if conversion fails.
    """
    if time_col not in df.columns:
        raise KeyError(f"Time column '{time_col}' not found. Columns: {list(df.columns)}")

    # Strip whitespace from textual cells to help conversion (works inplace)
    df = df.copy()
    df[time_col] = df[time_col].astype(str).str.strip()

    # Remove thousands separators like commas if present (e.g., "60,000")
    df[time_col] = df[time_col].str.replace(r'[,\s]+', '', regex=True)

    # Convert to numeric (ms). Non-convertible entries become NaN
    df[time_col] = pd.to_numeric(df[time_col], errors='coerce')

    # Detect invalid rows
    bad_mask = df[time_col].isna()
    if bad_mask.any():
        bad_rows = df.loc[bad_mask, time_col].index.tolist()
        sample = df.loc[bad_mask].head(10)  # show up to 10 bad rows
        raise ValueError(
            f"Could not convert {len(bad_rows)} rows in '{time_col}' to numbers. "
            f"Indices (up to 10 shown): {bad_rows[:10]}\n"
            f"Sample bad rows (first 10):\n{sample}\n"
            "Hint: open the Excel and look for headers/merged cells or non-numeric tokens in T1."
        )

    # Convert ms -> seconds
    df[time_col] = df[time_col].astype(float) / 1000.0
    return df

def parse_protocol_sheet(path, sheet_name):
    """
    Parse protocol excel sheet and extract intervals for:
      - CS+ (CS+ column 'on' ... '!on')
      - CS- (CS- column 'on' ... '!on')
      - Shock (LED2(1,3) column 'ON' ... '!on')
      - LED3 (LED3(1,4) column 'on' ... '!on') for protocol-wide on/off
    
    Returns dict with lists of (start_s, end_s) in seconds (relative to protocol start).
    """
    df = pd.read_excel(path, sheet_name=sheet_name, dtype=str).fillna("")
    # Normalize columns (strip whitespaces)
    df.columns = [str(c).strip() for c in df.columns]

    time_col = 'T1'
    csplus_col = 'CS+'
    csminus_col = 'CS-'
    shk_col = 'LED2(1,3)'
    led3_col = 'LED3(1,4)'
    
    # Sanitize time column robustly and convert to seconds
    df = _sanitize_time_column(df, time_col='T1')

    # ───────────────────────────────────────────────
    # Helper: extract ON/OFF intervals for a column
    # ───────────────────────────────────────────────
    def extract_intervals(df, signal_col):
        intervals = []
        current_start = None
        current_time = 0.0  # running protocol time
        
        for idx, row in df.iterrows():
            duration = row[time_col]
            cell = str(row.get(signal_col, "")).lower()

            if 'on' == cell:          # ON starts now
                if current_start is None:
                    current_start = current_time
                    
            elif '!on' in cell:        # OFF at start of this row
                if current_start is not None:
                    intervals.append((current_start, current_time))
                    current_start = None
            
            current_time += duration

        # Close if sheet ends while ON:
        if current_start is not None:
            intervals.append((current_start, current_time))
        
        return intervals

    # Extract LED3 → defines protocol start/stop
    led3 = extract_intervals(df, led3_col)

    # Extract CS+ / CS– from columns
    cs_plus = extract_intervals(df, csplus_col)
    cs_minus = extract_intervals(df, csminus_col)
    
    # Extract shock intervals
    shock = extract_intervals(df, shk_col)

    return {
        "CS+": cs_plus,
        "CS-": cs_minus,
        "Shock": shock,
        "LED3": led3
    }

def get_protocol_remapping(protocol_start_df, proto):
    """
    Compute effective slope and start by comparing LED3 duration in
    Doric time (from corrected Bonsai CSV) vs Imetronic time (from Excel).
    
    Returns: protocol_start (Doric), effective_slope
    """
    # Doric-corrected anchors from LED3.csv
    protocol_start = protocol_start_df['Time(s)'].values[0]
    protocol_stop  = protocol_start_df['Time(s)'].values[-1]
    doric_duration = protocol_stop - protocol_start

    # Imetronic anchors from Excel (relative seconds)
    led3_onset, led3_offset = proto["LED3"][0]
    imetronic_duration = led3_offset - led3_onset

    effective_slope = doric_duration / imetronic_duration
    #protocol_start=protocol_start-0.08
    print(f"Protocol start: {protocol_start}")

    print(f"Protocol duration — Imetronic : {imetronic_duration:.3f} s")
    print(f"Protocol duration — Doric     : {doric_duration:.3f} s")
    print(f"Effective slope               : {effective_slope:.8f}")
    print(f"Accumulated drift corrected   : {(doric_duration - imetronic_duration)*1000:.1f} ms")

    return protocol_start, effective_slope

def convert_to_absolute(intervals, protocol_start, effective_slope):
    """
    Remap relative Imetronic intervals to absolute Doric timestamps
    using empirically measured effective slope from LED3 anchors.
    """
    return [
        (protocol_start + (onset * effective_slope),
         protocol_start + (offset * effective_slope))
        for onset, offset in intervals
    ]

def add_interval_column(fp_df, intervals, colname, time_col='Time(s)'):
    """
    Mark fp_df[colname] = 1 for rows with time between any interval start/end (inclusive).
    """
    if colname in fp_df.columns:
        raise ValueError(f"Column {colname} already exists in fp_df.")
    fp_df[colname] = 0
    if len(intervals) == 0:
        return fp_df
    # convert to numpy arrays for speed
    times = fp_df[time_col].values
    mask_total = np.zeros_like(times, dtype=bool)
    for start, end in intervals:
        mask = (times >= start) & (times <= end)
        mask_total = mask_total | mask
    fp_df.loc[mask_total, colname] = 1
    return fp_df

def rms_sliding_window(speed, window_frames):
    """
    Compute sliding-window RMS for a 1D speed array.
    Uses efficient convolution via uniform_filter1d.
    """
    squared = speed ** 2
    mean_sq = uniform_filter1d(squared, size=window_frames, mode="nearest")
    rms = np.sqrt(mean_sq)
    return rms 

def filter_outside_arena(coordinates_df, arena_scale, bodyparts=['center','tail_base'],
                          margin_px=0, interpolate=True, max_interp_gap=10):
    """

    max_interp_gap: max consecutive NaN frames to interpolate across; longer gaps
                     are left as NaN rather than interpolated across blindly
    """
def filter_outside_arena(coordinates_df, arena_scale, bodyparts=None,
                          margin_px=0, interpolate=True, flag_gap_threshold=10):
    """
    NaNs-out x/y coordinates for selected bodyparts, in any frame, where the point
    falls outside the arena rectangle (catches reflections, which can have high DLC
    likelihood despite being spatially impossible). Optionally interpolates across
    short resulting gaps.

    coordinates_df: output of get_dlc_data (flat bodypart_x / bodypart_y columns)
    arena_scale: dict from get_scale_and_arena_rect, with 'Arena_rectangle_px'
                 as {'x1','y1','x2','y2'} corners
    bodyparts: list of bodypart names to filter (e.g. ['center', 'tail_base']).
               If None, filters all bodyparts found in coordinates_df.
               Exclude bodyparts here that can legitimately leave the arena
               footprint (e.g. 'nose' during rearing).
    margin_px: shrink (+) or grow (-) the valid region in pixels
    interpolate: if True, linearly interpolate across NaN gaps created by filtering
    additionally returns a per-bodypart boolean mask
    marking frames that came from a gap longer than flag_gap_threshold — these
    are 'low confidence' interpolations you may want to exclude from freezing
    bouts rather than trust as real movement/stillness.
    """
    df = coordinates_df.copy()

    rect = arena_scale["Arena_rectangle_px"]
    x_min = min(rect["x1"], rect["x2"]) + margin_px
    x_max = max(rect["x1"], rect["x2"]) - margin_px
    y_min = min(rect["y1"], rect["y2"]) + margin_px
    y_max = max(rect["y1"], rect["y2"]) - margin_px

    all_bodyparts = [c[:-2] for c in df.columns if c.endswith('_x')]
    if bodyparts is None:
        bodyparts = all_bodyparts

    long_gap_mask = pd.DataFrame(False, index=df.index, columns=bodyparts)
    n_flagged_total = 0

    for bp in bodyparts:
        x_col, y_col = f'{bp}_x', f'{bp}_y'
        x, y = df[x_col], df[y_col]
        outside = (x < x_min) | (x > x_max) | (y < y_min) | (y > y_max)
        n_flagged_total += outside.sum()
        df.loc[outside, [x_col, y_col]] = np.nan

        # mark which of the flagged frames sit in a run longer than flag_gap_threshold
        run_id = (outside != outside.shift()).cumsum()
        run_lengths = outside.groupby(run_id).transform('sum')
        long_gap_mask[bp] = outside & (run_lengths > flag_gap_threshold)

        if interpolate:
            # no limit: always fill, so downstream never sees NaN
            df[x_col] = df[x_col].interpolate(method='linear', limit_direction='both')
            df[y_col] = df[y_col].interpolate(method='linear', limit_direction='both')

    print(f'[filter_outside_arena] flagged {n_flagged_total} bodypart-frames outside '
          f'arena bounds across {len(bodyparts)} bodypart(s): {bodyparts}')
    return df, long_gap_mask

def detect_freezing(dlc_df, arena_scale, threshold=None,
                     rms_window_sec=0.4, min_bout_sec=1.0, max_gap_sec=0.1):
    """
    df: DLC dataframe with coordinate columns.
    arena_scale: dict with 'Scale_cm_per_px' and 'Video_fps'
    threshold: speed threshold in cm/s (interactive prompt if None)
    rms_window_sec: window for RMS-smoothing speed before thresholding (jitter robustness)
    min_bout_sec: minimum continuous freezing duration to count as a bout
    max_gap_sec: brief above-threshold gaps within a bout shorter than this are
                 tolerated (bridged) rather than splitting the bout — this is what
                 makes the duration check robust to single noisy frames
    """
    dist_scaling = arena_scale["Scale_cm_per_px"]
    fps = arena_scale['Video_fps']

    # --- Compute local speed for each body part ---
    s_nose = mp.compute_speed(dlc_df, dist_scale=dist_scaling, frame_rate=fps, bodypart='nose')
    s_center = mp.compute_speed(dlc_df, dist_scale=dist_scaling, frame_rate=fps, bodypart='center')
    s_tail = mp.compute_speed(dlc_df, dist_scale=dist_scaling, frame_rate=fps, bodypart='tail_base')

    # --- Smooth speeds with sliding-window RMS to absorb tracking jitter ---
    window_frames = max(1, int(round(rms_window_sec * fps)))
    rms_nose = rms_sliding_window(s_nose["Speed"].to_numpy(), window_frames)
    rms_center = rms_sliding_window(s_center["Speed"].to_numpy(), window_frames)
    rms_tail = rms_sliding_window(s_tail["Speed"].to_numpy(), window_frames)

    if threshold is None:
        plt.plot(rms_nose, label="nose (RMS)")
        plt.plot(rms_center, label="center (RMS)")
        plt.plot(rms_tail, label="tail (RMS)")
        plt.legend()
        plt.show()
        threshold = float(input("Enter speed threshold (cm/s): "))
        plt.close()

    # --- Freeze = sustained low movement, on smoothed speed ---
    freeze = (
        (rms_center <= threshold) &
        (rms_nose   <= threshold * 2) &
        (rms_tail   <= threshold)
    ).astype(int)

    # --- Bout detection: find runs of freeze==1, bridge short gaps, keep long-enough runs ---
    samples_1s = int(round(min_bout_sec * fps))
    max_gap_frames = int(round(max_gap_sec * fps))

    freeze_bouts = _bridge_and_filter_runs(freeze, min_run=samples_1s, max_gap=max_gap_frames)

    freezing_df = pd.DataFrame({'Freezing': freeze_bouts})
    total_speed = s_nose["Speed"] + s_center["Speed"] + s_tail["Speed"]
    total_speed_df = pd.DataFrame({'Speed': total_speed})
    rms_df = pd.DataFrame({'Speed_RMS_nose': rms_nose,
                            'Speed_RMS_center': rms_center,
                            'Speed_RMS_tail': rms_tail})

    behav_df = pd.concat([dlc_df, freezing_df, total_speed_df, rms_df], axis=1)
    return behav_df

def _bridge_and_filter_runs(binary_arr, min_run, max_gap):
    """
    Given a 0/1 array:
    1. Bridge (fill in) gaps of 0s shorter than max_gap that sit between two 1-runs,
       so a single noisy frame doesn't split one long bout into two short ones.
    2. Discard remaining runs of 1s shorter than min_run.
    Returns a 0/1 array of the same length.
    """
    arr = binary_arr.copy().astype(int)
    n = len(arr)

    # --- Step 1: bridge short gaps ---
    i = 0
    while i < n:
        if arr[i] == 0:
            start = i
            while i < n and arr[i] == 0:
                i += 1
            gap_len = i - start
            # bridge only if gap has a 1-run on both sides (not at the very edges)
            if gap_len <= max_gap and start > 0 and i < n:
                arr[start:i] = 1
        else:
            i += 1

    # --- Step 2: remove runs shorter than min_run ---
    out = np.zeros(n, dtype=int)
    i = 0
    while i < n:
        if arr[i] == 1:
            start = i
            while i < n and arr[i] == 1:
                i += 1
            run_len = i - start
            if run_len >= min_run:
                out[start:i] = 1
        else:
            i += 1

    return out

