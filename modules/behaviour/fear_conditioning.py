import pandas as pd
import numpy as np
from scipy.ndimage import uniform_filter1d
import json
import os
import importlib
import matplotlib.pyplot as plt

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

def convert_to_absolute(intervals, protocol_start_time):
    """
    intervals: list of (start_s, end_s) relative to protocol start
    protocol_start_time: bonsai timestamp (seconds) for protocol start
    returns list of (abs_start_s, abs_end_s)
    """
    return [(protocol_start_time + s, protocol_start_time + e) for s, e in intervals]

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

def detect_freezing_rms(speeds, fps=20, window_sec=1.0, threshold=1.0):
    """
    speeds: dict with speeds from different keypoints
            e.g. {"nose": s_nose, "center": s_center, "tail": s_tail}
    fps: sampling rate
    window_sec: RMS window length
    threshold: RMS threshold in cm/s
    """
    window_frames = int(window_sec * fps)

    # Compute RMS for each keypoint
    rms_dict = {}
    for key, v in speeds.items():
        rms_dict[key] = rms_sliding_window(v, window_frames)

    # Combine speeds: use max or mean
    # max = stricter, mean = smoother
    combined_rms = np.maximum.reduce(list(rms_dict.values()))

    # Freeze if RMS < threshold
    freezing = (combined_rms < threshold).astype(int)

    return freezing, combined_rms, rms_dict

def detect_freezing(dlc_df, scale_file, fps=20):
    """
    df: DLC dataframe with coordinate columns.
    video_scale in px/cm
    """

    # --- Load or calculate pixel-to-cm scaling ---
    if os.path.exists(scale_file):
        dist_scaling = json.load(open(scale_file))["Scale_cm_per_px"]
    else:
        print("❗dist.json not found. Scale set to 0.1 cm/px")
        dist_scaling = 0.1

    # --- Compute local speed for each body part ---
    s_nose = mp.compute_speed(dlc_df, dist_scale=dist_scaling, frame_rate=fps, bodypart='nose')
    s_center = mp.compute_speed(dlc_df, dist_scale=dist_scaling, frame_rate=fps, bodypart='center')
    s_tail = mp.compute_speed(dlc_df, dist_scale=dist_scaling, frame_rate=fps, bodypart='tail_base')

    speeds = {"nose": s_nose, "center": s_center, "tail": s_tail}

    # --- Load or set speed threshold ---
    plt.plot(s_nose["Speed"], label="nose")
    plt.plot(s_center["Speed"], label="center")
    plt.plot(s_tail["Speed"], label="tail")
    plt.legend()
    plt.show()
    threshold = float(input("Enter speed threshold (cm/s): "))
    plt.close()

    # --- Freeze = sustained low movement ---
    freeze = (
        (s_center["Speed"] <= threshold) &
        (s_nose["Speed"]   <= threshold) &
        (s_tail["Speed"]   <= threshold)
    ).astype(int)

    # --- Detect freezing bouts ---
    diff_f = np.diff(freeze)
    freeze_bouts = np.zeros_like(freeze)

    # Require >= 1 sec continuous freezing
    samples_1s = fps

    for i in range(samples_1s, len(freeze) - 2 * samples_1s):
        if diff_f[i] == 1 and np.mean(diff_f[i+1:i+samples_1s]) == 0:
            # find end of freezing
            end = np.where(diff_f[i:] < 0)[0]
            if len(end) > 0:
                end = end[0] + i
            else:
                end = len(freeze) - 1
            freeze_bouts[i:end] = 1
    freezing_df = pd.DataFrame({'Freezing':freeze_bouts})
    print(freezing_df)
    total_speed = s_nose["Speed"] + s_center["Speed"] + s_tail["Speed"]
    print(s_nose, s_nose["Speed"])
    total_speed_df = pd.DataFrame({'Speed': total_speed})
    print(total_speed_df)
    behav_df = pd.concat([dlc_df, freezing_df, total_speed_df], axis=0)

    return behav_df
 