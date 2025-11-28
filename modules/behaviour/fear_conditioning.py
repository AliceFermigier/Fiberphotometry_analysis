import pandas as pd
import numpy as np
from scipy.signal import savgol_filter
from scipy.ndimage import uniform_filter1d
import json
import os
import importlib

import modules.behaviour.mouse_position as mp
importlib.reload(mp)

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
    
    # Convert ms → seconds
    df[time_col] = df[time_col] / 1000.0

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

def speeds(df, fps=20, dist_file="dist.json", threshold_file="threshold.json"):
    """
    Reimplementation of the MATLAB freezing detection function.
    df: DLC dataframe with coordinate columns.
    """

    # --- Extract coordinates ---
    nose_x = df["nose_x"].values
    nose_y = df["nose_y"].values
    center_x = df["center_x"].values
    center_y = df["center_y"].values
    tail_x = df["tail_base_x"].values
    tail_y = df["tail_base_y"].values

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
 