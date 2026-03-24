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
    """
    Computes mapping from Doric time → Bonsai time using linear regression
    across all 10Hz TTL sync pulses, correcting for clock drift.

    Args:
        ttl_sync_df : DataFrame with 'Time(s)' column — onset times in Doric time
        led_df      : DataFrame with 'Time(s)' column — reception times in Bonsai time

    Returns:
        slope, intercept : such that bonsai_time ≈ slope * doric_time + intercept
    """
    ttl_times_doric  = ttl_sync_df['Time(s)'].values
    ttl_times_bonsai = led_df['Time(s)'].values

    # Trim to the same number of pulses in case of mismatches
    n = min(len(ttl_times_doric), len(ttl_times_bonsai))
    ttl_times_doric  = ttl_times_doric[:n]
    ttl_times_bonsai = ttl_times_bonsai[:n]

    # Linear regression: bonsai_time = slope * doric_time + intercept
    slope, intercept, r_value, _, _ = scipy.stats.linregress(ttl_times_doric, ttl_times_bonsai)

    session_duration = ttl_times_doric[-1] - ttl_times_doric[0]
    drift_ms = (slope - 1.0) * session_duration * 1000

    print(f"Pulses used       : {n}")
    print(f"R²                : {r_value**2:.8f}")
    print(f"Slope             : {slope:.6f}")
    print(f"Intercept         : {intercept:.4f} s")
    print(f"Accumulated drift : {drift_ms:.1f} ms")
    if abs(drift_ms) > 50:
        warnings.warn(f"Large drift detected: {drift_ms:.1f} ms. Check for missed/extra pulses.")

    return slope, intercept

def correct_behav_timestamps(behaviour_timestamps_df, slope, intercept, time_col='Time(s)'):
    behaviour_timestamps_df = behaviour_timestamps_df.copy()
    behaviour_timestamps_df[time_col] = (behaviour_timestamps_df[time_col] - intercept) / slope
    time = behaviour_timestamps_df[time_col].values
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

def align_camera_flashes(coordinates_df, frame_times_df):
    """
    Adds a real timestamp to each DLC / Boris frame.
    
    coordinates_df: DLC / Boris coordinates after filtering
    frame_times_df: df of times for each camera frame (len = n_frames)
    """

    frame_times = frame_times_df.values
    n_dlc = len(coordinates_df)
    n_cam = len(frame_times)

    if n_dlc != n_cam:
        print(f"[!] Truncating: DLC has {n_dlc}, camera flashes {n_cam}")
        min_len = min(n_dlc, n_cam)
        # truncate both so they match
        coordinates_df = coordinates_df.iloc[:min_len].copy()
        frame_times = frame_times[:min_len]

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

