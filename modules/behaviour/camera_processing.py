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

import modules.common.genplot as gp

#%%
###################
#DEFINED FUNCTIONS#
###################

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

def get_camera_flashes_from_csv(file_path):
    time_df = pd.read_csv(file_path, header=None)
    timestamps = time_df.loc[time_df[1]==0][0].values
    return pd.DataFrame({'Time(s)': timestamps})

def get_led_flashes_from_csv(file_path):
    timestamps_df = pd.read_csv(file_path, header=None)
    timestamps = timestamps_df[0].values
    return pd.DataFrame({'Time(s)': timestamps})

def align_fiber_with_led_flashes(deinterleaved_df, led_df):
    time_led = led_df['Time(s)']
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

def align_camera_flashes(behav_df, camera_df):
    if camera_df.empty or 'Time(s)' not in camera_df:
        raise ValueError("camera_df is empty or invalid — cannot align.")

    if len(behav_df) == len(camera_df['Time(s)']):
        time_df = camera_df
    else:
        start, stop = camera_df['Time(s)'].iloc[0], camera_df['Time(s)'].iloc[-1]
        time = np.linspace(start, stop, len(behav_df))
        time_df = pd.DataFrame({'Time(s)': time})

    return pd.concat([time_df, behav_df], axis=1)