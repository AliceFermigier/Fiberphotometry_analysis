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

def load_camera_df_doric(file_path):
    with h5py.File(file_path, 'r') as f:
        base = "DataAcquisition/FPConsole/Signals/Series0001/"
        
        camera = f[base + "DigitalIO/DIO03"][:]
        time = f[base + "DigitalIO/Time"][:]
        
    camera_df = pd.DataFrame({
        'Time(s)': time,
        'Camera flashes': camera,
    })

    return camera_df

def get_camera_flashes(file_path):
    camera_df = load_camera_df_doric(file_path)
    camera_diff = camera_df['Camera flashes'].diff()
    starts = np.where(camera_diff==1)[0].tolist()
    stops = np.where(camera_diff==-1)[0].tolist()
    timestamps=[]
    flash_indexes=[]
    for start,stop in zip(starts,stops):
        flash_index = round((start+stop)/2)
        flash_indexes.append(round(flash_index))
        timestamps.append(camera_df.loc[flash_index]['Time(s)'])
    
    camera_flashes_df = pd.DataFrame({
        'Time(s)': timestamps,
    })

    return camera_flashes_df

def align_camera_flashes(behav_df, camera_df):
    if len(behav_df)==len(camera_df['Time(s)']):
        time_df = camera_df
    else:
        [start, stop] = [camera_df['Time(s)'].iloc[0], camera_df['Time(s)'].iloc[-1]]
        time = np.linspace(start, stop, len(behav_df))
        time_df = pd.DataFrame({
            'Time(s)':time
        })

    return pd.concat([time_df, behav_df], axis=1)