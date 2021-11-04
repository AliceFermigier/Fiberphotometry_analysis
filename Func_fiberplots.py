#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  4 11:53:48 2021

Useful functions for fiberpho analysis

@author: alice
"""

#%%
##########
#IMPORTED#
##########

import pandas as pd
import numpy as np
import math
import sys

#Custom
#put path to directory where python files are stored
# sys.path.append('/Users/alice/Documents/GitHub/Fiberphotometry_analysis')
# sys.path.append('C:\\Users\\afermigier\\Documents\\GitHub\\Fiberphotometry_analysis')

#%%
###################
#DEFINED FUNCTIONS#
###################

def session_code(session):
    """
    Generate session code in file name
    """
    if session in ['Habituation','Training','S1']:
        code = '0'
    if session in ['S2','Test 1h','Test']:
        code = '1'
    elif session in ['S3','Test 24h']:
        code = '2'
        
    return(code)

def truncate(n, decimals=0):
    multiplier = 10 ** decimals
    return int(n * multiplier) / multiplier

def time_vector(fiberpho, SAMPLERATE) :
    """
    Creates timevector on which to plot the data, in pd format
    --> Parameters :
        fiberpho = float, duration of trial in secs
        samplerate = int, in Sps (for processed fiberpho data in Doric Neuroscience Studio, samplerate = 10Sps)
    --> Returns :
        timevector = pd series
    """
    #denoised_fiberpho = fiberpho['Analog In. | Ch.1 470 nm (Deinterleaved)_dF/F0_LowPass' + 
    #                             '-Analog In. | Ch.1 405 nm (Deinterleaved)_dF/F0_LowPass'].dropna()
    #--> if better timevector of the exact same lenght as fiberpho data
    
    duration =  math.ceil(fiberpho.at[len(fiberpho)-2,'Time(s)'])
    return pd.Series(np.linspace(0.0, duration, num = int(duration*SAMPLERATE)+1))

def timestamp_camera(camera) :
    """
    Function to extract the timestamps where the camera starts and stops
    --> Parameters
        camera : pd dataframe, camera I/O with sample rate = 12kSps
    --> Returns
        (camera_start, camera_stop) = timestamp when camera starts and stops in seconds (truncated to 0,1s) #camera_stop à enlever si pas besoin
    """
    ind_list = np.where(camera['Digital I/O | Ch.3 DI/O-3'] == 1)
    ind_list = ind_list[0].tolist()
    (ind_start, ind_stop) = (ind_list[0],ind_list[len(ind_list)-1])
    return (truncate(camera.at[ind_start, 'Time(s)'], 1),
            truncate(camera.at[ind_stop, 'Time(s)'], 1))

def timestamp_camera_fromraw(rawdata_df) :
    """
    Function to extract the timestamps where the camera starts and stops
    --> Parameters
        camera : pd dataframe, camera I/O with sample rate = 12kSps
    --> Returns
        (camera_start, camera_stop) = timestamp when camera starts and stops in seconds (truncated to 0,1s) #camera_stop à enlever si pas besoin
    """
    ind_list = np.where(rawdata_df['DI/O-3'] == 1)
    ind_list = ind_list[0].tolist()
    (ind_start, ind_stop) = (ind_list[0],ind_list[len(ind_list)-1])
    return (truncate(rawdata_df.at[ind_start, 'Time(s)'], 1),
            truncate(rawdata_df.at[ind_stop, 'Time(s)'], 1))