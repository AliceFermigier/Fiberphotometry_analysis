# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 14:02:26 2023

General functions for plotting

@author: alice fermigier
"""

#%%
##########
#IMPORTED#
##########

import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt

#%%
###################
#DEFINED FUNCTIONS#
###################

def session_code(session,exp):
    """
    Generate session code in file name
    """
    if exp=='NewContext':
        code = '0'
    else:
        if session in ['Habituation','Training','S1','Conditioning']:
            code = '0'
        elif session in ['S2','Test 1h','Test','Test1']:
            code = '1'
        elif session in ['S3','Test 24h','Test2']:
            code = '2'
        
    return(code)

def plot_rawdata(rawdata_df,exp,session,mouse):
    """
    Plots raw isosbestic and GCaMP traces
    """
    
    fig = plt.figure(figsize=(10,6))
    ax7 = fig.add_subplot(211)

    p1, = ax7.plot('Time(s)', '470 Deinterleaved', 
                   linewidth=1, color='deepskyblue', label='GCaMP', data = rawdata_df[100:]) 
    p2, = ax7.plot('Time(s)', '405 Deinterleaved', 
                   linewidth=1, color='blueviolet', label='ISOS', data = rawdata_df[100:])
    
    ax7.set_ylabel('V')
    ax7.set_xlabel('Time(s)')
    ax7.legend(handles=[p1,p2], loc='upper right')
    ax7.margins(0.01,0.3)
    ax7.set_title(f'GCaMP and Isosbestic raw traces - {exp} {session} {mouse}')
    
    return(fig)

def truncate(n, decimals=0):
    multiplier = 10 ** decimals
    return int(n * multiplier) / multiplier

def time_vector(fiberpho, samplerate) :
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
    return pd.Series(np.linspace(0.0, duration, num = int(duration*samplerate)+1))

def timestamp_camera(rawdata_df) :
    """
    Function to extract the timestamps where the camera starts and stops
    --> Parameters
        camera : pd dataframe, camera I/O with sample rate = 12kSps
    --> Returns
        (camera_start, camera_stop) = timestamp when camera starts and stops in seconds (truncated to 0,1s) #camera_stop à enlever si pas besoin
    """
    ind_list = np.where(rawdata_df['DI/O-3'] == 1)[0].tolist()
    (ind_start, ind_stop) = (ind_list[0],ind_list[len(ind_list)-1])
    return (truncate(rawdata_df.at[ind_start, 'Time(s)'], 1),
            truncate(rawdata_df.at[ind_stop, 'Time(s)'], 1))

def plot_fiberpho(fiberbehav_df,exp,session,mouse,method):
    """
    Plots isosbestic and Ca dependent deltaF/F or fitted
    """
       
    fig1 = plt.figure(figsize=(20,6))
    ax0 = fig1.add_subplot(111)
    
    p1, = ax0.plot('Time(s)', '470 dFF', linewidth=1, color='deepskyblue', 
                   label='GCaMP', data = fiberbehav_df[fiberbehav_df['Time(s)'] > 40]) 
    p2, = ax0.plot('Time(s)', '405 dFF', linewidth=1, color='blueviolet', 
                   label='ISOS', data = fiberbehav_df[fiberbehav_df['Time(s)'] > 40])
    
    ax0.set_ylabel(r'$\Delta$F/F')
    ax0.set_xlabel('Time(s)')
    ax0.legend(handles=[p1,p2], loc='upper right')
    ax0.margins(0.01,0.2)
    ax0.set_title(f'GCaMP and Isosbestic dFF - {exp} {session} {mouse} - {method}')
    
    return fig1