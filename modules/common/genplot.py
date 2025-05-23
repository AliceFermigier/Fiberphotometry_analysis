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

def session_code(session): #deprecated
    """
    Generate session code in file name.
    
    Parameters:
    session (str): Name of the session (e.g., 'S1', 'Test1', etc.).
    exp (str): Name of the experiment (e.g., 'NewContext', 'OtherExp').
    
    Returns:
    str: A session code ('0', '1', '2', etc.) based on the session and experiment.
    """
    
    session_codes = {
        '0': {'Habituation', 'Training', 'S1', 'Conditioning'},
        '1': {'S2', 'Test 1h', 'Test', 'Test1'},
        '2': {'S3', 'Test 24h', 'Test2'}
    }

    for code, sessions in session_codes.items():
        if session in sessions:
            return code
        else:
            print('WARNING : session name unknown. Please give in protocol.xlsx a session name referenced in gp.session_code.')
    
    # Default return value if the session does not match any known category
    return 'Unknown'

def plot_rawdata(rawdata_df, exp, mouse, crop=0):
    """
    Plots raw isosbestic (405 nm) and GCaMP (465 nm) traces.
    
    Parameters:
    rawdata_df (pd.DataFrame): DataFrame containing 'Time(s)', '465 Deinterleaved', and '405 Deinterleaved' columns.
    exp (str): Experiment name.
    session (str): Session name.
    mouse (str): Mouse identifier.
    
    Returns:
    matplotlib.figure.Figure: The generated plot figure.
    """
    # Slice data
    rawdata_subset = rawdata_df.iloc[crop:]  # Discard first ? rows
    
    # Create figure and axis
    fig, ax7 = plt.subplots(figsize=(10, 6))
    
    # Plot GCaMP (465) and Isosbestic (405) traces
    ax7.plot(rawdata_subset['Time(s)'], rawdata_subset['465 Deinterleaved'], 
             linewidth=1, color='deepskyblue', label='GCaMP')
    ax7.plot(rawdata_subset['Time(s)'], rawdata_subset['405 Deinterleaved'], 
             linewidth=1, color='blueviolet', label='ISOS')
    
    # Customize axis
    ax7.set_xlabel('Time (s)')
    ax7.set_ylabel('Voltage (V)')
    ax7.set_title(f'GCaMP and Isosbestic Raw Traces - {exp} {mouse}')
    ax7.legend(loc='upper right')
    ax7.margins(0, 0.3)  # Small margins around the data
    
    return fig

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

def plot_fiberpho(fiberbehav_df, exp, mouse, method):
    """
    Plots isosbestic and Ca dependent deltaF/F (dFF) and separate dFF plot
    """
    fig = plt.figure(figsize=(20, 10))  # increased height for two plots
    
    # First subplot: GCaMP and ISOS
    ax0 = fig.add_subplot(211)
    p1, = ax0.plot('Time(s)', '465 dFF', linewidth=1, color='deepskyblue', label='GCaMP', data=fiberbehav_df)
    p2, = ax0.plot('Time(s)', '405 dFF', linewidth=1, color='blueviolet', label='ISOS', data=fiberbehav_df)
    ax0.set_ylabel(r'$\Delta$F/F')
    ax0.set_xlabel('Time(s)')
    ax0.legend(handles=[p1, p2], loc='upper right')
    ax0.margins(0, 0.2)
    ax0.set_title(f'GCaMP and Isosbestic - {exp} {mouse} - {method}')
    
    # Second subplot: just 465 dFF (or any other dFF of interest)
    ax1 = fig.add_subplot(212)
    p3, = ax1.plot('Time(s)', 'Denoised dFF', linewidth=1, color='black', label='dFF', data=fiberbehav_df)
    ax1.set_ylabel(r'$\Delta$F/F')
    ax1.set_xlabel('Time(s)')
    ax1.legend(loc='upper right')
    ax1.margins(0, 0.2)
    ax1.set_title('dFF')
    
    plt.tight_layout()
    return fig