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
import plotly.express as px

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
    if '560 Deinterleaved' in rawdata_df.columns:
        ax7.plot(rawdata_subset['Time(s)'], rawdata_subset['560 Deinterleaved'], 
             linewidth=1, color='orange', label='rGECO')
    
    # Customize axis
    ax7.set_xlabel('Time (s)')
    ax7.set_ylabel('Voltage (V)')
    if '560 Deinterleaved' in rawdata_df.columns:
        title = f'GCaMP, rGECO and Isosbestic Raw Traces - {exp} {mouse}'
    else:
        title = f'GCaMP and Isosbestic Raw Traces - {exp} {mouse}'
    ax7.set_title(title)
    ax7.legend(loc='upper right')
    ax7.margins(0, 0.3)
    
    return fig

def plot_rawdata_interactive(rawdata_df, exp, mouse, crop=0):
    """
    Interactive plot of raw isosbestic (405 nm) and GCaMP (465 nm) traces, optionally rGECO (560 nm).

    Parameters:
    rawdata_df (pd.DataFrame): Must contain 'Time(s)', '465 Deinterleaved', '405 Deinterleaved'. Optional: '560 Deinterleaved'.
    exp (str): Experiment name.
    mouse (str): Mouse identifier.
    crop (int): Rows to skip from the beginning.

    Returns:
    plotly.graph_objects.Figure: Interactive figure.
    """
    # Slice data
    df = rawdata_df.iloc[crop:].copy()

    # Melt dataframe for px.line
    channels = ['465 Deinterleaved', '405 Deinterleaved']
    if '560 Deinterleaved' in df.columns:
        channels.append('560 Deinterleaved')

    df_long = df.melt(id_vars='Time(s)', value_vars=channels,
                      var_name='Signal', value_name='Voltage')

    # Map colors similar to original
    color_map = {
        '465 Deinterleaved': 'deepskyblue',
        '405 Deinterleaved': 'blueviolet',
        '560 Deinterleaved': 'orange'
    }

    # Create interactive line plot
    fig = px.line(df_long, x='Time(s)', y='Voltage', color='Signal',
                  color_discrete_map=color_map,
                  title=f"{'GCaMP, rGECO and ' if '560 Deinterleaved' in df.columns else ''}Isosbestic Raw Traces - {exp} {mouse}")

    # Add layout tweaks
    fig.update_layout(
        xaxis_title="Time (s)",
        yaxis_title="Voltage (V)",
        legend_title="Signal",
        margin=dict(l=50, r=50, t=50, b=50),
        hovermode="x unified"
    )

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

def plot_fiberpho(fiber_df, exp, mouse, method):
    """
    Plots isosbestic and Ca dependent deltaF/F (dFF) and separate dFF plot
    """

    fiber_df = fiber_df.iloc[20:-20]  # skip first and last 20 points
    fig = plt.figure(figsize=(20, 10))  # increased height for two plots
    
    # First subplot: GCaMP and ISOS
    ax0 = fig.add_subplot(211)
    p1, = ax0.plot('Time(s)', '465 Fitted', linewidth=1, color='deepskyblue', label='GCaMP', data=fiber_df)
    p2, = ax0.plot('Time(s)', '405 Fitted', linewidth=1, color='blueviolet', label='ISOS', data=fiber_df)
    ax0.set_ylabel(r'$\Delta$F/F')
    ax0.set_xlabel('Time(s)')
    ax0.legend(handles=[p1, p2], loc='upper right')
    ax0.margins(0, 0.2)
    ax0.set_title(f'GCaMP and Isosbestic - {exp} {mouse} - {method}')
    
    # Second subplot: just 465 dFF (or any other dFF of interest)
    ax1 = fig.add_subplot(212)
    p3, = ax1.plot('Time(s)', 'dFF', linewidth=1, color='black', label='dFF', data=fiber_df)
    ax1.set_ylabel(r'$\Delta$F/F')
    ax1.set_xlabel('Time(s)')
    ax1.legend(loc='upper right')
    ax1.margins(0, 0.2)
    ax1.set_title('dFF')
    
    plt.tight_layout()
    return fig

def plot_fiberpho_dualcolor(fiber_df, exp, mouse, method):
    """
    Plots isosbestic and Ca dependent deltaF/F (dFF) and separate dFF plot
    """
    fiber_df = fiber_df.iloc[20:-20]  # skip first and last 20 points
    if len(fiber_df.columns) == 5:
        fig = plt.figure(figsize=(20, 15))
        
        # First subplot: GCaMP and ISOS
        ax0 = fig.add_subplot(311)
        p1, = ax0.plot('Time(s)', '465 Fitted', linewidth=1, color='deepskyblue', label='GCaMP', data=fiber_df)
        p2, = ax0.plot('Time(s)', '405 Fitted', linewidth=1, color='blueviolet', label='fitted ISOS', data=fiber_df)
        ax0.set_ylabel(r'$\Delta$F/F')
        ax0.set_xlabel('Time(s)')
        ax0.legend(handles=[p1, p2], loc='upper right')
        ax0.margins(0, 0.2)
        ax0.set_title(f'Preprocessed data - {exp} {mouse} - {method}')
        
        # Second subplot: denoised 465 dFF
        ax1 = fig.add_subplot(312)
        p3, = ax1.plot('Time(s)', 'dFF', linewidth=1, color='deepskyblue', label='dFF', data=fiber_df)
        ax1.set_ylabel(r'$\Delta$F/F')
        ax1.set_xlabel('Time(s)')
        ax1.legend(loc='upper right')
        ax1.margins(0, 0.2)

        # Third subplot: 560 dFF
        ax2 = fig.add_subplot(313)
        p4, = ax2.plot('Time(s)', '560 dFF', linewidth=1, color='orange', label='560 dFF', data=fiber_df)
        ax2.set_ylabel(r'$\Delta$F/F')
        ax2.set_xlabel('Time(s)')
        ax2.legend(loc='upper right')
        ax2.margins(0, 0.2)
        
    elif len(fiber_df.columns) == 7:
        fig = plt.figure(figsize=(20, 20))
        
        # First subplot: GCaMP and ISOS
        ax0 = fig.add_subplot(411)
        p1, = ax0.plot('Time(s)', '465 Fitted', linewidth=1, color='deepskyblue', label='GCaMP', data=fiber_df)
        p2, = ax0.plot('Time(s)', '405 Fitted', linewidth=1, color='blueviolet', label='fitted ISOS', data=fiber_df)
        ax0.set_ylabel(r'$\Delta$F/F')
        ax0.set_xlabel('Time(s)')
        ax0.legend(handles=[p1, p2], loc='upper right')
        ax0.margins(0, 0.2)
        ax0.set_title(f'GCaMP and Isosbestic - {exp} {mouse} - {method}')
        
        # Second subplot: denoised 465 dFF
        ax1 = fig.add_subplot(412)
        p3, = ax1.plot('Time(s)', 'dFF', linewidth=1, color='deepskyblue', label='dFF', data=fiber_df)
        ax1.set_ylabel(r'$\Delta$F/F')
        ax1.set_xlabel('Time(s)')
        ax1.legend(loc='upper right')
        ax1.margins(0, 0.2)

        # Third subplot: 560 dFF and fitted 405
        ax2 = fig.add_subplot(413)
        p4, = ax2.plot('Time(s)', '560 Fitted', linewidth=1, color='orange', label='560', data=fiber_df)
        p5, = ax2.plot('Time(s)', '405 Fitted 560', linewidth=1, color='blueviolet', label='fitted ISOS', data=fiber_df)
        ax2.set_ylabel(r'$\Delta$F/F')
        ax2.set_xlabel('Time(s)')
        ax2.legend(handles=[p4, p5], loc='upper right')
        ax2.margins(0, 0.2)

        # Second subplot: denoised 465 dFF
        ax3 = fig.add_subplot(414)
        p6, = ax3.plot('Time(s)', '560 dFF', linewidth=1, color='orange', label='560 dFF', data=fiber_df)
        ax3.set_ylabel(r'$\Delta$F/F')
        ax3.set_xlabel('Time(s)')
        ax3.legend(loc='upper right')
        ax3.margins(0, 0.2)

    plt.tight_layout()
    return fig


def plot_denoised_photometry(denoised_dFFdata):
    """
    Create interactive Plotly Express visualization of photometry data.
    """
    
    # Reshape data for px (needs long format)
    data_long = pd.concat([
        denoised_dFFdata[['Time(s)']].assign(Signal='dFF raw', Value=denoised_dFFdata['Denoised dFF']),
        denoised_dFFdata[['Time(s)']].assign(Signal='dFF lowpass', Value=denoised_dFFdata['Denoised lowpass dFF'])
    ])
    
    fig = px.line(
        data_long,
        x='Time(s)',
        y='Value',
        color='Signal',
        color_discrete_map={'dFF raw': 'green', 'dFF lowpass': 'green'},
        labels={'Time(s)': 'Time (seconds)', 'Value': 'GCaMP Signal (V)'},
        title='Denoised signals'
    )
    
    # Adjust opacity for raw signal
    fig.data[0].update(opacity=0.3, line=dict(width=1))
    fig.data[1].update(line=dict(width=2))
    
    fig.update_layout(hovermode='x unified', template='plotly_white', height=600)
    fig.show()
    return fig