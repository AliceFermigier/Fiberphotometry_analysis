# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 11:30:12 2023

Functions for plotting with plethysmography data

@author: alice fermigier
"""

#%%
##########
#IMPORTED#
##########

#Libraries
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from ast import literal_eval
# citing pyabf : 
# Harden, SW (2020). pyABF 2.2.3. [Online]. Available: https://pypi.org/project/pyabf/

#%%
###################
#DEFINED FUNCTIONS#
###################

def align_sniffs(fiberpho_df, plethys_df, sniffs_df, sr, mouse):
    """
    Aligns fiberpho data with plethysmograph and sniffs data for a specific mouse.

    --> Parameters:
        fiberpho_df : pd.DataFrame 
            DataFrame containing fiber photometry data with 'Time(s)' and 'Denoised dFF' columns.
        plethys_df : pd.DataFrame 
            DataFrame containing plethysmograph data.
        sniffs_df : pd.DataFrame 
            DataFrame containing sniffing information, with 'Subject', 'Odor', 'Stim' and 'Start_Stop' columns.
        sr : float 
            Sampling rate of the system.
        mouse : str 
            The mouse identifier to align the data for.
    
    --> Returns:
        fibersniff_df : pd.DataFrame 
            DataFrame with fiber photometry, plethysmograph, and sniff data for the specified mouse.
    """
    
    # Filter sniffs data for the specific mouse
    sniffmouse_df = sniffs_df.loc[sniffs_df['Subject'] == mouse]
    
    # Create a base DataFrame with fiberpho 'Time(s)' and 'Denoised dFF' columns
    fibersniff_df = pd.DataFrame({
        'Time(s)': fiberpho_df['Time(s)'], 
        'Denoised dFF': fiberpho_df['Denoised dFF']
    })
    
    # Downsample plethysmography data to match the length of fiberpho data
    downsampled_indices = [round(i) for i in np.linspace(0, len(plethys_df) - 1, len(fiberpho_df))]
    plethys_df = plethys_df.loc[downsampled_indices]
    
    # Add the plethysmography data to the fibersniff DataFrame
    fibersniff_df['Plethysmograph'] = plethys_df['AIn-4'].values
    
    # Loop over each odor and align stim and sniff events
    for odor in set(sniffs_df['Odor']):
        
        # Process each stim associated with the odor
        for count, stim in enumerate(sniffmouse_df.loc[sniffmouse_df['Odor'] == odor, 'Stim'].values):
            stim_column_name = f'Stim {odor} {count}'
            fibersniff_df[stim_column_name] = 0  # Initialize a new column with zeros
            
            x_start, x_stop = literal_eval(stim)
            x_start_index, x_stop_index = round(x_start * sr), round(x_stop * sr)
            
            # Ensure the indices are within the bounds of the DataFrame
            x_start_index = max(0, min(x_start_index, len(fibersniff_df) - 1))
            x_stop_index = max(0, min(x_stop_index, len(fibersniff_df) - 1))
            
            if x_stop_index > x_start_index:  # Valid range
                fibersniff_df.loc[x_start_index:x_stop_index, stim_column_name] = 1

        # Process each sniff associated with the odor
        for count, list_sniffs in enumerate(sniffmouse_df.loc[sniffmouse_df['Odor'] == odor, 'Start_Stop'].values):
            sniff_column_name = f'Sniff {odor} {count}'
            fibersniff_df[sniff_column_name] = 0  # Initialize a new column with zeros
            
            for x_start, x_stop in literal_eval(list_sniffs):
                if [x_start, x_stop] != [0, 0]:  # Skip null sniff entries
                    x_start_index, x_stop_index = round(x_start * sr), round(x_stop * sr)
                    
                    # Ensure the indices are within the bounds of the DataFrame
                    x_start_index = max(0, min(x_start_index, len(fibersniff_df) - 1))
                    x_stop_index = max(0, min(x_stop_index, len(fibersniff_df) - 1))
                    
                    if x_stop_index > x_start_index:  # Valid range
                        fibersniff_df.loc[x_start_index:x_stop_index, sniff_column_name] = 1
    
    return fibersniff_df

def process_fibersniff(fibersniff_df, EVENT_TIME_THRESHOLD, THRESH_S, sr):
    """
    Fuse sniffs that are too close together and remove sniffs that are too short.
    
    --> Parameters:
        fibersniff_df : pd.DataFrame 
            DataFrame containing fiber photometry, plethysmograph, and sniff data.
        EVENT_TIME_THRESHOLD : float 
            Minimum duration (in seconds) for a sniff event to be valid.
        THRESH_S : float 
            Maximum gap (in seconds) between sniff events for them to be fused.
        sr : int 
            Sampling rate of the system.
    
    --> Returns:
        fibersniff_df : pd.DataFrame 
            Processed DataFrame with filtered and fused sniff data.
    """
    sr = round(sr)  # Ensure sr is an integer
    
    # --- 1. Fuse sniffs that are too close ---
    sniff_columns = fibersniff_df.filter(like='Sniff').columns  # Filter sniff-related columns
    for col in sniff_columns:
        prev_value = 0  # Store previous row value
        consecutive_count = 0  # Track the number of consecutive zeros
        for idx, current_value in zip(fibersniff_df.index, fibersniff_df[col]):
            if current_value == 0:
                if prev_value == 0:  # If previous and current values are 0, increase count
                    consecutive_count += 1
                else:  # Transition from 1 to 0, start counting
                    consecutive_count = 1
            else:  # When current_value == 1
                if consecutive_count <= THRESH_S * sr:  # If the gap between sniffs is short, fill the gap
                    fibersniff_df.loc[idx - consecutive_count:idx - 1, col] = 1
                consecutive_count = 0  # Reset the count when inside a sniff
            prev_value = current_value

    # --- 2. Remove sniffs that are too short ---
    for col in sniff_columns:
        prev_value = 0  # Store previous row value
        consecutive_count = 0  # Track the number of consecutive ones
        for idx, current_value in zip(fibersniff_df.index, fibersniff_df[col]):
            if current_value == 1:
                if prev_value == 1:  # If previous and current values are 1, increase count
                    consecutive_count += 1
                else:  # Transition from 0 to 1, start counting
                    consecutive_count = 1
            else:  # When current_value == 0
                if consecutive_count <= EVENT_TIME_THRESHOLD * sr:  # If the sniff is too short, remove it
                    fibersniff_df.loc[idx - consecutive_count:idx - 1, col] = 0
                consecutive_count = 0  # Reset the count when outside a sniff
            prev_value = current_value

    return fibersniff_df

def derive(fibersniff_df):
    """
    Calculate the derivative of behaviors of interest (like Sniffs) and put the result in a new DataFrame.
    This highlights when a behavior starts (1) and when it stops (-1).
    
    --> Parameters:
        fibersniff_df : pd.DataFrame 
            DataFrame containing binary behavior columns (0/1) for Sniffs and other behaviors.
    
    --> Returns:
        derived_df : pd.DataFrame 
            A new DataFrame where each behavior column contains 1 at the start of a behavior, 
            -1 at the end, and 0 elsewhere.
    """
    # Create a copy to avoid modifying the original DataFrame
    derived_df = fibersniff_df.copy()
    
    # Select only the columns of interest (starting from column 3) or those with 'Sniff' in the name
    behavior_columns = derived_df.filter(like='Sniff').columns
    
    # Calculate the derivative for each behavior column
    for col in behavior_columns:
        derived_df[col] = derived_df[col].diff().fillna(0)  # Derivative and fill NaNs with 0
    
    return derived_df

def plethyfiber_plot_raw(fiberpho_df, plethys_df, mouse):
    """
    Plots raw whole-body plethysmography (WBP) and fiberphotometry (GCaMP and ISOS) data.

    --> Parameters:
        fiberpho_df : pd.DataFrame 
            DataFrame containing the fiberphotometry data with the following columns:
            - 'Time(s)' : Time vector
            - 'Denoised dFF' : Denoised delta F/F
            - '470 dFF' : GCaMP signal
            - '405 dFF' : Isosbestic control signal
            
        plethys_df : pd.DataFrame 
            DataFrame containing the plethysmography data with the following columns:
            - 'Time(s)' : Time vector
            - 'AIn-4' : Plethysmograph signal (e.g., breathing) 
            
        mouse : str 
            Identifier for the mouse being analyzed (used for the title)
    
    --> Returns:
        fig : matplotlib.figure.Figure 
            The figure object containing all subplots.
    """
    
    # Create a figure with 3 vertically stacked subplots
    fig = plt.figure(figsize=(20, 9))
    
    ### PLOT 1: Whole Body Plethysmography (WBP) ###
    ax3 = fig.add_subplot(311)
    p1, = ax3.plot('Time(s)', 'AIn-4', linewidth=0.5, color='black', data=plethys_df, label='WBP')
    ax3.set_ylabel('WBP (A.U.)')
    ax3.set_title(f'Whole Body Plethysmography and Fiberphotometry - {mouse}', fontsize=14)
    ax3.margins(0, 0.1)
    ax3.legend(handles=[p1], loc='upper right', fontsize='small')
    
    ### PLOT 2: Denoised Fiberphotometry Data ###
    ax4 = fig.add_subplot(312, sharex=ax3)
    p2, = ax4.plot('Time(s)', 'Denoised dFF', linewidth=0.6, color='black', label='Denoised GCaMP-ISOS')
    ax4.set_ylabel(r'$\Delta$F/F')
    ax4.margins(0, 0.1)
    ax4.legend(handles=[p2], loc='upper right', fontsize='small')
    
    ### PLOT 3: GCaMP and Isosbestic Signals ###
    ax5 = fig.add_subplot(313, sharex=ax3)
    p3, = ax5.plot('Time(s)', '470 dFF', linewidth=0.5, color='deepskyblue', label='GCaMP (470 nm)')
    p4, = ax5.plot('Time(s)', '405 dFF', linewidth=0.5, color='blueviolet', label='ISOS (405 nm)')
    ax5.set_ylabel(r'$\Delta$F/F')
    ax5.set_xlabel('Time (s)')
    ax5.margins(0, 0.1)
    ax5.legend(handles=[p3, p4], loc='upper right', fontsize='small')
    
    # Return the figure object
    return fig

def plethyfiber_plot_sniffs(dfibersniff_df, sniffs_df, mouse):
    """
    Plots whole-body plethysmography (WBP) and fiberphotometry (Denoised dFF) data, highlighting sniff events 
    for different odors (Clean, HC, Novel) and marking stimulus start and end points with vertical lines.
    
    --> Parameters:
        dfibersniff_df : pd.DataFrame 
            DataFrame containing fiberphotometry and plethysmography data with the following columns:
            - 'Time(s)' : Time vector
            - 'Plethysmograph' : Plethysmograph signal (e.g., breathing)
            - 'Denoised dFF' : Denoised delta F/F
        
        sniffs_df : pd.DataFrame 
            DataFrame containing sniff event data with columns:
            - 'Subject' : Identifier for the subject (mouse)
            - 'Odor' : Odor type (e.g., 'Clean', 'HC', 'Novel')
            - 'Count' : Unique sniff event count
            - 'Stim' : Tuple [start, end] for stimulation times
            - 'Sniff' : List of tuples [start, end] for sniff times
        
        mouse : str 
            Identifier for the mouse being analyzed (used in the title)
    
    --> Returns:
        fig : matplotlib.figure.Figure 
            The figure object containing all subplots.
    """
    
    # Create figure and subplots
    fig = plt.figure(figsize=(20, 7))
    
    ### PLOT 1: Whole Body Plethysmography (WBP) ###
    ax8 = fig.add_subplot(211)
    p1, = ax8.plot('Time(s)', 'Plethysmograph', linewidth=0.5, color='black', data=dfibersniff_df, label='WBP')
    ax8.set_ylabel('WBP (A.U.)')
    ax8.set_title(f'Whole Body Plethysmography and Fiberphotometry - {mouse}', fontsize=14)
    ax8.margins(0, 0.1)
    ax8.legend(handles=[p1], loc='upper right', fontsize='small')
    
    ### PLOT 2: Denoised Fiberphotometry Data ###
    ax9 = fig.add_subplot(212, sharex=ax8)
    p2, = ax9.plot('Time(s)', 'Denoised dFF', linewidth=0.6, color='black', label='Denoised GCaMP-ISOS')
    ax9.set_ylabel(r'$\Delta$F/F')
    ax9.set_xlabel('Time (s)')
    ax9.margins(0, 0.1)
    ax9.legend(handles=[p2], loc='upper right', fontsize='small')
    
    ### FUNCTION TO HIGHLIGHT SNIFF EVENTS ###
    def plot_sniff_areas(df, ax, odor, color, odor_label_prefix='Sniff'):
        """
        Plots shaded regions on the plot for each sniff event of a given odor type.
        """
        for count in set(sniffs_df.loc[sniffs_df['Odor'] == odor, 'Count']):
            column_name = f'{odor_label_prefix} {odor} {count-1}'
            if column_name in df.columns:
                x_start = 0  # Initialize start point for sniff
                for x, y in zip(df['Time(s)'], df[column_name]):
                    if y == 1:
                        x_start = x
                    if y == -1 and x_start != 0:
                        ax.axvspan(x_start, x, facecolor=color, alpha=0.5, label=f'{odor_label_prefix} {odor}' if count == 1 else None)
                        x_start = 0

    # Call the function for each odor type
    sniff_colors = {'Clean': 'grey', 'HC': 'gold', 'Novel': 'purple'}
    for odor, color in sniff_colors.items():
        plot_sniff_areas(dfibersniff_df, ax9, odor, color)
    
    ### FUNCTION TO PLOT VERTICAL LINES FOR STIM EVENTS ###
    def plot_stim_lines(ax, df, mouse):
        """
        Plots vertical dashed lines to indicate stimulation start and end times.
        """
        for stim in df.loc[df['Subject'] == mouse, 'Stim']:
            x_start, x_stop = literal_eval(stim)
            ax.axvline(x_start, color='lightsteelblue', linestyle='--')
            ax.axvline(x_stop, color='slategrey', linestyle='--')
    
    # Plot stim lines for both plots
    plot_stim_lines(ax8, sniffs_df, mouse)
    plot_stim_lines(ax9, sniffs_df, mouse)
    
    # Return the figure object
    return fig

def PETH_sniff(dfibersniff_df, odor, event, timewindow, mouse, sr, PRE_EVENT_TIME=0, count=0, baselinewindow=False):
    """
    Creates a Peri-Event Time Histogram (PETH) of fiberphotometry data centered on sniff events for a specified odor.
    
    --> Parameters:
        dfibersniff_df : pd.DataFrame 
            DataFrame containing fiberphotometry and behavioral data with columns such as 'Denoised dFF' and 'Sniff {odor} {count}'
            
        odor : str 
            Odor type (e.g., 'Clean', 'HC', 'Novel')
            
        count : int 
            If 0, uses all sniff events for the given odor; otherwise, filters for a specific sniff event count.
            
        event : str 
            Either 'onset' or 'withdrawal', defining which moment to center the PETH on.
            
        timewindow : list 
            Two-element list specifying the [PRE_TIME, POST_TIME] window around the event (in seconds).
            
        mouse : str 
            Identifier for the subject
            
        sr : float 
            Sampling rate (samples per second) of the recorded data.
            
        PRE_EVENT_TIME : float (optional)
            Amount of time (in seconds) prior to the event to use for baseline normalization.
            
        baselinewindow : Bool
            Tells if standard deviation for z-score calculation is on a timewindow before behavior onset (True) or on whole trace (False)
            
    --> Returns:
        PETH_array : np.ndarray 
            2D array where each row corresponds to the z-scored fiberphotometry signal for a specific event, 
            centered on the event and spanning the timewindow [PRE_TIME, POST_TIME].
    """
    
    # Round sampling rate to ensure it's an integer
    sr = round(sr)
    
    # Set time window relative to event
    PRE_TIME = float(timewindow[0])
    POST_TIME = float(timewindow[1])
    WINDOW_BEFORE_ONSET = 1  # Window (in seconds) before event onset to identify the minimum dFF
    
    # Step 1: Identify indices for event onsets and withdrawals
    if count == 0:
        list_ind_event_o = [i for col in dfibersniff_df.columns[3:] 
                            if col.startswith(f'Sniff {odor}') 
                            for i in np.where(dfibersniff_df[col] == 1)[0]]

        list_ind_event_w = [j for col in dfibersniff_df.columns[3:] 
                            if col.startswith(f'Sniff {odor}') 
                            for j in np.where(dfibersniff_df[col] == -1)[0]]
    else:
        list_ind_event_o = np.where(dfibersniff_df[f'Sniff {odor} {count}'] == 1)[0].tolist()
        list_ind_event_w = np.where(dfibersniff_df[f'Sniff {odor} {count}'] == -1)[0].tolist()
    
    # Step 2: Select the event type (onset or withdrawal)
    list_event = list_ind_event_o if event == 'onset' else list_ind_event_w
    
    # Adjust for event 'onset' — find the minimum dFF in the window preceding the event
    if event == 'onset':
        for i, ind_onset in enumerate(list_event):
            if ind_onset - WINDOW_BEFORE_ONSET * sr >= 0:
                min_index = dfibersniff_df.loc[ind_onset - WINDOW_BEFORE_ONSET * sr : ind_onset, 'Denoised dFF'].idxmin()
                list_event[i] = min_index
    
    # Initialize the PETH array
    event_duration_samples = int((POST_TIME + PRE_TIME) * sr + 1)
    PETH_array = np.zeros((len(list_event), event_duration_samples))
    
    # Initialize mean dFF and std on whole trace
    F0 = dfibersniff_df['Denoised dFF'].mean()
    std0 = dfibersniff_df['Denoised dFF'].std()
    
    # Step 3: Calculate z-scored dFF centered on each event
    for i, ind_event in enumerate(list_event):
        # Ensure indices are valid (check if pre-event baseline window fits in the DataFrame)
        if ind_event - PRE_TIME * sr < 0 or ind_event + POST_TIME * sr >= len(dfibersniff_df):
            continue  # Skip this event if the window goes out of bounds
        
        # Calculate baseline F0 and standard deviation (std) from pre-event window
        baseline_start = int(ind_event - PRE_TIME * sr)
        baseline_end = int(ind_event - PRE_EVENT_TIME * sr)
        
        if baselinewindow:
            dFF_baseline = dfibersniff_df.loc[baseline_start:baseline_end, 'Denoised dFF']
            F0 = dFF_baseline.mean()
            std0 = dFF_baseline.std()
        
        # Extract the fiberphotometry signal for the full window around the event
        event_start = int(ind_event - PRE_TIME * sr)
        event_end = int(ind_event + POST_TIME * sr)
        dFF_signal = dfibersniff_df.loc[event_start:event_end, 'Denoised dFF']
        
        # Z-score the signal: z = (dFF - F0) / std0
        if std0 > 0:  # Avoid division by zero
            PETH_array[i] = (dFF_signal - F0) / std0
        else:
            PETH_array[i] = dFF_signal  # If std0 == 0, use the raw signal
    
    return PETH_array

def PETH_stim(dfibersniff_df, odor, event, timewindow, mouse, sr, PRE_EVENT_TIME=0, count=0, baselinewindow=False):
    """
    Creates a Peri-Event Time Histogram (PETH) of fiberphotometry data centered on stimulus events.
    
    --> Parameters:
        dfibersniff_df : pd.DataFrame 
            DataFrame containing fiberphotometry and behavioral data with columns such as 'Denoised dFF' and 'Stim {odor} {count}'
            
        odor : str 
            Odor type (e.g., 'Clean', 'HC', 'Novel')
            
        count : int (optional, default=0)
            If 0, uses all stimulus events for the given odor; otherwise, filters for a specific stimulus count.
            
        event : str 
            Either 'onset' or 'withdrawal', defining which moment to center the PETH on.
            
        timewindow : list 
            Two-element list specifying the [PRE_TIME, POST_TIME] window around the event (in seconds).
            
        mouse : str 
            Identifier for the subject (e.g., 'Mouse_A').
            
        sr : float 
            Sampling rate (samples per second) of the recorded data.
            
        PRE_EVENT_TIME : float (optional, default=0)
            Amount of time (in seconds) prior to the event to use for baseline normalization.
            
        baselinewindow : Bool (optional, default=False)
            Tells if standard deviation for z-score calculation is on a timewindow before behavior onset (True) or on whole trace (False)
            
    --> Returns:
        PETH_array : np.ndarray 
            2D array where each row corresponds to the z-scored fiberphotometry signal for a specific event, 
            centered on the event and spanning the timewindow [PRE_TIME, POST_TIME].
    """
    
    # Round sampling rate to ensure it's an integer
    sr = round(sr)
    
    # Set time window relative to event
    PRE_TIME = timewindow[0]
    POST_TIME = timewindow[1]
    BASELINE_DURATION = 10  # Duration in seconds for baseline calculation before PRE_EVENT_TIME
    
    # Step 1: Identify indices for event onsets and withdrawals
    list_ind_event_o = []
    list_ind_event_w = []
    
    if count == 0:
        for col in dfibersniff_df.columns[3:]:
            if col.startswith(f'Stim {odor}'):
                list_ind_event_o.extend(np.where(dfibersniff_df[col] == 1)[0].tolist())
                list_ind_event_w.extend(np.where(dfibersniff_df[col] == -1)[0].tolist())
    else:
        list_ind_event_o = np.where(dfibersniff_df[f'Stim {odor} {count}'] == 1)[0].tolist()
        list_ind_event_w = np.where(dfibersniff_df[f'Stim {odor} {count}'] == -1)[0].tolist()
    
    # Step 2: Select the event type (onset or withdrawal)
    list_ind_event = list_ind_event_o if event == 'onset' else list_ind_event_w

    # Handle case where there are no events
    if len(list_ind_event) == 0:
        print(f"No {event} events found for odor {odor} and count {count}.")
        return np.array([])  # Return empty array if no events found
    
    # Calculate the size of the PETH array
    event_duration_samples = int((POST_TIME + PRE_TIME) * sr + 1)
    PETH_array = np.zeros((len(list_ind_event), event_duration_samples))
    
    # Initiates mean dFF and std on whole trace
    F0 = dfibersniff_df['Denoised dFF'].mean()
    std0 = dfibersniff_df['Denoised dFF'].std()
    
    # Step 3: Calculate z-scored dFF centered on each event
    for i, ind_event in enumerate(list_ind_event):
        # Ensure indices are within bounds of DataFrame
        baseline_start = int(ind_event - BASELINE_DURATION * sr)
        baseline_end = int(ind_event - PRE_EVENT_TIME * sr)
        event_start = int(ind_event - PRE_TIME * sr)
        event_end = int(ind_event + POST_TIME * sr)
        
        if baseline_start < 0 or event_end >= len(dfibersniff_df):
            continue  # Skip event if it exceeds DataFrame boundaries
        
        if baselinewindow:
            # Calculate baseline F0 and standard deviation (std) for dFF
            dFF_baseline = dfibersniff_df.loc[baseline_start:baseline_end, 'Denoised dFF']
            F0 = dFF_baseline.mean()
            std0 = dFF_baseline.std()
            
        # Extract the dFF signal for the event window
        dFF_signal = dfibersniff_df.loc[event_start:event_end, 'Denoised dFF']
        
        # Z-score the signal: z = (dFF - F0) / std0
        if std0 > 0:  # Avoid division by zero
            PETH_array[i] = (dFF_signal - F0) / std0
        else:
            PETH_array[i] = dFF_signal  # If std0 == 0, return raw signal
    
    return PETH_array

def plot_PETH(PETH_data, odor, event, timewindow, BOI, sr, mouse):
    """
    Plots the PETH (Peri-Event Time Histogram) as both an average plot and a heatmap.
    
    --> Parameters:
        PETH_data : np.ndarray 
            2D array where each row corresponds to a single trial, and each column corresponds to 
            a time point relative to the event (z-scored fluorescence data)
            
        odor : str 
            The odor associated with the event (e.g., 'Clean', 'HC', 'Novel')
            
        event : str 
            The event type ('onset' or 'withdrawal')
            
        timewindow : list 
            A two-element list specifying the [PRE_TIME, POST_TIME] relative to the event (in seconds)
            
        BOI : str 
            Behavior of interest ('Sniff' or 'Stim')
            
        sr : float 
            Sampling rate (samples per second) of the recorded data
            
        mouse : str 
            Identifier for the subject (e.g., 'Mouse_A')
            
    --> Returns:
        fig : matplotlib.figure.Figure 
            The figure object containing the PETH plot.
    """
    
    # Extract pre and post event times
    PRE_TIME = timewindow[0]
    POST_TIME = timewindow[1]
    sr = round(sr)  # Ensure sampling rate is an integer
    
    # Handle case where PETH_data is empty
    if PETH_data.size == 0:
        print(f"No data available for {odor} {event} for mouse {mouse}. Returning None.")
        return None
    
    # Create time vector for the x-axis
    peri_time = np.arange(-PRE_TIME, POST_TIME + 1/sr, 1/sr)
    
    if PETH_data.shape[1] != len(peri_time):
        print(f"Warning: Data length ({PETH_data.shape[1]}) does not match peri_time length ({len(peri_time)}).")
    
    # Calculate mean and standard deviation for each time point
    mean_dFF_snips = np.mean(PETH_data, axis=0)
    std_dFF_snips = np.std(PETH_data, axis=0)
    
    # Create figure and subplots
    fig, (ax6, ax5) = plt.subplots(2, 1, figsize=(6, 10), gridspec_kw={'height_ratios': [1, 2]})
    
    # ----- Plot 1: Heatmap -----
    cax = ax6.imshow(
        PETH_data, 
        cmap='magma', 
        aspect='auto', 
        interpolation='none', 
        extent=[-PRE_TIME, POST_TIME, len(PETH_data), 0], 
        vmin=-6, 
        vmax=9
    )
    ax6.set_ylabel('Bout Number')
    ax6.set_yticks(np.arange(0.5, len(PETH_data), 2))
    ax6.set_yticklabels(np.arange(0, len(PETH_data), 2))
    ax6.axvline(x=0, linewidth=2, color='black', linestyle='--')
    ax6.set_title(f'{odor} {event} - Plethysmo {mouse}')
    
    # Colorbar for the heatmap
    fig.subplots_adjust(right=0.8, hspace=0.1)
    cbar_ax = fig.add_axes([0.85, 0.54, 0.02, 0.34])
    fig.colorbar(cax, cax=cbar_ax)
    
    # ----- Plot 2: Average Plot with Individual Traces -----
    for snip in PETH_data:
        ax5.plot(peri_time, snip, linewidth=0.5, color=[0.7, 0.7, 0.7], label='Individual trials' if snip is PETH_data[0] else None)
    
    ax5.plot(peri_time, mean_dFF_snips, linewidth=2, color='green', label='Mean response')
    
    ax5.fill_between(
        peri_time, 
        mean_dFF_snips + std_dFF_snips, 
        mean_dFF_snips - std_dFF_snips, 
        facecolor='green', 
        alpha=0.2, 
        label='±1 Std Dev'
    )
    
    ax5.axvline(x=0, linewidth=2, color='slategray', linestyle='--', label=f'{BOI} {event}')
    ax5.set_xlabel('Time (s)')
    ax5.set_ylabel('z-scored ΔF/F')
    ax5.marg
    
    return fig

def plot_PETH_pooled(included_groups, colorscheme, PETHarray_list, BOI, event, timewindow, exp, session):
    """
    Plots a pooled PETH (Peri-Event Time Histogram) averaged over multiple groups.
    
    --> Parameters:
        included_groups : list of str
            Names of the groups included in the plot (e.g., ['CD', 'HFD'])
        
        colorscheme : list of str
            List of colors for each group (e.g., ['cornflowerblue', 'coral'])
        
        PETHarray_list : list of np.ndarray
            Each element of the list is a 2D array for one group's PETH data 
            (trials x timepoints) relative to the event.
        
        BOI : str 
            Behavior of interest (e.g., 'Sniff' or 'Stim')
        
        event : str 
            Event type ('onset' or 'withdrawal')
        
        timewindow : list of int 
            Two-element list specifying the [PRE_TIME, POST_TIME] relative to the event (in seconds)
        
        exp : str 
            Name of the experiment (e.g., 'Experiment 1')
        
        session : str 
            Name of the session (e.g., 'Session A')
        
    --> Returns:
        fig : matplotlib.figure.Figure 
            The figure object containing the pooled PETH plot.
    """
    
    # Extract pre and post event times
    PRE_TIME, POST_TIME = timewindow
    peri_time = np.arange(-PRE_TIME, POST_TIME + 0.1, 0.1)  # Create peri-event time vector
    
    # Error check: Make sure the number of groups matches the number of PETH arrays
    if len(included_groups) != len(PETHarray_list):
        raise ValueError(f"Mismatch: {len(included_groups)} groups but {len(PETHarray_list)} PETH arrays.")
    
    if len(colorscheme) != len(included_groups):
        raise ValueError(f"Mismatch: {len(colorscheme)} colors but {len(included_groups)} groups.")
    
    # Create figure and axis
    fig, ax = plt.subplots(figsize=(6, 4))
    
    list_mean_dFF_snips = []
    list_sem_dFF_snips = []
    
    for i, (group, PETH_data) in enumerate(zip(included_groups, PETHarray_list)):
        
        if PETH_data.size == 0:
            print(f"Warning: No data available for group '{group}'. Skipping this group.")
            continue
        
        # Calculate mean and SEM for this group's PETH data
        mean_dFF_snips = np.mean(PETH_data, axis=0)
        sem_dFF_snips = np.std(PETH_data, axis=0) / np.sqrt(len(PETH_data))
        
        list_mean_dFF_snips.append(mean_dFF_snips)
        list_sem_dFF_snips.append(sem_dFF_snips)
        
        # Plot individual traces for this group
        for snip in PETH_data:
            ax.plot(peri_time, snip, linewidth=0.5, color=colorscheme[i], alpha=0.3)
        
        # Plot mean response for this group
        ax.plot(peri_time, mean_dFF_snips, linewidth=2, color=colorscheme[i], label=group)
        
        # Plot shaded region for standard error of the mean (SEM)
        ax.fill_between(
            peri_time, 
            mean_dFF_snips + sem_dFF_snips, 
            mean_dFF_snips - sem_dFF_snips, 
            color=colorscheme[i], 
            alpha=0.2
        )
    
    # Add vertical line for event (0s time point)
    ax.axvline(x=0, linewidth=2, color='slategray', linestyle='--', label=f'{BOI} {event}')
    
    # Customize axes
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('z-scored ΔF/F')
    ax.legend(loc='upper left', fontsize='small')
    ax.margins(0, 0.1)
    ax.set_title(f'{BOI} - {exp} {session} ({", ".join(included_groups)})')
    
    return fig