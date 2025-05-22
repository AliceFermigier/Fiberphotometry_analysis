 # -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 16:44:33 2023

Functions for plotting with behavioural data

@author: alice fermigier
"""

#%%
##########
#IMPORTED#
##########

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys

import modules.preprocess as pp

#%%
###################
#DEFINED FUNCTIONS#
###################

def align_behav(behav10Sps, fiberpho, timevector, timestart_camera, exp):
    """
    Aligns fiber photometry data with behavioral data from Boris on a time vector.
    
    Parameters 
    ----------
    behav10Sps : pd.DataFrame
        Binary behavioral data from Boris, with a sample rate of 10 samples per second (Sps).
        Each column corresponds to a specific behavior.
    fiberpho : pd.DataFrame
        Pre-processed fiber photometry data, including deltaF/F for 470nm and 405nm channels.
    timevector : pd.Series or list
        Time vector corresponding to the timestamps of the experiment.
    timestart_camera : float
        The time (in seconds) corresponding to when the camera starts relative to the time vector.
    exp : str
        Experiment name
    
    Returns
    -------
    pd.DataFrame
        Dataframe containing aligned time, denoised dFF, 405nm dFF, 470nm dFF, and behavioral event columns.
    """
    if behav10Sps.empty or fiberpho.empty or len(timevector) == 0:
        raise ValueError("Input data (behav10Sps, fiberpho, or timevector) is empty.")
    
    # Get list of behaviors from Boris
    list_behav = behav10Sps.columns[1:].tolist()
    
    # Identify where the camera starts, handle the case where timestart_camera is not in timevector
    try:
        indstart = np.where(np.round(timevector, 1) == timestart_camera)[0][0]
    except IndexError:
        raise ValueError(f"Camera start time ({timestart_camera}s) not found in time vector.")
    
    # Align behavior data
    behav_comp = [
        [0] * indstart + behav10Sps[behav].tolist()
        for behav in list_behav
    ]
    
    # Prepare fiber photometry data and time list
    denoised_fiberpho = fiberpho['Denoised dFF'].dropna().tolist()
    dff_405nm_list = fiberpho['405 dFF'].dropna().tolist()
    dff_470nm_list = fiberpho['470 dFF'].dropna().tolist()
    timelist = timevector.tolist()
    
    # Ensure all lists are of equal length
    min_length = min(len(timelist), len(denoised_fiberpho), len(behav_comp[0]), len(dff_405nm_list), len(dff_470nm_list))
    timelist = timelist[:min_length]
    denoised_fiberpho = denoised_fiberpho[:min_length]
    dff_405nm_list = dff_405nm_list[:min_length]
    dff_470nm_list = dff_470nm_list[:min_length]
    behav_crop = [behav[:min_length] for behav in behav_comp]
    
    # Create the final DataFrame
    fiberbehav_df = pd.DataFrame({
        'Time(s)': timelist, 
        'Denoised dFF': denoised_fiberpho,
        '405nm dFF': dff_405nm_list, 
        '470nm dFF': dff_470nm_list
    })
    
    # Insert behavioral columns
    for behav, data in zip(list_behav, behav_crop):
        fiberbehav_df[behav] = data
    
    return fiberbehav_df

def behav_process(fiberbehav_df, list_BOI, THRESH_S, EVENT_TIME_THRESHOLD, sr):
    """
    Processes behavioral events in the fiber photometry dataframe.

    Parameters
    ----------
    fiberbehav_df : pd.DataFrame
        Dataframe containing behavioral event columns (binary 0s and 1s) alongside fiber photometry data.
    list_BOI : list of str
        List of behavior columns (behaviors of interest, BOIs) to process.
    THRESH_S : float
        Maximum allowed time (in seconds) between two consecutive exploration events 
        to consider them as part of the same event.
    EVENT_TIME_THRESHOLD : float
        Minimum required duration (in seconds) for an event to be valid. Shorter events are removed.
    sr : float or int
        Sampling rate (in Hz) of the data.

    Returns
    -------
    pd.DataFrame
        The modified dataframe with updated behavioral event columns.
    """
    
    # Input validation
    if not isinstance(fiberbehav_df, pd.DataFrame):
        raise TypeError("fiberbehav_df must be a Pandas DataFrame.")
    if not isinstance(list_BOI, list) or not all(isinstance(boi, str) for boi in list_BOI):
        raise TypeError("list_BOI must be a list of strings.")
    if not all(boi in fiberbehav_df.columns for boi in list_BOI):
        missing_cols = [boi for boi in list_BOI if boi not in fiberbehav_df.columns]
        raise ValueError(f"The following behaviors of interest (BOI) are not in the DataFrame: {missing_cols}")
    if fiberbehav_df.empty:
        return fiberbehav_df

    sr = round(sr)  # Round sampling rate to ensure it is an integer
    
    for BOI in list_BOI:
        # 1. Fuse exploration events that are too close
        previous_state = 0
        consecutive_count = 0
        for index, current_state in zip(fiberbehav_df.index, fiberbehav_df[BOI]):
            if current_state == previous_state == 0:
                consecutive_count += 1
            elif current_state != previous_state and current_state == 0:
                consecutive_count = 1
            elif current_state != previous_state and current_state == 1:
                if 1 < consecutive_count <= THRESH_S * sr:  # Fuse if gap is small
                    fiberbehav_df.loc[index - consecutive_count:index - 1, BOI] = 1
            previous_state = current_state
        
        # 2. Remove short exploration events
        previous_state = 0
        consecutive_count = 0
        for index, current_state in zip(fiberbehav_df.index, fiberbehav_df[BOI]):
            if current_state == previous_state == 1:
                consecutive_count += 1
            elif current_state != previous_state and current_state == 1:
                consecutive_count = 1
            elif current_state != previous_state and current_state == 0:
                if 1 < consecutive_count <= EVENT_TIME_THRESHOLD * sr:  # Remove if duration is too short
                    fiberbehav_df.loc[index - consecutive_count:index - 1, BOI] = 0
            previous_state = current_state

    return fiberbehav_df

def derive(fiberbehav_df):
    """
    calculate the derivative of behav of interest and put in new df
    that way, it will show 1 when behaviour starts and -1 when it stops
    """
    for col in fiberbehav_df.columns[4:]:
        fiberbehav_df[col] = fiberbehav_df[col].diff()

    return fiberbehav_df 

def highlight_behavior_areas(ax, df, behavior_name, facecolor='grey', alpha=0.3, label_prefix=''):
    """
    Highlights areas on a plot corresponding to specific behavioral periods.
    
    Parameters:
    - ax: The matplotlib axis where the highlight will be drawn.
    - df: The DataFrame containing the behavior columns and time.
    - behavior_name: The name of the behavior column in the DataFrame.
    - facecolor: The color to use for the highlight.
    - alpha: The transparency level for the highlight.
    - label_prefix: A string prefix to make labels unique (to avoid duplicate legend entries).
    """
    i = 0
    x_start = 0
    for (x, y) in zip(df['Time(s)'].tolist(), df[behavior_name].tolist()):
        if y == 1:
            x_start = x
        if y == -1 and x_start != 0:
            ax.axvspan(x_start, x, facecolor=facecolor, alpha=alpha, label='_' * i + label_prefix + behavior_name)
            x_start = 0
            i += 1

def plot_fiberpho_behav(behavprocess_df, list_BOI, exp, mouse, THRESH_S, EVENT_TIME_THRESHOLD, batch):
    """
    Plots denoised deltaF/F aligned with behaviour (includes baseline)
    """
    
    # Adapt to crop plot to your liking
    behavprocesssnip_df = behavprocess_df[behavprocess_df['Time(s)'] > 0]
    
    fig2 = plt.figure(figsize=(20,5))
    ax1 = fig2.add_subplot(111)
    
    # Plot fiberpho trace for all cases
    p1, = ax1.plot('Time(s)', 'Denoised dFF', linewidth=0.6 if session not in ['Test', 'Test 1h', 'Test 24h', 'S3'] else 1, 
                   color='black', label='_GCaMP', data=behavprocesssnip_df)
    
    if exp == 'Fear' and session == 'Conditioning':
            x_1 = behavprocess_df.at[int(np.where(behavprocess_df['Shock'] == 1)[0][0]), 'Time(s)']
            x_2 = behavprocess_df.at[int(np.where(behavprocess_df['Shock'] == 1)[0][1]), 'Time(s)']
            ax1.axvline(x_1, color='yellow', ls='-', lw=2, label='Shock')
            ax1.axvline(x_2, color='yellow', ls='-', lw=2, label='Shock')
    
    behaviors_to_plot = {
        'Water consumption': ('cornflowerblue', 0.5),
        'Water ?': ('cornflowerblue', 0.3),
        'Saccharine consumption': ('gold', 0.5),
        'Saccharine ?': ('gold', 0.3),
        'Ethanol sniffing': ('purple', 0.3),
        'Tail suspension': ('red', 0.3),
        'Homecage': ('gold', 0.3),
        'Fear cage': ('blue', 0.3),
        'New context': ('darkturquoise', 0.3),
        'Exploration fam': ('gold', 0.3),
        'Exploration new': ('purple', 0.3),
        'Climbing': ('cornflowerblue', 0.3),
        'Rearing': ('mediumpurple', 0.3),
        'Exploration left': ('orange', 0.3),
        'Exploration right': ('darkturquoise', 0.3),
        'Exploration non social': ('grey', 0.3),
        'Exploration social': ('mediumvioletred', 0.3),
        'Center': ('mediumvioletred', 0.3),
        'Open arm': ('red', 0.3),
        'Closed arm': ('grey', 0.3)
    }
    
    # Highlight all behaviors in the session
    for behavior, (color, alpha) in behaviors_to_plot.items():
        if behavior in list_BOI and behavior in behavprocesssnip_df.columns:
            highlight_behavior_areas(ax1, behavprocesssnip_df, behavior, color, alpha) 

    # Vertical lines for specific events
    if 'Gate opens' in list_BOI and 'Gate opens' in behavprocesssnip_df.columns:
        x_entry = behavprocess_df.at[int(np.where(behavprocess_df['Gate opens'] == 1)[0][0]), 'Time(s)']
        ax1.axvline(x_entry, color='lightsteelblue', ls='--', label='Gate opens')
            
    if 'Entry in arena' in list_BOI and 'Entry in arena' in behavprocesssnip_df.columns:
        x_entry = behavprocess_df.at[int(np.where(behavprocess_df['Entry in arena'] == 1)[0][0]), 'Time(s)']
        ax1.axvline(x_entry, color='slategrey', ls='--', label='Entry in arena')

    ax1.set_ylabel(r'$\Delta$F/F')
    ax1.set_xlabel('Time(s)')
    ax1.legend(loc='upper left')
    ax1.margins(0.01, 0.03)
    ax1.set_title(f'dFF with Behavioural Scoring - {exp} {session} {mouse} {batch}- interbout {THRESH_S} - cut {EVENT_TIME_THRESHOLD}')
    
    return fig2

def PETH(behavprocess_df, BOI, event, timewindow, EVENT_TIME_THRESHOLD, PRE_EVENT_TIME=0, maxboutsnumber=None, baselinewindow=False):
    """
    Creates dataframe of fiberpho data centered on bout event for BOI.
    
    Parameters
    - behavprocess_df : pd.DataFrame 
            Aligned fiberpho and behavioral data for one mouse.
    - BOI : str 
            Behavior of interest (must match the column name in behavprocess_df).
    - event : str 
            'onset' or 'withdrawal' (event type to center on).
    - timewindow : list 
            Time before and after the event, [PRE_TIME, POST_TIME].
    - EVENT_TIME_THRESHOLD : float 
            Minimum time (in seconds) a bout must last to be included.
    - PRE_EVENT_TIME : float 
            Time window before event for baseline calculations.
    - maxboutsnumber : int or None 
            Maximum number of bouts to include.
    - baselinewindow : Bool 
        Tells if standard deviation for z-score calculation is on a timewindow before behavior onset (True) or on whole trace (False)
    
    Returns
    - PETH_array : np.ndarray 
            Z-scored fiberpho data centered on event with shape (num_bouts, timepoints).
    """
    # Set time window relative to event
    PRE_TIME, POST_TIME = timewindow

    # Ensure sampling rate is an integer
    sr = round(pp.samplerate(behavprocess_df))

    # Identify onset and withdrawal indices for the behavior of interest (BOI)
    list_ind_event_o = np.where(behavprocess_df[BOI] == 1)[0].tolist()
    list_ind_event_w = np.where(behavprocess_df[BOI] == -1)[0].tolist()

    # Limit to maximum number of bouts if maxboutsnumber is set
    if maxboutsnumber is not None:
        list_ind_event_o = list_ind_event_o[:maxboutsnumber]
        list_ind_event_w = list_ind_event_w[:maxboutsnumber]

    # Remove bouts that are too short
    min_event_duration = EVENT_TIME_THRESHOLD * sr
    valid_bouts = [(start, end) for start, end in zip(list_ind_event_o, list_ind_event_w) if 1 < end - start >= min_event_duration]
    list_ind_event_o, list_ind_event_w = zip(*valid_bouts) if valid_bouts else ([], [])

    # Choose the relevant event indices to align on (either onset or withdrawal)
    list_ind_event = list_ind_event_o if event == 'onset' else list_ind_event_w

    # Check if the event happens too late in the dataframe to process
    list_ind_event = [idx for idx in list_ind_event if idx + POST_TIME * sr < len(behavprocess_df)]

    # For onset events, adjust the start index to the minimum dFF within 1 second before and after the event
    if event == 'onset':
        list_ind_event = [
            behavprocess_df.loc[idx - 1 * sr : idx + 1 * sr, 'Denoised dFF'].idxmin() 
            for idx in list_ind_event
        ]

    # Preallocate the PETH array to store the z-scored traces
    n_bouts = len(list_ind_event)
    n_timepoints = (POST_TIME + PRE_TIME) * sr + 1
    PETH_array = np.zeros((n_bouts, n_timepoints))
    
    # Initialize mean and std on whole trace
    F0 = behavprocess_df['Denoised dFF'].mean()
    std0 = behavprocess_df['Denoised dFF'].std()

    # Loop through each event and extract the fiberpho trace centered on the event
    for i, ind_event in enumerate(list_ind_event):
        try: 
            if baselinewindow:
                # Calculate baseline mean (F0) and standard deviation (std0) for the time window before the event
                dFF_baseline = behavprocess_df.loc[ind_event - PRE_TIME * sr : ind_event - PRE_EVENT_TIME * sr, 'Denoised dFF']
                F0 = dFF_baseline.mean() 
                std0 = dFF_baseline.std()

            # Extract the fiberpho trace for the time window around the event
            event_window = behavprocess_df.loc[ind_event - PRE_TIME * sr : ind_event + POST_TIME * sr, 'Denoised dFF']
            
            # Ensure the event window has the correct length to avoid shape mismatch
            if len(event_window) == n_timepoints:
                PETH_array[i] = (event_window - F0) / std0
        except Exception as e:
            print(f"Error processing event at index {ind_event}: {e}")

    return PETH_array

def plot_PETH(PETH_data, BOI, event, timewindow, exp, session, mouse, group, 
              trace_color='black', fill_alpha=0.2, trace_linewidth=2, heatmap_cmap='magma'):
    """
    Plots PETH average and heatmap.
    
    --> Parameters:
        PETH_data : np.ndarray 
            Array of z-scored dFF data centered on events (shape: num_bouts x timepoints).
        BOI : str 
            Behavior of interest.
        event : str 
            Event type ('onset' or 'withdrawal').
        timewindow : list 
            Time window relative to event [PRE_TIME, POST_TIME].
        exp : str 
            Experiment identifier.
        session : str 
            Session identifier.
        mouse : str 
            Mouse identifier.
        group : str 
            Group identifier.
        trace_color : str (default: 'green') 
            Color for mean response line and fill.
        fill_alpha : float (default: 0.2) 
            Alpha (opacity) for the standard error fill.
        trace_linewidth : float (default: 2) 
            Line width for the mean response line.
        heatmap_cmap : str (default: 'magma') 
            Color map for the heatmap.
    
    --> Returns:
        fig : matplotlib.figure.Figure 
            The figure containing the plots.
    """
    
    # Unpack time window
    PRE_TIME, POST_TIME = timewindow
    
    # Check if PETH_data is valid
    if PETH_data is None or len(PETH_data) == 0:
        raise ValueError("PETH_data is empty or None. Please provide valid PETH data.")
    
    # Create figure and axes
    fig, (ax_heatmap, ax_trace) = plt.subplots(2, 1, figsize=(15, 10), gridspec_kw={'height_ratios': [1, 2]})
    
    # Create time vector for the x-axis
    peri_time = np.arange(-PRE_TIME, POST_TIME + 0.1, 0.1)
    
    # Calculate mean and standard error of the z-scored ΔF/F traces
    mean_dFF_snips = np.mean(PETH_data, axis=0)
    std_dFF_snips = np.std(PETH_data, axis=0) / np.sqrt(len(PETH_data))
    
    ## ----------------- Heatmap Plot ----------------- ##
    im = ax_heatmap.imshow(
        PETH_data, 
        cmap=heatmap_cmap, 
        aspect='auto', 
        interpolation='none', 
        extent=[-PRE_TIME, POST_TIME, len(PETH_data), 0],
        vmin=-3, 
        vmax=4
    )
    ax_heatmap.axvline(x=0, linewidth=2, color='black', linestyle='--', label=f'{event.capitalize()} event')
    ax_heatmap.set_ylabel('Bout Number')
    ax_heatmap.set_yticks(np.arange(0.5, len(PETH_data), 2))
    ax_heatmap.set_yticklabels(np.arange(0, len(PETH_data), 2))
    ax_heatmap.set_title(f'{BOI} {event.capitalize()} - {exp}, {session}, Mouse: {mouse}, Group: {group}')
    
    # Add colorbar
    cbar_ax = fig.add_axes([0.85, 0.54, 0.02, 0.34])  # Custom position for colorbar
    fig.colorbar(im, cax=cbar_ax, label='Z-scored ΔF/F')
    
    ## ----------------- Trace Plot ----------------- ##
    for trial_snip in PETH_data:
        ax_trace.plot(peri_time, trial_snip, linewidth=0.5, color=[0.7, 0.7, 0.7], label='Individual trials')
    
    ax_trace.plot(
        peri_time, 
        mean_dFF_snips, 
        linewidth=trace_linewidth, 
        color=trace_color, 
        label='Mean response'
    )
    
    ax_trace.fill_between(
        peri_time, 
        mean_dFF_snips + std_dFF_snips, 
        mean_dFF_snips - std_dFF_snips, 
        facecolor=trace_color, 
        alpha=fill_alpha, 
        label='Standard error'
    )
    
    ax_trace.axvline(x=0, linewidth=2, color='slategray', linestyle='--', label=f'{event.capitalize()} event')
    
    ax_trace.set_xlabel('Time (s)')
    ax_trace.set_ylabel('Z-scored ΔF/F')
    ax_trace.legend(loc='upper left', fontsize='small')
    ax_trace.margins(0, 0.01)
    
    # Finalize layout
    fig.subplots_adjust(right=0.8, hspace=0.1)
    
    return fig

def plot_PETH_pooled(PETH_array, BOI, event, timewindow, exp, session, group, 
                     trace_color='cornflowerblue', trace_alpha=0.3, fill_alpha=0.5,
                     line_width=1, fill=True):
    """
    Plots PETH averaged over 1 group

    --> Parameters:
        PETH_array : np.ndarray 
            Array of PETH traces (shape: num_bouts x timepoints).
        BOI : str 
            Behavior of interest.
        event : str 
            Event type ('onset' or 'withdrawal').
        timewindow : list 
            Time window relative to event [PRE_TIME, POST_TIME].
        exp : str 
            Experiment identifier.
        session : str 
            Session identifier.
        group : str 
            Group identifier.
        trace_color : str (default: 'cornflowerblue') 
            Color for individual trace lines.
        mean_color : str (default: 'slategrey') 
            Color for the mean response line and fill.
        fill_alpha : float (default: 0.5) 
            Alpha (opacity) for the standard error fill.
        line_width : float (default: 1) 
            Line width for the mean response line.
        fill : bool (default: True) 
            Whether to fill the area around the mean trace with standard error.
    
    --> Returns:
        fig : matplotlib.figure.Figure 
            The figure containing the plot.
    """
    
    # Unpack time window
    PRE_TIME, POST_TIME = float(timewindow[0]), float(timewindow[1])
    
    # Check if PETH_array is valid
    if PETH_array is None or len(PETH_array) == 0:
        raise ValueError("PETH_array is empty or None. Please provide valid PETH data.")
    
    # Determine if the input is a list of traces or or a single trace
    is_single = isinstance(PETH_array[0], (np.float64, float))
    
    if not is_single:
        listmean_dFF_snips = np.mean(PETH_array, axis=0)
        listsem_dFF_snips = np.std(PETH_array, axis=0) / np.sqrt(len(PETH_array))
    else:
        listmean_dFF_snips = PETH_array
        listsem_dFF_snips = np.full(len(PETH_array), 0)  # Default SEM if pre-averaged = 0
        fill = False
    
    # Create figure and axes
    fig, ax = plt.subplots(figsize=(6, 4))
    
    # Create time vector for the x-axis
    peri_time = np.arange(-PRE_TIME, POST_TIME + 0.1, 0.1)
    
    # Ensure peri_time matches the dimensions of PETH_array
    if len(peri_time) != len(listmean_dFF_snips):
        raise ValueError(f"Time vector of length {len(peri_time)} does not match PETH array of length {len(listmean_dFF_snips)}")
    
    ## ----------------- Plot Individual Traces ----------------- ##
    for snip in PETH_array:
        ax.plot(peri_time, snip, linewidth=0.5, color=trace_color, alpha=trace_alpha)
    
    ## ----------------- Plot Mean Trace ----------------- ##
    ax.plot(
        peri_time, 
        listmean_dFF_snips, 
        linewidth=line_width, 
        color=trace_color, 
        label=group
    )
    
    ## ----------------- Fill Between (SEM) ----------------- ##
    if fill:
        ax.fill_between(
            peri_time, 
            listmean_dFF_snips + listsem_dFF_snips, 
            listmean_dFF_snips - listsem_dFF_snips, 
            facecolor=trace_color, 
            alpha=fill_alpha
        )
    
    ## ----------------- Event Line ----------------- ##
    ax.axvline(x=0, linewidth=2, color='slategray', linestyle='--', label=f'{event.capitalize()} Event')
    
    ## ----------------- Axis Labels and Limits ----------------- ##
    ax.set_xlabel('Seconds')
    ax.set_ylabel('z-scored $\Delta$F/F')
    ax.legend(loc='upper left', fontsize='small')
    ax.set_ylim(-2, 4)
    ax.margins(0, 0.1)
    ax.set_title(f'{BOI} - {exp} {session} {group}')
    
    return fig
