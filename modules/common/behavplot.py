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
from pathlib import Path
import sys
import json

import modules.common.preprocess as pp

from scripts.loader import project_root

#%%
###################
#DEFINED FUNCTIONS#
###################

def correct_time_behav(fiberpho_raw, behav_df):
    first_timestamp = fiberpho_raw['Time(s)'][0]
    behav_df['Time(s)']=(behav_df['time'].values)+first_timestamp

    return behav_df

def align_dlc_to_fiber(fiberpho_df, dlc_df, time_col="Time(s)"):
    """
    Interpolates DLC (x,y) coordinates onto the fiber photometry time axis.
    """
    fp_times = fiberpho_df[time_col].values
    dlc_times = dlc_df[time_col].values
    
    aligned = fiberpho_df.copy()

    # For each DLC coordinate column
    for col in dlc_df.columns:
        if col in ["Time(s)"]:
            continue  # skip timestamp
        series = pd.Series(dlc_df[col].values, index=dlc_times)
        aligned[col] = np.interp(fp_times, dlc_times, series)
    
    return aligned

def align_behav(behav_df, fiberpho, list_BOI):
    """
    Aligns fiber photometry data with behavioral data from Boris or DLC on a time vector.
    """

    [start,stop]=[behav_df['Time(s)'].values[0],  behav_df['Time(s)'].values[-1]]
    behav_time = fiberpho.loc[(fiberpho['Time(s)'] >= start) & (fiberpho['Time(s)'] <= stop), 'Time(s)']

    pad_begin = np.empty(len(fiberpho.loc[fiberpho['Time(s)'] < start]), dtype=float)
    pad_end = np.empty(len(fiberpho['Time(s)'])-(len(behav_time)+len(pad_begin)), dtype=float)

    for col in behav_df.columns[1:]:
        if col in list_BOI:
            pad_begin.fill(0.0)
            pad_end.fill(0.0)
        else:
            pad_end.fill(np.nan)
            pad_begin.fill(np.nan)
        fiberpho[col] = np.concatenate(
            [pad_begin,
            np.interp(behav_time.values, behav_df['Time(s)'].values, behav_df[col].values),
            pad_end]
            )
    
    return fiberpho

def behav_process(df, list_BOI, THRESH_S, EVENT_TIME_THRESHOLD):
    sr = pp.samplerate(df)

    for BOI in list_BOI:
        if BOI not in df.columns:
            print(f"[!] BOI '{BOI}' not in dataframe")
            continue
        x = df[BOI].round().values.astype(int)

        # --- 1. Detect starts and ends of bouts ---
        diff = np.diff(np.r_[0, x, 0])
        starts = np.where(diff == 1)[0]
        ends   = np.where(diff == -1)[0]

        bouts = list(zip(starts, ends))  
        if len(bouts) == 0:
            print(f"[i] No bouts detected for {BOI}")
            continue

        # --- 2. Merge bouts separated by < THRESH_S seconds ---
        merged = []
        prev_start, prev_end = bouts[0]

        for start, end in bouts[1:]:
            gap = (start - prev_end) / sr
            if gap <= THRESH_S:
                # merge with previous
                prev_end = end
            else:
                merged.append((prev_start, prev_end))
                prev_start, prev_end = start, end

        merged.append((prev_start, prev_end))

        # --- 3. Remove short bouts ---
        cleaned = [
            (s, e) for (s, e) in merged 
            if (e - s) / sr >= EVENT_TIME_THRESHOLD
        ]

        # --- 4. Rewrite the BOI column ---
        new_x = np.zeros_like(x)
        for s, e in cleaned:
            new_x[s:e] = 1

        df[BOI] = new_x

    return df

def derive(fiberbehav_df, list_BOI):
    """
    Calculate the derivative of behaviors of interest and store in the same DataFrame.
    The result will show 1 when behavior starts and -1 when it stops.
    """
    for col in list_BOI:
        # Ensure values are only 0 or 1
        fiberbehav_df[col] = fiberbehav_df[col].apply(lambda x: 1 if x == 1 else 0)
        
        # Compute difference to detect transitions
        fiberbehav_df[col] = fiberbehav_df[col].diff().fillna(0)

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

def plot_fiberpho_behav(behavprocess_df, list_BOI, exp, mouse, THRESH_S, EVENT_TIME_THRESHOLD, batch, scaled=True):
    """
    Plots denoised deltaF/F aligned with behaviour (includes baseline). Adds Speed subplot only if present.
    """
    behavprocesssnip_df = behavprocess_df[behavprocess_df['Time(s)'] > 2]
    has_speed = 'Speed' in behavprocesssnip_df.columns
    has_560 = 'Denoised 560 dFF' in behavprocesssnip_df.columns

    if has_speed and has_560:
        fig = plt.figure(figsize=(20, 15))
        ax1 = fig.add_subplot(311)
    
    elif has_speed:
        fig = plt.figure(figsize=(20, 10))
        ax1 = fig.add_subplot(211)

    else:
        fig = plt.figure(figsize=(20, 5))
        ax1 = fig.add_subplot(111)       

    # Plot dFF trace
    ax1.plot('Time(s)', 'Denoised dFF', linewidth=1, color='black', label='_GCaMP', data=behavprocesssnip_df)
    if has_speed and has_560:
        ax2 = fig.add_subplot(312)
        ax2.plot('Time(s)', 'Denoised 560 dFF', linewidth=1, color='black', label='560 dFF', data=behavprocesssnip_df)

    # Highlight behaviors
    behavior_colors_path = Path(project_root) / "modules/behaviour/behaviour_colors.json"
    with open(behavior_colors_path, "r") as f:
        behaviors_to_plot = json.load(f)
    
    for behavior in list_BOI:
        if behavior in behavprocesssnip_df.columns:
            color, alpha = behaviors_to_plot.get(behavior, ('grey',0.05))
            highlight_behavior_areas(ax1, behavprocesssnip_df, behavior, color, alpha)
            if has_speed and has_560:
                highlight_behavior_areas(ax2, behavprocesssnip_df, behavior, color, alpha)

    # Add event lines
    for event, color, label in [('Gate opens', 'lightsteelblue', 'Gate opens'),
                                ('Entry in arena', 'slategrey', 'Entry in arena'),
                                ('Airpuffs', 'lime', 'Airpuffs')]:
        if event in list_BOI and event in behavprocesssnip_df.columns:
            event_indices = np.where(behavprocess_df[event] == 1)[0]
            if len(event_indices) > 0:
                for i, event_index in enumerate(event_indices):
                    x = behavprocess_df.at[event_index, 'Time(s)']
                    ax1.axvline(x, color=color, ls='--', label='_'*i + label)

    # Labels and formatting
    fs_mult = 4
    ax1.set_ylabel(r'$\Delta$F/F', fontsize=5 * fs_mult)
    ax1.set_xlabel('Time(s)', fontsize=5 * fs_mult)
    ax1.set_title(f'dFF with Behavioural Scoring - {exp} {mouse} {batch} - interbout {THRESH_S} - cut {EVENT_TIME_THRESHOLD}',
                  fontsize=5 * fs_mult)
    ax1.tick_params(axis='both', labelsize=4 * fs_mult)
    ax1.legend(loc='upper right', fontsize=4 * fs_mult)
    ax1.margins(0, 0.2)
    if scaled:
        ax1.set_ylim([-0.27, 0.75])
    
    if has_speed and has_560:
    # Labels and formatting
        fs_mult = 4
        ax2.set_ylabel(r'$\Delta$F/F', fontsize=5 * fs_mult)
        ax2.set_xlabel('Time(s)', fontsize=5 * fs_mult)
        ax2.tick_params(axis='both', labelsize=4 * fs_mult)
        ax2.legend(loc='upper right', fontsize=4 * fs_mult)
        ax2.margins(0, 0.2)
        if scaled:
            ax2.set_ylim([-0.27, 0.75])

    # Plot speed if available
    if has_speed and has_560:
        ax3 = fig.add_subplot(313)
        ax3.plot('Time(s)', 'Speed', linewidth=1, color='black', label='Speed', data=behavprocesssnip_df)
        ax3.set_ylabel('Speed (cm/s)', fontsize=5 * fs_mult)
        ax3.set_xlabel('Time(s)', fontsize=5 * fs_mult)
        ax3.tick_params(axis='both', labelsize=4 * fs_mult)
        ax3.margins(0, 0.2)
        if scaled:
            ax3.set_ylim([-1, 50])
    elif has_speed:
        ax2 = fig.add_subplot(212)
        ax2.plot('Time(s)', 'Speed', linewidth=1, color='black', label='Speed', data=behavprocesssnip_df)
        ax2.set_ylabel('Speed (cm/s)', fontsize=5 * fs_mult)
        ax2.set_xlabel('Time(s)', fontsize=5 * fs_mult)
        ax2.tick_params(axis='both', labelsize=4 * fs_mult)
        ax2.margins(0, 0.2)
        if scaled:
            ax2.set_ylim([-1, 50])

    plt.tight_layout()
    return fig

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
    #min_event_duration = EVENT_TIME_THRESHOLD * sr
    #valid_bouts = [(start, end) for start, end in zip(list_ind_event_o, list_ind_event_w) if 1 < end - start >= min_event_duration]
    #list_ind_event_o, list_ind_event_w = zip(*valid_bouts) if valid_bouts else ([], [])

    # Choose the relevant event indices to align on (either onset or withdrawal)
    list_ind_event = list_ind_event_o if event == 'onset' else list_ind_event_w

    # Check if the event happens too late in the dataframe to process
    list_ind_event = [idx for idx in list_ind_event if idx + POST_TIME * sr < len(behavprocess_df)]

    # For onset events, adjust the start index to the minimum dFF within 1 second before and after the event
    #if event == 'onset':
       # list_ind_event = [
  #          behavprocess_df.loc[idx - 1 * sr : idx + 1 * sr, 'Denoised dFF'].idxmin() 
        #    for idx in list_ind_event
      #  ]

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
                dFF_baseline = behavprocess_df.loc[ind_event - 1 * sr : ind_event - PRE_EVENT_TIME * sr, 'Denoised dFF']
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

def plot_PETH(PETH_data, BOI, event, timewindow, exp, mouse, group, 
              trace_color='black', fill_alpha=0.2, trace_linewidth=2, heatmap_cmap='magma'):
    import matplotlib.pyplot as plt
    import numpy as np

    # Unpack time window
    PRE_TIME, POST_TIME = timewindow

    # Check if PETH_data is valid
    if PETH_data is None or len(PETH_data) == 0:
        raise ValueError("PETH_data is empty or None. Please provide valid PETH data.")

    # Text size multiplier
    text_size = 20  # 4x typical 12 pt font

    # Create figure and axes
    fig, (ax_heatmap, ax_trace) = plt.subplots(2, 1, figsize=(15, 10), gridspec_kw={'height_ratios': [1, 2]})

    # Create time vector for the x-axis
    peri_time = np.linspace(-PRE_TIME, POST_TIME, PETH_data.shape[1])

    # Calculate mean and standard error of the z-scored ΔF/F traces
    mean_dFF_snips = np.mean(PETH_data, axis=0)
    std_dFF_snips = np.std(PETH_data, axis=0) / np.sqrt(len(PETH_data))

    ## ----------------- Heatmap Plot ----------------- ##
    vmin = np.min(PETH_data)
    vmax = np.max(PETH_data)
    im = ax_heatmap.imshow(
        PETH_data,
        cmap=heatmap_cmap,
        aspect='auto',
        interpolation='none',
        extent=[-PRE_TIME, POST_TIME, len(PETH_data), 0],
        vmin=vmin ,
        vmax=vmax
    )
    ax_heatmap.axvline(x=0, linewidth=2, color='black', linestyle='--', label=f'{event.capitalize()} event')
    ax_heatmap.set_ylabel('Bout Number', fontsize=text_size)
    ax_heatmap.set_yticks(np.arange(0.5, len(PETH_data), 2))
    ax_heatmap.set_yticklabels(np.arange(0, len(PETH_data), 2), fontsize=text_size * 0.8)
    ax_heatmap.set_title(f'{BOI} {event.capitalize()} - {exp}, Mouse: {mouse}, Group: {group}', fontsize=text_size)
    ax_heatmap.set_xticks([])
    ax_heatmap.set_xticklabels([])
    ax_heatmap.set_xlabel('')

    # Add colorbar
    cbar_ax = fig.add_axes([0.85, 0.54, 0.02, 0.34])  # Custom position for colorbar
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label('Z-scored ΔF/F', fontsize=text_size)
    cbar.ax.tick_params(labelsize=text_size * 0.8)

    ## ----------------- Trace Plot ----------------- ##
    for trial_snip in PETH_data:
        ax_trace.plot(peri_time, trial_snip, linewidth=0.5, color=[0.7, 0.7, 0.7])
    
    # Add single legend entry for individual trials
    ax_trace.plot([], [], linewidth=0.5, color=[0.7, 0.7, 0.7], label='Individual trials')

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

    ax_trace.set_xlabel('Time (s)', fontsize=text_size)
    ax_trace.set_ylabel('Z-scored ΔF/F', fontsize=text_size)
    ax_trace.tick_params(labelsize=text_size * 0.8)
    ax_trace.legend(loc='upper left', fontsize=text_size * 0.6)
    ax_trace.margins(0, 0.01)

    # Finalize layout
    fig.subplots_adjust(right=0.8, hspace=0.1)

    return fig

def plot_PETH_pooled(PETH_array, BOI, event, timewindow, exp, group, 
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
    peri_time = np.linspace(-PRE_TIME, POST_TIME, len(listmean_dFF_snips))
    
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
    ax.set_ylabel(r'z-scored $\Delta$F/F')
    ax.legend(loc='upper left', fontsize='small')
    ax.set_ylim(-1, 8)
    ax.margins(0, 0.1)
    ax.set_title(f'{BOI} - {exp} {group}')
    
    return fig
