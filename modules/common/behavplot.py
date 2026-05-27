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
    
    # Trim fiberpho to DLC time window
    mask = (fp_times >= dlc_times[0]) & (fp_times <= dlc_times[-1])
    aligned = fiberpho_df[mask].copy()
    fp_times = aligned[time_col].values

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
    Binary BOI columns use nearest-neighbor interpolation to preserve clean 0/1 values.
    Continuous columns (dFF, speed, etc.) use linear interpolation.
    """
    start, stop = behav_df['Time(s)'].values[0], behav_df['Time(s)'].values[-1]
    behav_time  = fiberpho.loc[
        (fiberpho['Time(s)'] >= start) & (fiberpho['Time(s)'] <= stop), 'Time(s)'
    ]

    n_before = len(fiberpho.loc[fiberpho['Time(s)'] < start])
    n_after  = len(fiberpho) - (len(behav_time) + n_before)

    pad_begin = np.empty(n_before, dtype=float)
    pad_end   = np.empty(n_after,  dtype=float)

    behav_times_arr = behav_df['Time(s)'].values
    fiber_times_arr = behav_time.values

    for col in behav_df.columns[1:]:
        if col in list_BOI:
            pad_begin.fill(0.0)
            pad_end.fill(0.0)

            # ── Nearest-neighbor: find the closest behav frame for each fiber frame ──
            indices      = np.searchsorted(behav_times_arr, fiber_times_arr, side='left')
            indices      = np.clip(indices, 0, len(behav_df) - 1)
            interpolated = behav_df[col].values[indices].astype(float)

        else:
            pad_begin.fill(np.nan)
            pad_end.fill(np.nan)

            # ── Linear interpolation for continuous signals ────────────────────────
            interpolated = np.interp(fiber_times_arr, behav_times_arr, behav_df[col].values)

        fiberpho[col] = np.concatenate([pad_begin, interpolated, pad_end])

    return fiberpho

def behav_process(df, list_BOI, THRESH_S, EVENT_TIME_THRESHOLD):
    sr = pp.samplerate(df)

    for BOI in list_BOI:
        if BOI not in df.columns:
            print(f"  [!] BOI '{BOI}' not in dataframe")
            continue

        x = df[BOI].round().values.astype(int)

        # --- 1. Detect starts and ends of bouts ---
        diff   = np.diff(np.r_[0, x, 0])
        starts = np.where(diff ==  1)[0]
        ends   = np.where(diff == -1)[0]
        bouts  = list(zip(starts, ends))
        print(f"{BOI} — raw bouts : {len(bouts)} ")

        if len(bouts) == 0:
            print(f"  [i] No bouts detected for {BOI}")
            continue

        # --- 2. Merge bouts separated by < THRESH_S seconds ---
        merged = []
        prev_start, prev_end = bouts[0]
        for start, end in bouts[1:]:
            gap = (start - prev_end) / sr
            if gap <= THRESH_S:
                prev_end = end
            else:
                merged.append((prev_start, prev_end))
                prev_start, prev_end = start, end
        merged.append((prev_start, prev_end))
        print(f"{BOI} — after merging (THRESH_S={THRESH_S}s): {len(merged)} bouts")

        # --- 3. Remove short bouts ---
        cleaned = [
            (s, e) for (s, e) in merged
            if (e - s) / sr >= EVENT_TIME_THRESHOLD
        ]
        print(f"{BOI} — after filtering (EVENT_TIME_THRESHOLD={EVENT_TIME_THRESHOLD}s): "
              f"{len(cleaned)} bouts")

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
    derived_df = fiberbehav_df.copy()  
    for col in list_BOI:
        # Ensure values are only 0 or 1
        derived_df[col] = derived_df[col].apply(lambda x: 1 if x == 1 else 0)
        
        # Compute difference to detect transitions
        derived_df[col] = derived_df[col].diff().fillna(0)

    return derived_df

def highlight_behavior_areas(ax, df, behavior_name, facecolor='grey', alpha=0.3, label_prefix=''):
    i = 0
    x_start = None
    ones  = (df[behavior_name] == 1).sum()
    minus = (df[behavior_name] == -1).sum()
    print(f"  highlight_behavior_areas: {behavior_name} → {ones} starts, {minus} ends")  # ← debug
    for x, y in zip(df['Time(s)'].tolist(), df[behavior_name].tolist()):
        if y == 1:
            x_start = x
        if y == -1 and x_start is not None:
            ax.axvspan(x_start, x, facecolor=facecolor, alpha=alpha,
                       label='_' * i + label_prefix + behavior_name)
            x_start = None
            i += 1
    print(f"  → {i} spans drawn") 

def plot_fiberpho_behav(behavprocess_df, list_BOI, exp, mouse, THRESH_S, EVENT_TIME_THRESHOLD, batch, scaled=True):
    behavprocesssnip_df = behavprocess_df.dropna()
    has_speed = 'Speed' in behavprocesssnip_df.columns
    has_560   = '560 dFF' in behavprocesssnip_df.columns

    # ── 1. Create ALL axes upfront ───────────────────────────────────────────
    if has_speed and has_560:
        print('Plotting 465 dFF, 560 dFF and speed')
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(20, 15))
    elif has_560 or has_speed:
        nrows = 2
        label = '465 dFF and ' + ('560 dFF' if has_560 else 'speed')
        print(f'Plotting {label}')
        fig, (ax1, ax2) = plt.subplots(nrows, 1, figsize=(20, 10))
        ax3 = None
    else:
        print('Plotting 465 dFF')
        fig, ax1 = plt.subplots(1, 1, figsize=(20, 5))
        ax2 = ax3 = None

    # ── 2. Plot traces ───────────────────────────────────────────────────────
    ax1.plot('Time(s)', 'dFF', linewidth=1, color='black', label='465 dFF',
             data=behavprocesssnip_df)

    if has_560:
        ax2.plot('Time(s)', '560 dFF', linewidth=1, color='black', label='560 dFF',
                 data=behavprocesssnip_df)

    speed_ax = None
    if has_speed and has_560:
        speed_ax = ax3
    elif has_speed:
        speed_ax = ax2
    if speed_ax is not None:
        speed_ax.plot('Time(s)', 'Speed', linewidth=1, color='black', label='Speed',
                      data=behavprocesssnip_df)

    # ── 3. Highlight behaviors on ALL axes ───────────────────────────────────
    behavior_colors_path = Path(project_root) / "modules/behaviour/behaviour_colors.json"
    with open(behavior_colors_path, "r") as f:
        behaviors_to_plot = json.load(f)

    all_axes = [ax for ax in [ax1, ax2, ax3] if ax is not None]

    for behavior in list_BOI:
        if behavior in behavprocesssnip_df.columns:
            color, alpha = behaviors_to_plot.get(behavior, ('grey', 0.05))
            for ax in all_axes:
                highlight_behavior_areas(ax, behavprocesssnip_df, behavior, color, alpha)
        else:
            print(f'{behavior} not found in data')

    # ── 4. Event lines on ax1 ────────────────────────────────────────────────
    for event, color, label in [('Gate opens',     'lightsteelblue', 'Gate opens'),
                                 ('Entry in arena', 'slategrey',      'Entry in arena'),
                                 ('Airpuffs',       'lime',           'Airpuffs')]:
        if event in list_BOI and event in behavprocesssnip_df.columns:
            event_indices = np.where(behavprocess_df[event] == 1)[0]
            for i, idx in enumerate(event_indices):
                x = behavprocess_df.at[idx, 'Time(s)']
                ax1.axvline(x, color=color, ls='--', label='_' * i + label)

    # ── 5. Formatting ────────────────────────────────────────────────────────
    fs_mult = 4
    ax1.set_ylabel(r'$\Delta$F/F', fontsize=5 * fs_mult)
    ax1.set_xlabel('Time(s)', fontsize=5 * fs_mult)
    ax1.set_title(
        f'dFF with Behavioural Scoring - {exp} {mouse} {batch}'
        f' - interbout {THRESH_S} - cut {EVENT_TIME_THRESHOLD}',
        fontsize=5 * fs_mult)
    ax1.tick_params(axis='both', labelsize=4 * fs_mult)
    ax1.legend(loc='upper right', fontsize=4 * fs_mult)
    ax1.margins(0, 0.2)
    if scaled:
        ax1.set_ylim([-0.27, 0.75])

    if has_560 and ax2 is not None:
        ax2.set_ylabel(r'$\Delta$F/F', fontsize=5 * fs_mult)
        ax2.set_xlabel('Time(s)', fontsize=5 * fs_mult)
        ax2.tick_params(axis='both', labelsize=4 * fs_mult)
        ax2.legend(loc='upper right', fontsize=4 * fs_mult)
        ax2.margins(0, 0.2)
        if scaled:
            ax2.set_ylim([-0.27, 0.75])

    if speed_ax is not None:
        speed_ax.set_ylabel('Speed (cm/s)', fontsize=5 * fs_mult)
        speed_ax.set_xlabel('Time(s)', fontsize=5 * fs_mult)
        speed_ax.tick_params(axis='both', labelsize=4 * fs_mult)
        speed_ax.legend(loc='upper right', fontsize=4 * fs_mult)
        speed_ax.margins(0, 0.2)
        if scaled:
            speed_ax.set_ylim([-1, 50])

    plt.tight_layout()
    return fig

def PETH(behavprocess_df, BOI, event, timewindow, EVENT_TIME_THRESHOLD, 
         PRE_EVENT_TIME=0, maxboutsnumber=None, baselinewindow=False,
         dFF_column = 'dFF'):
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

    # Preallocate the PETH array to store the z-scored traces
    n_bouts = len(list_ind_event)
    n_timepoints = (POST_TIME + PRE_TIME) * sr + 1
    PETH_array = np.zeros((n_bouts, n_timepoints))
    
    # Initialize mean and std on whole trace
    F0 = behavprocess_df[dFF_column].mean()
    std0 = behavprocess_df[dFF_column].std()

    # Loop through each event and extract the fiberpho trace centered on the event
    for i, ind_event in enumerate(list_ind_event):
        try: 
            if baselinewindow:
                # Calculate baseline mean (F0) and standard deviation (std0) for the time window before the event
                dFF_baseline = behavprocess_df.loc[ind_event - 1 * sr : ind_event - PRE_EVENT_TIME * sr, dFF_column]
                F0 = dFF_baseline.mean() 
                std0 = dFF_baseline.std()

            # Extract the fiberpho trace for the time window around the event
            event_window = behavprocess_df.loc[ind_event - PRE_TIME * sr : ind_event + POST_TIME * sr, dFF_column]
            
            # Ensure the event window has the correct length to avoid shape mismatch
            if len(event_window) == n_timepoints:
                PETH_array[i] = (event_window - F0) / std0
        except Exception as e:
            print(f"Error processing event at index {ind_event}: {e}")

    return PETH_array

def plot_PETH(PETH_data, BOI, event, timewindow, exp, batch, mouse, group, ylim = None,
              trace_color='black', fill_alpha=0.2, trace_linewidth=2, heatmap_cmap='RdBu_r',
              dff_column = '465'):
    
    # Unpack time window
    PRE_TIME, POST_TIME = timewindow

    # Check if PETH_data is valid
    if PETH_data is None or len(PETH_data) == 0:
        raise ValueError("PETH_data is empty or None. Please provide valid PETH data.")

    # Text size multiplier
    text_size = 38  # 4x typical 12 pt font

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
    ax_heatmap.set_ylabel('Bout #', fontsize=text_size)
    ax_heatmap.set_yticks([0.5, len(PETH_data) - 0.5])
    ax_heatmap.set_yticklabels([1, len(PETH_data)], fontsize=text_size * 0.9)
    ax_heatmap.set_title(f'{BOI} {event.capitalize()} - {exp}, Mouse: {mouse}, Batch: {batch}, Group: {group}', fontsize=text_size*0.6)
    ax_heatmap.set_xticks([])
    ax_heatmap.set_xticklabels([])
    ax_heatmap.set_xlabel('')

    # Add colorbar
    cbar_ax = fig.add_axes([0.85, 0.54, 0.02, 0.34])  # Custom position for colorbar
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label(f'Z-scored {dff_column} ΔF/F', fontsize=text_size)
    cbar.ax.tick_params(labelsize=text_size * 0.9)

    ## ----------------- Trace Plot ----------------- ##
    for trial_snip in PETH_data:
        ax_trace.plot(peri_time, trial_snip, linewidth=0.5, color=[0.7, 0.7, 0.7])
    
    # Add single legend entry for individual trials
    ax_trace.plot([], [], linewidth=0.5, color=[0.7, 0.7, 0.7], label='_Individual trials')

    ax_trace.plot(
        peri_time,
        mean_dFF_snips,
        linewidth=trace_linewidth,
        color=trace_color,
        label='_Mean response'
    )

    ax_trace.fill_between(
        peri_time,
        mean_dFF_snips + std_dFF_snips,
        mean_dFF_snips - std_dFF_snips,
        facecolor=trace_color,
        alpha=fill_alpha,
        label='_Standard error'
    )

    ax_trace.axvline(x=0, linewidth=2, color='slategray', linestyle='--', label=f'{event.capitalize()} {BOI}')

    ax_trace.set_xlabel('Time (s)', fontsize=text_size)
    ax_trace.set_ylabel(f'Z-scored {dff_column} ΔF/F', fontsize=text_size)
    ax_trace.tick_params(labelsize=text_size * 0.9)
    ax_trace.legend(loc='upper left', fontsize=text_size * 0.9)
    ax_trace.margins(0, 0.01)
    if ylim != None:
        ax_trace.set_ylim(ylim[0],ylim[1])

    # Finalize layout
    fig.subplots_adjust(right=0.8, hspace=0.1)

    return fig

def plot_PETH_pooled(PETH_array, BOI, event, timewindow, exp, group, ylim=None,
                     trace_color='cornflowerblue', trace_alpha=0.3, fill_alpha=0.5,
                     line_width=1, fill=True, dff_column='465'):
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
    ax.set_xlabel('Time(s)')
    ax.set_ylabel(f'z-scored {dff_column} ΔF/F')
    ax.legend(loc='upper right', fontsize='medium')
    if ylim != None:
        ax.set_ylim(ylim[0],ylim[1])
    ax.margins(0, 0.1)
    ax.set_title(f'{BOI} - {exp} {group}')
    
    return fig

def remove_first_bout(dfiberbehav_df, behavior):
    """Zero out the first bout of a behavior to exclude it from PETH."""
    dfiberbehav_df = dfiberbehav_df.copy()
    diff = dfiberbehav_df[behavior]
    onsets  = diff[diff == 1].index
    offsets = diff[diff == -1].index

    if len(onsets) == 0:
        return dfiberbehav_df  # No bout found, return unchanged

    first_onset = onsets[0]
    # Find the first offset that comes after the first onset
    subsequent_offsets = offsets[offsets > first_onset]
    first_offset = subsequent_offsets[0] if len(subsequent_offsets) > 0 else dfiberbehav_df.index[-1]

    dfiberbehav_df.loc[first_onset:first_offset, behavior] = 0
    return dfiberbehav_df

def PETH_by_bout(PETH_list, max_bouts=None, min_mice=1, step=1):
    """
    Reorganise per-mouse PETH arrays by bout number and compute
    group mean ± SEM for each bout position across mice.

    Parameters
    ----------
    PETH_list : list of np.ndarray, each shape (n_bouts_i, timepoints)
    max_bouts : int, optional
    min_mice  : int
        Minimum number of mice required for a position to be kept.
    step : int
        Number of consecutive bouts to group together.
        step=1 keeps individual bouts (default behaviour).
        If max_bouts is not a multiple of step, trailing bouts are
        dropped and a warning is raised.

    Returns
    -------
    bout_means : np.ndarray, shape (n_positions, timepoints)
    bout_sems  : np.ndarray, shape (n_positions, timepoints)
    bout_n     : np.ndarray, shape (n_positions,)
    """
    import warnings

    if max_bouts is None:
        max_bouts = max(len(p) for p in PETH_list)

    timepoints = PETH_list[0].shape[1]

    # ── Handle step truncation ────────────────────────────────────────────────
    remainder = max_bouts % step
    if remainder != 0:
        print(
            f"max_bouts={max_bouts} is not a multiple of step={step}. "
            f"The last {remainder} bout(s) will be ignored "
            f"({max_bouts // step} groups of {step} instead of {max_bouts})."
        )
    n_positions = max_bouts // step

    bout_means = np.full((n_positions, timepoints), np.nan)
    bout_sems  = np.full((n_positions, timepoints), np.nan)
    bout_n     = np.zeros(n_positions, dtype=int)

    for pos_idx in range(n_positions):
        bout_indices = range(pos_idx * step, (pos_idx + 1) * step)

        # Per-mouse: average across the bouts in this group that the mouse has
        mouse_means = []
        for peth in PETH_list:
            traces = [peth[b] for b in bout_indices if len(peth) > b]
            if traces:
                mouse_means.append(np.mean(traces, axis=0))

        n = len(mouse_means)
        bout_n[pos_idx] = n
        if n >= min_mice:
            arr = np.array(mouse_means)          # (n_mice, timepoints)
            bout_means[pos_idx] = arr.mean(axis=0)
            bout_sems[pos_idx]  = arr.std(axis=0) / np.sqrt(n)

    return bout_means, bout_sems, bout_n


def plot_PETH_by_bout(bout_means, bout_sems, bout_n,
                      BOI, event, timewindow, exp, group,
                      max_bouts_to_show=None, min_mice=1, step=1,
                      cmap_name='plasma', ylim=None, fill_alpha=0.15,
                      dff_column='465'):

    PRE_TIME, POST_TIME = float(timewindow[0]), float(timewindow[1])
    peri_time = np.linspace(-PRE_TIME, POST_TIME, bout_means.shape[1])

    valid = np.where(bout_n >= min_mice)[0]
    if max_bouts_to_show is not None:
        valid = valid[:max_bouts_to_show]
    if len(valid) == 0:
        raise ValueError(f"No bout positions have >= {min_mice} mice.")

    n_valid = len(valid)
    cmap    = plt.cm.get_cmap(cmap_name, n_valid)
    colors  = [cmap(i) for i in range(n_valid)]

    def _label(rank):
        i = valid[rank]
        n = bout_n[i]
        if step == 1:
            return f'Bout #{i + 1}  (n={n})'
        else:
            return f'Bouts {i*step + 1}–{(i+1)*step}  (n={n})'

    # ── Layout: 2 rows × 2 cols; col 1 is narrow colorbar slot ──────────────
    fig = plt.figure(figsize=(13, 8))
    gs  = fig.add_gridspec(
        2, 2,
        width_ratios  = [20, 1],
        height_ratios = [1,  2],
        hspace = 0.05,          # tight vertical gap → x-axes visually touch
        wspace = 0.05,
    )
    ax_hm   = fig.add_subplot(gs[0, 0])
    ax_cbar = fig.add_subplot(gs[0, 1])
    ax_tr   = fig.add_subplot(gs[1, 0], sharex=ax_hm)   # ← shared x-axis
    fig.add_subplot(gs[1, 1]).set_visible(False)         # empty corner

    # ── Heatmap ───────────────────────────────────────────────────────────────
    heatmap_data = bout_means[valid]
    abs_max = np.nanmax(np.abs(heatmap_data))

    im = ax_hm.imshow(
        heatmap_data,
        cmap='RdBu_r', aspect='auto', interpolation='none',
        extent=[-PRE_TIME, POST_TIME, n_valid + 0.5, 0.5],
        vmin=-abs_max, vmax=abs_max
    )
    ax_hm.axvline(x=0, linewidth=1.5, color='black', linestyle='--')
    ax_hm.set_ylabel('Bout group' if step > 1 else 'Bout #', fontsize=14)
    ax_hm.set_yticks(np.arange(1, n_valid + 1))
    ax_hm.set_yticklabels([_label(r) for r in range(n_valid)], fontsize=9)
    ax_hm.set_title(
        f'{BOI} {event.capitalize()} by bout — {exp}, {group}'
        + (f'  [step={step}]' if step > 1 else ''),
        fontsize=14
    )
    plt.setp(ax_hm.get_xticklabels(), visible=False)    # hidden: shared with ax_tr

    cbar = fig.colorbar(im, cax=ax_cbar)                # colorbar in its own column
    cbar.set_label(f'Z-scored {dff_column} ΔF/F', fontsize=11)
    cbar.ax.tick_params(labelsize=9)

    # ── Traces ────────────────────────────────────────────────────────────────
    for rank, bout_idx in enumerate(valid):
        mean  = bout_means[bout_idx]
        sem   = bout_sems[bout_idx]
        color = colors[rank]

        ax_tr.plot(peri_time, mean, color=color, linewidth=1.5,
                   label=_label(rank))
        ax_tr.fill_between(peri_time, mean - sem, mean + sem,
                           color=color, alpha=fill_alpha)

    ax_tr.axvline(x=0, linewidth=1.5, color='slategray', linestyle='--',
                  label=f'{event.capitalize()} {BOI}')
    ax_tr.axhline(y=0, linewidth=0.8, color='black', linestyle=':')
    ax_tr.set_xlabel('Time (s)', fontsize=14)
    ax_tr.set_ylabel(f'Z-scored {dff_column} ΔF/F', fontsize=14)
    ax_tr.tick_params(labelsize=11)
    ax_tr.legend(loc='upper right', fontsize=9, ncol=2)
    ax_tr.margins(0, 0.05)
    if ylim is not None:
        ax_tr.set_ylim(ylim)

    return fig