# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 17:16:16 2023

@author: alice fermigier
"""

#%%
##########
#IMPORTED#
##########

import pandas as pd
import numpy as np

import modules.common.transients as tr
import modules.common.preprocess as pp

#%%
###################
#DEFINED FUNCTIONS#
###################

def meandFF_behav(list_BOI, fiberbehav_df, exp, session, mouse, group, batch):
    """
    Calculates mean dFF during each behavior of interest (BOI).
    Outputs a DataFrame with 'Subject', 'Group', 'Mean dFF', 'Baseline', 'Post_baseline', and mean dFF for each behavior.
    
    Args:
    - list_BOI (list): List of behaviors of interest to analyze.
    - fiberbehav_df (DataFrame): Dataframe containing behavior onset/offsets and 'dFF' signal.
    - exp (str): The experiment type (e.g., 'Fear', 'NewContext', etc.).
    - session (str): The session type (e.g., 'Test', 'Conditioning', etc.).
    - mouse (str): The mouse identifier.
    - group (str): The group identifier.
    
    Returns:
    - DataFrame: A dataframe containing Subject, Group, Mean dFF, Baseline, Post_baseline, and individual behavior dFFs.
    """
    
    # Filter the DataFrame to ignore the first 15 seconds
    fiberbehavsnip_df = fiberbehav_df[fiberbehav_df['Time(s)'] > 15]
    
    # Get the index of when the trial begins (gate opens or equivalent)
    try:
        ind_start_trial = fiberbehavsnip_df.index[fiberbehavsnip_df['Entry in arena'] == 1].tolist()[0]
    except IndexError:
        print(f"Warning: Start trial not found for Mouse {mouse} in Experiment {exp} - Session {session}. Defaulting to 0.")
        ind_start_trial = 0

    # Initialize the DataFrame to hold the calculated mean dFFs
    results = {
        'Batch': batch,
        'Subject': mouse,
        'Group': group
    }

    # --- Calculate mean dFF for each behavior ---
    list_behav_analyzed = []
    for behav in list_BOI:
        if behav in fiberbehav_df.columns[2:].tolist():
            if fiberbehavsnip_df[behav].sum() > 2:  # Check if behavior has enough activity
                try:
                    # Calculate mean dFF for the behavior
                    meandFF_behav_df = fiberbehavsnip_df.groupby(behav, as_index=False)['dFF'].mean()
                    mean_dFF_value = meandFF_behav_df.loc[meandFF_behav_df[behav] == 1, 'dFF'].values
                    if mean_dFF_value.size > 0:
                        results[behav] = mean_dFF_value[0]
                        list_behav_analyzed.append(behav)
                    else:
                        results[behav] = np.nan
                except Exception as e:
                    print(f"Error processing behavior '{behav}' for Mouse {mouse}: {e}")
                    results[behav] = np.nan
            else:
                results[behav] = np.nan
        else:
            results[behav] = np.nan

    # --- Calculate baseline dFF (before trial starts) ---
    try:
        meandFF_baseline = fiberbehavsnip_df.loc[:ind_start_trial, 'dFF'].mean()
    except Exception as e:
        print(f"Error calculating baseline dFF for Mouse {mouse}: {e}")
        meandFF_baseline = np.nan

    results['Baseline'] = meandFF_baseline

    # --- Calculate post-baseline dFF (after gate opens, excluding exploration) ---
    try:
        meandFF_postbaseline_df = fiberbehavsnip_df.loc[ind_start_trial:]
        for behav in list_behav_analyzed:
            meandFF_postbaseline_df = meandFF_postbaseline_df.loc[meandFF_postbaseline_df[behav] == 0]

        meandFF_postbaseline = meandFF_postbaseline_df['dFF'].mean()
    except Exception as e:
        print(f"Error calculating post-baseline dFF for Mouse {mouse}: {e}")
        meandFF_postbaseline = np.nan

    results['Post_baseline'] = meandFF_postbaseline

    # --- Calculate mean dFF when no exploration happens during the total trial ---
    try:
        meandFF_df = fiberbehavsnip_df
        for behav in list_behav_analyzed:
            meandFF_df = meandFF_df.loc[meandFF_df[behav] == 0]

        meandFF = meandFF_df['dFF'].mean()
    except Exception as e:
        print(f"Error calculating mean dFF for Mouse {mouse}: {e}")
        meandFF = np.nan

    results['Mean dFF'] = meandFF

    # --- Create a DataFrame with the results ---
    meandFFs_df = pd.DataFrame([results])
    
    return meandFFs_df

def diffmeanmaxdFF_behav(behavprocess_df, list_BOI, mouse, group, batch):
    """
    Calculates mean dFF difference between the beginning and end of bouts.
    Also calculates mean mean and max dFF across all bouts (will give a different value than meandFF_behav : here 1 bout = 1 value)
    
    Args:
    - behavprocess_df (DataFrame): Dataframe containing behavior onset/offsets and 'dFF' signal.
    - list_BOI (list): List of behaviors of interest to analyze.
    - mouse (str): The mouse identifier.
    - group (str): The group identifier.
    
    Returns:
    - DataFrame: A dataframe containing Subject, Group, Behavior, Mean dFF, Max dFF, and Delta dFF for each behavior.
    """
    
    # Filter out specific behaviors that shouldn't be analyzed
    list_behav_analyzed = [behav for behav in list_BOI if behav not in ['Entry in arena', 'Gate opens', 'Shock', 'Tail suspension']]
    
    # Initialize result lists
    results = {'Batch': [],
               'Subject': [], 
               'Group': [], 
               'Behaviour': [], 
               'Mean dFF': [], 
               'Max dFF': [], 
               'Delta dFF': []}
    
    for behav in list_behav_analyzed:
        
        # Check if behavior column exists in dataframe
        if behav in behavprocess_df.columns[2:]:
            
            # Identify start and stop indices for the current behavior
            list_starts = np.where(behavprocess_df[behav] == 1)[0].tolist()
            list_stops = np.where(behavprocess_df[behav] == -1)[0].tolist()
            
            # Check for matching pairs of starts and stops
            if len(list_starts) != len(list_stops):
                print(f'Warning: Mismatched start/stop counts for behavior "{behav}" (Mouse: {mouse})')
            
            list_deltadFF_behav = []
            list_meandFF_behav = []
            list_maxdFF_behav = []
            
            for (start, stop) in zip(list_starts, list_stops):
                try:
                    # Handle potential index out-of-bound issues
                    start_window = max(0, start - 5)
                    stop_window = min(len(behavprocess_df) - 1, stop + 5)
                    
                    # Calculate means at the start and stop of the bout (±5 indices around start/stop)
                    mean_start = behavprocess_df.loc[start_window:start + 5, 'dFF'].mean()
                    mean_stop = behavprocess_df.loc[stop - 5:stop_window, 'dFF'].mean()
                    
                    # Calculate delta dFF
                    delta_dFF = mean_stop - mean_start
                    
                    # Calculate mean and max dFF within the bout
                    mean_dFF = behavprocess_df.loc[start:stop, 'dFF'].mean()
                    max_dFF = behavprocess_df.loc[start:stop, 'dFF'].max()
                    
                    list_deltadFF_behav.append(delta_dFF)
                    list_meandFF_behav.append(mean_dFF)
                    list_maxdFF_behav.append(max_dFF)
                
                except Exception as e:
                    print(f"Error processing bout for behavior '{behav}' in Mouse {mouse}: {e}")
            
            # Aggregate mean values for the behavior
            mean_deltadFF = np.nanmean(list_deltadFF_behav) if list_deltadFF_behav else np.nan
            mean_meandFF = np.nanmean(list_meandFF_behav) if list_meandFF_behav else np.nan
            mean_maxdFF = np.nanmean(list_maxdFF_behav) if list_maxdFF_behav else np.nan
            
        else:
            # If behavior column doesn't exist, set NaN values
            mean_deltadFF = np.nan
            mean_meandFF = np.nan
            mean_maxdFF = np.nan
        
        # Store results for this behavior
        results['Batch'].append(batch)
        results['Subject'].append(mouse)
        results['Group'].append(group)
        results['Behaviour'].append(behav)
        results['Mean dFF'].append(mean_meandFF)
        results['Max dFF'].append(mean_maxdFF)
        results['Delta dFF'].append(mean_deltadFF)
    
    # Create a DataFrame from the results dictionary
    diffdFF_df = pd.DataFrame(results)
    
    return diffdFF_df

def diffmeanmaxdFF_behav_perbout(behavprocess_df, list_BOI, mouse, group, batch):
    """
    Same as diffmeanmaxdFF_behav but with detailed data for each bout

    Args:
    - behavprocess_df (DataFrame): DataFrame containing bout onsets/offsets and 'dFF'.
    - list_BOI (list): List of behaviors of interest to analyze.
    - mouse (str): The mouse identifier.
    - group (str): The group identifier.

    Returns:
    - DataFrame: A dataframe containing Subject, Group, Behaviour, Bout number, Mean dFF, Max dFF, and Delta dFF for each bout.
    """
    
    # Filter out behaviors that should not be analyzed
    list_behav_analyzed = [behav for behav in list_BOI if behav not in ['Entry in arena', 'Gate opens', 'Shock']]
    
    # Initialize lists to store bout-related dFF data
    list_deltadFF = []
    list_meandFF = []
    list_maxdFF = []
    listind_behav = []
    listind_bouts = []
    
    for behav in list_behav_analyzed:
        # Identify start and stop indices of the behavior
        list_starts = np.where(behavprocess_df[behav] == 1)[0].tolist()
        list_stops = np.where(behavprocess_df[behav] == -1)[0].tolist()
        
        if len(list_starts) == 0 or len(list_stops) == 0:
            print(f"Warning: No starts or stops found for behavior '{behav}' in Mouse {mouse}. Skipping this behavior.")
            continue
        
        # Ensure starts and stops are paired properly
        if len(list_starts) != len(list_stops):
            print(f"Warning: Mismatched start/stop pairs for behavior '{behav}' in Mouse {mouse}.")
        
        bout_n = 1  # Counter for each bout of the behavior
        for start_idx, stop_idx in zip(list_starts, list_stops):
            try:
                # Calculate mean dFF around start and stop (±5 frames) while handling index overflow
                start_window = behavprocess_df.loc[max(0, start_idx-5):min(start_idx+5, len(behavprocess_df)-1), 'dFF']
                stop_window = behavprocess_df.loc[max(0, stop_idx-5):min(stop_idx+5, len(behavprocess_df)-1), 'dFF']
                
                mean_start = start_window.mean() if not start_window.empty else np.nan
                mean_stop = stop_window.mean() if not stop_window.empty else np.nan
                
                # Calculate delta dFF (difference from start to stop)
                delta_dFF = mean_stop - mean_start if not (np.isnan(mean_start) or np.isnan(mean_stop)) else np.nan
                
                # Calculate mean and max dFF during the behavior bout
                if start_idx < stop_idx:  # To avoid inverted intervals
                    bout_window = behavprocess_df.loc[start_idx:stop_idx, 'dFF']
                    mean_dFF = bout_window.mean() if not bout_window.empty else np.nan
                    max_dFF = bout_window.max() if not bout_window.empty else np.nan
                else:
                    print(f"Warning: Start index ({start_idx}) is greater than Stop index ({stop_idx}) for behavior '{behav}' in Mouse {mouse}.")
                    mean_dFF = np.nan
                    max_dFF = np.nan
                
                # Store the results
                list_deltadFF.append(delta_dFF)
                list_meandFF.append(mean_dFF)
                list_maxdFF.append(max_dFF)
                listind_behav.append(behav)
                listind_bouts.append(bout_n)
                
                bout_n += 1  # Increment the bout number

            except Exception as e:
                print(f"Error processing bout for behavior '{behav}' (Mouse {mouse}, Bout {bout_n}): {e}")
                list_deltadFF.append(np.nan)
                list_meandFF.append(np.nan)
                list_maxdFF.append(np.nan)
                listind_behav.append(behav)
                listind_bouts.append(bout_n)
                bout_n += 1

    # Create DataFrame from collected lists
    diffdFF_df = pd.DataFrame({
        'Batch': [batch] * len(listind_behav),
        'Subject': [mouse] * len(listind_behav),
        'Group': [group] * len(listind_behav),
        'Behaviour': listind_behav,
        'Bout': listind_bouts,
        'Mean dFF': list_meandFF,
        'Max dFF': list_maxdFF,
        'Delta dFF': list_deltadFF
    })
    
    return diffdFF_df

def variance_transients(fiberbehav_df, list_BOI, mouse, group, exp, batch, threshold):
    """
    Calculates variance, transient frequency, and amplitude during whole trace and during behaviours
    
    Output:
    - A dataframe with variance, transient frequency and amplitude.
    """
    
    # Calculate variance during whole trace, baseline and post-baseline periods
    variance = np.var(fiberbehav_df['dFF'])
    
    # Calculate transients for whole trace, baseline and post-baseline periods
    peaks_df, peak_frequency, peak_amplitude, transients_fig = tr.transients(fiberbehav_df, threshold)

    # Store the results in a dataframe
    results_df = pd.DataFrame({
        'Batch': batch,
        'Subject': mouse,
        'Group': group,
        'Variance': variance,
        'Transients Frequency': peak_frequency,
        'Transients Amplitude': peak_amplitude,
    }, index=[0])

    if list_BOI != []:
        for behavior in list_BOI:
            if behavior not in fiberbehav_df.columns:
                print(f"Warning: Behavior '{behavior}' not found in DataFrame.")
                continue

            # Find rows where behavior is active
            behavior_mask = fiberbehav_df[behavior] == 1

            # Peaks during that behavior
            behavior_peaks = peaks_df[(peaks_df['Peaks'] == 1) & behavior_mask]

            # Duration of behavior period in seconds
            behavior_time = fiberbehav_df.loc[behavior_mask, 'Time(s)']
            if behavior_time.empty:
                duration = np.nan
            else:
                duration = behavior_time.iloc[-1] - behavior_time.iloc[0]

            # Calculate frequency and mean amplitude during behavior
            if duration and len(behavior_peaks) > 0:
                freq = len(behavior_peaks) / duration
                amp = behavior_peaks['Filtered dFF'].mean()
            else:
                freq = np.nan
                amp = np.nan

            # Add to results
            results_df[f'{behavior} Transients Frequency'] = freq
            results_df[f'{behavior} Transients Amplitude'] = amp
    
    return results_df, transients_fig

def process_event(fiberbehav_df, ind_event, TIME_MEANMAX):
    """Calculate mean and max dFF before and after a behavioral event."""
    try:
        sr=pp.samplerate(fiberbehav_df)
        # Calculate mean and max dFF after the event
        dFF_after_event = fiberbehav_df.loc[ind_event : ind_event + TIME_MEANMAX * sr, 'dFF']
        mean_dFF_after = dFF_after_event.mean()
        max_dFF_after = dFF_after_event.max()
        
        # Calculate mean and max dFF before the event
        dFF_before_event = fiberbehav_df.loc[ind_event - TIME_MEANMAX * sr : ind_event, 'dFF']
        mean_dFF_before = dFF_before_event.mean()
        max_dFF_before = dFF_before_event.max()
        
        return mean_dFF_after, max_dFF_after, mean_dFF_before, max_dFF_before
    except Exception as e:
        print(f"Error processing event at index {ind_event}: {e}")
        return None, None, None, None
    
def process_mouse_meanmax_beforeafter(fiberbehav_df, mouse, group, BOI, batch, MAX_BOUTS_NUMBER=None):
    """Process the fiberbehav file for a single mouse, calculating mean and max dFF for events."""
    results = pd.DataFrame(columns=['Subject', 'Group', 'Mean dFF before entry', 
                                    'Mean dFF entry', 'Max dFF before entry', 'Max dFF entry'])
    
    try:
        sr = pp.samplerate(fiberbehav_df)  # Sample rate
        
        if BOI not in fiberbehav_df.columns[2:]:
            print(f'{BOI} not in data for mouse {mouse}')
            return results  # Return empty DataFrame
        
        ind_event_list = np.where(fiberbehav_df[BOI] == 1)[0].tolist()
        print(f"Events for mouse {mouse}: {ind_event_list}")
        
        if MAX_BOUTS_NUMBER is not None:
            ind_event_list = ind_event_list[:MAX_BOUTS_NUMBER]  # Limit the number of bouts if specified

        for ind_event in ind_event_list:
            mean_dFF_entry, max_dFF_entry, mean_dFF_before_entry, max_dFF_before_entry = process_event(fiberbehav_df, ind_event, sr)
            
            if mean_dFF_entry is not None:
                results = pd.DataFrame({
                    'Batch': batch,
                    'Subject': mouse,
                    'Group': group,
                    f'Mean dFF before {BOI}': mean_dFF_before_entry,
                    f'Mean dFF after {BOI}': mean_dFF_entry,
                    f'Max dFF before {BOI}': max_dFF_before_entry,
                    f'Max dFF after {BOI}': max_dFF_entry
                }, index=[0])
        
        return results
    except Exception as e:
        print(f"Error processing mouse {mouse}: {e}")
        return results

def meandFF_sniffs(fiberbehav_df, exp, session, mouse, group, batch, joined=True):
    """
    Calculates mean dFF (delta F/F) during each behavior for a given mouse and session.

    Parameters:
    -----------
    fiberbehav_df : DataFrame
        DataFrame with fiber photometry and behavioral data.
    exp : str
        Experiment identifier.
    session : str
        Session identifier.
    mouse : str
        Subject/mouse identifier.
    group : str
        Group identifier (e.g., control or experimental).
    joined : bool, optional (default=True)
        If True, combines 'stim' and 'sniff' columns with the same odor name 
        (removes the last 2 characters to group them together).
    
    Returns:
    --------
    meandFFs_df : DataFrame
        DataFrame with columns ['Subject', 'Group', 'Mean dFF for behavior 1', 'Mean dFF for behavior 2', ...]
        Contains the mean dFF for each behavior.
    """
    
    # List to store mean dFF for each behavior
    mean_dFF_list = []
    
    # Exclude the first 15 seconds from the data (likely to avoid transient artifacts)
    fiberbehavsnip_df = fiberbehav_df[fiberbehav_df['Time(s)'] > 15]
    
    # If `joined=True`, group 'stim' and 'sniff' with the same odor together
    if joined:
        fiberbehavsnip_df = fiberbehavsnip_df.rename(
            columns=lambda col: col[:-2] if col not in ['Time(s)', 'dFF'] else col
        )
        fiberbehavsnip_df = fiberbehavsnip_df.groupby(level=0, axis=1, sort=False).sum()
        print(f'Joined columns: {fiberbehavsnip_df.columns}')  # Debugging print to verify joined columns
    
    # Calculate the mean dFF for each behavior
    for behavior in fiberbehavsnip_df.columns[3:]:  # Skip first 3 columns: ['Time(s)', 'dFF', ...]
        # Group by the binary presence of the behavior and compute the mean
        mean_dFF_behav_df = fiberbehavsnip_df.groupby([behavior], as_index=False).mean()
        
        # Extract the mean dFF where the behavior is present (behavior == 1)
        mean_dFF = mean_dFF_behav_df.loc[mean_dFF_behav_df[behavior] == 1, 'dFF']
        
        if not mean_dFF.empty:
            mean_dFF_list.append(mean_dFF.values[0])
        else:
            mean_dFF_list.append(np.nan)  # Handle case where no behavior == 1 occurs
            
    # Create the final output DataFrame
    dFF_values = [batch, mouse, group] + mean_dFF_list
    dFF_columns = ['Batch', 'Subject', 'Group'] + list(fiberbehavsnip_df.columns[3:])
    
    meandFFs_df = pd.DataFrame(data=[dFF_values], columns=dFF_columns)
    
    return meandFFs_df

def meanmax_dFF_stims(behavprocess_df, list_BOI, mouse, group, TIME_MEANMAX, sr, batch):
    """
    Calculates the difference in dFF between the beginning and end of behavioral bouts.
    Also calculates the mean and max dFF during each bout.

    Parameters:
    -----------
    behavprocess_df : DataFrame
        DataFrame containing behavioral information and dFF signal.
    list_BOI : list
        List of behaviors of interest (columns to analyze).
    mouse : str
        Subject/mouse identifier.
    group : str
        Group identifier (e.g., control or experimental).
    TIME_MEANMAX : int
        Time window (in seconds) used for calculating mean and max dFF.
    sr : int
        Sampling rate (in Hz) of the signal.
    
    Returns:
    --------
    meanmaxdFF_df : DataFrame
        DataFrame with the following columns:
        ['Subject', 'Group', 'Behaviour', 'Mean dFF before', 'Mean dFF during', 
         'Max dFF before', 'Max dFF during']
    """
    
    # Initialize lists to store analysis results
    analyzed_behaviors = []  # Tracks which behaviors are analyzed
    mean_dFF_before = []     # Mean dFF before the behavior onset
    mean_dFF_after = []      # Mean dFF during the behavior
    max_dFF_before = []      # Maximum dFF before the behavior onset
    max_dFF_after = []       # Maximum dFF during the behavior
    
    # Loop through each behavior in the list of Behaviors of Interest (BOI)
    for behavior in list_BOI:
        if behavior not in behavprocess_df.columns:
            print(f"Warning: Behavior '{behavior}' not found in DataFrame columns. Skipping...")
            continue
        
        # Identify the start index of the behavior (where behavior == 1)
        behavior_starts = np.where(behavprocess_df[behavior] == 1)[0]
        
        if len(behavior_starts) == 0:
            print(f"No occurrences of behavior '{behavior}' found. Skipping...")
            continue
        
        start_idx = behavior_starts[0]  # Take the first occurrence of the behavior
        
        # Calculate valid index ranges, ensuring no out-of-bounds errors
        start_before = max(0, start_idx - TIME_MEANMAX * sr)
        start_after = start_idx
        end_after = min(len(behavprocess_df), start_idx + TIME_MEANMAX * sr)
        
        # Extract dFF data for before and after behavior onset
        dFF_before = behavprocess_df.iloc[start_before:start_after]['dFF']
        dFF_after = behavprocess_df.iloc[start_after:end_after]['dFF']
        
        # Calculate statistics for before and after
        mean_dFF_before.append(dFF_before.mean() if not dFF_before.empty else np.nan)
        max_dFF_before.append(dFF_before.max() if not dFF_before.empty else np.nan)
        mean_dFF_after.append(dFF_after.mean() if not dFF_after.empty else np.nan)
        max_dFF_after.append(dFF_after.max() if not dFF_after.empty else np.nan)
        
        # Track the analyzed behavior
        analyzed_behaviors.append(behavior)
    
    # Create a DataFrame to store results
    meanmaxdFF_df = pd.DataFrame({
        'Batch': [batch] * len(analyzed_behaviors),
        'Subject': [mouse] * len(analyzed_behaviors),
        'Group': [group] * len(analyzed_behaviors),
        'Behaviour': analyzed_behaviors,
        'Mean dFF before': mean_dFF_before,
        'Mean dFF during': mean_dFF_after,
        'Max dFF before': max_dFF_before,
        'Max dFF during': max_dFF_after
    })
    
    return meanmaxdFF_df
