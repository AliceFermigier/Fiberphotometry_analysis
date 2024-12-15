# -*- coding: utf-8 -*-
"""
Created on Sun Dec 15 16:53:24 2024

To run fiberphotometry analysis with behaviour or plethysmography data
2 - ANALYSIS WITH BEHAVIOUR BORIS FILE

@author: alice
"""

#%%IMPORTED
###########

import pandas as pd
import numpy as np
from pathlib import Path
import os
import warnings

#import functions
import preprocess as pp
import genplot as gp
import behavplot as bp
import statcalc as sc
import transients as tr

from script_preprocess import experiment_path, exp, ORDER, CUT_FREQ

#%%LOADER
##########

analysis_path = experiment_path / 'Analysis' 
data_path = experiment_path / 'Data'
os.chdir(experiment_path)
os.getcwd() 

#import ID and groups of all mice
subjects_df = pd.read_excel(experiment_path / 'subjects.xlsx', sheet_name='Included')
#import tasks in protocol
proto_df = pd.read_excel(experiment_path / 'protocol.xlsx')

############
#PARAMETERS#
############

#threshold to fuse behaviour if bouts are too close, in secs
THRESH_S = 1
#threshold for PETH : if events are too short do not plot them and do not include them in PETH, in seconds
EVENT_TIME_THRESHOLD = 1
    
#%% 2 - ANALYSIS - BEHAVIOUR
############################

exp_path = analysis_path / exp
data_path_exp = data_path / proto_df.loc[proto_df['Task']==exp, 'Data_path'].values[0]
pp_path = data_path_exp / 'Preprocessing'

# 2.1 - Align with behaviour, create corresponding excel, plot fiberpho data with behaviour

for session_path in [Path(f.path) for f in os.scandir(exp_path) if f.is_dir()]:
    session = session_path.name  # Extract session name
    print('##########################################')
    print(f'EXPERIMENT : {exp} - SESSION : {session}')
    print('##########################################')
    
    # Generate session code for the current session
    code = gp.session_code(session, exp)
    
    # Create repository path where fiberbehav data will be stored
    repo_path = session_path / f'length{EVENT_TIME_THRESHOLD}_interbout{THRESH_S}_o{ORDER}f{CUT_FREQ}'
    repo_path.mkdir(exist_ok=True)  # Create directory if it doesn't exist
    
    for mouse in subjects_df['Subject']:
        print("--------------")
        print(f'MOUSE : {mouse}')
        print("--------------")
        
        # Create common file prefix for files related to the current mouse and session
        file_prefix = f'{mouse}_{code}'
        
        # Define paths for raw, behavioral, and fiberphotometry data
        rawdata_path = data_path_exp / f'{file_prefix}.csv'
        behav_path = data_path_exp / f'behav_{code}_{mouse}.csv'
        fiberpho_path = pp_path / f'{file_prefix}_dFFfilt.csv'
        
        # Paths for output files to be checked
        fiberbehav_path = repo_path / f'{file_prefix}_fiberbehav.csv'
        fiberbehav_notderived_path = repo_path / f'{file_prefix}_fiberbehavnotderived.csv'
        fiberbehav_plot_pdf_path = repo_path / f'{file_prefix}_fiberbehav.pdf'
        fiberbehav_plot_png_path = repo_path / f'{file_prefix}_fiberbehav.png'

        # Determine if we should process this mouse's data
        ready = behav_path.exists()
        done = fiberbehav_path.exists()

        print(f'Ready for analysis? {ready}')
        print(f'Analysis already done? {done}')

        # Proceed only if behavior data exists and analysis hasn't been done
        if ready and not done:
            print(f'Processing behavior alignment for {mouse}...')
            
            # 1 Load raw and behavioral data
            rawdata_cam_df = pd.read_csv(
                rawdata_path, 
                skiprows=1, 
                usecols=['Time(s)', 'DI/O-3'], 
                encoding="ISO-8859-1"
            )
            behav10Sps = pd.read_csv(behav_path)
            fiberpho = pd.read_csv(fiberpho_path, index_col=0)
            
            # Get the sampling rate for the fiberphotometry data
            sr = pp.samplerate(fiberpho)
            
            print(f'Experiment: {exp}, Session: {session}, Mouse: {mouse}')
            
            # 2 List of behaviors of interest (BOI)
            list_BOI = behav10Sps.columns[1:].tolist()
            print(f'Behaviours to analyze: {list_BOI}')
                
            # 3 Align behavior and fiberphotometry data
            timevector = gp.time_vector(fiberpho, sr)
            timestart_camera = gp.timestamp_camera(rawdata_cam_df)[0]
            print(f'Start camera timestamp: {timestart_camera}')
            
            fiberbehav_df = bp.align_behav(behav10Sps, fiberpho, timevector, timestart_camera, exp)
            fiberbehav_df = bp.behav_process(fiberbehav_df, list_BOI, THRESH_S, EVENT_TIME_THRESHOLD, sr)
            
            # 4 Save the intermediate "not derived" data
            fiberbehav_df.to_csv(fiberbehav_notderived_path, index=False)
            
            # 5 Derive and save final processed data
            dfiberbehav_df = bp.derive(fiberbehav_df)
            dfiberbehav_df.to_csv(fiberbehav_path, index=False)
            
            # 6 Plot fiberphotometry data with behavior
            fig_fiberbehav = bp.plot_fiberpho_behav(
                dfiberbehav_df, 
                list_BOI, 
                exp, 
                session, 
                mouse, 
                THRESH_S, 
                EVENT_TIME_THRESHOLD
            )
            
            # Save plots in both PDF and PNG formats
            fig_fiberbehav.savefig(fiberbehav_plot_pdf_path)
            fig_fiberbehav.savefig(fiberbehav_plot_png_path)
            
            print(f'Analysis for mouse {mouse} complete. Data saved in {repo_path}')
                        
#%% 2.2 - Calculate mean, max, and delta dFF within behavioural states (for state behaviours)

for session_path in [Path(f.path) for f in os.scandir(exp_path) if f.is_dir()]:
    session = session_path.name  # Extract session name
    print('##########################################')
    print(f'EXPERIMENT: {exp} - SESSION: {session}')
    print('##########################################')
    
    # Generate session code
    code = gp.session_code(session, exp)
    
    # Create necessary paths
    repo_path = session_path / f'length{EVENT_TIME_THRESHOLD}_interbout{THRESH_S}_o{ORDER}f{CUT_FREQ}'
    groupanalysis_path = repo_path / 'Group_analysis'
    groupanalysis_path.mkdir(exist_ok=True)  # Create directory if it doesn't exist
    
    # Create lists to store data for export
    mean_dFFs_list = []
    diffmeanmaxdFF_list = []
    diffmeanmaxdFF_perbout_list = []
    
    for mouse in subjects_df['Subject']:
        file_prefix = f'{mouse}_{code}'
        
        # Paths to the relevant fiberbehav files
        fiberbehav_path = repo_path / f'{file_prefix}_fiberbehav.csv'
        fiberbehav_notderived_path = repo_path / f'{file_prefix}_fiberbehavnotderived.csv'
        
        # Check if required data exists for this mouse
        if fiberbehav_path.exists():
            print("--------------")
            print(f'PROCESSING MOUSE: {mouse}')
            print("--------------")
            group = subjects_df.loc[subjects_df['Subject'] == mouse, 'Group'].values[0]
            
            # 1 Load data from CSV files
            dfiberbehav_df = pd.read_csv(fiberbehav_path, index_col=0)
            fiberbehav_df = pd.read_csv(fiberbehav_notderived_path, index_col=0)
            
            # 2 Calculate metrics (mean, max, and delta dFF)
            mean_dFF_result = sc.meandFF_behav(list_BOI, dfiberbehav_df, exp, session, mouse, group)
            diffmeanmax_dFF_result = sc.diffmeanmaxdFF_behav(fiberbehav_df, list_BOI, mouse, group)
            diffmeanmax_dFF_perbout_result = sc.diffmeanmaxdFF_behav_perbout(fiberbehav_df, list_BOI, mouse, group)
            
            # 3 Add the results to the list for later export
            mean_dFFs_list.append(mean_dFF_result)
            diffmeanmaxdFF_list.append(diffmeanmax_dFF_result)
    
    # 4 Concatenate results and export to Excel
    if mean_dFFs_list:
        meandFFs_allmice = pd.concat(mean_dFFs_list, ignore_index=True)
        meandFFs_allmice_path = groupanalysis_path / f'{exp}_{session}_globmeandFFs.xlsx'
        meandFFs_allmice.to_excel(meandFFs_allmice_path, index=False)
        print(f'Saved global mean dFFs to {meandFFs_allmice_path}')
    else:
        print('No mean dFF data to export.')

    if diffmeanmaxdFF_list:
        diffdFFs_allmice = pd.concat(diffmeanmaxdFF_list, ignore_index=True)
        diffdFFs_allmice_path = groupanalysis_path / f'{exp}_{session}_diffmaxmeandFFs.xlsx'
        diffdFFs_allmice.to_excel(diffdFFs_allmice_path, index=False)
        print(f'Saved delta max mean dFFs to {diffdFFs_allmice_path}')
    else:
        print('No delta max mean dFF data to export.')
        
#%% 2.3 - Calculate mean and max before and after behaviour onset for whole group
# ----------------------------- #
BOI = 'Shock'  # Behavior of Interest
TIME_MEANMAX = 5  # Time window for mean and max calculation, before and after behaviour (seconds)
TIME_BEFORE = 1  # Time before event to analyze (seconds)
MAX_BOUTS_NUMBER = None  # Limit the number of bouts to analyze, if specified
# ----------------------------- #

for session_path in [Path(f.path) for f in os.scandir(exp_path) if f.is_dir()]:
    session = session_path.name  # Get the session name from the path
    print('##########################################')
    print(f'EXPERIMENT: {exp} - SESSION: {session}')
    print('##########################################')
    
    code = gp.session_code(session, exp)
    repo_path = session_path / f'length{EVENT_TIME_THRESHOLD}_interbout{THRESH_S}_o{ORDER}f{CUT_FREQ}'
    groupanalysis_path = repo_path / 'Group_analysis'
    
    if len(os.listdir(repo_path)) > 1:
        if not groupanalysis_path.exists():
            groupanalysis_path.mkdir()
        
        # Initialize lists to store data
        meandFF_list = []
        meandFF_before_list = []
        subject_list = []
        group_list = []
        maxdFF_list = []
        maxdFF_before_list = []
        
        for mouse in subjects_df['Subject']:
            print("--------------")
            print(f'MOUSE: {mouse}')
            print("--------------")
            
            fiberbehav_file = repo_path / f'{mouse}_{code}_fiberbehav.csv'
            
            if not fiberbehav_file.exists():
                print(f"File not found: {fiberbehav_file}")
                continue
            
            try:
                # Read in the CSV file
                fiberbehav_df = pd.read_csv(fiberbehav_file, index_col=0)
                group = subjects_df.loc[subjects_df['Subject'] == mouse, 'Group'].values[0]
                sr = pp.samplerate(fiberbehav_df)  # Calculate the sample rate
                
                if BOI in fiberbehav_df.columns[2:].tolist():
                    ind_event_list = np.where(fiberbehav_df[BOI] == 1)[0].tolist()
                    
                    if MAX_BOUTS_NUMBER is not None:
                        ind_event_list = ind_event_list[:MAX_BOUTS_NUMBER]  # Limit the number of bouts if specified

                    for ind_event in ind_event_list:
                        try:
                            # Identify the time of the event (index of minimum dFF before onset)
                            ind_event_adjusted = fiberbehav_df.loc[ind_event - TIME_BEFORE * sr : ind_event, 'Denoised dFF'].idxmin()
                            
                            # Calculate mean and max dFF after the event
                            dFF_after_event = fiberbehav_df.loc[ind_event_adjusted : ind_event_adjusted + TIME_MEANMAX * sr, 'Denoised dFF']
                            mean_dFF_entry = dFF_after_event.mean() if not dFF_after_event.empty else np.nan
                            max_dFF_entry = dFF_after_event.max() if not dFF_after_event.empty else np.nan
                            
                            # Calculate mean and max dFF before the event
                            dFF_before_event = fiberbehav_df.loc[ind_event_adjusted - TIME_MEANMAX * sr : ind_event_adjusted, 'Denoised dFF']
                            mean_dFF_before_entry = dFF_before_event.mean() if not dFF_before_event.empty else np.nan
                            max_dFF_before_entry = dFF_before_event.max() if not dFF_before_event.empty else np.nan
                            
                            # Append results to lists
                            subject_list.append(mouse)
                            group_list.append(group)
                            meandFF_list.append(mean_dFF_entry)
                            meandFF_before_list.append(mean_dFF_before_entry)
                            maxdFF_list.append(max_dFF_entry)
                            maxdFF_before_list.append(max_dFF_before_entry)
                        
                        except Exception as e:
                            print(f"Error processing event at index {ind_event} for mouse {mouse}: {e}")
                            continue
            
            except Exception as e:
                print(f"Error processing file {fiberbehav_file}: {e}")
                continue
        
        # Export data to Excel
        output_file = groupanalysis_path / f'{exp}_{session}_{TIME_MEANMAX}sbeforeafter_maxmeandFF.xlsx'
        
        if not output_file.exists():
            try:
                meanmaxdFFs_df = pd.DataFrame({
                    'Subject': subject_list,
                    'Group': group_list,
                    'Mean dFF before entry': meandFF_before_list,
                    'Mean dFF entry': meandFF_list,
                    'Max dFF before entry': maxdFF_before_list,
                    'Max dFF entry': maxdFF_list
                })
                meanmaxdFFs_df.to_excel(output_file, index=False)
                print(f"Data successfully exported to {output_file}")
            except Exception as e:
                print(f"Error exporting data to {output_file}: {e}")

#%% 2.4 - Plot PETH for each mouse

# PETH parameters 
EVENT_LIST = ['onset']  # Event triggers, e.g., onset, withdrawal
TIME_WINDOWS = [[5, 30]]  # Time window for PETH calculation (pre, post)

# Loop over each session folder in the experiment path
for session_path in [Path(f.path) for f in os.scandir(exp_path) if f.is_dir()]:
    session = session_path.name
    print('##########################################')
    print(f'EXPERIMENT : {exp} - SESSION : {session}')
    print('##########################################')

    # Generate the session code
    code = gp.session_code(session, exp)
    
    # Create the repository and PETH paths
    repo_path = session_path / f'length{EVENT_TIME_THRESHOLD}_interbout{THRESH_S}_o{ORDER}f{CUT_FREQ}'
    peth_path = repo_path / 'PETH'
    peth_path.mkdir(parents=True, exist_ok=True)  # Create directory if it doesn't exist
    
    # Loop over each mouse in the subjects DataFrame
    for mouse in subjects_df['Subject']:
        group = subjects_df.loc[subjects_df['Subject'] == mouse, 'Group'].values[0]
        fiberbehav_path = repo_path / f'{mouse}_{code}_fiberbehav.csv'

        if fiberbehav_path.exists():  # Check if fiber behavior file exists for this mouse
            print("--------------")
            print(f'MOUSE : {mouse}')
            print("--------------")

            try:
                # Read the fiber behavior file
                dfiberbehav_df = pd.read_csv(fiberbehav_path, index_col=0)
            except Exception as e:
                warnings.warn(f"Failed to read file {fiberbehav_path}: {e}")
                continue

            # List all behaviors of interest (BOI) by excluding specific behaviors
            behaviors_of_interest = [col for col in dfiberbehav_df.columns[4:] 
                                     if col not in ['Entry in arena', 'Gate opens', 'Tail suspension']]
            
            for behavior in behaviors_of_interest:
                for event, time_window in zip(EVENT_LIST, TIME_WINDOWS):  
                    try:
                        # Generate the PETH data for the current behavior, event, and time window
                        peth_data = bp.PETH(dfiberbehav_df, behavior, event, time_window, EVENT_TIME_THRESHOLD)
                        
                        # Create a DataFrame from the PETH data
                        time_index = np.arange(-time_window[0], time_window[1] + 0.1, 0.1)
                        peth_df = pd.DataFrame(np.transpose(peth_data), index=time_index)
                        
                        # Plot the PETH and save the figure 
                        peth_plot = bp.plot_PETH(peth_data, behavior, event, time_window, exp, session, mouse, group)
                        plot_filename = f'{mouse}_{code}_{behavior}_{event[0]}{time_window[0] - time_window[1]}_PETH.png'
                        peth_plot_path = peth_path / plot_filename
                        peth_plot.savefig(peth_plot_path)
                        peth_plot.close()
                        
                    except Exception as e:
                        warnings.warn(f"Error while processing PETH for {mouse}, {behavior}, {event}: {e}")
                        
#%% 2.5 Plot PETH for each group and extract mean and max Z-scored data

# ----------------------------- #
# Parameters
BOI = 'Shock'
TIME_WINDOW = [5, 30]  # In seconds
MAXBOUTSNUMBER = 1
# ----------------------------- #

for session_path in [Path(f.path) for f in os.scandir(exp_path) if f.is_dir()]:
    session = session_path.name  # Extract session name
    print('##########################################')
    print(f'EXPERIMENT: {exp} - SESSION: {session}')
    print('##########################################')
    
    code = gp.session_code(session, exp)
    repo_path = session_path / f'length{EVENT_TIME_THRESHOLD}_interbout{THRESH_S}_o{ORDER}f{CUT_FREQ}'
    peth_path = repo_path / 'PETH'
    peth_path.mkdir(parents=True, exist_ok=True)  # Create PETH directory if it doesn't exist

    # Initialize data storage lists
    subject_list = []
    group_list = []
    PETH_array = None
    PETH_mean_list = []
    PETH_max_list = []

    # Loop over each subject (mouse)
    for mouse in subjects_df['Subject']:
        print("--------------")
        print(f'MOUSE: {mouse}')
        print("--------------")
        
        fiberbehav_file = repo_path / f'{mouse}_{code}_fiberbehav.csv'
        
        if not fiberbehav_file.exists():
            print(f"File not found: {fiberbehav_file}")
            continue
        
        try:
            fiberbehav_df = pd.read_csv(fiberbehav_file, index_col=0)
            group = subjects_df.loc[subjects_df['Subject'] == mouse, 'Group'].values[0]
            sr = pp.samplerate(fiberbehav_df)
            
            if BOI in fiberbehav_df.columns[2:].tolist():
                subject_list.append(mouse)
                group_list.append(group)
                print(f'PETH {BOI} for {mouse}')
                
                # Calculate PETH for the current mouse
                PETH_mouse = bp.PETH(fiberbehav_df, BOI, 'onset', TIME_WINDOW, EVENT_TIME_THRESHOLD, maxboutsnumber=MAXBOUTSNUMBER)
                
                if PETH_array is None:
                    PETH_array = PETH_mouse
                    print('Initialized PETH_array successfully')
                else:
                    PETH_array = np.concatenate((PETH_array, PETH_mouse))  # Stack new data
                
                # Calculate mean and max dFF before and after the event (PETH)
                mean_before = np.mean(PETH_mouse[:TIME_WINDOW[0]])  # Mean before event
                mean_after = np.mean(PETH_mouse[TIME_WINDOW[0]:])   # Mean after event
                max_before = np.max(PETH_mouse[:TIME_WINDOW[0]])    # Max before event
                max_after = np.max(PETH_mouse[TIME_WINDOW[0]:])     # Max after event
                
                PETH_mean_list.append((mean_before, mean_after))
                PETH_max_list.append((max_before, max_after))

        except Exception as e:
            print(f"Error processing {fiberbehav_file}: {e}")
            continue

    # Export mean/max PETH data to Excel
    meanmaxPETH_df = pd.DataFrame({
        'Subject': subject_list,
        'Group': group_list,
        f'Mean dFF before {BOI}': [x[0] for x in PETH_mean_list],
        f'Mean dFF after {BOI}': [x[1] for x in PETH_mean_list],
        f'Max dFF before {BOI}': [x[0] for x in PETH_max_list],
        f'Max dFF after {BOI}': [x[1] for x in PETH_max_list]
    })
    meanmaxPETH_df.to_excel(peth_path / f'{BOI}_{TIME_WINDOW[0]}_{TIME_WINDOW[1]}_PETHmeanmax.xlsx')

    # Plot PETH for each group
    included_groups = ['CD', 'HFD']
    for group in included_groups:
        # Filter PETH data for the current group
        group_indices = [i for i, g in enumerate(group_list) if g == group]
        PETH_array_group = PETH_array[group_indices]
        
        print(f"Group {group} PETH data size: {PETH_array_group.shape}")

        # Plot pooled PETH for the group
        fig_PETHpooled = bp.plot_PETH_pooled(PETH_array_group, BOI, 'onset', TIME_WINDOW, exp, session, group)
        fig_PETHpooled.savefig(peth_path / f'{group}_{BOI}_{TIME_WINDOW[1]}_PETH.pdf')
        fig_PETHpooled.savefig(peth_path / f'{group}_{BOI}_{TIME_WINDOW[1]}_PETH.png')

#%% 2.6 - Compute variance and transients on whole trace and pre/post baseline
# ----------------------------- #
# Parameters
TIME_BEGIN = 15  # In seconds
lowcut = 0.0001  # Lowcut frequency for bandpass filter (Hz)
highcut = 0.5    # Highcut frequency for bandpass filter (Hz)
# ----------------------------- #

var_transients_list = []

# Loop over each subject (mouse)
for mouse in subjects_df['Subject']:
    print("--------------")
    print(f'MOUSE: {mouse}')
    print("--------------")
    
    # Set file path for the fiber behavior CSV
    fiber_file = repo_path / f'{mouse}_{code}_fiberbehav.csv'
    
    # Check if the file exists
    if not fiber_file.exists():
        print(f"File not found: {fiber_file}")
        continue

    try:
        # Load data
        dfiber_df = pd.read_csv(fiber_file)
        
        # Crop time at the beginning to avoid artifacts
        dfiber_df = dfiber_df[TIME_BEGIN*sr:].reset_index(drop=True)
        
        # Apply bandpass filter to the dFF signal
        dfiber_df['Filtered dFF'] = tr.bandpass_filter(dfiber_df['Denoised dFF'], lowcut, highcut, sr)
        
        # Plot signal and spectrum
        tr.plot_signal_and_spectrum(dfiber_df['Time(s)'], dfiber_df['Denoised dFF'], dfiber_df['Filtered dFF'], sr)
        
        # Retrieve group information
        group = subjects_df.loc[subjects_df['Subject'] == mouse, 'Group'].values[0]
        
        #Calculate variance and transients characteristics
        mouse_df = sc.variance_transients(dfiber_df, mouse, group, exp, session)
        var_transients_list.append(mouse_df)

    except Exception as e:
        print(f"Error processing {fiber_file}: {e}")
        continue
    
# Concatenate results and export to Excel
variability_df = pd.concat(var_transients_list, ignore_index=True)
output_file = groupanalysis_path / f'Variability_1o{ORDER}f{lowcut}_{highcut}.xlsx'
variability_df.to_excel(output_file, index=False)

print(f"Variability data saved to {output_file}")