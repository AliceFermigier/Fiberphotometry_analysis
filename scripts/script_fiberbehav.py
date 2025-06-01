# -*- coding: utf-8 -*-
"""
Created on Sun Dec 15 16:53:24 2024

To run fiberphotometry analysis with behaviour or plethysmography data
2 - ANALYSIS WITH BEHAVIOUR BORIS FILE

@author: alice fermigier
"""

#%%IMPORTED
###########

import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import os
import warnings
import json
import importlib

#import functions
import modules.common.preprocess as pp
importlib.reload(pp)
import modules.common.genplot as gp
importlib.reload(gp)
import modules.common.behavplot as bp
importlib.reload(bp)
import modules.common.statcalc as sc
importlib.reload(sc)
import modules.common.transients as tr
importlib.reload(tr)
import modules.common.nomenclature as nom
importlib.reload(nom)
import modules.behaviour.mouse_position as mp
importlib.reload(mp)
import modules.behaviour.epm as epm
importlib.reload(epm)
import modules.behaviour.camera_processing as cp
importlib.reload(cp)
import modules.common.clean_signal as cs
importlib.reload(cs)

from scripts.loader import analysis_path, data_path, exp, ORDER, CUT_FREQ, proto_df, subjects_df, THRESH_S, EVENT_TIME_THRESHOLD, batches

#%% 2 - ANALYSIS - BEHAVIOUR
############################
automated_alignment = True
arena_analysis = False
dlc_data = False
boris = True

exp = 'OF_1uL_Ctrl_2'
exp_path = analysis_path / exp
datapath_exp_dict = nom.get_experiment_data_path(batches, proto_df, data_path, exp)

# 2.1 - Align with behaviour, create corresponding excel, plot fiberpho data with behaviour
print('###################')
print(f'EXPERIMENT : {exp}')
print('###################')

# Create repository path where fiberbehav data will be stored
repo_path = exp_path / f'length{EVENT_TIME_THRESHOLD}_interbout{THRESH_S}_o{ORDER}f{CUT_FREQ}'
repo_path.mkdir(exist_ok=True)  # Create directory if it doesn't exist

# Loop through each mouse in the subject DataFrame
for mouse, batch in zip(subjects_df['Subject'], subjects_df['Batch']):
    print("-----------------------------") 
    print(f'BATCH : {batch}, MOUSE : {mouse}')
    print("-----------------------------")
    try:
        data_path_exp = datapath_exp_dict[batch]
        pp_path = data_path_exp / 'Preprocessing'
        behav_path_exp = data_path_exp / 'Behaviour'
        
        # Define paths for raw, behavioral, and fiberphotometry data
        rawdata_path = data_path_exp / f'{mouse}.doric'
        deinterleaved_raw_path = pp_path / f'{mouse}_deinterleaved.csv'
        dlc_path = behav_path_exp / f'{mouse}_reducedDLC_resnet50_FiberMEC_EPMMay14shuffle1_100000_filtered.csv'
        fiberpho_path = pp_path / f'{mouse}_dFFfilt.csv'
        
        # Path to arena boundary
        try:
            arena_coordinates_file = [coordinates for coordinates in os.listdir(behav_path_exp) if coordinates.endswith('.json')][0]
            print(f'Arena boundaries : {arena_coordinates_file}')
            arena_coordinates_path = behav_path_exp / arena_coordinates_file
        except:
            print(f'Arena boundaries not found for mouse {mouse}. Create json file via get_epm_coordinates module.')
        
        # Paths for output files to be checked
        behav_path = behav_path_exp / f'behav_{mouse}.csv' 

        fiberbehav_path = repo_path / f'{batch}_{mouse}_fiberbehav.csv'
        fiberbehav_notderived_path = repo_path / f'{batch}_{mouse}_fiberbehavnotderived.csv'

        fiberbehav_plot_pdf_path = repo_path / f'{batch}_{mouse}_fiberbehav.pdf'
        fiberbehav_plot_png_path = repo_path / f'{batch}_{mouse}_fiberbehav.png'


        print(f'Behavior and coordinates alignment for {mouse}...')

        # Get and filter DLC data ; analyse coordinates depending on test
        if arena_analysis:
            print('Open json file')
            with open(arena_coordinates_path, 'r') as file:
                arena_coordinates = json.load(file)

        if dlc_data:
            print('Get dlc data')
            coordinates_df = mp.get_dlc_data(dlc_path, threshold=0.99)

        if not automated_alignment:
            print('Get camera flashes')
            camera_df = cp.get_camera_flashes(rawdata_path)

            print('Aligning time with coordinates')
            coordinates_df = cp.align_camera_flashes(coordinates_df, camera_df)

        print(f'Analyze mouse position for {exp}')
        if 'EPM' in exp:
            list_BOI = ['Open arm', 'Closed arm', 'Center']
            behav_df = epm.analyze_mouse_position(coordinates_df, arena_coordinates, bodypart='nose')
            behav_df.to_csv(behav_path)
            print(f'Behaviour file exported to {behav_path}')
        
        elif boris:
            behav_df = pd.read_csv(behav_path)
            deinterleaved_df = pd.read_csv(deinterleaved_raw_path)
            list_BOI = [col for col in behav_df.columns if col not in ['time']]
            behav_df = bp.correct_time_behav(deinterleaved_df, behav_df)

        # Load fiberphotometry data
        fiberpho = pd.read_csv(fiberpho_path)
        fiberpho = pp.butterfilt(fiberpho, ORDER, CUT_FREQ)
            
        # Align behavior and fiberphotometry data
        print('Aligning fiberphotometry and behaviour data')
        fiberbehav_df = bp.align_behav(behav_df, fiberpho, list_BOI)
        fiberbehav_df = bp.behav_process(fiberbehav_df, list_BOI, THRESH_S, EVENT_TIME_THRESHOLD)
        
        # 4 Save the intermediate "not derived" data
        fiberbehav_df.to_csv(fiberbehav_notderived_path, index=False)
        
        # 5 Derive and save final processed data
        dfiberbehav_df = bp.derive(fiberbehav_df, list_BOI)
        dfiberbehav_df.to_csv(fiberbehav_path, index=False)
        
        # 6 Plot fiberphotometry data with behavior
        fig_fiberbehav = bp.plot_fiberpho_behav(
            dfiberbehav_df, 
            list_BOI, 
            exp, 
            mouse, 
            THRESH_S, 
            EVENT_TIME_THRESHOLD,
            batch
        )
        
        # Save plots in both PDF and PNG formats
        fig_fiberbehav.savefig(fiberbehav_plot_pdf_path)
        fig_fiberbehav.savefig(fiberbehav_plot_png_path)
        plt.close(fig_fiberbehav)

    except Exception as e:
        print(f'Error while processing mouse {mouse} : {e}')

print(f'Analysis for {exp} complete. Data saved in {repo_path}')
                        
 #%% 2.2 - Calculate mean, max, and delta dFF within behavioural states (for state behaviours)

for session_path in [Path(f.path) for f in os.scandir(exp_path) if f.is_dir()]:
    session = session_path.name  # Extract session name
    print('###################')
    print(f'EXPERIMENT: {exp}')
    print('###################')
    
    # Create necessary paths
    repo_path = session_path / f'length{EVENT_TIME_THRESHOLD}_interbout{THRESH_S}_o{ORDER}f{CUT_FREQ}'
    groupanalysis_path = repo_path / 'Group_analysis'
    groupanalysis_path.mkdir(exist_ok=True)  # Create directory if it doesn't exist
    
    # Create lists to store data for export 
    mean_dFFs_list = []
    diffmeanmaxdFF_list = []
    diffmeanmaxdFF_perbout_list = []
    
    for mouse, batch, group in zip(subjects_df['Subject'], subjects_df['Batch'], subjects_df['Group']):
        
        # Paths to the relevant fiberbehav files
        fiberbehav_path = repo_path / f'{mouse}_fiberbehav.csv'
        fiberbehav_notderived_path = repo_path / f'{mouse}_fiberbehavnotderived.csv'
        
        # Check if required data exists for this mouse
        if fiberbehav_path.exists():
            print("-------------------")
            print(f'PROCESSING MOUSE: {mouse} batch {batch}')
            print("-------------------")
            
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
        
        # Initialize result Dataframe
        results_df = pd.DataFrame(columns=['Subject',  
                                           'Group', 
                                           f'Mean dFF before {BOI}', 
                                           f'Mean dFF after {BOI}', 
                                           f'Max dFF before {BOI}', 
                                           f'Max dFF after {BOI}'])
        
        for mouse, batch, group in zip(subjects_df['Subject'], subjects_df['Batch'], subjects_df['Group']):
            file_prefix = f'{batch}_{mouse}_{code}'
            print("------------------")
            print(f'MOUSE: {mouse} {batch}')
            print("------------------")
            
            fiberbehav_file = repo_path / f'{file_prefix}_fiberbehav.csv'
            if not fiberbehav_file.exists():
                print(f"File not found: {fiberbehav_file}")
                continue
            
            try:
                # Read in the CSV file, get mouse group
                fiberbehav_df = pd.read_csv(fiberbehav_file)
                group = subjects_df.loc[subjects_df['Subject'] == mouse, 'Group'].values[0]
                # Calculate results for mouse
                results_mouse_df = sc.process_mouse_meanmax_beforeafter(fiberbehav_df, mouse, group, BOI, MAX_BOUTS_NUMBER)
                # Put results in general Dataframe
                results_df = pd.concat([results_df, results_mouse_df], ignore_index=True)
                
            except Exception as e:
                print(f"Error processing file {fiberbehav_file}: {e}")
                continue
        
        # Export data to Excel
        output_file = groupanalysis_path / f'{exp}_{session}_{TIME_MEANMAX}sbeforeafter{BOI}_maxmeandFF.xlsx'
        
        try:
            results_df.to_excel(output_file, index=False)
            print(f"Data successfully exported to {output_file}")
        except Exception as e:
            print(f"Error exporting data to {output_file}: {e}")

#%% 2.4 - Plot PETH for each mouse

# PETH parameters 
EVENT_LIST = ['onset']  # Event triggers, e.g., onset, withdrawal
TIME_WINDOWS = [[5, 30]]  # Time window for PETH calculation (pre, post)

# Loop over each session folder in the experiment path
print('##########################################')
print(f'EXPERIMENT : {exp}')
print('##########################################')

# Create the repository and PETH paths
repo_path = exp_path / f'length{EVENT_TIME_THRESHOLD}_interbout{THRESH_S}_o{ORDER}f{CUT_FREQ}'
peth_path = repo_path / 'PETH'
peth_path.mkdir(parents=True, exist_ok=True)  # Create directory if it doesn't exist

# Loop over each mouse in the subjects DataFrame
for mouse, batch, group in zip(subjects_df['Subject'], subjects_df['Batch'], subjects_df['Group']):
    fiberbehav_path = repo_path / f'{batch}_{mouse}_{code}_fiberbehav.csv'

    if fiberbehav_path.exists():  # Check if fiber behavior file exists for this mouse
        print("--------------")
        print(f'MOUSE : {mouse} {batch}')
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
MAXBOUTSNUMBER = None
# ----------------------------- #

print('##########################################')
print(f'EXPERIMENT: {exp}')
print('##########################################')

repo_path = exp_path / f'length{EVENT_TIME_THRESHOLD}_interbout{THRESH_S}_o{ORDER}f{CUT_FREQ}'
peth_path = repo_path / 'PETH'
peth_path.mkdir(parents=True, exist_ok=True)  # Create PETH directory if it doesn't exist

# Initialize data storage lists
subject_list = []
group_list = []
PETH_array = None
PETH_mean_list = []
PETH_max_list = []

# Loop over each subject (mouse)
for mouse, batch, group in zip(subjects_df['Subject'], subjects_df['Batch'], subjects_df['Group']):
    print("--------------")
    print(f'MOUSE: {mouse} {batch}')
    print("--------------")
    
    fiberbehav_file = repo_path / f'{batch}_{mouse}_{code}_fiberbehav.csv'
    
    if not fiberbehav_file.exists():
        print(f"File not found: {fiberbehav_file}")
        continue
    
    try:
        fiberbehav_df = pd.read_csv(fiberbehav_file, index_col=0)
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
exp = 'EPM_3'
if 'EPM' in exp:
    list_BOI = ['Open arm', 'Closed arm', 'Center']
lowcut = 0.1  # Lowcut frequency for bandpass filter (Hz)
highcut = 3.0    # Highcut frequency for bandpass filter (Hz)
# ----------------------------- #

exp_path = analysis_path / exp
repo_path = exp_path / f'length{EVENT_TIME_THRESHOLD}_interbout{THRESH_S}_o{ORDER}f{CUT_FREQ}'
groupanalysis_path = repo_path / 'Group_analysis'
groupanalysis_path.mkdir(parents=True, exist_ok=True)

var_transients_list = []

# Loop over each subject (mouse)
for mouse, batch in zip(subjects_df['Subject'], subjects_df['Batch']):
    print("--------------")
    print(f'MOUSE: {mouse}, BATCH: {batch}')
    print("--------------")
    
    # Set file path for the fiber behavior CSV
    fiber_file = repo_path / f'{batch}_{mouse}_fiberbehavnotderived.csv'
    
    # Check if the file exists
    if not fiber_file.exists():
        print(f"File not found: {fiber_file}")
        continue

    # Load data
    dfiber_df = pd.read_csv(fiber_file)
    
    # Apply bandpass filter to the dFF signal
    dfiber_df['Filtered dFF'] = tr.bandpass_filter(dfiber_df, lowcut, highcut)
    
    # Plot signal and spectrum
    tr.plot_signal_and_spectrum(dfiber_df)
    
    # Retrieve group information
    group = subjects_df.loc[subjects_df['Subject'] == mouse, 'Group'].values[0]
    
    #Calculate variance and transients characteristics
    mouse_df = sc.variance_transients(dfiber_df, list_BOI, mouse, group, exp, batch)
    var_transients_list.append(mouse_df)

    
# Concatenate results and export to Excel
variability_df = pd.concat(var_transients_list, ignore_index=True)
output_file = groupanalysis_path / f'Variability_1o{ORDER}f{lowcut}_{highcut}.xlsx'
variability_df.to_excel(output_file, index=False)

print(f"Variability data saved to {output_file}")
 # %%
