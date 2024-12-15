# -*- coding: utf-8 -*-
"""
Created on Sun Dec 15 19:41:58 2024

To run fiberphotometry analysis with behaviour or plethysmography data
3 - ANALYSIS WITH PLETHYSMOGRAPHY DATA

@author: alice
"""

#%%
##########
#IMPORTED# 
##########

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'figure.max_open_warning' : 0})
from pathlib import Path
import os
from dash import Dash, dcc, html
import plotly.express as px

#import functions
import preprocess as pp
import genplot as gp
import behavplot as bp
import plethyplot as plp
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

exp_path = analysis_path / exp
data_path_exp = data_path / proto_df.loc[proto_df['Task']==exp, 'Data_path'].values[0]
pp_path = data_path_exp / 'Preprocessing'


############
#PARAMETERS#
############

#threshold to fuse behaviour if bouts are too close, in secs
THRESH_S = 1
#threshold for PETH : if events are too short do not plot them and do not include them in PETH, in seconds
EVENT_TIME_THRESHOLD = 1

#%% 3 - ANALYSIS - PLETHYSMOGRAPH
#################################

# 3.1 - Visualize and score manually sniffs and stims in excel file

#------------------#
session = 'Test'
mouse = 'B2f'
#------------------#

code = gp.session_code(session,exp)
rawdata_path = data_path_exp / f'{mouse}_{code}.csv'
plethys_df = pd.read_csv(rawdata_path, skiprows=1, usecols=['Time(s)','AIn-4'])

app = Dash(__name__)
data_df = plethys_df.loc[[i for i in range(0,len(plethys_df),600)]] #downsample plethysmo data
fig = px.line(data_df, x='Time(s)', y='AIn-4')
app.layout = html.Div([
    html.H4(f'{exp} {session} {mouse}'),
    dcc.Graph(
        id=f'{exp}',
        figure=fig
    )
])
if __name__ == '__main__':
    app.run_server(debug=False, use_reloader=False)

#%% 3.2 - Align with sniffs, create corresponding csv, plot fiberpho data with sniffs and stims

#import sniffs file 
sniffs_df = pd.read_excel(data_path_exp / 'Sniffs.xlsx')
session_path = exp_path / 'Test'
session = str(session_path).split('\\')[-1]

# Print session details
print('##########################################')
print(f'EXPERIMENT : {exp} - SESSION : {session}')
print('##########################################')

# Generate session code
code = gp.session_code(session, exp)

# Define and create necessary directories
repo_path = session_path / f'length{EVENT_TIME_THRESHOLD}_interbout{THRESH_S}_o{ORDER}f{CUT_FREQ}'
raw_path = repo_path / 'Raw'

# Create directories if they do not exist
os.makedirs(repo_path, exist_ok=True)
os.makedirs(raw_path, exist_ok=True)

# Loop over each mouse in the subject list
for mouse in subjects_df['Subject']:
    print("--------------")
    print(f'MOUSE : {mouse}')
    print("--------------")
    
    if mouse in set(sniffs_df['Subject']):
        # Load fiber photometry data
        fiberpho_file = pp_path / f'{mouse}_{code}_dFFfilt.csv'
        fiberpho_df = pd.read_csv(fiberpho_file, index_col=0)
        
        # Check if fibersniff file already exists, skip processing if it does
        fibersniff_file = repo_path / f'{mouse}_{code}_fibersniff.csv'
        if not fibersniff_file.is_file():
            # Load plethysmograph data
            rawdata_path = data_path_exp / f'{mouse}_{code}.csv'
            sr = pp.samplerate(fiberpho_df)
            
            try:
                plethys_df = pd.read_csv(rawdata_path, skiprows=1, usecols=['Time(s)', 'AIn-4'])
            except:
                print('Columns are not defined as usual, trying option without skiprows')
                plethys_df = pd.read_csv(rawdata_path, skiprows=1, usecols=['Time(s)', 'AIn-4'])
                continue

            # Align sniffs and process fibersniff data
            fibersniff_df = plp.align_sniffs(fiberpho_df, plethys_df, sniffs_df, sr, mouse)
            fibersniff_df = plp.process_fibersniff(fibersniff_df, EVENT_TIME_THRESHOLD, THRESH_S, sr)
            dfibersniff_df = plp.derive(fibersniff_df)
            
            # Save fibersniff files
            fibersniff_not_derived_file = repo_path / f'{mouse}_{code}_fibersniffnotderived.csv'
            fibersniff_derived_file = repo_path / f'{mouse}_{code}_fibersniff.csv'
            
            fibersniff_df.to_csv(fibersniff_not_derived_file)
            dfibersniff_df.to_csv(fibersniff_derived_file)
        
        # Plot and save raw plethysmograph and fiber photometry signal
        raw_plot_file = raw_path / f'{mouse}_WBPfiberpho_raw.png'
        if not raw_plot_file.is_file():
            fig_raw = plp.plethyfiber_plot_raw(fiberpho_df, plethys_df, mouse)
            fig_raw.savefig(raw_plot_file)
        
        # Plot and save sniff-aligned signals (if not already saved)
        sniff_plot_pdf = repo_path / f'{mouse}_WBPfiberpho_sniffs.pdf'
        sniff_plot_png = repo_path / f'{mouse}_WBPfiberpho_sniffs.png'
        
        if not sniff_plot_pdf.is_file() or not sniff_plot_png.is_file():
            fig_sniffs = plp.plethyfiber_plot_sniffs(dfibersniff_df, sniffs_df, mouse)
            fig_sniffs.savefig(sniff_plot_png)
            fig_sniffs.savefig(sniff_plot_pdf)
        else:
            print('Sniff-aligned plot already exists for', mouse)
            
    # Close all open plot figures to free up memory
    plt.close('all')
                
#%% 3.3 - Plot PETH for each mouse, sniffs, and stim   

# Import sniffs file 
sniffs_df = pd.read_excel(data_path_exp / 'Sniffs.xlsx')

# Set up session path and details
session_path = exp_path / 'Test'
session = str(session_path).split('\\')[-1]   
EVENT_LIST = ['onset']  # Event triggers, e.g., onset, withdrawal
TIME_WINDOWS = [[5, 30]]  # Time window for PETH calculation (pre, post)        

# Print session details
print('##########################################')
print(f'EXPERIMENT : {exp} - SESSION : {session}')
print('##########################################')

# Generate session code and paths
code = gp.session_code(session, exp)
repo_path = session_path / f'length{EVENT_TIME_THRESHOLD}_interbout{THRESH_S}_o{ORDER}f{CUT_FREQ}'
PETH_path = repo_path / 'PETH'

# Create necessary directories if they don't exist
os.makedirs(repo_path, exist_ok=True)
os.makedirs(PETH_path, exist_ok=True)

# Loop over each subject (mouse) in the experiment
for mouse in subjects_df['Subject']:
    print("--------------")
    print(f'MOUSE : {mouse}')
    print("--------------")
    
    # Only process if mouse is in the sniffs DataFrame
    if mouse in set(sniffs_df['Subject']): 
        # Load the derived fibersniff data for the current mouse
        dfibersniff_file = repo_path / f'{mouse}_{code}_fibersniff.csv'
        dfibersniff_df = pd.read_csv(dfibersniff_file, index_col=0)
        
        # Get the sampling rate
        sr = pp.samplerate(dfibersniff_df)
        
        # Loop over each odor in the sniffs DataFrame
        for odor in set(sniffs_df['Odor']):
            print(f'Processing odor: {odor}')
            
            # Loop through event types and time windows
            for event, timewindow in zip(EVENT_LIST, TIME_WINDOWS):
                print(f'Processing event: {event}')
                
                # Generate paths for the PETH output files
                sniff_png_path = PETH_path / f'{mouse}{odor}_PETHsniff{event[0]}.png'
                sniff_pdf_path = PETH_path / f'{mouse}{odor}_PETHsniff{event[0]}.pdf'
                stim_png_path = PETH_path / f'{mouse}{odor}_PETHstim{event[0]}.png'
                stim_pdf_path = PETH_path / f'{mouse}{odor}_PETHstim{event[0]}.pdf'
                
                # Only process if PETH files do not exist
                if not sniff_pdf_path.is_file() or not stim_pdf_path.is_file():
                    
                    # --- PETH for sniff ---
                    print(f'Generating PETH for sniff: {mouse} - {odor} - {event}')
                    PETHsniff_data = plp.PETH_sniff(dfibersniff_df, odor, event, timewindow, mouse, sr)
                    fig_PETHsniff = plp.plot_PETH(PETHsniff_data, odor, event, timewindow, 'Sniff', sr, mouse)
                    
                    # Save PETH plots for sniff
                    fig_PETHsniff.savefig(sniff_png_path)
                    fig_PETHsniff.savefig(sniff_pdf_path)
                    
                    # --- PETH for stim ---
                    print(f'Generating PETH for stim: {mouse} - {odor} - {event}')
                    PETHstim_data = plp.PETH_stim(dfibersniff_df, odor, event, timewindow, mouse, sr)
                    fig_PETHstim = plp.plot_PETH(PETHstim_data, odor, event, timewindow, 'Stim', sr, mouse)
                    
                    # Save PETH plots for stim
                    fig_PETHstim.savefig(stim_png_path)
                    fig_PETHstim.savefig(stim_pdf_path)
        
    # Close all open plot figures to free up memory
    plt.close('all')

        
#%% 3.4.1 - Calculate mean, max before and during stims for whole group

#------------------------------#
TIME_MEANMAX = 15  # Time window in seconds
#------------------------------#

# Print session and experiment details
print('##########################################')
print(f'EXPERIMENT: {exp} - SESSION: {session}')
print('##########################################')

# Generate paths for analysis
code = gp.session_code(session, exp)
repo_path = session_path / f'length{EVENT_TIME_THRESHOLD}_interbout{THRESH_S}_o{ORDER}f{CUT_FREQ}'
groupanalysis_path = repo_path / 'Group analysis'

# Create directory if it doesn't exist
os.makedirs(groupanalysis_path, exist_ok=True)

# Initialize list to store mean/max dFF for each mouse
meanmaxdFF_list = []

# Set 'Subject' as the index for faster lookup
subjects_df.set_index('Subject', inplace=True)

# Loop through each mouse in the experiment
for mouse in subjects_df.index:
    print("--------------")
    print(f'MOUSE: {mouse}')
    print("--------------")
    
    if mouse in sniffs_df['Subject'].values:
        try:
            # Get the group for this mouse (e.g., control, experimental)
            group = subjects_df.at[mouse, 'Group']
            
            # Path to the dfibersniff CSV file
            fiber_sniff_path = repo_path / f'{mouse}_{code}_fibersniff.csv'
            
            # Check if file exists before attempting to read
            if not fiber_sniff_path.is_file():
                print(f"File not found for mouse '{mouse}': {fiber_sniff_path}")
                continue
            
            # Load dfibersniff DataFrame
            dfibersniff_df = pd.read_csv(fiber_sniff_path, index_col=0)
            
            # Calculate sample rate
            sr = pp.samplerate(dfibersniff_df)
            
            # List of behaviors of interest (BOI) - all 'Stim' columns
            list_BOI = [col for col in dfibersniff_df.columns[3:] if 'Stim' in col]
            
            # Calculate mean and max dFF before and after for each bout
            meanmaxdFF = sc.meanmax_dFF_stims(dfibersniff_df, list_BOI, mouse, group, TIME_MEANMAX, sr)
            meanmaxdFF_list.append(meanmaxdFF)
        
        except Exception as e:
            print(f"Error processing mouse '{mouse}': {e}")
    else:
        print(f"Mouse '{mouse}' not found in sniffs DataFrame. Skipping...")

# Combine all results into a single DataFrame
if meanmaxdFF_list:
    meanmaxdFF_allmice = pd.concat(meanmaxdFF_list, ignore_index=True)
    
    # Export combined results to Excel
    excel_path = groupanalysis_path / f'{exp}_{session}_{TIME_MEANMAX}s_maxmeandFFstims.xlsx'
    meanmaxdFF_allmice.to_excel(excel_path, index=False)
    print(f"Exported mean/max dFF results for all mice to {excel_path}")
else:
    print("No dFF data was processed for this experiment.")
        
#%% 3.4.2 - Calculate mean, max and diff dFF for sniffs

print('##########################################')
print(f'EXPERIMENT: {exp} - SESSION: {session}')
print('##########################################')

# Generate paths for analysis
code = gp.session_code(session, exp)
repo_path = session_path / f'length{EVENT_TIME_THRESHOLD}_interbout{THRESH_S}_o{ORDER}f{CUT_FREQ}'
groupanalysis_path = repo_path / 'Group analysis'

# Create directory if it doesn't exist
os.makedirs(groupanalysis_path, exist_ok=True)

# Initialize lists to store results for all mice
mean_dFFs_list = []       # Stores the mean dFF for each behavior for each mouse
diffmeanmaxdFF_list = []  # Stores the dFF difference between beginning and end of bouts for each mouse

# Set 'Subject' as index for faster lookup
subjects_df.set_index('Subject', inplace=True)

# Loop through each mouse in the experiment
for mouse in subjects_df.index:
    print("--------------")
    print(f'MOUSE: {mouse}')
    print("--------------")
    
    if mouse in sniffs_df['Subject'].values:
        try:
            # Get group for the current mouse
            group = subjects_df.at[mouse, 'Group']
            
            # Paths to the CSV files for this mouse
            fiber_sniff_path = repo_path / f'{mouse}_{code}_fibersniff.csv'
            fiber_sniff_notderived_path = repo_path / f'{mouse}_{code}_fibersniffnotderived.csv'
            
            # Check if the required files exist before processing
            if not fiber_sniff_path.is_file():
                print(f"File not found for mouse '{mouse}': {fiber_sniff_path}")
                continue
            
            if not fiber_sniff_notderived_path.is_file():
                print(f"File not found for mouse '{mouse}': {fiber_sniff_notderived_path}")
                continue
            
            # Load the fiber sniff CSV files
            dfibersniff_df = pd.read_csv(fiber_sniff_path, index_col=0)
            fibersniff_df = pd.read_csv(fiber_sniff_notderived_path, index_col=0)
            
            # Calculate the sampling rate
            sr = pp.samplerate(dfibersniff_df)
            
            # Get the list of behaviors of interest (BOI)
            list_BOI = dfibersniff_df.columns[3:]
            
            # Calculate the difference between the mean and max dFF for each bout
            diff_dFF = sc.diffmeanmax_dFF(dfibersniff_df, list_BOI, mouse, group)
            diffmeanmaxdFF_list.append(diff_dFF)
            
            # Calculate the mean dFF during each behavior
            mean_dFF = sc.meandFF_sniffs(fibersniff_df, exp, session, mouse, group, joined=True)
            mean_dFFs_list.append(mean_dFF)
        
        except Exception as e:
            print(f"Error processing mouse '{mouse}': {e}")
    else:
        print(f"Mouse '{mouse}' not found in sniffs DataFrame. Skipping...")

# Combine and export the results
try:
    if diffmeanmaxdFF_list:
        diffdFFs_allmice = pd.concat(diffmeanmaxdFF_list, ignore_index=True)
        diff_dFFs_path = groupanalysis_path / f'{exp}_{session}_diffmaxmeandFFs.xlsx'
        diffdFFs_allmice.to_excel(diff_dFFs_path, index=False)
        print(f"Exported dFF difference results for all mice to {diff_dFFs_path}")
    else:
        print("No dFF difference data was processed for this experiment.")
except Exception as e:
    print(f"Error exporting diffmaxmeandFFs.xlsx: {e}")

try:
    if mean_dFFs_list:
        mean_dFFs_allmice = pd.concat(mean_dFFs_list, ignore_index=True)
        mean_dFFs_path = groupanalysis_path / f'{exp}_{session}_meandFFs.xlsx'
        mean_dFFs_allmice.to_excel(mean_dFFs_path, index=False)
        print(f"Exported mean dFF results for all mice to {mean_dFFs_path}")
    else:
        print("No mean dFF data was processed for this experiment.")
except Exception as e:
    print(f"Error exporting meandFFs.xlsx: {e}")

    
#%% 3.5 - Plot PETH for all mice

#------------------------------#
# Settings for PETH plotting
included_groups = ['CD', 'HFD']  # Groups to be included in the PETH plot
event = 'onset'  # or 'withdrawal'
timewindow = [-6, 15]  # Time window for the PETH plot
count = 1  # Count for PETH function
#------------------------------#

print('##########################################')
print(f'EXPERIMENT: {exp} - SESSION: {session}')
print('##########################################')

# Generate paths for analysis
code = gp.session_code(session, exp)
repo_path = session_path / f'length{EVENT_TIME_THRESHOLD}_interbout{THRESH_S}_o{ORDER}f{CUT_FREQ}'
groupanalysis_path = repo_path / 'Group analysis'

# Create directory if it doesn't exist
os.makedirs(groupanalysis_path, exist_ok=True)

# Prepare list of all mouse subjects for quick lookup
valid_mouse_subjects = set(sniffs_df['Subject'].values)
subjects_df.set_index('Subject', inplace=True)

# Main loop for processing PETH for each odor
for type_event in ['Stim']:
    for odor in sniffs_df['Odor'].unique():
        print(f'Processing odor: {odor}')
        
        # Initialize lists to hold PETH data
        group_list = []  # List of group labels for each mouse
        PETH_list = []  # List to store PETH data for each mouse
        
        for mouse in subjects_df.index:
            if mouse in valid_mouse_subjects:
                print(f'Processing mouse: {mouse}')
                
                try:
                    # Get group for the current mouse
                    group = subjects_df.at[mouse, 'Group']
                    group_list.append(group)
                    
                    # Paths to the CSV file for this mouse
                    fiber_sniff_path = repo_path / f'{mouse}_{code}_fibersniff.csv'
                    
                    # Check if the file exists
                    if not fiber_sniff_path.is_file():
                        print(f"File not found for mouse '{mouse}': {fiber_sniff_path}")
                        continue
                    
                    # Load the fiber sniff CSV file for the mouse
                    dfibersniff_df = pd.read_csv(fiber_sniff_path, index_col=0)
                    sr = pp.samplerate(dfibersniff_df)
                    
                    # Calculate PETH data for the current mouse
                    if type_event == 'Stim':
                        PETH_data = plp.PETH_stim(dfibersniff_df, odor, event, timewindow, mouse, sr, count=count)
                    elif type_event == 'Sniff':
                        PETH_data = plp.PETH_sniff(dfibersniff_df, odor, event, timewindow, mouse, sr, count=count)
                    
                    # Calculate mean PETH for the current mouse and add to the list
                    PETH_mean = np.mean(PETH_data, axis=0)
                    PETH_list.append(PETH_mean)
                
                except Exception as e:
                    print(f"Error processing mouse '{mouse}': {e}")
            else:
                print(f"Mouse '{mouse}' not found in sniffs DataFrame. Skipping...")
        
        # Combine PETH data and create plot
        try:
            PETH_array = np.stack(PETH_list)
            BOI = f'{type_event} {odor}'  # Behavior of interest
            
            PETHarray_list = []  # List of PETH arrays for each group
            for group in included_groups:
                PETH_array_group = PETH_array.copy()  # Copy the PETH array
                list_todelete = [i for i, group_mouse in enumerate(group_list) if group not in group_mouse]
                PETH_group_array = np.delete(PETH_array_group, list_todelete, axis=0)
                PETHarray_list.append(PETH_group_array)
            
            # Plot the PETH for pooled data of all groups
            fig_PETHpooled = bp.plot_PETH_pooled(included_groups, PETHarray_list, BOI, event, timewindow, exp, session)
            
            # Save the plot as PDF and PNG
            filename_base = f'{included_groups[0]}{included_groups[1]}_{type_event}_{odor}_{event[0]}_{timewindow[1]}_PETH'
            pdf_path = groupanalysis_path / f'{filename_base}.pdf'
            png_path = groupanalysis_path / f'{filename_base}.png'
            
            fig_PETHpooled.savefig(pdf_path)
            fig_PETHpooled.savefig(png_path)
            
            print(f"Saved PETH plot as PDF: {pdf_path}")
            print(f"Saved PETH plot as PNG: {png_path}")
            
            # Close all open plot figures to free up memory
            plt.close('all')
        
        except Exception as e:
            print(f"Error creating PETH plot for odor '{odor}': {e}")

#%% 3.6 - Compute variance and transients

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
# Concatenate results and export to Excel
variability_df = pd.concat(var_transients_list, ignore_index=True)
output_file = groupanalysis_path / f'Variability_1o{ORDER}f{lowcut}_{highcut}.xlsx'
variability_df.to_excel(output_file, index=False)

print(f"Variability data saved to {output_file}")