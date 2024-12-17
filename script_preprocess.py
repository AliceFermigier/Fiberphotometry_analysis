# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 13:39:29 2023

To run fiberphotometry analysis with behaviour or plethysmography data
1 - PREPROCESSING 

@author: alice fermigier
"""

#%%
##########
#IMPORTED# 
##########

import pandas as pd
from pathlib import Path
import os
from dash import Dash, dcc, html, Input, Output, State
import plotly.express as px
import sys  
import matplotlib.pyplot as plt

#path to other scripts in sys.path
path_to_gitrepo=r'C:\Users\alice\Documents\GitHub\Fiberphotometry_analysis'
if path_to_gitrepo not in sys.path:
    sys.path.append(path_to_gitrepo)

#import functions
import preprocess as pp
import genplot as gp
import nomenclature as nom

#%%
########
#LOADER#
########

experiment_path = Path(r'C:\Users\alice\Desktop\Data_Fiber')
analysis_path = experiment_path / 'Analysis' 
data_path = experiment_path / 'Data'
os.chdir(experiment_path)
os.getcwd()

#import ID and groups of all mice
subjects_df = pd.read_excel(experiment_path / 'subjects.xlsx', sheet_name='Included')
#import tasks in protocol
proto_df = pd.read_excel(experiment_path / 'protocol.xlsx')
#create artifacts file if not existent
artifact_file = experiment_path / 'artifacts.xlsx' # File to store artifact timestamps
nom.create_or_load_artifacts_file(artifact_file, option='create_only')

############
#PARAMETERS#
############

#time to crop at the beginning of the trial for , in seconds
TIME_BEGIN = 60
#filter characteristics
ORDER = 4
CUT_FREQ = 1 #in Hz

#%% 1 - PREPROCESSING
#####################

#------------------#
exp = 'OdDis1'
#------------------#

# Step 1: Get the list of sessions for the experiment
session_names = nom.get_experiment_sessions(proto_df, exp)
print(f"Sessions for experiment '{exp}': {session_names}")

# Step 2: Create main experiment folder and session subfolders
exp_path = nom.setup_experiment_directory(analysis_path, exp, session_names)
print(f"Experiment directory created at: {exp_path}")

# Step 3: Get the path to the raw data folder for the experiment
data_path_exp = nom.get_experiment_data_path(proto_df, data_path, exp)
print(f"Raw data path for experiment '{exp}': {data_path_exp}")

# Step 4: Create preprocessing directory inside the raw data path
pp_path = nom.setup_preprocessing_directory(data_path_exp)
print(f"Preprocessing directory created at: {pp_path}")

#%% 1.1 - Deinterleave data and save in separate file

# Loop through each session directory
for session_path in [Path(f.path) for f in os.scandir(exp_path) if f.is_dir()]:
    session = session_path.name  # Extract session name
    print('##########################################')
    print(f'EXPERIMENT : {exp} - SESSION : {session}')
    print('##########################################')
    
    # Generate session code for the current session
    code = gp.session_code(session, exp)
    
    # Loop through each mouse in the subject DataFrame
    for mouse in subjects_df['Subject']:
        print("--------------") 
        print(f'MOUSE : {mouse}')
        print("--------------")
        
        # Create a common file prefix
        file_prefix = f'{mouse}_{code}'
        
        # Paths for input raw data and output deinterleaved and plot files
        raw_data_path = data_path_exp / f'{file_prefix}.csv'
        deinterleaved_path = pp_path / f'{file_prefix}_deinterleaved.csv'
        raw_plot_path = pp_path / f'{file_prefix}_rawdata.png'
        
        # Check if raw data exists and deinterleaved data does not exist
        if raw_data_path.exists() and not deinterleaved_path.exists():
            
            #1 Load raw data (only necessary columns) and deinterleave
            rawdata_df = pd.read_csv(
                raw_data_path, 
                skiprows=1, 
                usecols=['Time(s)', 'AIn-1', 'DI/O-1', 'DI/O-2'], 
                encoding="ISO-8859-1"
            )
            
            #2 Deinterleave the data and save to CSV
            deinterleaved_df = pp.deinterleave(rawdata_df)
            deinterleaved_df.to_csv(deinterleaved_path, index=False)
            
            #3 Plot raw data and save as PNG
            fig_raw = gp.plot_rawdata(deinterleaved_df, exp, session, mouse)
            fig_raw.savefig(raw_plot_path)
            plt.close(fig_raw)

#%% 1.2 - Open artifacted data and score artifacts

#------------------#
session = 'Test'
mouse = 'A2m'
filecode = f'{exp}_{session}_{mouse}'
#------------------#

# in excel 'Filecode', put '{exp}_{session}_{mouse}'
code = gp.session_code(session,exp)
deinterleaved_df = pd.read_csv(pp_path/f'{mouse}_{code}_deinterleaved.csv')

# Create the Dash app
app = Dash(__name__)

# Create the figure
fig = px.line(deinterleaved_df[TIME_BEGIN:], x='Time(s)', y='405 Deinterleaved')

# App layout
app.layout = html.Div([
    html.H4(f'{exp} {session} {mouse}'),
    
    dcc.Graph(
        id='plot',
        figure=fig,
        config={'displayModeBar': True}  # Add buttons for zooming, panning, etc.
    ),
    
    html.Div(id='artifact-message', style={'color': 'green', 'fontWeight': 'bold'}),
    
    html.Button("Save Artifacts", id="save-button", n_clicks=0),
    
    dcc.Store(id='artifact-storage', data=[]),  # Store artifact tuples (start, end)
    dcc.Store(id='click-tracker', data=None)  # Keep track of first/second click
])

# Callback to handle user clicks and record artifact intervals
@app.callback(
    [Output('artifact-storage', 'data'),
     Output('artifact-message', 'children'),
     Output('click-tracker', 'data')],
    Input('plot', 'clickData'),
    [State('artifact-storage', 'data'),
     State('click-tracker', 'data')]
)
def capture_artifact(click_data, artifact_intervals, click_state):
    """
    Handles clicks on the plot. 
    On the first click, the start of the artifact is captured. 
    On the second click, the end of the artifact is captured, and the interval is saved.
    """
    if click_data:
        time_clicked = click_data['points'][0]['x']
        
        if click_state is None:  # First click (start of the artifact)
            click_state = time_clicked
            message = f'Artifact start marked at {time_clicked:.2f} seconds. Now click the end point.'
            print(message)
        else:  # Second click (end of the artifact)
            start = min(click_state, time_clicked)
            end = max(click_state, time_clicked)
            artifact_intervals.append((start, end))
            message = f'Artifact interval ({start:.2f}s, {end:.2f}s) saved. Click to start a new interval.'
            print(message)
            click_state = None  # Reset click state for next pair of clicks

    else:
        message = 'Click on the graph to mark the start of an artifact.'

    return artifact_intervals, message, click_state

# Callback to save artifact intervals to an Excel file
@app.callback(
    Output('save-button', 'children'),
    Input('save-button', 'n_clicks'),
    State('artifact-storage', 'data')
)
def save_artifacts_to_excel(n_clicks, artifact_intervals):
    """
    Saves the artifact intervals to an Excel file when the save button is pressed.
    Each row in the Excel file contains the start and end times of each artifact.
    """
    if n_clicks > 0:
        if len(artifact_intervals) > 0:
            print(f"\n--- Processing filecode: {filecode} ---")
            print(f"Artifacts to store: {artifact_intervals}")
            pp.update_artifacts_file(artifact_file, filecode, artifact_intervals)
            print(f"Saved {len(artifact_intervals)} artifact intervals to {artifact_file}")
            return f'Saved {len(artifact_intervals)} Artifacts'
        else:
            print("No artifacts to save.")
            return "No artifacts to save"

    return "Save Artifacts"

# Run the server
if __name__ == '__main__':
    app.run_server(debug=False, use_reloader=False)

# Dash is running on http://127.0.0.1:8050/
# You can change port if 8050 already taken (8051, etc)

#%% 1.4 - Artifact correction and dFF calculation

#import artifacts boundaries
artifacts_df = pd.read_excel(experiment_path / 'artifacts.xlsx')
method = 'fit'

for session_path in [Path(f.path) for f in os.scandir(exp_path) if f.is_dir()]:
    session = str(session_path).split('\\')[-1]
    print('##########################################')
    print(f'EXPERIMENT : {exp} - SESSION : {session}')
    print('##########################################')
    code = gp.session_code(session,exp)
    for mouse in subjects_df['Subject']:
        if os.path.exists(pp_path/f'{mouse}_{code}_deinterleaved.csv'):
            try:
                print("--------------")
                print(f'MOUSE : {mouse}')
                print("--------------")
                deinterleaved_df = pd.read_csv(pp_path/f'{mouse}_{code}_deinterleaved.csv')
                filecode = f'{exp}_{session}_{mouse}'
                
                # calculate dFF with artifacts removal, then interpolate missing data
                dFFdata_df = pp.dFF(deinterleaved_df,artifacts_df,filecode,method)
                interpdFFdata_df = pp.interpolate_dFFdata(dFFdata_df, method='linear')
                #sometimes 1st timestamps=Nan instead of 0, raises an error
                interpdFFdata_df['Time(s)'] = interpdFFdata_df['Time(s)'].fillna(0) 
                
                # smooth data with butterworth filter or simple moving average (SMA)
                #smoothdFF_df=pp.smoothing_SMA(interpdFFdata_df,win_size=7)
                smoothdFF_df=pp.butterfilt(interpdFFdata_df, ORDER, CUT_FREQ)
                smoothdFF_df.to_csv(pp_path/f'{mouse}_{code}_dFFfilt.csv')
                
                #plotted GCaMP and isosbestic curves after dFF or fitting
                fig_dFF = gp.plot_fiberpho(smoothdFF_df,exp,session,mouse,method)
                fig_dFF.savefig(pp_path/f'{mouse}_{code}_{method}dFF.png')
                plt.close(fig_dFF) 
            except Exception as e:
                    print(f'Problem in processing mouse {mouse} : {e}')
            