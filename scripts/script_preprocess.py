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
import importlib

#path to other scripts in sys.path
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

#import functions
import modules.preprocess as pp
importlib.reload(pp)
import modules.genplot as gp
importlib.reload(gp)
import modules.nomenclature as nom
importlib.reload(nom)
import modules.clean_signal as cs
importlib.reload(cs)

from scripts.loader import experiment_path, analysis_path, data_path, exp, ORDER, CUT_FREQ, proto_df, subjects_df, artifact_file, TIME_BEGIN, batches

#%% 
# 1 - PREPROCESSING
#####################

# Step 1: Create main experiment folder and session subfolders
exp_path = nom.setup_experiment_directory(analysis_path, exp)
print(f"Experiment directory created at: {exp_path}")

# Step 2: Get the path to the raw data folder for the experiment
datapath_exp_dict = nom.get_experiment_data_path(batches, proto_df, data_path, exp)
print(f"Raw data paths for experiment '{exp}': {datapath_exp_dict}")

#%% 
# 1.1 - Deinterleave and clean data and save in separate file
# Loop through each session directory
print('######################')
print(f'EXPERIMENT : {exp}')
print('######################')

# Loop through each mouse in the subject DataFrame
for mouse, batch in zip(subjects_df['Subject'], subjects_df['Batch']):
    print("-----------------------------") 
    print(f'BATCH : {batch}, MOUSE : {mouse}')
    print("-----------------------------")
    
    data_path_exp = datapath_exp_dict[batch]
    
    # Create preprocessing directory inside the raw data path
    pp_path = nom.setup_preprocessing_directory(data_path_exp)

    # Find raw data
    raw_data_path = data_path_exp / f'{mouse}.doric'
    
    # Paths for output deinterleaved and plot files
    deinterleaved_path = pp_path / f'{mouse}_deinterleaved.csv'
    cleaned_path = pp_path / f'{mouse}_deinterleaved_cleaned.csv'
    raw_plot_path = pp_path / f'{mouse}_rawdata.png'
    cleaned_plot_path = pp_path / f'{mouse}_cleaned.png'
    
    # Check if raw data exists and deinterleaved data does not exist
    if raw_data_path.exists() and not deinterleaved_path.exists():
        
        #1 Load deinterleaved raw data and clean data
        deinterleaved_df = pp.load_deinterleaved_doric(raw_data_path)
        cleaned_df = cs.clean_signal(deinterleaved_df, detrending=False, apply_hampel=True, apply_filter=False)
        
        #2 Save to CSV
        deinterleaved_df.to_csv(deinterleaved_path, index=False)
        cleaned_df.to_csv(cleaned_path, index=False)
        
        #3 Plot raw data and cleaned data and save as PNG
        fig_raw = gp.plot_rawdata(deinterleaved_df, exp, mouse)
        fig_cleaned = gp.plot_rawdata(cleaned_df, exp, mouse)
        fig_raw.savefig(raw_plot_path)
        fig_cleaned.savefig(cleaned_plot_path)
        plt.close(fig_raw)
        plt.close(fig_cleaned)

#%% 
# 1.3 - Open artifacted data and score artifacts (when big artifacts due to patch cord disconnection)

#------------------#
session = 'Conditioning'
mouse = 'mCD2l'
batch = 'B10'
filecode = f'{exp}_{mouse}'
#------------------#

# in excel 'Filecode', put '{exp}_{mouse}'
code = gp.session_code(session,exp)
pp_path = datapath_exp_dict[batch] / 'Preprocessing'
deinterleaved_df = pd.read_csv(pp_path/f'{mouse}_deinterleaved.csv')

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
    
    html.Div(id='artifact-message', style={'color': 'black', 'fontWeight': 'bold'}),
    
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


#%% 
# 1.4 - Artifact correction and dFF calculation

#import artifacts boundaries
artifacts_df = pd.read_excel(experiment_path / 'artifacts.xlsx')
method = 'fit'

print('#####################')
print(f'EXPERIMENT : {exp}')
print('#####################')
for mouse, batch in zip(subjects_df['Subject'], subjects_df['Batch']):
    pp_path = datapath_exp_dict[batch] / 'Preprocessing'
    if os.path.exists(pp_path/f'{mouse}_deinterleaved.csv'):
        try:
            print("-----------------------------") 
            print(f'BATCH : {batch}, MOUSE : {mouse}')
            print("-----------------------------")
            cleaned_df = pd.read_csv(pp_path/f'{mouse}_deinterleaved_cleaned.csv')
            filecode = f'{exp}_{mouse}'
            
            # calculate dFF with artifacts removal, then interpolate missing data
            dFFdata_df = pp.dFF(cleaned_df,artifacts_df,filecode,method)
            interpdFFdata_df = pp.interpolate_dFFdata(dFFdata_df, method='linear')
            #sometimes 1st timestamps=Nan instead of 0, raises an error
            interpdFFdata_df['Time(s)'] = interpdFFdata_df['Time(s)'].fillna(0) 
            
            # smooth data with butterworth filter or simple moving average (SMA)
            #smoothdFF_df=pp.smoothing_SMA(interpdFFdata_df,win_size=7)
            smoothdFF_df=pp.butterfilt(interpdFFdata_df, ORDER, CUT_FREQ)
            smoothdFF_df.to_csv(pp_path/f'{mouse}_dFFfilt.csv')
            
            #plotted GCaMP and isosbestic curves after dFF or fitting
            fig_dFF = gp.plot_fiberpho(smoothdFF_df,exp,mouse,method)
            fig_dFF.savefig(pp_path/f'{mouse}_{method}dFF.png')
            plt.close(fig_dFF) 
        except Exception as e:
                print(f'Problem in processing mouse {mouse} : {e}')
