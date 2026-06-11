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
import modules.common.preprocess as pp
importlib.reload(pp)
import modules.common.genplot as gp
importlib.reload(gp)
import modules.common.nomenclature as nom
importlib.reload(nom)
import modules.common.clean_signal as cs
importlib.reload(cs)
import modules.common.median_filtering as mf
importlib.reload(mf)

from scripts.loader import experiment_path, analysis_path, data_path, proto_df, subjects_df, artifact_file, TIME_BEGIN, batches, exclusion_file

#%% 
# 1 - PREPROCESSING
#####################

exp = 'Reward_Hab1'
dual_color = False
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
    raw_data_path = data_path_exp / f'{mouse}_0000.doric'
    
    # Paths for output deinterleaved and plot files
    deinterleaved_path = pp_path / f'{mouse}_deinterleaved.csv'
    cleaned_path = pp_path / f'{mouse}_deinterleaved_cleaned.csv'
    raw_plot_path = pp_path / f'{mouse}_rawdata.png'
    cleaned_plot_path = pp_path / f'{mouse}_cleaned.png'
    
    # Check if raw data exists and deinterleaved data does not exist
    if raw_data_path.exists() and not deinterleaved_path.exists():
        if dual_color:
            # Load deinterleaved raw data and clean data.
            deinterleaved_df = pp.load_lockin_dualcolor_doric(raw_data_path)
            downsampled_df = pp.downsample(deinterleaved_df, target_frequency=40)
            cleaned_df = cs.remove_high_artifacts_dualcolor(downsampled_df)

            # Save to CSV
            deinterleaved_df.to_csv(deinterleaved_path, index=False)
            cleaned_df.to_csv(cleaned_path, index=False)

            # Plot raw data and cleaned data and save as PNG
            fig_raw = gp.plot_rawdata(deinterleaved_df, exp, mouse)
            fig_cleaned = gp.plot_rawdata(cleaned_df, exp, mouse)
            fig_raw.savefig(raw_plot_path)
            fig_cleaned.savefig(cleaned_plot_path)
            plt.close(fig_raw)
            plt.close(fig_cleaned)
        
        else:
            #1 Load deinterleaved raw data and clean data
            deinterleaved_df = pp.load_deinterleaved_doric(raw_data_path)
            cleaned_df = cs.remove_high_artifacts(deinterleaved_df)
            
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
mouse = '904'
batch = 1
filecode = f'{exp}_{mouse}'
#------------------# 

# in excel 'Filecode', put '{exp}_{mouse}'
pp_path = datapath_exp_dict[batch] / 'Preprocessing'
deinterleaved_df = pd.read_csv(pp_path/f'{mouse}_deinterleaved.csv')
downsampled_df = pp.downsample(deinterleaved_df, target_frequency=40)

# Create the Dash app
app = Dash(__name__)

# Create the figure
fig = px.line(downsampled_df[TIME_BEGIN:], x='Time(s)', y='405 Deinterleaved')

# App layout
app.layout = html.Div([
    html.H4(f'{exp} {mouse}'),
    
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
    app.run(debug=False, use_reloader=False)

# Dash is running on http://127.0.0.1:8050/
# You can change port if 8050 already taken (8051, etc)


#%% 
# 1.4 - Artifact correction, high-pass filtering and dFF calculation

#import artifacts boundaries
artifacts_df = pd.read_excel(experiment_path / 'artifacts.xlsx')
method = 'fit'
correct_photobleach_method = 'highpass'

print('#####################')
print(f'EXPERIMENT : {exp}')
print('#####################')
for mouse, batch in zip(subjects_df['Subject'], subjects_df['Batch']):
    pp_path = datapath_exp_dict[batch] / 'Preprocessing'
    if os.path.exists(pp_path/f'{mouse}_deinterleaved.csv'):
        print("-----------------------------") 
        print(f'BATCH : {batch}, MOUSE : {mouse}')
        print("-----------------------------")
        cleaned_df = pd.read_csv(pp_path/f'{mouse}_deinterleaved_cleaned.csv')
        filecode = f'{exp}_{mouse}'
        
        # calculate dFF with artifacts removal, then interpolate missing data
        if dual_color:
            dFFdata_df = pp.dFF_dualcolor(cleaned_df, artifacts_df, filecode, method_560='lowess', apply_median_filter=True)
        else:
            dFFdata_df = pp.dFF(cleaned_df,
                                artifacts_df,
                                filecode,
                                method,
                                apply_median_filter=True,
                                fit_model_name='huber')
            #_, _, median_fig = mf.iterative_median_filter(cleaned_df, '465 Deinterleaved')
            #median_fig.savefig(pp_path/f'{mouse}_median_filtering.png')

        # interpolate missing data
        interpdFFdata_df = pp.interpolate_dFFdata(dFFdata_df, method='linear')
        # sometimes 1st timestamps=Nan instead of 0, raises an error
        interpdFFdata_df['Time(s)'] = interpdFFdata_df['Time(s)'].fillna(0)

        if correct_photobleach_method == 'highpass':
            #high-pass filter to remove slow oscillations
            filtered_dFFdata_df = cs.highpass_filter_dff(interpdFFdata_df, dual_color, cutoff_freq = 0.006)

        elif correct_photobleach_method == 'exponential':
            #exponential detrend to remove slow oscillations
            filtered_dFFdata_df = cs.exponential_detrend(interpdFFdata_df, dual_color)

        # Save output to csv
        filtered_dFFdata_df.to_csv(pp_path/f'{mouse}_dFF_corrected.csv', index=False)

        #plotted GCaMP and isosbestic curves after dFF and photobleach correction
        if dual_color:
            fig_dFF = gp.plot_fiberpho_dualcolor(filtered_dFFdata_df,exp,mouse,method)
            fig_dFF.savefig(pp_path/f'{mouse}_{method}dFF_corrected.png')
            plt.close(fig_dFF) 

        else:
            fig_dFF = gp.plot_fiberpho(filtered_dFFdata_df,exp,mouse,method)
            fig_dFF.savefig(pp_path/f'{mouse}_{method}dFF_corrected.png')
            plt.close(fig_dFF) 

# %%
# 1.5 - Manually remove corrupted data if some are left

#------------------#
mouse = '913'
batch = 1
filecode = f'{exp}_{mouse}'
#------------------# 

# in excel 'Filecode', put '{exp}_{mouse}'
pp_path = datapath_exp_dict[batch] / 'Preprocessing'
filtered_dFFdata_df = pd.read_csv(pp_path/f'{mouse}_dFF_corrected.csv')
downsampled_df = pp.downsample(filtered_dFFdata_df, target_frequency=40)

# Create the Dash app
app = Dash(__name__)

# Create the figure
fig = px.line(downsampled_df[TIME_BEGIN:], x='Time(s)', y='dFF')

# App layout
app.layout = html.Div([
    html.H4(f'{exp} {mouse}'),
    
    dcc.Graph(
        id='plot',
        figure=fig,
        config={'displayModeBar': True}  # Add buttons for zooming, panning, etc.
    ),
    
    html.Div(id='artifact-message', style={'color': 'black', 'fontWeight': 'bold'}),
    
    html.Button("Save Excluded Regions", id="save-button", n_clicks=0),
    
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
def capture_exclusions(click_data, exclusion_intervals, click_state):
    """
    Handles clicks on the plot. 
    On the first click, the start of the artifact is captured. 
    On the second click, the end of the artifact is captured, and the interval is saved.
    """
    if click_data:
        time_clicked = click_data['points'][0]['x']
        
        if click_state is None:  # First click (start of the artifact)
            click_state = time_clicked
            message = f'Excluded Region start marked at {time_clicked:.2f} seconds. Now click the end point.'
            print(message)
        else:  # Second click (end of the artifact)
            start = min(click_state, time_clicked)
            end = max(click_state, time_clicked)
            exclusion_intervals.append((start, end))
            message = f'Excluded Region interval ({start:.2f}s, {end:.2f}s) saved. Click to start a new interval.'
            print(message)
            click_state = None  # Reset click state for next pair of clicks

    else:
        message = 'Click on the graph to mark the start of an excluded region.'

    return exclusion_intervals, message, click_state

# Callback to save artifact intervals to an Excel file
@app.callback(
    Output('save-button', 'children'),
    Input('save-button', 'n_clicks'),
    State('artifact-storage', 'data')
)
def save_exclusions_to_excel(n_clicks, exclusion_intervals):
    """
    Saves the artifact intervals to an Excel file when the save button is pressed.
    Each row in the Excel file contains the start and end times of each artifact.
    """
    if n_clicks > 0:
        if len(exclusion_intervals) > 0:
            print(f"\n--- Processing filecode: {filecode} ---")
            print(f"Excluded Regions to store: {exclusion_intervals}")
            pp.update_artifacts_file(exclusion_file, filecode, exclusion_intervals)
            print(f"Saved {len(exclusion_intervals)} artifact intervals to {exclusion_file}")
            return f'Saved {len(exclusion_intervals)} Excluded Regions'
        else:
            print("No excluded regions to save.")
            return "No excluded regions to save"

    return "Save Excluded Regions"

# Run the server
if __name__ == '__main__':
    app.run(debug=False, use_reloader=False)

 # %%
# 1.6 - Replace corrupted data with Nan, then interpolate
exclusion_df = pd.read_excel(experiment_path / 'exclusions.xlsx')

for mouse, batch in zip(subjects_df['Subject'], subjects_df['Batch']):
    pp_path = datapath_exp_dict[batch] / 'Preprocessing'
    if os.path.exists(pp_path/f'{mouse}_deinterleaved.csv'):
        print("-----------------------------") 
        print(f'BATCH : {batch}, MOUSE : {mouse}')
        print("-----------------------------")
        filecode = f'{exp}_{mouse}'
        filtered_dFFdata_df = pd.read_csv(pp_path/f'{mouse}_dFF_corrected.csv')

        final_dFF = filtered_dFFdata_df.copy()
        final_dFF = pp.apply_excluded_regions(final_dFF,'dFF',exclusion_df,filecode)

        if dual_color:
            final_dFF = pp.apply_excluded_regions(final_dFF,'560 dFF',exclusion_df,filecode)

        interpdFFdata_df = pp.interpolate_dFFdata(final_dFF, method='linear')

        interpdFFdata_df.to_csv(pp_path/f'{mouse}_dFF_corrected_final.csv', index=False)

        #plotted curves after corrupted zones removal
        if dual_color:
            fig_dFF = gp.plot_fiberpho_dualcolor(interpdFFdata_df,exp,mouse,method)
            fig_dFF.savefig(pp_path/f'{mouse}_{method}dFF_corrected_final.png')
            plt.close(fig_dFF) 

        else:
            fig_dFF = gp.plot_fiberpho(interpdFFdata_df,exp,mouse,method)
            fig_dFF.savefig(pp_path/f'{mouse}_{method}dFF_corrected_final.png')
            plt.close(fig_dFF) 

# %%
