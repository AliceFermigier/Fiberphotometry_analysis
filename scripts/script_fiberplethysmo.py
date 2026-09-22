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
import os
from dash import Dash, dcc, html, Input, Output, State
import plotly.express as px
import importlib

#import functions
import modules.common.preprocess as pp
importlib.reload(pp)
import modules.common.genplot as gp
importlib.reload(gp)
import modules.common.behavplot as bp
importlib.reload(bp)
import modules.common.plethyplot as plp
importlib.reload(plp)
import modules.common.statcalc as sc
importlib.reload(sc)
import modules.common.transients as tr
importlib.reload(tr)
import modules.common.nomenclature as nom
importlib.reload(nom)
import modules.common.quantification as quantif
importlib.reload(quantif)
import modules.common.clean_signal as cs
importlib.reload(cs)

from scripts.loader import experiment_path, analysis_path, data_path, proto_df, subjects_df, batches

#%%LOADER
##########

exp = 'Plethysmo'

exp_path = analysis_path / exp
data_path_exp = data_path / proto_df.loc[proto_df['Task']==exp, 'Data_path'].values[0]
pp_path = data_path_exp / 'Preprocessing'

############
#PARAMETERS#
############

dual_color = False

#filter characteristics
ORDER = 4
CUT_FREQ = None #in Hz

#threshold to fuse behaviour if bouts are too close, in secs
THRESH_S = 0
#threshold for PETH : if events are too short do not plot them and do not include them in PETH, in seconds
EVENT_TIME_THRESHOLD = 0

exp = 'Plethysmo'
#list_BOI = ['Stim Clean','Stim HC','Stim Novel','Sniff Clean','Sniff HC','Sniff Novel']

exp_path = analysis_path / exp
datapath_exp_dict = nom.get_experiment_data_path(batches, proto_df, data_path, exp)

#%% 3 - ANALYSIS - PLETHYSMOGRAPH
#################################

# 3.1 - Visualize and score manually sniffs and stims in excel file

#------------------#
mouse = 'B1f'
batch = 5
group = subjects_df.loc[subjects_df['Subject'] == mouse, 'Group'].values[0]
#------------------#

data_path_exp = datapath_exp_dict[batch]
behav_path_exp = data_path_exp / 'Behaviour'
rawdata_path = data_path_exp / f'{mouse}_1.csv'
plethys_df = pd.read_csv(rawdata_path, skiprows=1, usecols=['Time(s)','AIn-4'])

# Dropdown options for odors
odors = ['Clean 1', 'Clean 2', 'Clean 3', 'HC 1', 'HC 2', 'Novel 1', 'Novel 2']

# Paths and file setup
sniffs_stims_file = behav_path_exp / 'Stims_Sniffs.xlsx'
nom.create_or_load_sniffs_file(sniffs_stims_file)

app = Dash(__name__)
data_df = plethys_df.loc[[i for i in range(0,len(plethys_df),600)]] #downsample plethysmo data
fig = px.line(data_df, x='Time(s)', y='AIn-4')

# App layout
app.layout = html.Div([
    html.H4(f'{exp} {mouse}'),
    dcc.Dropdown(id='odor-dropdown', options=[{'label': odor, 'value': odor} for odor in odors],
                 placeholder='Select Odor'),
    dcc.Dropdown(id='stim_sniff-dropdown', options=[{'label': stim_sniffs, 'value': stim_sniffs} for stim_sniffs in ['Stim', 'Sniffs']],
                 placeholder='Select Interval'),
    dcc.Graph(
        id='plot',
        figure=fig,
        config={'displayModeBar': True}  # Add buttons for zooming, panning, etc.
    ),
    html.Div(id='div-message', style={'color': 'black', 'fontWeight': 'bold'}),
    html.Button("Save Stim and Sniffs", id="save-button", n_clicks=0),
    dcc.Store(id='intervals-store', data={'Stim': [], 'Sniffs': []}),
    dcc.Store(id='click-tracker', data=None)  # Keep track of first/second click
])

# Callback to handle user clicks and record artifact intervals
@app.callback(
    [Output('intervals-store', 'data'),
     Output('div-message', 'children'),
     Output('click-tracker', 'data')],
    Input('plot', 'clickData'),
    [State('intervals-store', 'data'),
     State('click-tracker', 'data'),
     State('odor-dropdown', 'value'),
     State('stim_sniff-dropdown', 'value'),]
    )

def capture_interval(click_data, intervals, click_state, odor, stim_sniff):
    """
    Handles clicks on the plot. 
    On the first click, the start of an interval is captured. 
    On the first click, the stop of an interval is captured, and the interval is saved.
    
    Parameters: (takes parameters and returns from @app.callback, in the defined order)
    - click_data : Input('plot', 'clickData')
    - intervals : State('intervals-store', 'data')
    - click_state : State('click-tracker', 'data')
    - odor : State('odor-dropdown', 'value')
    - stim_sniff : State('stim_sniff-dropdown', 'value')
    
    Returns:
    - intervals : Output('intervals-store', 'data')
    - message : Output('div-message', 'children')
    - click_state : Output('click-tracker', 'data')
    """
    if click_data:
        time_clicked = click_data['points'][0]['x']
        
        if click_state is None:  # First click (start of the artifact)
            click_state = time_clicked
            message = f'{stim_sniff} {odor} start marked at {time_clicked:.2f} seconds. Now click the end point.'
            print(message)
        else:  # Second click (end of the artifact)
            start = min(click_state, time_clicked)
            end = max(click_state, time_clicked)
            if stim_sniff == 'Stim':
                intervals['Stim']=[start, end]
                message = f'Stim {odor} ({start:.2f}s, {end:.2f}s) saved.'
            elif stim_sniff == 'Sniffs':
                intervals['Sniffs'].append([start, end])
                message = f'Sniff {odor} ({start:.2f}s, {end:.2f}s) saved. Click on new sniff to score it.'
            print(message)
            click_state = None  # Reset click state for next pair of clicks

    else: 
        message = 'Click on the graph to mark the start of your interval.'

    return intervals, message, click_state

# Callback to save stim and sniffs to an Excel file
@app.callback(
    Output('save-button', 'children'),
    Input('save-button', 'n_clicks'),
    [State('intervals-store', 'data'),
     State('odor-dropdown', 'value'),]
    )

def save_to_excel(n_clicks, intervals, odor):
    """
    Saves the artifact intervals to an Excel file when the save button is pressed.
    Each row in the Excel file contains the start and end times of each artifact.
    """
    if n_clicks > 0:
        if len(intervals['Stim']) > 0:
            print(f"\n--- Processing {mouse} {odor} ---")
            print(f"To store: {intervals}")
            odor_true = odor.split()[0]
            count = odor.split()[1]
            plp.update_sniffs_file(sniffs_stims_file, mouse, group, odor_true, count, intervals['Sniffs'], intervals['Stim'])
            print(f"Saved stim and sniffs intervals to {sniffs_stims_file}")
            return f'Successfully saved {len(intervals['Stim'])/2} stim and {len(intervals['Sniffs'])} sniffs'
        else:
            print("No stim to save.")
            return "No stim to save"

    return "Save Stim and Sniffs"

# Run App
if __name__ == '__main__':
    app.run_server(debug=False, use_reloader=False)
    
# Dash is running on http://127.0.0.1:8050/
# You can change port if 8050 already taken (8051, etc)

#%% 3.2 - Align with sniffs, create corresponding csv, plot fiberpho data with sniffs and stims

# Print session details
print('#####################')
print(f'EXPERIMENT : {exp}')
print('#####################')

# Define and create necessary directories
repo_path = exp_path / f'length{EVENT_TIME_THRESHOLD}_interbout{THRESH_S}_o{ORDER}f{CUT_FREQ}'
repo_path.mkdir(exist_ok=True)

# Loop over each mouse in the subject list
for mouse, batch in zip(subjects_df['Subject'],subjects_df['Batch']):
    print("--------------")
    print(f'MOUSE : {mouse}')
    print("--------------")

    #import sniffs file 
    behav_path_exp = data_path_exp / 'Behaviour'
    sniffs_df = pd.read_excel(behav_path_exp / 'Stims_Sniffs.xlsx')
    
    if mouse in set(sniffs_df['Subject']):
        # Load fiber photometry data and filter is specified
        print('Loading fiberphotometry data')
        fiberpho_file = pp_path / f'{mouse}_dFF_corrected_final.csv'
        fiberpho_df = pd.read_csv(fiberpho_file)
        if CUT_FREQ is not None:
            fiberpho_df = cs.lowpass_dFF(fiberpho_df, dual_color, order = ORDER, cut_freq = CUT_FREQ)

        # Check if fibersniff file already exists, skip processing if it does
        fibersniff_derived_file = repo_path / f'{batch}_{mouse}_fiberbehav.csv'

        if not fibersniff_derived_file.is_file():
            # Load plethysmograph data
            print('Loading plethysmography data')
            rawdata_path = data_path_exp / f'{mouse}_1.csv'
            sr = pp.samplerate(fiberpho_df)
            
            try:
                plethys_df = pd.read_csv(rawdata_path, skiprows=1, usecols=['Time(s)', 'AIn-4'])
            except:
                print('Columns are not defined as usual, trying option without skiprows')
                plethys_df = pd.read_csv(rawdata_path, usecols=['Time(s)', 'AIn-4'])
 
            # Align sniffs and process fibersniff data
            print('Aligning stims and sniffs with fiberphotometry data')
            fibersniff_df = plp.align_sniffs(fiberpho_df, plethys_df, sniffs_df, sr, mouse)
            fibersniff_df = plp.process_fibersniff(fibersniff_df, EVENT_TIME_THRESHOLD, THRESH_S, sr)
            dfibersniff_df = plp.derive(fibersniff_df)
            
            # Save fibersniff files
            fibersniff_not_derived_file = repo_path / f'{batch}_{mouse}_fiberbehavnotderived.csv'
            fibersniff_derived_file = repo_path / f'{batch}_{mouse}_fiberbehav.csv'
            
            fibersniff_df.to_csv(fibersniff_not_derived_file)
            dfibersniff_df.to_csv(fibersniff_derived_file)

            # Save fibersniff files concatenated
            fibersniffconcat_df = plp.concat_fibersniff_columns(fibersniff_df)
            dfibersniffconcat_df = plp.derive(fibersniffconcat_df)
            
            fibersniffconcat_not_derived_file = repo_path / f'{batch}_{mouse}_fiberbehavconcatnotderived.csv'
            fibersniffconcat_derived_file = repo_path / f'{batch}_{mouse}_fiberbehavconcat.csv'
            
            fibersniffconcat_df.to_csv(fibersniffconcat_not_derived_file)
            dfibersniffconcat_df.to_csv(fibersniffconcat_derived_file)
        
        # Plot and save sniff-aligned signals
        print('Plotting')
        sniff_plot_pdf = repo_path / f'{batch}_{mouse}_WBPfiberpho_sniffs.pdf'
        sniff_plot_png = repo_path / f'{batch}_{mouse}_WBPfiberpho_sniffs.png'
        if not sniff_plot_pdf.is_file() or not sniff_plot_png.is_file():
            fig_sniffs = plp.plethyfiber_plot_sniffs(dfibersniff_df, sniffs_df, mouse, batch)
            fig_sniffs.savefig(sniff_plot_png)
            fig_sniffs.savefig(sniff_plot_pdf)
        else:
            print('Sniff-aligned plot already exists for', mouse)
            
    # Close all open plot figures to free up memory
    plt.close('all')

print(f'\n✅ Analysis for {exp} complete.\nData saved in: {repo_path}')

# %% 3.3 - Quantify dFF during stims and sniffs.

subjects_df['Group'] = subjects_df['Group'].fillna('')  # if group = NaN, replaces it with an empty string

# Define necessary directories
repo_path = exp_path / f'length{EVENT_TIME_THRESHOLD}_interbout{THRESH_S}_o{ORDER}f{CUT_FREQ}'

# Load excluded subjects
excluded_subjects_df = pd.read_excel(experiment_path / 'subjects.xlsx',
                                      sheet_name=f'Excluded_{exp}')

zone_cols  = ['Stim Novel 0', 'Stim Novel 1', 'Stim HC 0', 'Stim HC 1',
              'Stim Clean 0', 'Stim Clean 1', 'Stim Clean 2']
behav_cols = ['Sniff Novel 0', 'Sniff Novel 1', 'Sniff HC 0', 'Sniff HC 1',
              'Sniff Clean 0', 'Sniff Clean 1', 'Sniff Clean 2']

merge_into = {
    'Sniff Novel 0': 'Stim Novel 0', 'Sniff Novel 1': 'Stim Novel 1',
    'Sniff HC 0':    'Stim HC 0',    'Sniff HC 1':    'Stim HC 1',
    'Sniff Clean 0': 'Stim Clean 0', 'Sniff Clean 1': 'Stim Clean 1', 'Sniff Clean 2': 'Stim Clean 2',
}

# stim -> sniff (fixed direction), built from zone/behav cols matching by trailing number
stim_sniff_pairs = {
    'Stim Novel 0': 'Sniff Novel 0', 'Stim Novel 1': 'Sniff Novel 1',
    'Stim HC 0':    'Sniff HC 0',    'Stim HC 1':    'Sniff HC 1',
    'Stim Clean 0': 'Sniff Clean 0', 'Stim Clean 1': 'Sniff Clean 1', 'Stim Clean 2': 'Sniff Clean 2',
}

# One dict of lists per (channel, use_zscore) combo, each holding {'main': [...], 'baseline': [...]}
result_sets = {
    ('dFF',     False): {'main': [], 'baseline': []},
    ('dFF',     True):  {'main': [], 'baseline': []},
    ('560 dFF', False): {'main': [], 'baseline': []},
    ('560 dFF', True):  {'main': [], 'baseline': []},
}

excluded_ids = set(excluded_subjects_df['Subject'].astype(str))

for mouse, batch, group in zip(subjects_df['Subject'], subjects_df['Batch'], subjects_df['Group']):
    print(f"--- {mouse} {batch} {group} ---")
    fiberbehav_file = repo_path / f'{batch}_{mouse}_fiberbehavnotderived.csv'

    if str(mouse) in excluded_ids:
        print(f"Mouse {mouse} excluded")
        continue

    if not fiberbehav_file.exists():
        print(f"  File not found, skipping.")
        continue

    fiberbehav_df = pd.read_csv(fiberbehav_file, index_col=0)
    has_560 = '560 dFF' in fiberbehav_df.columns

    for (dff_col, use_zscore), storage in result_sets.items():
        if dff_col == '560 dFF' and not has_560:
            continue

        record, baseline_record = quantif.extract_dff_summary(
            fiberbehav_df    = fiberbehav_df,
            mouse            = mouse,
            batch            = batch,
            group            = group,
            zone_cols        = zone_cols,
            behav_cols       = behav_cols,
            baseline_col     = None,
            merge_into       = merge_into,
            dff_col          = dff_col,
            fps              = 10,
            use_zscore       = use_zscore,
            stim_sniff_pairs = stim_sniff_pairs,
        )
        storage['main'].append(record)
        storage['baseline'].append(baseline_record)

# ── Save to Excel (one file per raw/zscored x channel combo, main + baseline sheets) ──
file_map = {
    ('dFF',     False): 'dFF_summary_raw.xlsx',
    ('dFF',     True):  'dFF_summary_zscored.xlsx',
    ('560 dFF', False): 'dFF_summary_raw_560.xlsx',
    ('560 dFF', True):  'dFF_summary_zscored_560.xlsx',
}

for key, filename in file_map.items():
    storage = result_sets[key]
    if not storage['main']:
        continue  # e.g. no 560 channel in this dataset at all
    with pd.ExcelWriter(repo_path / filename) as writer:
        pd.DataFrame(storage['main']).to_excel(writer, sheet_name='Summary', index=False)
        pd.DataFrame(storage['baseline']).to_excel(writer, sheet_name='Baseline', index=False)

print(f"Saved {len(result_sets[('dFF', False)]['main'])} mice.\nData saved in: {repo_path}")

# %% 3.4 - Quantify dFF during stims and sniffs concatenated.

subjects_df['Group'] = subjects_df['Group'].fillna('') # if group = Nan, replaces it with an empty string

# Define necessary directories
repo_path = exp_path / f'length{EVENT_TIME_THRESHOLD}_interbout{THRESH_S}_o{ORDER}f{CUT_FREQ}'

#Load excluded subjects
excluded_subjects_df = pd.read_excel(experiment_path / 'subjects.xlsx', 
                                     sheet_name=f'Excluded_{exp}')

dFF_records_raw_concat = []
dFF_records_zscored_concat = []
dFF_records_raw_560_concat = []
dFF_records_zscored_560_concat = []

for mouse, batch, group in zip(subjects_df['Subject'], subjects_df['Batch'], subjects_df['Group']):
    print(f"--- {mouse} {batch} {group} ---")
    fiberbehav_file = repo_path / f'{batch}_{mouse}_fiberbehavconcatnotderived.csv'

    excluded_ids = set(excluded_subjects_df['Subject'].astype(str))
    if str(mouse) in excluded_ids:
        print(f"Mouse {mouse} excluded")
        continue
    # Getting mean and AUC of dFF during behaviours
    if not fiberbehav_file.exists():
        print(f"File not found, skipping.")
        continue

    fiberbehav_df = pd.read_csv(fiberbehav_file, index_col=0)

    for records, use_zscore in [(dFF_records_raw_concat, False), (dFF_records_zscored_concat, True)]:

        record = quantif.extract_dff_summary(
            fiberbehav_df = fiberbehav_df,
            mouse         = mouse,
            batch         = batch,
            group         = group,
            zone_cols     = ['Stim Novel', 'Stim HC', 'Stim Clean'],
            behav_cols    = ['Sniff Novel', 'Sniff HC', 'Sniff Clean'],
            baseline_col  = None,
            merge_into    = {'Sniff Novel': 'Stim Novel', 'Sniff HC':'Stim HC', 'Sniff Clean':'Stim Clean'},
            dff_col       = 'dFF',
            fps           = 10,
            use_zscore    = use_zscore,
            stim_sniff_pairs={'Sniff Novel': 'Stim Novel', 'Sniff HC':'Stim HC', 'Sniff Clean':'Stim Clean'}
        )
        records.append(record)

    if '560 dFF' in fiberbehav_df.columns:
        for records, use_zscore in [(dFF_records_raw_560_concat, False), (dFF_records_zscored_560_concat, True)]:
            record = quantif.extract_dff_summary(
                fiberbehav_df = fiberbehav_df,
                mouse         = mouse,
                batch         = batch,
                group         = group,
                zone_cols     = ['Stim Novel', 'Stim HC', 'Stim Clean'],
                behav_cols    = ['Sniff Novel', 'Sniff HC', 'Sniff Clean'],
                baseline_col  = None,
                merge_into    = {'Sniff Novel': 'Stim Novel', 'Sniff HC':'Stim HC', 'Sniff Clean':'Stim Clean'},
                dff_col       = '560 dFF',
                fps           = 10,
                use_zscore    = use_zscore,
                stim_sniff_pairs={'Sniff Novel': 'Stim Novel', 'Sniff HC':'Stim HC', 'Sniff Clean':'Stim Clean'}
            )
            records.append(record)


# ── Save to Excel ─────────────────────────────────────────────────────────────
pd.DataFrame(dFF_records_raw_concat).to_excel(
    repo_path / 'dFF_summary_raw_concat.xlsx', index=False)
pd.DataFrame(dFF_records_zscored_concat).to_excel(
    repo_path / 'dFF_summary_zscored_concat.xlsx', index=False)
if '560 dFF' in fiberbehav_df.columns:
    pd.DataFrame(dFF_records_raw_560_concat).to_excel(
        repo_path / 'dFF_summary_raw_560_concat.xlsx', index=False)
    pd.DataFrame(dFF_records_zscored_560_concat).to_excel(
        repo_path / 'dFF_summary_zscored_560_concat.xlsx', index=False)
print(f"Saved {len(dFF_records_raw_560_concat)} mice to dFF_summary_concat.xlsx and dFF_summary_zscored_concat.xlsx\nData saved in: {repo_path}")
# %%
