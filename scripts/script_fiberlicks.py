# -*- coding: utf-8 -*-
"""
2 - ANALYSIS WITH LICKS AND AIRPUFF FILES

@author: alice fermigier
"""

#%%IMPORTED
###########

import pandas as pd
import numpy as np
import os
from pathlib import Path
import matplotlib.pyplot as plt
from dash import Dash, dcc, html
from dash.dependencies import Input, Output
import plotly.express as px
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
import modules.behaviour.lick_detection as ld
importlib.reload(ld)

from scripts.loader import analysis_path, data_path, proto_df, subjects_df, batches

#%% 2 - ANALYSIS - BEHAVIOUR
############################

dlc_data = False

#filter characteristics
ORDER = 4
CUT_FREQ = None #in Hz

#threshold to fuse behaviour if bouts are too close, in secs
THRESH_S = 6
#threshold for PETH : if events are too short do not plot them and do not include them in PETH, in seconds
EVENT_TIME_THRESHOLD = 0

#capacitance threshold to extract licks
THRESH_LICKS = 200

exp = 'Reward_Airpuffs'
list_BOI = ['Licks', 'Airpuffs']
exp_path = analysis_path / exp
datapath_exp_dict = nom.get_experiment_data_path(batches, proto_df, data_path, exp)


#%% 2.1 - Detect licks in capacitance data

# Load all mice data
mouse_data = {}
for mouse, batch in zip(subjects_df['Subject'], subjects_df['Batch']):
    licks_df = ld.load_mouse_data(mouse, batch, datapath_exp_dict)
    plot_df = licks_df.groupby(licks_df.index // 10).mean()
    mouse_data[f'{batch} - {mouse}'] = plot_df

# Create Dash app
app = Dash(__name__)
app.title = f"Licks Viewer - {exp}"

app.layout = html.Div([
    html.H3("Capacitance Traces"),
    
    # Dropdown to select mouse
    dcc.Dropdown(
        id='mouse-dropdown',
        options=[{'label': name, 'value': name} for name in mouse_data.keys()],
        value=list(mouse_data.keys())[0],  # default to first mouse
        clearable=False
    ),
    
    # Plot area
    dcc.Graph(id='licks-plot')
])

@app.callback(
    Output('licks-plot', 'figure'),
    Input('mouse-dropdown', 'value')
)
def update_plot(selected_mouse):
    df = mouse_data[selected_mouse]
    fig = px.line(df, x='time(s)', y='capacitance',
                  title=f'{selected_mouse} - Capacitance Trace')
    fig.add_hline(
        y=THRESH_LICKS,
        line_dash="dash",
        line_color="red",
        annotation_text=f"Threshold = {THRESH_LICKS}",
        annotation_position="top right"
    )
    fig.update_layout(
        xaxis_title='Time (s)',
        yaxis_title='Capacitance (a.u.)',
        template='plotly_white',
        showlegend=False
    )
    return fig

if __name__ == '__main__':
    app.run(debug=False, use_reloader=False)

#%% 2.1 - Detect licks in capacitance data

for mouse, batch in zip(subjects_df['Subject'], subjects_df['Batch']):
    print("-----------------------------") 
    print(f'BATCH : {batch}, MOUSE : {mouse}')
    print("-----------------------------")
    
    data_path_exp = datapath_exp_dict[batch]
    licks_file = data_path_exp / f'licks_from_cp_{mouse}.csv' # File to store licks

    # ------------------ #
    # Create the Dash app
    app = Dash(__name__)
    app.title = f"Licks - {mouse}"

    # Load data relative to mouse
    licks_df = ld.load_mouse_data(mouse, batch, datapath_exp_dict)

    # Create initial figure
    plot_df = licks_df.groupby(licks_df.index // 10).mean()
    fig = px.line(plot_df, x='time(s)', y='capacitance')
    # here plot a line with the threshold at THRESH_LICKS = 200

    # ------------------ #
    # App layout
    app.layout = html.Div([
        html.H4(f'{exp} {mouse}'),

        dcc.Graph(
            id='plot',
            figure=fig,
            config={'displayModeBar': True}
        ),
    ])

    # Run the server
    if __name__ == '__main__':
        app.run(debug=False, use_reloader=False)

#%% 2.2 - Align with behaviour, create corresponding excel, plot fiberpho data with behaviour
print('###################')
print(f'EXPERIMENT : {exp}')
print('###################')

# Create repository path where fiberbehav data will be stored
repo_path = exp_path / f'length{EVENT_TIME_THRESHOLD}_interbout{THRESH_S}_o{ORDER}f{CUT_FREQ}'
repo_path.mkdir(exist_ok=True)

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
        led_flashes_path = data_path_exp / f'miniscope_sync_{mouse}.csv'
        camera_flashes_path = data_path_exp / f'camera_flashes_{mouse}.csv'
        licks_path = data_path_exp / f'licks_{mouse}.csv'
        airpuff_path = data_path_exp / f'airpuffs_{mouse}.csv'
        deinterleaved_raw_path = pp_path / f'{mouse}_deinterleaved.csv'
        fiberpho_path = pp_path / f'{mouse}_dFFfilt.csv'
        dlc_path = behav_path_exp / f'{mouse}DLC_resnet50_FiberMEC_EPMMay14shuffle1_100000_filtered.csv'

        led_df = cp.get_timestamps_from_bonsai_csv(led_flashes_path)
        deinterleaved_df = pd.read_csv(deinterleaved_raw_path)
        time_gap = cp.time_gap(deinterleaved_df, led_df)

        fiberpho_df = pd.read_csv(fiberpho_path)

        # Align licks and airpuff timestamps to dFF data
        licks_df = cp.get_timestamps_from_bonsai_csv(licks_path)
        licks_df = cp.correct_behav_timestamps(licks_df, time_gap)
        fiberbehav_df = cp.align_behav_timestamps(fiberpho_df, licks_df, "Licks")

        if airpuff_path.exists():
            airpuff_df = cp.get_timestamps_from_bonsai_csv(airpuff_path)
            airpuff_df = cp.correct_behav_timestamps(airpuff_df, time_gap)
            fiberbehav_df = cp.align_behav_timestamps(fiberbehav_df, airpuff_df, "Airpuffs")

        fiberbehav_df = bp.behav_process(fiberbehav_df, list_BOI, THRESH_S, EVENT_TIME_THRESHOLD)

        # DLC data
        coordinates_df = None
        if dlc_data:
            try:
                print('Get DLC data')
                coordinates_df = mp.get_dlc_data(dlc_path, threshold=0.99)
            except Exception as e:
                print(f'[!] DLC file error for {mouse}: {e}')

        # Save outputs
        fiberbehav_notderived_path = repo_path / f'{batch}_{mouse}_fiberbehavnotderived.csv'
        fiberbehav_path = repo_path / f'{batch}_{mouse}_fiberbehav.csv'
        fiberbehav_df.to_csv(fiberbehav_notderived_path, index=False)

        dfiberbehav_df = bp.derive(fiberbehav_df, list_BOI)
        print(f'n_licks {mouse} : {len(np.where(dfiberbehav_df["Licks"]==1)[0])}')
        print(f'n_airpuffs {mouse} : {len(np.where(dfiberbehav_df["Airpuffs"]==1)[0])}')
        dfiberbehav_df.to_csv(fiberbehav_path, index=False)

        # Plotting
        fig = bp.plot_fiberpho_behav(
            dfiberbehav_df, list_BOI, exp, mouse,
            THRESH_S, EVENT_TIME_THRESHOLD, batch,
            scaled = False
        )
        fig.savefig(repo_path / f'{batch}_{mouse}_fiberbehav.pdf')
        fig.savefig(repo_path / f'{batch}_{mouse}_fiberbehav.png')
        plt.close(fig)

    except Exception as e:
        print(f'[!] Error while processing mouse {mouse}: {e}')

print(f'\n✅ Analysis for {exp} complete.\nData saved in: {repo_path}')
# %%
