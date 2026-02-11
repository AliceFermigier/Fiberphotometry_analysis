# -*- coding: utf-8 -*-
"""
Created on Sun Dec 15 16:53:24 2024

To run fiberphotometry analysis with EPM data
2 - ANALYSIS WITH EPM DeepLabCut data

@author: alice fermigier
"""

#%%IMPORTED
###########

import pandas as pd
import numpy as np
from pathlib import Path
import os
import warnings
import json
import importlib
import modules.common.switch_matplotlib_backends as smb
importlib.reload(smb)
plt = smb.with_agg() #imports matplotlib.pyplot with Agg backend
import json

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
import modules.behaviour.get_epm_coordinates as getepm
importlib.reload(getepm)
import modules.behaviour.get_video_scale as getvid
importlib.reload(getvid)

from scripts.loader import analysis_path, data_path, proto_df, subjects_df, batches

#%% 2 - ANALYSIS - BEHAVIOUR
############################

# If data was recorded via Bonsai, with a camera_flashes csv file, set to True (else camera flashes were directly recorded in the Doric console)
bonsai_setup = True

#filter characteristics
ORDER = 4
CUT_FREQ = None #in Hz

#threshold to fuse behaviour if bouts are too close, in secs
THRESH_S = 0
#threshold for PETH : if events are too short do not plot them and do not include them in PETH, in seconds
EVENT_TIME_THRESHOLD = 0

exp = 'EPM'
exp_path = analysis_path / exp
datapath_exp_dict = nom.get_experiment_data_path(batches, proto_df, data_path, exp)

# 2.1 - Align with behaviour, create corresponding excel, plot fiberpho data with behaviour
print('###################')
print(f'EXPERIMENT : {exp}')
print('###################')

# Create repository path where fiberbehav data will be stored
repo_path = exp_path / f'length{EVENT_TIME_THRESHOLD}_interbout{THRESH_S}_o{ORDER}f{CUT_FREQ}'
repo_path.mkdir(exist_ok=True)

#%% 2.1.1 - Get coordinates of EPM arena

for mouse, batch in zip(subjects_df['Subject'], subjects_df['Batch']):
    print("-----------------------------") 
    print(f'BATCH : {batch}, MOUSE : {mouse}')
    print("-----------------------------")
    
    data_path_exp = datapath_exp_dict[batch]
    behav_path_exp = data_path_exp / 'Behaviour'
    video_path = data_path_exp / f"{mouse}.avi"
    output_json = behav_path_exp / f"{mouse}_epm_coordinates.json"

    #Indicate the emplacement of the lickport and the two airpuff ports
    if output_json.is_file():
        print('EPM json already exists')
    else:
        print('Get EPM coordinates')
        plt = smb.with_qt5agg() 
        arena = getepm.define_epm_boundaries(video_path)
        getepm.save_boundaries_to_json(arena, output_json)
plt = smb.with_agg()

#%% 2.1.2 - Get scale and area coordinates for each video

for mouse, batch in zip(subjects_df['Subject'], subjects_df['Batch']):
    print("-----------------------------") 
    print(f'BATCH : {batch}, MOUSE : {mouse}')
    print("-----------------------------")
    
    data_path_exp = datapath_exp_dict[batch]
    behav_path_exp = data_path_exp / 'Behaviour'
    video_path = data_path_exp / f"{mouse}.avi"
    scale_json = behav_path_exp / f"{mouse}_scale_and_arena_coordinates.json"
    real_world_distance_cm=20
    real_world_distance_name="2 open arms lenght"

    #Indicate the arena boundaries and the coordinates of the known distance
    if scale_json.is_file():
        print('Scale json already exists')
    else:
        print('Get arena coordinates')
        plt = smb.with_qt5agg()
        scale_and_coords = getvid.get_scale_and_arena_rect(video_path, real_world_distance_cm, real_world_distance_name)
        getvid.save_to_json(scale_and_coords, scale_json)
plt = smb.with_agg()

#%% 2.2 - Analyze fiberpho data alongside EPM data

for mouse, batch in zip(subjects_df['Subject'], subjects_df['Batch']):
    print("-----------------------------") 
    print(f'BATCH : {batch}, MOUSE : {mouse}')
    print("-----------------------------")
    try:
        data_path_exp = datapath_exp_dict[batch]
        pp_path = data_path_exp / 'Preprocessing'
        behav_path_exp = data_path_exp / 'Behaviour'

        # Define paths for raw, behavioral, and fiberphotometry data
        camera_flashes_path = data_path_exp / f'camera_flashes_{mouse}.csv'
        rawdata_path = data_path_exp / f'{mouse}_0000.doric'
        led_flashes_path = data_path_exp / f'miniscope_sync_{mouse}.csv'
        deinterleaved_raw_path = pp_path / f'{mouse}_deinterleaved.csv'
        dlc_path = behav_path_exp / f'{mouse}DLC_Resnet101_RewardBox_FiberMECNov12shuffle5_snapshot_090_filtered.csv'
        fiberpho_path = pp_path / f'{mouse}_dFFfilt.csv'
        
        # Arena boundaries
        try:
            arena_json = behav_path_exp / f"{mouse}_epm_coordinates.json"
            with open(arena_json, 'r') as f:
                arena_coordinates = json.load(f)
        except Exception as e:
            print(f'[!] Arena boundary file missing or invalid for {mouse}: {e}')
            arena_coordinates = {}

        # Arena scale
        try:
            scale_json = behav_path_exp / f"{mouse}_scale_and_arena_coordinates.json"
            with open(scale_json, 'r') as f:
                arena_scale = json.load(f)
        except Exception as e:
            print(f'[!] Scale file missing or invalid for {mouse}: {e}')
            arena_scale = {}

        # DLC data
        coordinates_df = None
        try:
            print('Get DLC data')
            coordinates_df = mp.get_dlc_data(dlc_path, threshold=0.95)
        except Exception as e:
            print(f'[!] DLC file error for {mouse}: {e}')

        # Time alignment
        behav_df = None
        if bonsai_setup:
                led_df = cp.get_timestamps_from_bonsai_csv(led_flashes_path) # gets led flashes from Bonsai files
                deinterleaved_df = pd.read_csv(deinterleaved_raw_path)
                time_gap = cp.time_gap(deinterleaved_df, led_df)

                fiberpho_df = pd.read_csv(fiberpho_path)

        else :
            print('Get camera flashes')
            camera_df = cp.get_camera_flashes(rawdata_path)
            print(camera_df)

            if not camera_df.empty:
                if coordinates_df is not None:
                    print('Aligning time with coordinates')
                    coordinates_df = cp.align_camera_flashes(coordinates_df, camera_df)
                else:
                    raise ValueError("No valid data to align with camera.")
            else:
                print('[!] No camera flashes found, falling back to automated alignment.')
                automated_alignment = True

        # Analyze EPM position
        if 'EPM' in exp and coordinates_df is not None:
            print(f'Analyzing mouse position for {mouse}')
            list_BOI = ['Open arm', 'Closed arm', 'Center']
            behav_df = epm.analyze_mouse_position(coordinates_df, arena_coordinates, bodypart='nose')

        if behav_df is None:
            raise ValueError("Behavior dataframe is missing or invalid.")

        behav_path = behav_path_exp / f'behav_{mouse}.csv'
        behav_df.to_csv(behav_path)
        print(f'Behaviour file exported to {behav_path}')

        # Fiber photometry
        fiberpho = pd.read_csv(fiberpho_path)
        if CUT_FREQ is not None:
            fiberpho = pp.butterfilt(fiberpho, ORDER, CUT_FREQ)

        # Align behavior and fiber data
        print('Aligning fiberphotometry and behaviour data')
        fiberbehav_df = bp.align_behav(behav_df, fiberpho, list_BOI)
        fiberbehav_df = bp.behav_process(fiberbehav_df, list_BOI, THRESH_S, EVENT_TIME_THRESHOLD)

        # Save outputs
        fiberbehav_notderived_path = repo_path / f'{batch}_{mouse}_fiberbehavnotderived.csv'
        fiberbehav_path = repo_path / f'{batch}_{mouse}_fiberbehav.csv'
        fiberbehav_df.to_csv(fiberbehav_notderived_path, index=False)

        dfiberbehav_df = bp.derive(fiberbehav_df, list_BOI)
        dfiberbehav_df.to_csv(fiberbehav_path, index=False)

        # Plotting
        fig = bp.plot_fiberpho_behav(
            dfiberbehav_df, list_BOI, exp, mouse,
            THRESH_S, EVENT_TIME_THRESHOLD, batch,
            scaled = True
        )
        fig.savefig(repo_path / f'{batch}_{mouse}_fiberbehav.pdf')
        fig.savefig(repo_path / f'{batch}_{mouse}_fiberbehav.png')
        plt.close(fig)

    except Exception as e:
        print(f'[!] Error while processing mouse {mouse}: {e}')

print(f'\n✅ Analysis for {exp} complete.\nData saved in: {repo_path}')
                    
# %%
