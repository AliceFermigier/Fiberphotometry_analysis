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

from scripts.loader import analysis_path, data_path, proto_df, subjects_df, batches

#%% 2 - ANALYSIS - BEHAVIOUR
############################
automated_alignment = True
arena_analysis = False
dlc_data = False
boris = True

#filter characteristics
ORDER = 4
CUT_FREQ = None #in Hz

#threshold to fuse behaviour if bouts are too close, in secs
THRESH_S = 0
#threshold for PETH : if events are too short do not plot them and do not include them in PETH, in seconds
EVENT_TIME_THRESHOLD = 0

exp = 'Social_Interaction'
exp_path = analysis_path / exp
datapath_exp_dict = nom.get_experiment_data_path(batches, proto_df, data_path, exp)

# 2.1 - Align with behaviour, create corresponding excel, plot fiberpho data with behaviour
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
        rawdata_path = data_path_exp / f'{mouse}.doric'
        led_flashes_path = data_path_exp / f'{mouse}_fiber.csv'
        deinterleaved_raw_path = pp_path / f'{mouse}_deinterleaved.csv'
        dlc_path = behav_path_exp / f'{mouse}DLC_resnet50_FiberMEC_EPMMay14shuffle1_100000_filtered.csv'
        fiberpho_path = pp_path / f'{mouse}_dFFfilt.csv'
        boris_path = behav_path_exp / f'behav_boris_{mouse}.csv' 
        
        # Arena boundaries
        try:
            arena_coordinates_file = [f for f in os.listdir(behav_path_exp) if f.endswith('.json')][0]
            print(f'Arena boundaries : {arena_coordinates_file}')
            arena_coordinates_path = behav_path_exp / arena_coordinates_file
            with open(arena_coordinates_path, 'r') as f:
                arena_coordinates = json.load(f)
        except Exception as e:
            print(f'[!] Arena boundary file missing or invalid for {mouse}: {e}')
            arena_coordinates = {}

        # DLC data
        coordinates_df = None
        if arena_analysis and dlc_data:
            try:
                print('Get DLC data')
                coordinates_df = mp.get_dlc_data(dlc_path, threshold=0.99)
            except Exception as e:
                print(f'[!] DLC file error for {mouse}: {e}')

        # BORIS data
        boris_df = None
        list_BOI = ['Decoy']
        if boris:
            print('Get BORIS data')
            try:
                boris_df = pd.read_csv(boris_path)
                list_BOI = [col for col in boris_df.columns if col != 'time']
            except Exception as e:
                print(f'[!] BORIS file not found, creating decoy for {mouse}: {e}')
                deinterleaved_df = pd.read_csv(deinterleaved_raw_path)
                boris_df = pd.DataFrame({
                    'Time(s)': np.arange(deinterleaved_df['Time(s)'].iloc[0],
                                         deinterleaved_df['Time(s)'].iloc[-1],
                                         0.01),
                    'Decoy': 0
                })
        
        # Time alignment
        behav_df = None
        if not automated_alignment:
            print('Get camera flashes')
            camera_df = cp.get_camera_flashes(rawdata_path)
            print(camera_df)

            if not camera_df.empty:
                if coordinates_df is not None:
                    print('Aligning time with coordinates')
                    coordinates_df = cp.align_camera_flashes(coordinates_df, camera_df)
                elif boris_df is not None:
                    print('Aligning time with BORIS scoring')
                    behav_df = cp.align_camera_flashes(boris_df, camera_df)
                else:
                    raise ValueError("No valid data to align with camera.")
            else:
                print('[!] No camera flashes found, falling back to automated alignment.')
                automated_alignment = True

        if automated_alignment:
            print('Automated alignment')
            led_df = cp.get_led_flashes_from_csv(led_flashes_path)
            deinterleaved_df = pd.read_csv(deinterleaved_raw_path)
            aligned_deinterleaved_df = cp.align_fiber_with_led_flashes(deinterleaved_df, led_df)
            if boris_df is not None:
                behav_df = bp.correct_time_behav(deinterleaved_df, boris_df)
            else:
                raise ValueError("BORIS data required for automated alignment.")

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
