# -*- coding: utf-8 -*-
"""
2 - ANALYSIS WITH LICKS AND AIRPUFF FILES

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
THRESH_S = 0
#threshold for PETH : if events are too short do not plot them and do not include them in PETH, in seconds
EVENT_TIME_THRESHOLD = 0

exp = 'Reward_Hab'
exp_path = analysis_path / exp
datapath_exp_dict = nom.get_experiment_data_path(batches, proto_df, data_path, exp)

#%% 2.1 - Detect licks in capacitance data

for mouse, batch in zip(subjects_df['Subject'], subjects_df['Batch']):
    print("-----------------------------") 
    print(f'BATCH : {batch}, MOUSE : {mouse}')
    print("-----------------------------")
    try:
        data_path_exp = datapath_exp_dict[batch]
        behav_path_exp = data_path_exp / 'Behaviour'
        capacitance_txt_path = data_path_exp / f'{mouse}.txt'
        licks_df = ld.txt_to_df(capacitance_txt_path)

        # Script to plot licks_df with Plotly, define threshold for licks, mark artifacts and remove them from final binary file 

    except Exception as e:
        print(f'[!] Error while processing capacitance data from mouse {mouse}: {e}')

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
        deinterleaved_raw_path = pp_path / f'{mouse}_deinterleaved.csv'
        fiberpho_path = pp_path / f'{mouse}_dFFfilt.csv'
        dlc_path = behav_path_exp / f'{mouse}DLC_resnet50_FiberMEC_EPMMay14shuffle1_100000_filtered.csv'

        led_df = cp.get_timestamps_from_bonsai_csv(led_flashes_path)
        camera_flashes_df = cp.get_timestamps_from_bonsai_csv(camera_flashes_path)
        deinterleaved_df = pd.read_csv(deinterleaved_raw_path)
        aligned_deinterleaved_df = cp.align_fiber_with_led_flashes(deinterleaved_df, led_df)

        # DLC data
        coordinates_df = None
        if dlc_data:
            try:
                print('Get DLC data')
                coordinates_df = mp.get_dlc_data(dlc_path, threshold=0.99)
            except Exception as e:
                print(f'[!] DLC file error for {mouse}: {e}')

    except Exception as e:
        print(f'[!] Error while processing mouse {mouse}: {e}')

print(f'\n✅ Analysis for {exp} complete.\nData saved in: {repo_path}')