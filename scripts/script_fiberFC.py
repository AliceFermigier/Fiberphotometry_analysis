# -*- coding: utf-8 -*-
"""
2 - ANALYSIS WITH LICKS AND AIRPUFF FILES

@author: alice fermigier
"""

#%%IMPORTED
###########

import pandas as pd
import numpy as np
import importlib
import os
from pathlib import Path
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
import modules.behaviour.get_lick_and_airpuff_ports_coordinates as getlap
importlib.reload(getlap)
import modules.behaviour.camera_processing as cp
importlib.reload(cp)
import modules.common.clean_signal as cs
importlib.reload(cs)
import modules.behaviour.lick_detection as ld
importlib.reload(ld)
import modules.behaviour.behaviour_metrics as bm
importlib.reload(bm)
import modules.behaviour.fear_conditioning as fc
importlib.reload(fc)

from scripts.loader import analysis_path, experiment_path, data_path, proto_df, subjects_df, batches

#%% 2 - ANALYSIS - BEHAVIOUR
############################

#filter characteristics
ORDER = 4
CUT_FREQ = None #in Hz

#threshold to fuse behaviour if bouts are too close, in secs
THRESH_S = 0
#threshold for PETH : if events are too short do not plot them and do not include them in PETH, in seconds
EVENT_TIME_THRESHOLD = 0

exp = 'Fear_Conditioning'
if 'Conditioning' in exp:
    list_BOI = ['Shock','Freezing','CS+','CS-']
    dlc_suffix = 'DLC_resnet50_FearConditioning18shuffle1_100000'
    sheet = 'Conditioning'
else:
    list_BOI = ['Freezing','CS+','CS-']
    dlc_suffix = 'DLC_resnet50_FearHab18shuffle1_100000'
    if 'Habituation' in exp:
        sheet = 'Habituation'
    else:
        sheet = 'Retrieval'
exp_path = analysis_path / exp
datapath_exp_dict = nom.get_experiment_data_path(batches, proto_df, data_path, exp)

#%% 2.2 - Align with behaviour, create corresponding excel, plot fiberpho data with behaviour
print('###################')
print(f'EXPERIMENT : {exp}')
print('###################')

dlc_data = True

# Create repository path where fiberbehav data will be stored
repo_path = exp_path / f'length{EVENT_TIME_THRESHOLD}_interbout{THRESH_S}_o{ORDER}f{CUT_FREQ}'
repo_path.mkdir(exist_ok=True)

# Loop through each mouse in the subject DataFrame
for mouse, batch in zip(subjects_df['Subject'], subjects_df['Batch']):
    print("-----------------------------") 
    print(f'BATCH : {batch}, MOUSE : {mouse}')
    print("-----------------------------")

    data_path_exp = datapath_exp_dict[batch]
    pp_path = data_path_exp / 'Preprocessing'
    behav_path_exp = data_path_exp / 'Behaviour'

    # Define paths for raw, behavioral, and fiberphotometry data
    protocol_file = experiment_path / "fear_protocol.xlsx"
    led_flashes_path = data_path_exp / f'miniscope_sync_{mouse}.csv'
    camera_flashes_path = data_path_exp / f'camera_flashes_{mouse}.csv'
    LED3_path = data_path_exp / f'LED3_{mouse}.csv'
    deinterleaved_raw_path = pp_path / f'{mouse}_deinterleaved.csv'
    fiberpho_path = pp_path / f'{mouse}_dFFfilt.csv'
    dlc_path = behav_path_exp / f'{mouse}{dlc_suffix}.csv'

    led_df = cp.get_timestamps_from_bonsai_csv(led_flashes_path) # gets led flashes from Bonsai files
    deinterleaved_df = pd.read_csv(deinterleaved_raw_path)
    time_gap = cp.time_gap(deinterleaved_df, led_df)

    fiberpho_df = pd.read_csv(fiberpho_path)

    # Parse protocol
    print('Aligning Imetronic data')
    proto = fc.parse_protocol_sheet(protocol_file, sheet)

    # Get absolute protocol start from bonsai
    protocol_start_df = cp.get_timestamps_from_bonsai_csv(LED3_path)
    protocol_start_df = cp.correct_behav_timestamps(protocol_start_df, time_gap)
    protocol_start = protocol_start_df['Time(s)'].values[0]

    # Convert relative intervals to absolute Bonsai time
    cs_plus_abs = fc.convert_to_absolute(proto["CS+"], protocol_start)
    cs_minus_abs = fc.convert_to_absolute(proto["CS-"], protocol_start)
    led3_abs = fc.convert_to_absolute(proto["LED3"], protocol_start)
    if 'Shock' in list_BOI:
        shock_abs = fc.convert_to_absolute(proto["Shock"], protocol_start)

    # Add interval columns to your fiberphotometry data
    fiberpho_df = fc.add_interval_column(fiberpho_df, cs_plus_abs, "CS_plus")
    fiberpho_df = fc.add_interval_column(fiberpho_df, cs_minus_abs, "CS_minus")
    fiberpho_df = fc.add_interval_column(fiberpho_df, led3_abs, "Protocol_start")
    if 'Shock' in list_BOI:
        fiberpho_df = fc.add_interval_column(fiberpho_df, shock_abs, "Shock")
    print('Imetronic data aligned')

    # DLC data
    frame_times_df = cp.get_timestamps_from_bonsai_csv(camera_flashes_path)
    frame_times_df = cp.correct_behav_timestamps(frame_times_df, time_gap)
    coordinates_df = None
    if dlc_data:
        try:
            print('Getting DLC data')
            coordinates_df = mp.get_dlc_data(dlc_path, threshold=0.95)
            coordinates_df = cp.align_camera_flashes(coordinates_df, frame_times_df)
        except Exception as e:
            print(f'[!] DLC file error for {mouse}: {e}')

    # Align DLC and fiber data
    print('Aligning fiberphotometry and behaviour data')
    fiberbehav_df = bp.align_dlc_to_fiber(fiberbehav_df, coordinates_df)

    # Compute freezing bouts using DLC data
    print('Computing freezing bouts')

    # Post-process data (fuse behaviours that are too close and delete the ones that are too short)
    fiberbehav_df = bp.behav_process(fiberbehav_df, list_BOI, THRESH_S, EVENT_TIME_THRESHOLD)

    # Save outputs
    fiberbehav_notderived_path = repo_path / f'{batch}_{mouse}_fiberbehavnotderived.csv'
    fiberbehav_path = repo_path / f'{batch}_{mouse}_fiberbehav.csv'
    fiberbehav_df.to_csv(fiberbehav_notderived_path, index=False)

    dfiberbehav_df = bp.derive(fiberbehav_df, list_BOI)