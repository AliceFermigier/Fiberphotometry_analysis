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
import matplotlib.pyplot as plt
from pathlib import Path
import modules.common.switch_matplotlib_backends as smb
importlib.reload(smb)
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
import modules.behaviour.get_video_scale as getvid
importlib.reload(getvid)
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
import modules.common.quantification as quantif
importlib.reload(quantif)

from scripts.loader import analysis_path, experiment_path, data_path, proto_df, subjects_df, batches

#%% 2 - ANALYSIS - BEHAVIOUR
############################

#filter characteristics
ORDER = 4
CUT_FREQ = 20 #in Hz

#threshold to fuse behaviour if bouts are too close, in secs
THRESH_S = 2
#threshold for PETH : if events are too short do not plot them and do not include them in PETH, in seconds
EVENT_TIME_THRESHOLD = 0

exp = 'FearConditioning'
if 'Cond' in exp:
    list_BOI = ['Freezing','Shock','CS+','CS-']
    dlc_suffix = 'DLC_Resnet50_Fear_conditioningMar2shuffle1_snapshot_110_filtered'
    sheet = 'Conditioning'
else:
    list_BOI = ['Freezing','CS+','CS-']
    if 'Hab' in exp:
        sheet = 'Habituation'
        dlc_suffix = 'DLC_Resnet50_20260220_Fear_hab_and_retFeb20shuffle1_snapshot_090_filtered'
    else:
        sheet = 'Retrieval'
        dlc_suffix = 'DLC_Resnet50_20260220_Fear_hab_and_retFeb20shuffle1_snapshot_best-90_filtered'
exp_path = analysis_path / exp
datapath_exp_dict = nom.get_experiment_data_path(batches, proto_df, data_path, exp)

#%% 2.1 - Get scale and area coordinates for each video

for mouse, batch in zip(subjects_df['Subject'], subjects_df['Batch']):
    print("-----------------------------") 
    print(f'BATCH : {batch}, MOUSE : {mouse}')
    print("-----------------------------")
    
    data_path_exp = datapath_exp_dict[batch]
    behav_path_exp = data_path_exp / 'Behaviour'
    video_path = behav_path_exp / f"{mouse}.avi"
    arena_json = behav_path_exp / f"{mouse}_arena_coordinates.json"
    real_world_distance_cm=30

    #Indicate the arena boundaries and the coordinates of the known distance
    if arena_json.is_file():
        print('Arena json already exists')
    else:
        print('Get arena coordinates')
        plt = smb.with_qt5agg()
        scale_and_coords = getvid.get_scale_and_arena_rect(video_path, real_world_distance_cm)
        getvid.save_to_json(scale_and_coords, arena_json)
plt = smb.with_agg()

#%% 2.2 - Align with behaviour, create corresponding excel, plot fiberpho data with behaviour
print('###################')
print(f'EXPERIMENT : {exp}')
print('###################')

dlc_data = True
dual_color = True

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
    save_dir_QC = repo_path / 'Camera_alignment_QC'

    # Define paths for raw, behavioral, and fiberphotometry data
    raw_doric_path = data_path_exp / f'{mouse}_0000.doric'
    protocol_file = experiment_path / "fear_protocol.xlsx"
    led_flashes_path = data_path_exp / f'miniscope_sync_{mouse}.csv'
    camera_flashes_path = data_path_exp / f'camera_flashes_{mouse}.csv'
    LED3_path = data_path_exp / f'LED3_{mouse}.csv'
    deinterleaved_raw_path = pp_path / f'{mouse}_deinterleaved.csv'
    fiberpho_path = pp_path / f'{mouse}_dFF_corrected_final.csv'
    dlc_path = behav_path_exp / f'{mouse}{dlc_suffix}.csv'
    scale_json = behav_path_exp / f"{mouse}_arena_coordinates.json"

    if not os.path.isfile(fiberpho_path):
        print(f"Fiberphotometry data {fiberpho_path} does not exist. Skipping.")
        continue

    fiberpho_df = pd.read_csv(fiberpho_path)
    if CUT_FREQ != None:
        print(f"Filtering dFF data : order = {ORDER}; cutting frequency = {CUT_FREQ}")
        fiberpho_df = cs.lowpass_dFF(fiberpho_df, dual_color, order = 2, cut_freq = 6)

    # Parse protocol
    print('Aligning Imetronic data')
    proto = fc.parse_protocol_sheet(protocol_file, sheet)

    # Extract sync channel from Doric raw data and Bonsai corresponding sync data
    ttl_sync_df = cp.extract_sync_channel(raw_doric_path, sync_channel = "DIO04")
    led_flashes_path = data_path_exp / f'miniscope_sync_{mouse}.csv'
    led_df = cp.get_timestamps_from_bonsai_csv(led_flashes_path) 
    deinterleaved_df = pd.read_csv(deinterleaved_raw_path)
    # Compute linear regression to correct for differences between clocks
    print(f"Syncing Doric and Bonsai clocks")
    slope, intercept = cp.time_mapping(ttl_sync_df, led_df)

    # Get absolute protocol start from bonsai
    try:
        protocol_start_df = cp.get_timestamps_from_bonsai_csv(LED3_path)
        protocol_start_df = cp.correct_behav_timestamps(protocol_start_df, slope, intercept)
        protocol_start, effective_slope = fc.get_protocol_remapping(protocol_start_df, proto)
    except Exception as e:
        print(f"[!] LED3 unavailable ({e}), falling back to manual shock scoring")
        score_path = behav_path_exp / f'{batch}_{mouse}_manual_shock_scoring.json'
        with open(score_path, 'r') as f:
            score_json = json.load(f)
        protocol_start = score_json['protocol_start']
        effective_slope = score_json['effective_slope']

    # Convert relative intervals to absolute Bonsai time
    cs_plus_abs  = fc.convert_to_absolute(proto["CS+"],   protocol_start, effective_slope)
    cs_minus_abs = fc.convert_to_absolute(proto["CS-"],   protocol_start, effective_slope)
    led3_abs     = fc.convert_to_absolute(proto["LED3"],  protocol_start, effective_slope)
    if 'Shock' in list_BOI:
        shock_abs    = fc.convert_to_absolute(proto["Shock"], protocol_start, effective_slope)

    # Add interval columns to fiberphotometry data
    fiberbehav_df = fc.add_interval_column(fiberpho_df, cs_plus_abs, "CS+")
    fiberbehav_df = fc.add_interval_column(fiberbehav_df, cs_minus_abs, "CS-")
    fiberbehav_df = fc.add_interval_column(fiberbehav_df, led3_abs, "Protocol_start")
    if 'Shock' in list_BOI:
        fiberbehav_df = fc.add_interval_column(fiberbehav_df, shock_abs, "Shock")
    print('Imetronic data aligned')

    # DLC data
    frame_times_df = cp.get_timestamps_from_bonsai_csv(camera_flashes_path)
    frame_times_df = cp.correct_behav_timestamps(frame_times_df, slope, intercept)
    coordinates_df = None

    # Arena scale
    try:
        with open(scale_json, 'r') as f:
            arena_scale = json.load(f)
    except Exception as e:
        print(f'[!] Scale file missing or invalid for {mouse}: {e}')
        arena_scale = {}

    if dlc_data:

        print('Getting DLC data')
        coordinates_df = mp.get_dlc_data(dlc_path, threshold=0.6)
        coordinates_df = cp.align_camera_flashes(coordinates_df, 
                                                    frame_times_df, 
                                                    mouse=mouse,  
                                                    batch=batch, 
                                                    save_dir_QC=save_dir_QC, 
                                                    expected_fps=arena_scale["Video_fps"])
        coordinates_df, _ = fc.filter_outside_arena(coordinates_df,
                                                arena_scale, 
                                                bodyparts=['center','tail_base'], 
                                                margin_px=0)


        # Compute freezing bouts using DLC data
        print('Computing freezing bouts')
        behav_df = fc.detect_freezing(coordinates_df, arena_scale, threshold=0.5)

        # Align DLC and fiber data
        print('Aligning fiberphotometry and behaviour data')
        fiberbehav_df = bp.align_dlc_to_fiber(fiberbehav_df, behav_df)

    # Post-process data (fuse behaviours that are too close and delete the ones that are too short)
    fiberbehav_df = bp.behav_process(fiberbehav_df, list_BOI, THRESH_S, EVENT_TIME_THRESHOLD)

    # Save outputs
    fiberbehav_notderived_path = repo_path / f'{batch}_{mouse}_fiberbehavnotderived.csv'
    fiberbehav_path = repo_path / f'{batch}_{mouse}_fiberbehav.csv'
    fiberbehav_df.to_csv(fiberbehav_notderived_path)

    dfiberbehav_df = bp.derive(fiberbehav_df, list_BOI)
    dfiberbehav_df.to_csv(fiberbehav_path)

    # Plotting
    fig = bp.plot_fiberpho_behav(
        dfiberbehav_df, list_BOI, exp, mouse,
        THRESH_S, EVENT_TIME_THRESHOLD, batch,
        scaled = False, speed_ylim = [-1, 80]
    )
    fig.savefig(repo_path / f'{batch}_{mouse}_fiberbehav.pdf')
    fig.savefig(repo_path / f'{batch}_{mouse}_fiberbehav.png')
    plt.close()


#%% 2.3 - Plot behavioural metrics

print('###################')
print(f'EXPERIMENT : {exp}')
print('###################')

behaviors_to_plot = list_BOI

# Create repository path where data will be stored
exp_path = analysis_path / exp
repo_path = exp_path / f'length0_interbout0_o{ORDER}f{CUT_FREQ}'
behavioural_analysis_path = exp_path / 'Behavioural_analysis'
behavioural_analysis_path.mkdir(exist_ok=True)
fig_dir = behavioural_analysis_path / 'Figures'

all_metrics = {}

# Loop through each mouse in the subject DataFrame
for mouse, batch in zip(subjects_df['Subject'], subjects_df['Batch']):
    print("-----------------------------") 
    print(f'BATCH : {batch}, MOUSE : {mouse}')
    print("-----------------------------")

    fiberbehav_notderived_path = repo_path / f'{batch}_{mouse}_fiberbehavnotderived.csv'

    if not os.path.isfile(fiberbehav_notderived_path):
        print(f"Fiber behav file {fiberbehav_notderived_path} not found. Skipping.")
        continue

    fiberbehav_notderived_df = pd.read_csv(fiberbehav_notderived_path)

    metrics = bm.compute_behavior_metrics_FC(fiberbehav_notderived_df)
    all_metrics[mouse] = metrics

    # Plot raster
    bm.plot_behavior_raster(fiberbehav_notderived_df, 
                            mouse, batch, behaviors=behaviors_to_plot, 
                            save_dir=fig_dir)

    print(f"\n=== Analysis complete for mouse {batch}_{mouse}. Plots stored in {behavioural_analysis_path}. ===")

# Export behavioral metrics to Excel
print("\nExporting all behavioral metrics to Excel...")

try:
    # Concatenate all metrics into a single DataFrame
    metrics_list = []
    for mouse, df in all_metrics.items():
        df = df.copy()
        df["Mouse"] = mouse
        metrics_list.append(df)

    all_metrics_df = pd.concat(metrics_list, ignore_index=True)

    # Define output path
    excel_path = behavioural_analysis_path / "behavioral_metrics.xlsx"

    # Export to Excel
    all_metrics_df.to_excel(excel_path, index=False)

    print(f"✔ Behavioral metrics successfully exported to:\n    {excel_path}")

except Exception as e:
    print(f"[!] Error while exporting behavioral metrics: {e}")

# %% 2.3 - Quantify dFF during behaviours.

subjects_df['Group'] = subjects_df['Group'].fillna('') # if group = Nan, replaces it with an empty string

#Load excluded subjects
excluded_subjects_df = pd.read_excel(experiment_path / 'subjects.xlsx', 
                                     sheet_name=f'Excluded_{exp}')

dFF_records_raw = []
dFF_records_zscored = []
dFF_records_raw_560 = []
dFF_records_zscored_560 = []

for mouse, batch, group in zip(subjects_df['Subject'], subjects_df['Batch'], subjects_df['Group']):
    print(f"--- {mouse} {batch} {group} ---")
    fiberbehav_file = repo_path / f'{batch}_{mouse}_fiberbehavnotderived.csv'

    if int(mouse) in excluded_subjects_df['Subject'].values:
        print(f"Mouse {mouse} excluded")
        continue

    # Getting mean and AUC of dFF during behaviours
    if not fiberbehav_file.exists():
        print(f"  File not found, skipping.")
        continue

    fiberbehav_df = pd.read_csv(fiberbehav_file, index_col=0)

    for records, use_zscore in [(dFF_records_raw, False), (dFF_records_zscored, True)]:

        record = quantif.extract_dff_summary(
            fiberbehav_df = fiberbehav_df,
            mouse         = mouse,
            batch         = batch,
            group         = group,
            zone_cols     = [],
            behav_cols    = list_BOI,
            dff_col       = 'dFF',
            fps           = arena_scale['Video_fps'], 
            use_zscore    = use_zscore,
        )
        records.append(record)

    if '560 dFF' in fiberbehav_df.columns:
        for records, use_zscore in [(dFF_records_raw_560, False), (dFF_records_zscored_560, True)]:
            record = quantif.extract_dff_summary(
                fiberbehav_df = fiberbehav_df,
                mouse         = mouse,
                batch         = batch,
                group         = group,
                zone_cols     = [],
                behav_cols    = list_BOI,
                dff_col       = '560 dFF',
                fps           = arena_scale['Video_fps'],
                use_zscore    = use_zscore,
            )
            records.append(record)


# ── Save to Excel ─────────────────────────────────────────────────────────────
pd.DataFrame(dFF_records_raw).to_excel(
    repo_path / 'dFF_summary_raw.xlsx', index=False)
pd.DataFrame(dFF_records_zscored).to_excel(
    repo_path / 'dFF_summary_zscored.xlsx', index=False)
if '560 dFF' in fiberbehav_df.columns:
    pd.DataFrame(dFF_records_raw_560).to_excel(
        repo_path / 'dFF_summary_raw_560.xlsx', index=False)
    pd.DataFrame(dFF_records_zscored_560).to_excel(
        repo_path / 'dFF_summary_zscored_560.xlsx', index=False)
print(f"Saved {len(dFF_records_raw)} mice to dFF_summary.xlsx and dFF_summary_zscored.xlsx")
# %%
