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

from scripts.loader import analysis_path, data_path, proto_df, subjects_df, batches

#%% 2 - ANALYSIS - BEHAVIOUR
############################

#filter characteristics
ORDER = 4
CUT_FREQ = None #in Hz

#threshold to fuse behaviour if bouts are too close, in secs
THRESH_S = 5
#threshold for PETH : if events are too short do not plot them and do not include them in PETH, in seconds
EVENT_TIME_THRESHOLD = 0

exp = 'RewardAirpuff'
list_BOI = ['Licks', 'Licks_filtered', 'Nose_in_any_airport', 'Airpuffs']
#['Licks', 'Airpuffs']
exp_path = analysis_path / exp
datapath_exp_dict = nom.get_experiment_data_path(batches, proto_df, data_path, exp)

#%% 2.1.1 - Get coordinates of lickport and airpuff ports

for mouse, batch in zip(subjects_df['Subject'], subjects_df['Batch']):
    print("-----------------------------") 
    print(f'BATCH : {batch}, MOUSE : {mouse}')
    print("-----------------------------")
    
    data_path_exp = datapath_exp_dict[batch]
    behav_path_exp = data_path_exp / 'Behaviour'
    video_path = behav_path_exp / f"{mouse}.avi"
    output_json = behav_path_exp / f"{mouse}_ports_coordinates.json"

    #Indicate the emplacement of the lickport and the two airpuff ports
    if output_json.is_file():
        print('Ports json already exists')
    else:
        print('Get ports coordinates')
        plt = smb.with_qt5agg() 
        ports = getlap.define_ports(video_path)
        getlap.save_ports_to_json(ports, output_json)
plt.close('all') 
plt = smb.with_agg()

#%% 2.1.2 - Get scale and area coordinates for each video

for mouse, batch in zip(subjects_df['Subject'], subjects_df['Batch']):
    print("-----------------------------") 
    print(f'BATCH : {batch}, MOUSE : {mouse}')
    print("-----------------------------")
    
    data_path_exp = datapath_exp_dict[batch]
    behav_path_exp = data_path_exp / 'Behaviour'
    video_path = behav_path_exp / f"{mouse}.avi"
    arena_json = behav_path_exp / f"{mouse}_arena_coordinates.json"
    real_world_distance_cm=20

    #Indicate the arena boundaries and the coordinates of the known distance
    if arena_json.is_file():
        print('Arena json already exists')
    else:
        print('Get arena coordinates')
        plt = smb.with_qt5agg()
        scale_and_coords = getvid.get_scale_and_arena_rect(video_path, real_world_distance_cm)
        getvid.save_to_json(scale_and_coords, arena_json)
plt.close('all') 
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

    # Define paths for raw, behavioral, and fiberphotometry data
    camera_flashes_path = data_path_exp / f'camera_flashes_{mouse}.csv'
    licks_path = data_path_exp / f'licks_{mouse}.csv'
    airpuff_path = data_path_exp / f'airpuffs_{mouse}.csv'
    deinterleaved_raw_path = pp_path / f'{mouse}_deinterleaved.csv'
    fiberpho_path = pp_path / f'{mouse}_dFF_corrected.csv'
    dlc_path = behav_path_exp / f'{mouse}DLC_Resnet50_RewardAirpuff_20260217Feb17shuffle1_snapshot_100_filtered.csv'
    output_json = behav_path_exp / f"{mouse}_ports_coordinates.json"
    arena_json = behav_path_exp / f"{mouse}_arena_coordinates.json"

    led_flashes_path = data_path_exp / f'miniscope_sync_{mouse}.csv'
    led_df = cp.get_timestamps_from_bonsai_csv(led_flashes_path) 
    deinterleaved_df = pd.read_csv(deinterleaved_raw_path)
    slope, intercept = cp.time_gap(deinterleaved_df, led_df)

    fiberpho_df = pd.read_csv(fiberpho_path)

    # Align licks and airpuff timestamps to dFF data
    print("Aligning licks")
    licks_df = cp.get_timestamps_from_bonsai_csv(licks_path)
    licks_df = cp.correct_behav_timestamps(licks_df, slope, intercept)
    fiberbehav_df = cp.align_behav_timestamps(fiberpho_df, licks_df, "Licks")

    if airpuff_path.exists():
        print("Aligning airpuffs")
        airpuff_df = cp.get_timestamps_from_bonsai_csv(airpuff_path)
        airpuff_df = cp.correct_behav_timestamps(airpuff_df, slope, intercept)
        fiberbehav_df = cp.align_behav_timestamps(fiberbehav_df, airpuff_df, "Airpuffs")

    # DLC data
    print("Aligning camera frames")
    frame_times_df = cp.get_timestamps_from_bonsai_csv(camera_flashes_path)
    frame_times_df = cp.correct_behav_timestamps(frame_times_df, slope, intercept)
    coordinates_df = None
    if dlc_data:
        try:
            print('Get DLC data')
            coordinates_df = mp.get_dlc_data(dlc_path, threshold=0.85)
            coordinates_df = cp.align_camera_flashes(coordinates_df, frame_times_df)
        except Exception as e:
            print(f'[!] DLC file error for {mouse}: {e}')

        # Align DLC and fiber data
        # Also takes out any fiber data that wasn't taken in video
        print('Aligning fiberphotometry and behaviour data')
        fiberbehav_df = bp.align_dlc_to_fiber(fiberbehav_df, coordinates_df)

        # Clean licking data. Radius in cm.
        print('Cleaning licking data')
        ports = json.load(open(output_json, "r"))
        scale_and_coords = json.load(open(arena_json, "r"))
        fiberbehav_df = ld.filter_licking(fiberbehav_df, ports, scale_and_coords, lick_col="Licks", lick_radius_cm=0.5)

        # Scoring nose-in-airport time. Radius in cm.
        fiberbehav_df = ld.detect_airpuff_entry(fiberbehav_df, ports, scale_and_coords, radius_cm=3.0)

    # Post-process data (fuse behaviours that are too close and delete the ones that are too short)
    fiberbehav_df = bp.behav_process(fiberbehav_df, list_BOI, THRESH_S, EVENT_TIME_THRESHOLD)

    # Save outputs
    fiberbehav_notderived_path = repo_path / f'{batch}_{mouse}_fiberbehavnotderived.csv'
    fiberbehav_path = repo_path / f'{batch}_{mouse}_fiberbehav.csv'
    fiberbehav_df.to_csv(fiberbehav_notderived_path, index=False)

    dfiberbehav_df = bp.derive(fiberbehav_df, list_BOI)
    print(f'n_licks {mouse} : {len(np.where(dfiberbehav_df["Licks"]==1)[0])}')
    if dlc_data:
        print(f'n_licks_filtered {mouse} : {len(np.where(dfiberbehav_df["Licks_filtered"]==1)[0])}')
    if airpuff_path.exists():
        print(f'n_airpuffs {mouse} : {len(np.where(dfiberbehav_df["Airpuffs"]==1)[0])}')
    dfiberbehav_df.to_csv(fiberbehav_path, index=False)

    # Plotting
    fig = bp.plot_fiberpho_behav(
        dfiberbehav_df, list_BOI, exp, mouse,
        THRESH_S, EVENT_TIME_THRESHOLD, batch,
        scaled = False)
    
    fig.savefig(repo_path / f'{batch}_{mouse}_fiberbehav.pdf')
    fig.savefig(repo_path / f'{batch}_{mouse}_fiberbehav.png')
    plt.close(fig)

print(f'\n✅ Analysis for {exp} complete.\nData saved in: {repo_path}')

#%% 2.3 - Plot behavioural metrics

exp = 'RewardAirpuff'

print('###################')
print(f'EXPERIMENT : {exp}')
print('###################')

BIN_SIZE = 1   # seconds
N_TIME_BINS_HEATMAP = 3
HEATMAP_BINS = (50, 50)  # x, y bins

behaviors_to_plot = [
    "Licks_filtered",
    "Nose_in_any_airport",
    "Airpuffs"
]

# Create repository path where data will be stored
exp_path = analysis_path / exp
repo_path = exp_path / f'length0_interbout0_o{ORDER}f{CUT_FREQ}'
behavioural_analysis_path = exp_path / 'Behavioural_analysis'
behavioural_analysis_path.mkdir(exist_ok=True)
data_path_exp = datapath_exp_dict[batch]
behav_path_exp = data_path_exp / 'Behaviour'

all_metrics = {}

# Loop through each mouse in the subject DataFrame
for mouse, batch in zip(subjects_df['Subject'], subjects_df['Batch']):
    print("-----------------------------") 
    print(f'BATCH : {batch}, MOUSE : {mouse}')
    print("-----------------------------")


    fiberbehav_notderived_path = repo_path / f'{batch}_{mouse}_fiberbehavnotderived.csv'
    fiberbehav_path = repo_path / f'{batch}_{mouse}_fiberbehav.csv'
    fiberbehav_notderived_df = pd.read_csv(fiberbehav_notderived_path)
    fiberbehav_df = pd.read_csv(fiberbehav_path)
    arena_json = behav_path_exp / f"{mouse}_arena_coordinates.json"
    ports_json = behav_path_exp / f"{mouse}_ports_coordinates.json"

    # Compute behavioral metrics
    metrics = bm.compute_behavior_metrics(fiberbehav_notderived_df, BIN_SIZE)
    all_metrics[mouse] = metrics

    # Plot behavioral metrics
    mouse_fig_dir = behavioural_analysis_path / 'Figures' / f'batch {batch} mouse {mouse}'
    bm.plot_behavior_metrics(metrics, mouse, batch, BIN_SIZE, save_dir=mouse_fig_dir)

    # Plot raster
    bm.plot_behavior_raster(fiberbehav_notderived_df, 
                            mouse, batch, behaviors=behaviors_to_plot, 
                            save_dir=mouse_fig_dir)

    # Heatmap
    bm.compute_and_plot_heatmap(fiberbehav_notderived_df,
                                mouse, batch,
                                ports_json, arena_json,
                                bins=HEATMAP_BINS, 
                                n_bins=N_TIME_BINS_HEATMAP, 
                                save_dir=mouse_fig_dir)

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

# %%
