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
import json
import importlib
import matplotlib.pyplot as plt
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
import modules.common.quantification as quantif
importlib.reload(quantif)
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
import modules.behaviour.behaviour_metrics as bm
importlib.reload(bm)

from scripts.loader import analysis_path, data_path, proto_df, subjects_df, batches

#%% 2 - ANALYSIS - BEHAVIOUR
############################

# If data was recorded via Bonsai, with a camera_flashes csv file, set to True (else camera flashes were directly recorded in the Doric console)
bonsai_setup = True

#filter characteristics
ORDER = 4
CUT_FREQ = None #in Hz

#threshold to fuse behaviour if bouts are too close, in secs
THRESH_S = 2
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
    video_path = behav_path_exp / f"{mouse}.avi"
    output_json = behav_path_exp / f"{mouse}_epm_coordinates.json"

    #Indicate the emplacement of the lickport and the two airpuff ports
    try:
        if output_json.is_file():
            print('EPM json already exists')
        else:
            print('Get EPM coordinates')
            plt = smb.with_qt5agg() 
            arena = getepm.define_epm_boundaries(video_path)
            getepm.save_boundaries_to_json(arena, output_json)
    except Exception as e:
        print(f"Problem while processing mouse {mouse} : {e}")
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
    scale_json = behav_path_exp / f"{mouse}_scale_and_arena_coordinates.json"
    real_world_distance_cm=76.5
    real_world_distance_name="2 open arms lenght"

    #Indicate the arena boundaries and the coordinates of the known distance
    try:
        if scale_json.is_file():
            print('Scale json already exists')
        else:
            print('Get arena coordinates')
            plt = smb.with_qt5agg()
            scale_and_coords = getvid.get_scale_and_arena_rect(video_path, real_world_distance_cm, real_world_distance_name)
            getvid.save_to_json(scale_and_coords, scale_json)
    except Exception as e:
        print(f"Problem while processing mouse {mouse} : {e}")
plt.close('all') 
plt = smb.with_agg()

#%% 2.2 - Analyze fiberpho data alongside EPM data

for mouse, batch, group in zip(subjects_df['Subject'], subjects_df['Batch'], subjects_df['Group']):
    print("-----------------------------") 
    print(f'BATCH : {batch}, MOUSE : {mouse}')
    print("-----------------------------")

    data_path_exp = datapath_exp_dict[batch]
    pp_path = data_path_exp / 'Preprocessing'
    behav_path_exp = data_path_exp / 'Behaviour'

    # Define paths for raw, behavioral, and fiberphotometry data
    camera_flashes_path = data_path_exp / f'camera_flashes_{mouse}.csv'
    rawdata_path = data_path_exp / f'{mouse}_0000.doric'
    led_flashes_path = data_path_exp / f'miniscope_sync_{mouse}.csv'
    deinterleaved_raw_path = pp_path / f'{mouse}_deinterleaved.csv'
    dlc_path = behav_path_exp / f'{mouse}DLC_Resnet50_EPM_quality50Feb10shuffle1_snapshot_110_filtered.csv'
    fiberpho_path = pp_path / f'{mouse}_dFF_corrected.csv'
    
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
        coordinates_df = mp.get_dlc_data(dlc_path, threshold=0.90)
    except Exception as e:
        print(f'[!] DLC file error for {mouse}: {e}')

    # Time alignment
    try:
        if bonsai_setup:
                led_df = cp.get_timestamps_from_bonsai_csv(led_flashes_path) # gets led flashes from Bonsai files
                ttl_sync_df = cp.extract_sync_channel(rawdata_path, sync_channel = "DIO04")
                deinterleaved_df = pd.read_csv(deinterleaved_raw_path)
                slope, intercept = cp.time_mapping(ttl_sync_df, led_df)

                frame_times_df = cp.get_timestamps_from_bonsai_csv(camera_flashes_path)
                frame_times_df = cp.correct_behav_timestamps(frame_times_df, slope, intercept)
                coordinates_df = cp.align_camera_flashes(coordinates_df, frame_times_df)

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
    except Exception as e:
        print(f'[!] DLC file error for {mouse}: {e}')

    # Analyze EPM position
    behav_df = None
    try:
        if 'EPM' in exp and coordinates_df is not None:
            print(f'Analyzing mouse position for {mouse}')
            list_BOI = ['Open arm', 'Closed arm', 'Center', 'Head dipping']

            behav_df = epm.analyze_mouse_position(coordinates_df, arena_coordinates, arena_scale, bodypart='nose')

        if behav_df is None:
            raise ValueError("Behavior dataframe is missing or invalid.")

        behav_path = behav_path_exp / f'behav_{mouse}.csv'
        behav_df.to_csv(behav_path)
        print(f'Behaviour file exported to {behav_path}')

        # Load fiber photometry data and filter if specified
        fiberpho = pd.read_csv(fiberpho_path)
        if CUT_FREQ is not None:
            fiberpho = cs.lowpass_dFF(fiberpho, ORDER, CUT_FREQ)

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

        # Plotting fiberphotometry data along behavioural classification
        fig = bp.plot_fiberpho_behav(
            dfiberbehav_df, list_BOI, exp, mouse,
            THRESH_S, EVENT_TIME_THRESHOLD, batch,
            scaled = False)
        plt.show()
        fig.savefig(repo_path / f'{batch}_{mouse}_fiberbehav.pdf')
        fig.savefig(repo_path / f'{batch}_{mouse}_fiberbehav.png')
        plt.close('all')

    except Exception as e:
        print(f'[!] DLC file error for {mouse}: {e}')

print(f'\n✅ Analysis for {exp} complete.\nData saved in: {repo_path}')
                    
# %% 2.3 - Averaged heatmaps on all mice, grouped according to specified groups in subject file

###### TO SET ######
bodypart = 'center'
use_zscore = True

# ── Data collection ───────────────────────────────────────────────────────────
subjects_df['Group'] = subjects_df['Group'].fillna('')
included_groups = set(subjects_df['Group'])
subject_list = []
group_list   = []
x_list       = []
y_list       = []
dFF_list     = []

for mouse, batch, group in zip(subjects_df['Subject'], subjects_df['Batch'], subjects_df['Group']):
    print(f"--- {mouse} {batch} {group}---")
    dfiberbehav_file = repo_path / f'{batch}_{mouse}_fiberbehav.csv'
    if not dfiberbehav_file.exists():
        print(f"  File not found, skipping.")
        continue

    dfiberbehav_df = pd.read_csv(dfiberbehav_file, index_col=0)
    subject_list.append(mouse)
    group_list.append(group)
    x_list.append(dfiberbehav_df[f'{bodypart}_x'].values)
    y_list.append(dfiberbehav_df[f'{bodypart}_y'].values)
    dFF_list.append(dfiberbehav_df['dFF'].values)


# ── Grouped heatmap plotting ──────────────────────────────────────────────────
for group in included_groups:
    group_indices = [i for i, g in enumerate(group_list) if g == group]
    if not group_indices:
        print(f"No subjects found for group: {group}")
        continue

    fig = epm.plot_epm_dff_heatmap_grouped(
        x_list         = [x_list[i]   for i in group_indices],
        y_list         = [y_list[i]   for i in group_indices],
        dFF_list       = [dFF_list[i] for i in group_indices],
        subject_list   = [subject_list[i] for i in group_indices],
        epm_coordinates = arena_coordinates,
        group          = group,
        bodypart       = bodypart,
        bins           = (50, 50),
        show_individual= False,   # set False for group average only
        use_zscore = use_zscore,
        save_dir       = repo_path / 'dFF_heatmaps_grouped',
    )
    plt.show()

# %% 2.4 - Quantify dFF in open arm, closed arm and center. Plotting ang getting behavioural data.

subjects_df['Group'] = subjects_df['Group'].fillna('') # if group = Nan, replaces it with an empty string

behav_records = []
dFF_records_raw = []
dFF_records_zscored = []

for mouse, batch, group in zip(subjects_df['Subject'], subjects_df['Batch'], subjects_df['Group']):
    print(f"--- {mouse} {batch} {group} ---")
    behav_file = behav_path_exp / f'behav_{mouse}.csv'
    fiberbehav_file = repo_path / f'{batch}_{mouse}_fiberbehavnotderived.csv'

    if not behav_file.exists():
        print(f"Behav file not found, skipping.")
        continue

    behav_df = pd.read_csv(behav_file, index_col=0)

    # Plotting behavioural data for each mouse
    plt = smb.with_qt5agg()
    epm.plot_epm_behavior(behav_df, arena_coordinates, mouse, batch, bodypart='nose', 
                    n_bins=1, bins=(50, 50), save_dir=repo_path)
    plt.close('all')
    plt = smb.with_agg()

    # Getting behavioural metrics
    record = bm.extract_behav_summary(
        behav_df   = behav_df,
        mouse      = mouse,
        batch      = batch,
        group      = group,
        fps        = arena_scale['Video_fps'],
        behav_cols = ['Closed arm', 'Open arm', 'Center','Head dipping'],
        speed_col  = 'Speed',
        immobility_threshold=0.1
    )
    behav_records.append(record)

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
            zone_cols     = ['Closed arm', 'Open arm', 'Center'],
            behav_cols    = ['Head dipping'],
            # Head dipping also counts toward Open arm
            merge_into    = {'Head dipping': 'Open arm'},
            dff_col       = 'dFF',
            fps           = arena_scale['Video_fps'],
            use_zscore    = use_zscore,
        )
        records.append(record)

# ── Save to Excel ─────────────────────────────────────────────────────────────
pd.DataFrame(behav_records).to_excel(
    repo_path / 'behav_summary.xlsx', index=False)
print(f"Saved {len(behav_records)} mice to behav_summary.xlsx")

pd.DataFrame(dFF_records_raw).to_excel(
    repo_path / 'dFF_summary_raw.xlsx', index=False)
pd.DataFrame(dFF_records_zscored).to_excel(
    repo_path / 'dFF_summary_zscored.xlsx', index=False)
print(f"Saved {len(dFF_records_raw)} mice to dFF_summary.xlsx and dFF_summary_zscored.xlsx")

# %%
