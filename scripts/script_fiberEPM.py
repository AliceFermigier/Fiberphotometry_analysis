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

from scripts.loader import analysis_path, experiment_path, data_path, proto_df, subjects_df, batches

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
EVENT_TIME_THRESHOLD = 0.2

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
    fiberpho_path = pp_path / f'{mouse}_dFF_corrected_final.csv'
    
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
        coordinates_df = mp.get_dlc_data(dlc_path, threshold=0.60)
    except Exception as e:
        print(f'[!] DLC or fiberpho file error for {mouse}: {e}')

    # Time alignment
    try:
        if bonsai_setup:
                # Extract sync channel from Doric raw data and Bonsai corresponding sync data
                ttl_sync_df = cp.extract_sync_channel(rawdata_path, sync_channel = "DIO04")
                led_flashes_path = data_path_exp / f'miniscope_sync_{mouse}.csv'

                _, _, duplicated, missed = cp.diagnose_bonsai_sync(led_flashes_path, ttl_sync_df)

                if duplicated+missed > 0:
                    led_df = cp.get_and_clean_sync_timestamps_from_bonsai_csv(led_flashes_path, ttl_sync_df)
                else:
                    led_df = cp.get_timestamps_from_bonsai_csv(led_flashes_path) 

                # Compute linear regression to correct for differences between clocks
                print(f"Syncing Doric and Bonsai clocks")
                slope, intercept = cp.time_mapping(ttl_sync_df, led_df)
                print("Checking sync quality")
                cp.diagnose_sync(ttl_sync_df, led_df)
 
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
        print(f'[!] Behaviour time alignment error for {mouse}: {e}')

    # Analyze EPM position
    behav_df = None

    if 'EPM' in exp and coordinates_df is not None:
        print(f'Analyzing mouse position for {mouse}')
        list_BOI = ['Open arm', 'Closed arm', 'Closed arm to Center', 'Open arm to Center', 'Head dipping']

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



print(f'\n✅ Analysis for {exp} complete.\nData saved in: {repo_path}')
                    
# %% 2.3 - Averaged heatmaps on all mice, grouped according to specified groups in subject file

###### TO SET ######
bodypart = 'center'
use_zscore = True

# ── Data collection ───────────────────────────────────────────────────────────
subjects_df['Group'] = subjects_df['Group'].fillna('')
included_groups = set(subjects_df['Group'])
subject_list    = []
group_list      = []
x_list          = []
y_list          = []
dFF_list        = []
dFF_560_list    = []
closed_arm_list = []

#Load excluded subjects
excluded_subjects_df = pd.read_excel(experiment_path / 'subjects.xlsx', 
                                     sheet_name=f'Excluded_{exp}')

for mouse, batch, group in zip(subjects_df['Subject'], subjects_df['Batch'], subjects_df['Group']):
    print(f"--- {mouse} {batch} {group} ---")

    if int(mouse) in excluded_subjects_df['Subject'].values:
        print(f"Mouse {mouse} excluded")
        continue

    dfiberbehav_file = repo_path / f'{batch}_{mouse}_fiberbehav.csv'
    if not dfiberbehav_file.exists():
        print(f"File not found, skipping")
        continue

    dfiberbehav_df = pd.read_csv(dfiberbehav_file, index_col=0)
    subject_list.append(mouse)
    group_list.append(group)
    x_list.append(dfiberbehav_df[f'{bodypart}_x'].values)
    y_list.append(dfiberbehav_df[f'{bodypart}_y'].values)
    dFF_list.append(dfiberbehav_df['dFF'].values)
    dFF_560_list.append(                                           
        dfiberbehav_df['560 dFF'].values
        if '560 dFF' in dfiberbehav_df.columns else None) 
    closed_arm_list.append(dfiberbehav_df[f'Closed arm'].values)
    has_dual = any(d is not None for d in dFF_560_list)

print('Colllected data:')
print(f'Subjects:{subject_list}')
print(f'Groups:{group_list}')

# ── Grouped heatmap plotting ──────────────────────────────────────────────────
for group in included_groups:
    group_indices = [i for i, g in enumerate(group_list) if g == group]
    if not group_indices:
        continue

    shared_kwargs = dict(
        epm_coordinates = arena_coordinates,
        group           = group,
        bodypart        = bodypart,
        bins            = (50, 50),
        show_individual = False,
        use_zscore      = use_zscore,
        save_dir        = repo_path / 'dFF_heatmaps_grouped',
    )

    # 465 nm — always
    fig_465 = epm.plot_epm_dff_heatmap_grouped(
        x_list          = [x_list[i]   for i in group_indices],
        y_list          = [y_list[i]   for i in group_indices],
        dFF_list        = [dFF_list[i] for i in group_indices],
        closed_arm_list = [closed_arm_list[i] for i in group_indices],
        subject_list    = [subject_list[i] for i in group_indices],
        signal_name     = '465nm',
        vmin=-3, vmax=3,
        **shared_kwargs,
    )
    plt.show()

    # 560 nm — only for subjects that have it
    if has_dual:
        valid = [(i, dFF_560_list[i]) for i in group_indices if dFF_560_list[i] is not None]
        if valid:
            idx_560, dffs_560 = zip(*valid)
            fig_560 = epm.plot_epm_dff_heatmap_grouped(
                x_list       = [x_list[i]   for i in idx_560],
                y_list       = [y_list[i]   for i in idx_560],
                dFF_list     = list(dffs_560),
                closed_arm_list = [closed_arm_list[i] for i in group_indices],
                subject_list = [subject_list[i] for i in idx_560],
                signal_name  = '560nm',
                vmin=-3, vmax=3,
                **shared_kwargs,
            )
            plt.show()

# %% 2.4 - Quantify dFF in open arm, closed arm and center. Plotting and getting behavioural data.

subjects_df['Group'] = subjects_df['Group'].fillna('') # if group = Nan, replaces it with an empty string

#Load excluded subjects
excluded_subjects_df = pd.read_excel(experiment_path / 'subjects.xlsx', 
                                     sheet_name=f'Excluded_{exp}')

behavioural_analysis_path = repo_path / 'Behavioural_analysis'

behav_records = []
dFF_records_raw = []
dFF_records_zscored = []
dFF_records_raw_560 = []
dFF_records_zscored_560 = []

for mouse, batch, group in zip(subjects_df['Subject'], subjects_df['Batch'], subjects_df['Group']):
    print(f"--- {mouse} {batch} {group} ---")
    behav_file = behav_path_exp / f'behav_{mouse}.csv'
    fiberbehav_file = repo_path / f'{batch}_{mouse}_fiberbehavnotderived.csv'

    if int(mouse) in excluded_subjects_df['Subject'].values:
        print(f"Mouse {mouse} excluded")
        continue

    if not behav_file.exists():
        print(f"Behav file not found, skipping.")
        continue

    behav_df = pd.read_csv(behav_file, index_col=0)

    # Plotting behavioural data for each mouse
    plt = smb.with_qt5agg()
    epm.plot_epm_behavior(behav_df, arena_coordinates, mouse, batch, bodypart='nose', 
                    n_bins=1, bins=(50, 50), save_dir=behavioural_analysis_path)
    plt.close('all')
    plt = smb.with_agg()

    # Getting behavioural metrics
    record = bm.extract_behav_summary(
        behav_df   = behav_df,
        mouse      = mouse,
        batch      = batch,
        group      = group,
        fps        = arena_scale['Video_fps'],
        behav_cols = ['Closed arm', 'Open arm', 'Closed arm to Center','Open arm to Center','Head dipping'],
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
            zone_cols     = ['Closed arm', 'Open arm', 'Closed arm to Center', 'Open arm to Center'],
            behav_cols    = ['Head dipping'],
            baseline_col  = 'Closed arm',
            # Head dipping also counts as Open arm
            merge_into    = {'Head dipping': 'Open arm'},
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
                zone_cols     = ['Closed arm', 'Open arm', 'Closed arm to Center', 'Open arm to Center'],
                behav_cols    = ['Head dipping'],
                baseline_col  = 'Closed arm',
                # Head dipping also counts as Open arm
                merge_into    = {'Head dipping': 'Open arm'},
                dff_col       = '560 dFF',
                fps           = arena_scale['Video_fps'],
                use_zscore    = use_zscore,
            )
            records.append(record)


# ── Save to Excel ─────────────────────────────────────────────────────────────
pd.DataFrame(behav_records).to_excel(
    behavioural_analysis_path / 'behav_summary.xlsx', index=False)
print(f"Saved {len(behav_records)} mice to behav_summary.xlsx")

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

# %% 2.5 Plot grouped EPM heatmap and pie chart

N_TIME_BINS_HEATMAP = 1

repo_path = exp_path / f'length0_interbout0_o{ORDER}f{CUT_FREQ}'
behavioural_analysis_path = repo_path / 'Behavioural_analysis'

data_path_exp = datapath_exp_dict[batch]
behav_path_exp = data_path_exp / 'Behaviour'

# ── Pass 1: collect landmarks and compute shared reference ────────────────────
all_epm_landmarks = {}
all_epm_coordinates = {}
for mouse in subjects_df['Subject']:
    arena_json = behav_path_exp / f"{mouse}_epm_coordinates.json"
    with open(arena_json, 'r') as f:
        arena_coordinates = json.load(f)
    all_epm_landmarks[mouse] = epm.load_epm_landmarks(arena_coordinates)
    all_epm_coordinates[mouse] = arena_coordinates

ref_epm_landmarks = epm.compute_reference_epm(all_epm_landmarks)

# ── Pass 2: align coordinates, collect group data ─────────────────────────────
all_aligned_pos  = []
all_behav_dfs    = []
all_groups_list  = []
ref_epm_coords_list = []  # transformed epm_coordinates per mouse → average for outlines

for mouse, batch, group in zip(subjects_df['Subject'], subjects_df['Batch'], subjects_df['Group']):
    print(f"--- {batch} {mouse} {group} ---")
    behav_file = behav_path_exp / f'behav_{mouse}.csv'
    if int(mouse) in excluded_subjects_df['Subject'].values:
        print(f"Mouse {mouse} excluded")
        continue

    if not behav_file.exists():
        print("  Behav file not found, skipping.")
        continue

    behav_df = pd.read_csv(behav_file, index_col=0)
    mouse_epm_coords = all_epm_coordinates[mouse]  # per-mouse epm_coordinates dict

    M = bm.estimate_port_transform(all_epm_landmarks[mouse], ref_epm_landmarks)

    x_aligned, y_aligned = bm.apply_transform(
        behav_df['nose_x'].values,
        behav_df['nose_y'].values, M)

    all_aligned_pos.append((x_aligned, y_aligned))
    all_behav_dfs.append(behav_df)
    all_groups_list.append(group)
    ref_epm_coords_list.append(epm.get_aligned_epm_coordinates(mouse_epm_coords, M))

# Mean transformed EPM coordinates → use for zone outlines on group heatmap
ref_epm_coords_keys = ref_epm_coords_list[0].keys()
ref_epm_coordinates = {
    k: np.mean([d[k] for d in ref_epm_coords_list])
    for k in ref_epm_coords_keys
}

# ── Pass 3: plot per group ────────────────────────────────────────────────────
for group in subjects_df['Group'].unique():
    group_indices = [i for i, g in enumerate(all_groups_list) if g == group]
    save_dir = behavioural_analysis_path / 'Grouped figures' / f'Group_{group}'

    epm.plot_group_epm_heatmap(
        [all_aligned_pos[i] for i in group_indices],
        ref_epm_coordinates=ref_epm_coordinates,
        bins=(50, 50),
        n_bins=N_TIME_BINS_HEATMAP,
        label=f'Group {group}',
        cmap='jet',
        save_dir=save_dir,
    )

    epm.plot_group_epm_pie(
        [all_behav_dfs[i] for i in group_indices],
        all_groups=[all_groups_list[i] for i in group_indices],
        label=f'Group {group}',
        save_dir=save_dir,
    )
# %%
