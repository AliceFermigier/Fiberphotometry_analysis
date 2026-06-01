#%%IMPORTED
###########

import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import warnings
import importlib
from scipy.ndimage import gaussian_filter1d
from scipy.stats import combine_pvalues


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
import modules.common.quantification as quantif
importlib.reload(quantif)
import modules.common.correlation as corr
importlib.reload(corr)

from scripts.loader import experiment_path, analysis_path, data_path, proto_df, subjects_df, batches

#%%
dual_color = True

#filter characteristics
ORDER = 4
CUT_FREQ = 20 #in Hz
#threshold to fuse behaviour if bouts are too close, in secs
THRESH_S = 3
#threshold for PETH : if events are too short do not plot them and do not include them in PETH, in seconds
EVENT_TIME_THRESHOLD = 0

#%% Compute and plot joint PETHs
# ----------------------------- #
# PETH parameters
exp = 'RewardAirpuff'
BOI = 'Airpuffs'
baseline = False
MAXBOUTSNUMBER = 40
event = 'onset'

# Plot parameters
TIME_WINDOW = [2, 2]
HEATMAP_MINMAX = [-0.5,0.5]
Y_LIM_COINCIDENCE = [-0.2,0.5]
 
# PETH by bout number
MIN_MICE_PER_BOUT = 2
MAX_BOUTS_TO_SHOW = MAXBOUTSNUMBER

# Behaviours to exclude from baseline
behaviours_excluded_baseline_list = ['Airpuffs','Licks_filtered']

if baseline:
    tag = f"windowedbaseline_maxbouts{MAXBOUTSNUMBER}"
else:
    tag = f"wholetrace_maxbouts{MAXBOUTSNUMBER}"

# Set groups
subjects_df['Group'] = subjects_df['Group'].fillna('')
included_groups = set(subjects_df['Group'])
# ----------------------------- #

print('##########################################')
print(f'EXPERIMENT: {exp}')
print('##########################################')

excluded_subjects_df = pd.read_excel(experiment_path / 'subjects.xlsx', 
                                     sheet_name=f'Excluded_{exp}')

exp_path = analysis_path / exp
repo_path = exp_path / f'length{EVENT_TIME_THRESHOLD}_interbout{THRESH_S}_o{ORDER}f{CUT_FREQ}'
corr_path = repo_path / f'PETH_correlation_{tag}'
corr_path.mkdir(parents=True, exist_ok=True)

corr_path_indiv = corr_path / 'Individual plots'
corr_path_indiv.mkdir(parents=True, exist_ok=True)

# Initialize data storage lists
subject_list = []
group_list = []
PETH_list     = []
PETH_list_560 = []
dfiberbehav_dict = {}

# Loop over each subject (mouse)
for mouse, batch, group in zip(subjects_df['Subject'], subjects_df['Batch'], subjects_df['Group']):
    print("--------------")
    print(f'MOUSE: {mouse} {batch}')
    print("--------------")

    fiberbehav_file = repo_path / f'{batch}_{mouse}_fiberbehav.csv'

    if not fiberbehav_file.exists():
        print(f"File not found: {fiberbehav_file}")
        continue
    if int(mouse) in excluded_subjects_df['Subject'].values:
        print(f"Mouse {mouse} excluded")
        continue

    dfiberbehav_df = pd.read_csv(fiberbehav_file, index_col=0)
    if BOI == 'Airpuffs':
        dfiberbehav_clean = bp.remove_first_bout(dfiberbehav_df.reset_index(drop=True), BOI)
    else:
        dfiberbehav_clean = dfiberbehav_df.reset_index(drop=True)

    sr = pp.samplerate(dfiberbehav_clean)
    dfiberbehav_dict[mouse] = dfiberbehav_clean

    if BOI in dfiberbehav_df.columns[2:].tolist():
        subject_list.append(mouse)
        group_list.append(group)
        print(f'PETH {BOI} for {mouse}')

        ## Get PETHs
        # --- 465 channel ---
        PETH_mouse = bp.PETH(
            dfiberbehav_clean, BOI, event, TIME_WINDOW, EVENT_TIME_THRESHOLD,
            baselinewindow=baseline, maxboutsnumber=MAXBOUTSNUMBER
        )
        print(f"PETH shape : {PETH_mouse.shape}")
        PETH_list.append(PETH_mouse)

        # --- 560 channel ---
        PETH_mouse_560 = bp.PETH(
            dfiberbehav_clean, BOI, event, TIME_WINDOW, EVENT_TIME_THRESHOLD,
            baselinewindow=baseline, maxboutsnumber=MAXBOUTSNUMBER, dFF_column='560 dFF'
        )
        PETH_list_560.append(PETH_mouse_560)

all_coinc_records = []   # collects one row per mouse across all groups

for group in included_groups:
    group_indices = [i for i, g in enumerate(group_list) if g == group]

    per_mouse_jpsth = []
    per_mouse_coinc = []

    per_mouse_jpsth_raw = []
    per_mouse_coinc_raw = []

    per_mouse_jpsth_corrected = []
    per_mouse_coinc_corrected = []

    per_mouse_jpsth_zscore = []
    per_mouse_coinc_zscore = []

    # ── Per-mouse JPSTH ───────────────────────────────────────────────────────
    for i in group_indices:
        mouse          = subject_list[i]
        peth_465_mouse = PETH_list[i]
        peth_560_mouse = PETH_list_560[i]
        n_bouts_mouse  = len(peth_465_mouse)

        # Per-mouse baseline JPSTH
        jpsth_bl, coinc_bl = corr.compute_baseline_jpsth(
            dfiberbehav_dict[mouse],
            behaviours_excluded_baseline_list,
            sr, TIME_WINDOW,
            pad_s=2,
            exclusion_col='dFF ExclusionMask'
        )

        # Per-mouse JPSTH
        jpsth_m_raw, predictor_m, jpsth_m, coinc_m, coinc_m_raw, predictor_diag_m = corr.compute_joint_psth(peth_465_mouse, peth_560_mouse)
       
        per_mouse_jpsth.append(jpsth_m)
        per_mouse_coinc.append(coinc_m)

        per_mouse_jpsth_raw.append(jpsth_m_raw)
        per_mouse_coinc_raw.append(coinc_m_raw)

        # JPSTH correction with baseline
        jpsth_corrected_m = jpsth_m  - jpsth_bl
        coinc_corrected_m = coinc_m  - coinc_bl

        per_mouse_jpsth_corrected.append(jpsth_corrected_m)
        per_mouse_coinc_corrected.append(coinc_corrected_m)

        # Per-mouse figure
        fig_m = corr.plot_joint_psth(
            jpsth_m, coinc_m, TIME_WINDOW, BOI, event, exp, group,
            n_bouts=n_bouts_mouse, mouse=mouse 
        )
        fig_m.savefig(corr_path_indiv / f'{group}_{mouse}_{BOI}_-{TIME_WINDOW[0]}_{TIME_WINDOW[1]}_JPETH.pdf')
        fig_m.savefig(corr_path_indiv / f'{group}_{mouse}_{BOI}_-{TIME_WINDOW[0]}_{TIME_WINDOW[1]}_JPETH.png')
        plt.close(fig_m)

        # Per-mouse figure raw
        fig_m_raw = corr.plot_joint_psth(
            jpsth_m_raw, coinc_m_raw, TIME_WINDOW, BOI, event, exp, group,
            n_bouts=n_bouts_mouse, mouse=mouse 
        )
        fig_m_raw.savefig(corr_path_indiv / f'{group}_{mouse}_{BOI}_-{TIME_WINDOW[0]}_{TIME_WINDOW[1]}_JPETH_raw.pdf')
        fig_m_raw.savefig(corr_path_indiv / f'{group}_{mouse}_{BOI}_-{TIME_WINDOW[0]}_{TIME_WINDOW[1]}_JPETH_raw.png')
        plt.close(fig_m_raw)

        # Per-mouse predictor
        fig_m_predictor = corr.plot_joint_psth(
            predictor_m, predictor_diag_m, TIME_WINDOW, BOI, event, exp, group,
            n_bouts=n_bouts_mouse, mouse=mouse 
        )
        fig_m_predictor.savefig(corr_path_indiv / f'{group}_{mouse}_{BOI}_-{TIME_WINDOW[0]}_{TIME_WINDOW[1]}_JPETH_predictor.pdf')
        fig_m_predictor.savefig(corr_path_indiv / f'{group}_{mouse}_{BOI}_-{TIME_WINDOW[0]}_{TIME_WINDOW[1]}_JPETH_predictor.png')
        plt.close(fig_m_predictor)

        # Per-mouse figure baseline
        fig_m_bl = corr.plot_joint_psth(
            jpsth_bl, coinc_bl, TIME_WINDOW, BOI, event, exp, group,
            n_bouts=n_bouts_mouse, mouse=mouse 
        )
        fig_m_bl.savefig(corr_path_indiv / f'{group}_{mouse}_{BOI}_-{TIME_WINDOW[0]}_{TIME_WINDOW[1]}_JPETH_bl.pdf')
        fig_m_bl.savefig(corr_path_indiv / f'{group}_{mouse}_{BOI}_-{TIME_WINDOW[0]}_{TIME_WINDOW[1]}_JPETH_bl.png')
        plt.close(fig_m_bl)

        # Per-mouse figure corrected
        fig_m_corrected = corr.plot_joint_psth(
            jpsth_corrected_m, coinc_corrected_m, TIME_WINDOW, BOI, event, exp, group,
            n_bouts=n_bouts_mouse, mouse=mouse 
        )
        fig_m_corrected.savefig(corr_path_indiv / f'{group}_{mouse}_{BOI}_-{TIME_WINDOW[0]}_{TIME_WINDOW[1]}_JPETH_corrected.pdf')
        fig_m_corrected.savefig(corr_path_indiv / f'{group}_{mouse}_{BOI}_-{TIME_WINDOW[0]}_{TIME_WINDOW[1]}_JPETH_corrected.png')
        plt.close(fig_m_corrected)

        jpsth_bl, coinc_bl

        # ── Shuffle null and z-score correction ───────────────────────────────
        jpsth_shuf_mean_560, jpsth_shuf_std_560 = corr.compute_shuffle_jpsth(
            dfiberbehav_dict[mouse],         
            peth_465_mouse,
            BOI, event, TIME_WINDOW, EVENT_TIME_THRESHOLD, sr,
            n_shuffles=100,
            baseline=baseline,
            sig_to_shuffle = '560 dFF',
            maxboutsnumber=MAXBOUTSNUMBER,
        )
        coinc_shuf_mean_560 = np.diag(jpsth_shuf_mean_560)

        jpsth_shuf_mean_465, jpsth_shuf_std_465 = corr.compute_shuffle_jpsth(
            dfiberbehav_dict[mouse],         
            peth_465_mouse,
            BOI, event, TIME_WINDOW, EVENT_TIME_THRESHOLD, sr,
            n_shuffles=100,
            baseline=baseline,
            sig_to_shuffle = 'dFF',
            maxboutsnumber=MAXBOUTSNUMBER,
        )
        coinc_shuf_mean_465 = np.diag(jpsth_shuf_mean_465)

        # Per-mouse figure shuffle
        fig_shuf_m_560 = corr.plot_joint_psth(
            jpsth_shuf_mean_560, coinc_shuf_mean_560, TIME_WINDOW, BOI, event, exp, group,
            n_bouts=n_bouts_mouse, mouse=mouse 
        )
        fig_shuf_m_560.savefig(corr_path_indiv / f'{group}_{mouse}_{BOI}_-{TIME_WINDOW[0]}_{TIME_WINDOW[1]}_JPETH_shuf_560.pdf')
        fig_shuf_m_560.savefig(corr_path_indiv / f'{group}_{mouse}_{BOI}_-{TIME_WINDOW[0]}_{TIME_WINDOW[1]}_JPETH_shuf_560.png')
        plt.close(fig_shuf_m_560)

        fig_shuf_m_465 = corr.plot_joint_psth(
            jpsth_shuf_mean_465, coinc_shuf_mean_465, TIME_WINDOW, BOI, event, exp, group,
            n_bouts=n_bouts_mouse, mouse=mouse 
        )
        fig_shuf_m_465.savefig(corr_path_indiv / f'{group}_{mouse}_{BOI}_-{TIME_WINDOW[0]}_{TIME_WINDOW[1]}_JPETH_shuf_465.pdf')
        fig_shuf_m_465.savefig(corr_path_indiv / f'{group}_{mouse}_{BOI}_-{TIME_WINDOW[0]}_{TIME_WINDOW[1]}_JPETH_shuf_465.png')
        plt.close(fig_shuf_m_465)

        # Z-score: how many SDs above the circular-shift null is each cell?
        jpsth_z_m = (jpsth_m_raw - (jpsth_shuf_mean_560+jpsth_shuf_mean_465)) / (jpsth_shuf_std_560+jpsth_shuf_std_465)
        coinc_z_m = np.diag(jpsth_z_m)

        per_mouse_jpsth_zscore.append(jpsth_z_m)
        per_mouse_coinc_zscore.append(coinc_z_m)

        # Per-mouse z-score figure
        fig_m_z = corr.plot_joint_psth(
            jpsth_z_m, coinc_z_m, TIME_WINDOW, BOI, event, exp, group,
            n_bouts=n_bouts_mouse, mouse=mouse,
            cmap='RdBu_r'     
        )
        fig_m_z.savefig(corr_path_indiv / f'{group}_{mouse}_{BOI}_-{TIME_WINDOW[0]}_{TIME_WINDOW[1]}_JPETH_zscore.pdf')
        fig_m_z.savefig(corr_path_indiv / f'{group}_{mouse}_{BOI}_-{TIME_WINDOW[0]}_{TIME_WINDOW[1]}_JPETH_zscore.png')
        plt.close(fig_m_z)

        # Use z-scored coincidence for metrics (more interpretable than baseline-subtracted)
        metrics = corr.extract_coincidence_metrics(coinc_z_m, TIME_WINDOW, step_s=0.1)
        all_coinc_records.append({
            'Mouse'   : mouse,
            'Group'   : group,
            'n_bouts' : n_bouts_mouse,
            **metrics,
        })

    # ── Group average ───────────────────────────────
    jpsth_stack   = np.stack(per_mouse_jpsth)          # (n_mice, n_tp, n_tp)
    coinc_stack   = np.stack(per_mouse_coinc)          # (n_mice, n_tp)
    jpsth_group   = np.nanmean(jpsth_stack, axis=0)
    coinc_group   = np.nanmean(coinc_stack, axis=0)
    coinc_sem     = np.nanstd(coinc_stack, axis=0) / np.sqrt(len(group_indices))
    n_bouts_group = len(group_indices)

    # Group JPSTH figure
    fig_g = corr.plot_joint_psth(
        jpsth_group, coinc_group, TIME_WINDOW, BOI, event, exp, group,
        n_bouts=n_bouts_group, coincidence_sem=coinc_sem
    )
    fig_g.savefig(corr_path / f'{group}_{BOI}_-{TIME_WINDOW[0]}_{TIME_WINDOW[1]}_JPETH.pdf')
    fig_g.savefig(corr_path / f'{group}_{BOI}_-{TIME_WINDOW[0]}_{TIME_WINDOW[1]}_JPETH.png')
    plt.close(fig_g)

    # Standalone coincidence figure
    fig_coinc = corr.plot_coincidence(
        coinc_group, TIME_WINDOW, BOI, event, exp, group,
        n_bouts=n_bouts_group, coincidence_sem=coinc_sem
    )
    fig_coinc.savefig(corr_path / f'{group}_{BOI}_-{TIME_WINDOW[0]}_{TIME_WINDOW[1]}_coincidence.pdf')
    fig_coinc.savefig(corr_path / f'{group}_{BOI}_-{TIME_WINDOW[0]}_{TIME_WINDOW[1]}_coincidence.png')
    plt.close(fig_coinc)

    # ── Group average raw ───────────────────────────────
    jpsth_stack_raw   = np.stack(per_mouse_jpsth_raw)          # (n_mice, n_tp, n_tp)
    coinc_stack_raw   = np.stack(per_mouse_coinc_raw)          # (n_mice, n_tp)
    jpsth_group_raw   = np.nanmean(jpsth_stack_raw, axis=0)
    coinc_group_raw   = np.nanmean(coinc_stack_raw, axis=0)
    coinc_sem_raw     = np.nanstd(coinc_stack_raw, axis=0) / np.sqrt(len(group_indices))

    # Group JPSTH figure
    fig_g_raw = corr.plot_joint_psth(
        jpsth_group_raw, coinc_group_raw, TIME_WINDOW, BOI, event, exp, group,
        n_bouts=n_bouts_group, coincidence_sem=coinc_sem_raw
    )
    fig_g_raw.savefig(corr_path / f'{group}_{BOI}_-{TIME_WINDOW[0]}_{TIME_WINDOW[1]}_JPETH_raw.pdf')
    fig_g_raw.savefig(corr_path / f'{group}_{BOI}_-{TIME_WINDOW[0]}_{TIME_WINDOW[1]}_JPETH_raw.png')
    plt.close(fig_g_raw)

    # Standalone coincidence figure
    fig_coinc_raw = corr.plot_coincidence(
        coinc_group_raw, TIME_WINDOW, BOI, event, exp, group,
        n_bouts=n_bouts_group, coincidence_sem=coinc_sem_raw
    )
    fig_coinc_raw.savefig(corr_path / f'{group}_{BOI}_-{TIME_WINDOW[0]}_{TIME_WINDOW[1]}_coincidence_raw.pdf')
    fig_coinc_raw.savefig(corr_path / f'{group}_{BOI}_-{TIME_WINDOW[0]}_{TIME_WINDOW[1]}_coincidence_raw.png')
    plt.close(fig_coinc_raw)

    # ── Group average baseline-corrected ───────────────────────────────
    jpsth_stack_corrected   = np.stack(per_mouse_jpsth_corrected)          # (n_mice, n_tp, n_tp)
    coinc_stack_corrected   = np.stack(per_mouse_coinc_corrected)          # (n_mice, n_tp)
    jpsth_group_corrected   = np.nanmean(jpsth_stack_corrected, axis=0)
    coinc_group_corrected   = np.nanmean(coinc_stack_corrected, axis=0)
    coinc_sem_corrected     = np.nanstd(coinc_stack_corrected, axis=0) / np.sqrt(len(group_indices))

    # Group JPSTH figure
    fig_g_corrected = corr.plot_joint_psth(
        jpsth_group_corrected, coinc_group_corrected, TIME_WINDOW, BOI, event, exp, group,
        n_bouts=n_bouts_group, vmin=HEATMAP_MINMAX[0], vmax=HEATMAP_MINMAX[1], coincidence_sem=coinc_sem_corrected
    )
    fig_g_corrected.savefig(corr_path / f'{group}_{BOI}_-{TIME_WINDOW[0]}_{TIME_WINDOW[1]}_JPETH_corrected.pdf')
    fig_g_corrected.savefig(corr_path / f'{group}_{BOI}_-{TIME_WINDOW[0]}_{TIME_WINDOW[1]}_JPETH_corrected.png')
    plt.close(fig_g_corrected)

    # Standalone coincidence figure
    fig_coinc_corrected = corr.plot_coincidence(
        coinc_group_corrected, TIME_WINDOW, BOI, event, exp, group,
        n_bouts=n_bouts_group, ylim=Y_LIM_COINCIDENCE, coincidence_sem=coinc_sem_corrected
    )
    fig_coinc_corrected.savefig(corr_path / f'{group}_{BOI}_-{TIME_WINDOW[0]}_{TIME_WINDOW[1]}_coincidence_corrected.pdf')
    fig_coinc_corrected.savefig(corr_path / f'{group}_{BOI}_-{TIME_WINDOW[0]}_{TIME_WINDOW[1]}_coincidence_corrected.png')
    plt.close(fig_coinc_corrected)

    # ── Group average z-scored ────────────────────────────────────────────────
    jpsth_stack_z = np.stack(per_mouse_jpsth_zscore)
    coinc_stack_z = np.stack(per_mouse_coinc_zscore)
    jpsth_group_z = np.nanmean(jpsth_stack_z, axis=0)
    coinc_group_z = np.nanmean(coinc_stack_z, axis=0)
    coinc_sem_z   = np.nanstd(coinc_stack_z,  axis=0) / np.sqrt(len(group_indices))

    fig_g_z = corr.plot_joint_psth(
        jpsth_group_z, coinc_group_z, TIME_WINDOW, BOI, event, exp, group,
        n_bouts=n_bouts_group, vmin=-4, vmax=4,
        coincidence_sem=coinc_sem_z
    )
    fig_g_z.savefig(corr_path / f'{group}_{BOI}_-{TIME_WINDOW[0]}_{TIME_WINDOW[1]}_JPETH_zscore.pdf')
    fig_g_z.savefig(corr_path / f'{group}_{BOI}_-{TIME_WINDOW[0]}_{TIME_WINDOW[1]}_JPETH_zscore.png')
    plt.close(fig_g_z)

    fig_coinc_z = corr.plot_coincidence(
        coinc_group_z, TIME_WINDOW, BOI, event, exp, group,
        n_bouts=n_bouts_group, coincidence_sem=coinc_sem_z,
        ylim=[-0.5,5], zscore=True
    )
    fig_coinc_z.savefig(corr_path / f'{group}_{BOI}_-{TIME_WINDOW[0]}_{TIME_WINDOW[1]}_coincidence_zscore.pdf')
    fig_coinc_z.savefig(corr_path / f'{group}_{BOI}_-{TIME_WINDOW[0]}_{TIME_WINDOW[1]}_coincidence_zscore.png')
    plt.close(fig_coinc_z)

# ── Export all per-mouse metrics ──────────────────────────────────────────────
pd.DataFrame(all_coinc_records).to_excel(
    corr_path / f'{BOI}_coincidence_metrics.xlsx', index=False
)
print(f"✔ Coincidence metrics exported to: {corr_path}")

 #%% Compute and plot cross-correlation 
# ----------------------------- #
# PETH parameters
exp = 'FearRetrieval'
BOI = 'CS+'
baseline = False
MAXBOUTSNUMBER = None
event = 'onset'

# Plot parameters
TIME_WINDOW = [2, 2]
Y_LIM = [-2,2.5]
Y_LIM_DUAL = [-2,2.5]

# PETH by bout number
MIN_MICE_PER_BOUT = 2
MAX_BOUTS_TO_SHOW = MAXBOUTSNUMBER

# Granger causality max lag 
MAX_LAG_S = 0.15

if baseline:
    tag = f"windowedbaseline_maxbouts{MAXBOUTSNUMBER}"
else:
    tag = f"wholetrace_maxbouts{MAXBOUTSNUMBER}"

# Set groups
subjects_df['Group'] = subjects_df['Group'].fillna('')
included_groups = set(subjects_df['Group'])
# ----------------------------- #

print('##########################################')
print(f'EXPERIMENT: {exp}')
print('##########################################')

excluded_subjects_df = pd.read_excel(experiment_path / 'subjects.xlsx', 
                                     sheet_name=f'Excluded_{exp}')

exp_path = analysis_path / exp
repo_path = exp_path / f'length{EVENT_TIME_THRESHOLD}_interbout{THRESH_S}_o{ORDER}f{CUT_FREQ}'
corr_path = repo_path / f'PETH_correlation_{tag}'
corr_path.mkdir(parents=True, exist_ok=True)

# Initialize data storage lists
subject_list = []
group_list = []
PETH_list     = []
PETH_list_560 = []
dfiberbehav_dict = {}

# Loop over each subject (mouse)
for mouse, batch, group in zip(subjects_df['Subject'], subjects_df['Batch'], subjects_df['Group']):
    print("--------------")
    print(f'MOUSE: {mouse} {batch}')
    print("--------------")

    fiberbehav_file = repo_path / f'{batch}_{mouse}_fiberbehav.csv'

    if not fiberbehav_file.exists():
        print(f"File not found: {fiberbehav_file}")
        continue
    if int(mouse) in excluded_subjects_df['Subject'].values:
        print(f"Mouse {mouse} excluded")
        continue

    dfiberbehav_df = pd.read_csv(fiberbehav_file, index_col=0)
    if BOI == 'Airpuffs':
        dfiberbehav_clean = bp.remove_first_bout(dfiberbehav_df.reset_index(drop=True), BOI)
    else:
        dfiberbehav_clean = dfiberbehav_df.reset_index(drop=True)

    sr = pp.samplerate(dfiberbehav_clean)
    dfiberbehav_dict[mouse] = dfiberbehav_clean

    if BOI in dfiberbehav_df.columns[2:].tolist():
        subject_list.append(mouse)
        group_list.append(group)
        print(f'PETH {BOI} for {mouse}')

        ## Get PETHs
        # --- 465 channel ---
        PETH_mouse = bp.PETH(
            dfiberbehav_clean, BOI, event, TIME_WINDOW, EVENT_TIME_THRESHOLD,
            baselinewindow=baseline, maxboutsnumber=MAXBOUTSNUMBER
        )
        print(f"PETH shape : {PETH_mouse.shape}")
        PETH_list.append(PETH_mouse)

        # --- 560 channel ---
        PETH_mouse_560 = bp.PETH(
            dfiberbehav_clean, BOI, event, TIME_WINDOW, EVENT_TIME_THRESHOLD,
            baselinewindow=baseline, maxboutsnumber=MAXBOUTSNUMBER, dFF_column='560 dFF'
        )
        PETH_list_560.append(PETH_mouse_560)

## Compute correlation metrics per group
for group in included_groups:
    group_indices = [i for i, g in enumerate(group_list) if g == group]
    PETH_list_group   = [PETH_list[i] for i in group_indices]
    PETH_list_560_group = [PETH_list_560[i] for i in group_indices]

    ## Compute cross-correlation
    lags_s, mean_xcorr, sem_xcorr, peak_lag_s, _ = corr.compute_peth_crosscorr(
            PETH_list_group, PETH_list_560_group, sr
        )

    fig_corr = corr.plot_peth_crosscorr(lags_s, mean_xcorr, sem_xcorr, peak_lag_s,
                    BOI, exp, group, MAXBOUTSNUMBER,
                    color='cornflowerblue', fill_alpha=0.25)
    
    fig_corr.savefig(corr_path / f'{group}_{BOI}_-{TIME_WINDOW[0]}_{TIME_WINDOW[1]}_crosscorrelation.pdf')
    fig_corr.savefig(corr_path / f'{group}_{BOI}_-{TIME_WINDOW[0]}_{TIME_WINDOW[1]}_crosscorrelation.png')
    plt.close(fig_corr)

    for method_corrsig in ['phase_randomization','trial_permutation']:
        ## Compute cross-correlation significance
        lags_s_sig, mean_xcorr_sig, sem_xcorr_sig, peak_lag_s_sig, ci_low, ci_high, is_sig = corr.compute_crosscorr_significance(
            PETH_list_group, PETH_list_560_group, sr, method=method_corrsig, max_lag_s=TIME_WINDOW[0], n_shuffles=1000, ci=95)
        
        fig_corr_sig = corr.plot_crosscorr_with_significance(lags_s_sig, mean_xcorr_sig, sem_xcorr_sig,
                                        peak_lag_s_sig, ci_low, ci_high, is_sig,
                                        BOI, exp, group, MAXBOUTSNUMBER,
                                        color='cornflowerblue',
                                        sig_style='overlay')
        
        fig_corr_sig.savefig(corr_path / f'{group}_{BOI}_-{TIME_WINDOW[0]}_{TIME_WINDOW[0]}_{method_corrsig}_crosscorrelation_sig.pdf')
        fig_corr_sig.savefig(corr_path / f'{group}_{BOI}_-{TIME_WINDOW[0]}_{TIME_WINDOW[0]}_{method_corrsig}_crosscorrelation_sig.png')
        plt.close(fig_corr_sig)
    
    ## Compute deconvolved correlation (matches R-GECO signal to GRAB-ACh dynamics)
    peth_560_deconv_list = [
        np.array([corr.deconvolve_rgeco(row, sr) for row in peth])
        for peth in PETH_list_560_group
    ]
    lags_s_deconvolved, mean_xcorr_deconvolved, sem_xcorr_deconvolved, peak_lag_s_deconvolved, _ = corr.compute_peth_crosscorr(
        PETH_list_group, peth_560_deconv_list, sr
    )

    fig_corr_deconvolved = corr.plot_peth_crosscorr(lags_s_deconvolved, mean_xcorr_deconvolved, sem_xcorr_deconvolved, peak_lag_s_deconvolved,
                    BOI, exp, group, MAXBOUTSNUMBER,
                    color='cornflowerblue', fill_alpha=0.25)
    
    fig_corr_deconvolved.savefig(corr_path / f'{group}_{BOI}_-{TIME_WINDOW[0]}_{TIME_WINDOW[1]}_crosscorrelation_deconvolved.pdf')
    fig_corr_deconvolved.savefig(corr_path / f'{group}_{BOI}_-{TIME_WINDOW[0]}_{TIME_WINDOW[1]}_crosscorrelation_deconvolved.png')
    plt.close(fig_corr_deconvolved)

    # ── Granger causality: pool across all mice in this group ─────────────────
    all_lag_results  = {}   # lag → {'F': [...], 'p': [...]}
    total_accepted   = 0
    total_rejected   = 0

    for i in group_indices:
        mouse    = subject_list[i]
        df_mouse = dfiberbehav_dict.get(mouse)
        if df_mouse is None:
            print(f"  [!] No dataframe found for {mouse}, skipping Granger.")
            continue

        _, raw, acc, rej = corr.test_granger_causality(
            df_mouse, BOI,
            max_lag_s          = MAX_LAG_S,
            sr                 = None,
            alpha              = 0.05,
            return_raw         = True,
        )
        total_accepted += acc
        total_rejected += rej

        for lag, data in raw.items():
            if lag not in all_lag_results:
                all_lag_results[lag] = {'F': [], 'p': []}
            all_lag_results[lag]['F'].extend(data['F'])
            all_lag_results[lag]['p'].extend(data['p'])

    # ── Combine across mice + bouts with Fisher's method ─────────────────────
    sr_granger = round(pp.samplerate(next(iter(dfiberbehav_dict.values()))))
    group_output = []

    for lag in sorted(all_lag_results.keys()):
        F_vals = all_lag_results[lag]['F']
        p_vals = all_lag_results[lag]['p']
        if not F_vals:
            continue
        combined_p = combine_pvalues(p_vals, method='fisher')[1]
        group_output.append({
            'lag_samples'    : lag,
            'lag_s'          : lag / sr_granger,
            'mean_F'         : np.mean(F_vals),
            'std_F'          : np.std(F_vals),
            'combined_p'     : combined_p,
            'significant'    : combined_p < 0.05,
            'n_bouts_total'  : len(F_vals),
            'accepted_bouts' : total_accepted,
            'rejected_bouts' : total_rejected,
        })

    granger_df = pd.DataFrame(group_output)

    granger_df.to_excel(corr_path / f'{group}_{BOI}_maxlag{MAX_LAG_S}s_granger.xlsx',
                        index=False)
    fig_granger = corr.plot_granger_results(granger_df, BOI, exp, group, alpha=0.05)
    fig_granger.savefig(corr_path / f'{group}_{BOI}_maxlag{MAX_LAG_S}s_granger.pdf')
    fig_granger.savefig(corr_path / f'{group}_{BOI}_maxlag{MAX_LAG_S}s_granger.png')
    plt.close(fig_granger)

    print(f"✔ Correlation and causality results and plots exported to:{corr_path}")