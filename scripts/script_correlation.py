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
THRESH_S = 2
#threshold for PETH : if events are too short do not plot them and do not include them in PETH, in seconds
EVENT_TIME_THRESHOLD = 0

#%% Compute and plot cross-correlation 
# ----------------------------- #
# PETH parameters
exp = 'EPM'
BOI = 'Open arm'
baseline = False
MAXBOUTSNUMBER = 40
event = 'onset'

# MAX LAG
MAX_LAG_XCORR_S = 1

# Plot parameters
TIME_WINDOW = [2, 2]
Y_LIM = [-2,2.5]
Y_LIM_DUAL = [-2,2.5]
BASELINE_STARTSTOP = [TIME_WINDOW[0],0.5]

# Behaviours to exclude from baseline
behaviours_excluded_baseline_list = ['Open arm','Head dipping','Open arm to Center','Closed arm to Center']

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

    dfiberbehav_df = pd.read_csv(fiberbehav_file)
    if BOI == 'Airpuffs':
        dfiberbehav_clean = bp.remove_first_bout(dfiberbehav_df.reset_index(drop=True), BOI)
    else:
        dfiberbehav_clean = dfiberbehav_df#.reset_index(drop=True)

    sr = pp.samplerate(dfiberbehav_clean)
    dfiberbehav_dict[mouse] = dfiberbehav_clean

    if BOI in dfiberbehav_df.columns[2:].tolist():
        subject_list.append(mouse)
        group_list.append(group)
        print(f'PETH {BOI} for {mouse}')

        ## Get PETHs
        # --- 465 channel ---
        PETH_mouse = bp.PETH(
            dfiberbehav_clean, BOI, event, TIME_WINDOW,
            behaviours_excluded_baseline_list,
            baselinewindow=baseline, maxboutsnumber=MAXBOUTSNUMBER,
            baseline_start_stop_s=BASELINE_STARTSTOP,
            baseline_method='median'
        )

        print(f"PETH shape : {PETH_mouse.shape}")
        PETH_list.append(PETH_mouse)

        # --- 560 channel ---
        PETH_mouse_560 = bp.PETH(
            dfiberbehav_clean, BOI, event, TIME_WINDOW,
            behaviours_excluded_baseline_list,
            baselinewindow=baseline, maxboutsnumber=MAXBOUTSNUMBER, dFF_column='560 dFF',
            baseline_start_stop_s=BASELINE_STARTSTOP,
            baseline_method='median'
        )
        PETH_list_560.append(PETH_mouse_560)

## Compute correlation metrics per group
export_rows = []
for group in included_groups:
    group_indices = [i for i, g in enumerate(group_list) if g == group]
    PETH_list_group   = [PETH_list[i] for i in group_indices]
    PETH_list_560_group = [PETH_list_560[i] for i in group_indices]

    ## Compute cross-correlation
    lags_s, mean_xcorr, sem_xcorr, peak_lag_s, per_mouse_arr, valid_idx = corr.compute_peth_crosscorr(
        PETH_list_group, PETH_list_560_group, sr, max_lag_s=MAX_LAG_XCORR_S
    )
    # --- per-mouse behavior-relative peak & lag, correctly aligned via valid_idx ---
    behavior_peak_by_mouse = {}
    for row_i, orig_i in enumerate(valid_idx):
        mouse = subject_list[group_indices[orig_i]]
        curve = per_mouse_arr[row_i]
        behavior_peak_by_mouse[mouse] = {
            'Behavior_Peak_Xcorr': float(curve.max()),
            'Behavior_Peak_Lag_s': float(lags_s[np.argmax(curve)]),
        }

    # Compute baseline cross-correlation
    baseline_xcorrs = []
    baseline_mice   = []
    lags_bl = None
    for i in group_indices:
        mouse  = subject_list[i]
        result = corr.compute_baseline_crosscorr(
            dfiberbehav_dict[mouse],
            behaviours_excluded_baseline_list,
            sr, pad_s=2, max_lag_s=MAX_LAG_XCORR_S,
            exclusion_col='dFF ExclusionMask'
        )
        if result[0] is not None:
            lags_bl, xcorr_bl, _ = result
            baseline_xcorrs.append(xcorr_bl)
            baseline_mice.append(mouse)
            
    baseline_peak_by_mouse = {
        mouse: {
            'Baseline_Peak_Xcorr': float(curve.max()),
            'Baseline_Peak_Lag_s': float(lags_bl[np.argmax(curve)]),
        }
        for mouse, curve in zip(baseline_mice, baseline_xcorrs)
    }

    mean_xcorr_bl = np.mean(baseline_xcorrs, axis=0) if baseline_xcorrs else None
    sem_xcorr_bl  = (np.std(baseline_xcorrs, axis=0) / np.sqrt(len(baseline_xcorrs))
                     if baseline_xcorrs else None)

    # --- merge behavior + baseline peaks per mouse into export_rows ---
    all_mice_in_group = set(behavior_peak_by_mouse) | set(baseline_peak_by_mouse)
    for mouse in all_mice_in_group:
        row = {'Mouse': mouse, 'Group': group}
        row.update(behavior_peak_by_mouse.get(mouse, {
            'Behavior_Peak_Xcorr': np.nan, 'Behavior_Peak_Lag_s': np.nan}))
        row.update(baseline_peak_by_mouse.get(mouse, {
            'Baseline_Peak_Xcorr': np.nan, 'Baseline_Peak_Lag_s': np.nan}))
        export_rows.append(row)

    for method_corrsig in ['phase_randomization','trial_permutation']:
        ## Compute cross-correlation significance
        mean_xcorr, sem_xcorr, peak_lag_s, ci_low, ci_high, is_sig = \
            corr.compute_crosscorr_significance(
                per_mouse_arr, lags_s,
                PETH_list_group, PETH_list_560_group, sr,
                method=method_corrsig, n_shuffles=1000
            )
        
        fig_corr_sig = corr.plot_crosscorr_with_significance(
            lags_s, mean_xcorr, sem_xcorr,
            peak_lag_s, ci_low, ci_high, is_sig,
            BOI, exp, group, len(group_indices),
            color='cornflowerblue', sig_style='overlay',
            baseline_xcorr=mean_xcorr_bl,
            baseline_sem=sem_xcorr_bl,
        )
        
        fig_corr_sig.savefig(corr_path / f'{group}_{BOI}_-{TIME_WINDOW[0]}_{TIME_WINDOW[0]}_{method_corrsig}_crosscorrelation_sig.pdf')
        fig_corr_sig.savefig(corr_path / f'{group}_{BOI}_-{TIME_WINDOW[0]}_{TIME_WINDOW[0]}_{method_corrsig}_crosscorrelation_sig.png')
        plt.close(fig_corr_sig)

print(f"✔ Cross-correlation plots exported to:{corr_path}")

summary_df = pd.DataFrame(export_rows)[
    ['Mouse', 'Group', 'Behavior_Peak_Xcorr', 'Behavior_Peak_Lag_s',
     'Baseline_Peak_Xcorr', 'Baseline_Peak_Lag_s']
].sort_values(['Group', 'Mouse'])

out_file = corr_path / f'{BOI}_peak_xcorr_summary.xlsx'
summary_df.to_excel(out_file, index=False)
print(f"✔ Peak cross-correlation summary exported to: {out_file}")

 #%% Compute and plot Granger causality test results

# Parameters
############################## 
exp = 'FearRetrieval'
BOI = 'CS+'
event = 'onset'

# Granger causality max lag 
MAX_LAG_S = 0.15
##############################

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

for group in included_groups:
    group_indices = [i for i, g in enumerate(group_list) if g == group]

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

    print(f"✔ Granger causality results and plots exported to:{corr_path}")