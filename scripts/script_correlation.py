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
CUT_FREQ = None #in Hz
#threshold to fuse behaviour if bouts are too close, in secs
THRESH_S = 3
#threshold for PETH : if events are too short do not plot them and do not include them in PETH, in seconds
EVENT_TIME_THRESHOLD = 0

#%% Compute and plot cross-correlation 

# ----------------------------- #
# PETH parameters
exp = 'RewardAirpuffs'
BOI = 'Licks_filtered'
baseline = False
MAXBOUTSNUMBER = 30
event = 'onset'

# Plot parameters
TIME_WINDOW = [2, 2]
Y_LIM = [-2,2.5]
Y_LIM_DUAL = [-2,2.5]

# PETH by bout number
MIN_MICE_PER_BOUT = 3
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

    ## Plot joint PETHs
    # Stack all bouts from all mice in this group
    peth_465_all = np.concatenate(PETH_list_group, axis=0)
    peth_560_all = np.concatenate(PETH_list_560_group, axis=0)

    _, _, jpsth_corr, coincidence = corr.compute_joint_psth(peth_465_all, peth_560_all)

    fig_jpsth = corr.plot_joint_psth(
        jpsth_corr, coincidence,
        TIME_WINDOW, BOI, event, exp, group,
        n_bouts=len(peth_465_all)
    )
    fig_jpsth.savefig(corr_path / f'{group}_{BOI}_-{TIME_WINDOW[0]}_{TIME_WINDOW[1]}_JPETH.pdf')
    fig_jpsth.savefig(corr_path / f'{group}_{BOI}_-{TIME_WINDOW[0]}_{TIME_WINDOW[1]}_JPETH.png')
    plt.close(fig_jpsth)

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

    ## Compute cross-correlation significance
    lags_s_sig, mean_xcorr_sig, sem_xcorr_sig, peak_lag_s_sig, ci_low, ci_high, is_sig = corr.compute_crosscorr_significance(
        PETH_list_group, PETH_list_560_group, sr, max_lag_s=TIME_WINDOW[0], n_shuffles=1000, ci=95)
    
    fig_corr_sig = corr.plot_crosscorr_with_significance(lags_s_sig, mean_xcorr_sig, sem_xcorr_sig,
                                      peak_lag_s_sig, ci_low, ci_high, is_sig,
                                      BOI, exp, group, MAXBOUTSNUMBER,
                                      color='cornflowerblue',
                                      sig_style='overlay')
    
    fig_corr_sig.savefig(corr_path / f'{group}_{BOI}_-{TIME_WINDOW[0]}_{TIME_WINDOW[0]}_crosscorrelation_sig.pdf')
    fig_corr_sig.savefig(corr_path / f'{group}_{BOI}_-{TIME_WINDOW[0]}_{TIME_WINDOW[0]}_crosscorrelation_sig.png')
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

# %%
