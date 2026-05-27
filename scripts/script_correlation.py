#%%IMPORTED
###########

import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import warnings
import importlib
from scipy.ndimage import gaussian_filter1d


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

from scripts.loader import analysis_path, data_path, proto_df, subjects_df, batches, exp_path

#%%
dual_color = True

#filter characteristics
ORDER = 4
CUT_FREQ = None #in Hz
#threshold to fuse behaviour if bouts are too close, in secs
THRESH_S = 3
#threshold for PETH : if events are too short do not plot them and do not include them in PETH, in seconds
EVENT_TIME_THRESHOLD = 0.5

#%% Compute and plot cross-correlation 

# ----------------------------- #
# PETH parameters
exp = 'RewardHab'
BOI = 'Licks_filtered'
baseline = False
MAXBOUTSNUMBER = 20
event = 'onset'

# Plot parameters
TIME_WINDOW = [5, 5]
Y_LIM = [-2,2.5]
Y_LIM_DUAL = [-2,2.5]

# ── PETH by bout number
MIN_MICE_PER_BOUT = 3    # hide bout positions covered by fewer mice
MAX_BOUTS_TO_SHOW = MAXBOUTSNUMBER
STEP = 5

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

repo_path = exp_path / f'length{EVENT_TIME_THRESHOLD}_interbout{THRESH_S}_o{ORDER}f{CUT_FREQ}'
corr_path = repo_path / f'PETH_correlation_{tag}'
corr_path.mkdir(parents=True, exist_ok=True)

# Initialize data storage lists
subject_list = []
group_list = []
PETH_list     = []
PETH_list_560 = []

# Loop over each subject (mouse)
for mouse, batch, group in zip(subjects_df['Subject'], subjects_df['Batch'], subjects_df['Group']):
    print("--------------")
    print(f'MOUSE: {mouse} {batch}')
    print("--------------")

    fiberbehav_file = repo_path / f'{batch}_{mouse}_fiberbehav.csv'

    if not fiberbehav_file.exists():
        print(f"File not found: {fiberbehav_file}")
        continue

    dfiberbehav_df = pd.read_csv(fiberbehav_file, index_col=0)
    if BOI == 'Airpuffs':
        dfiberbehav_clean = bp.remove_first_bout(dfiberbehav_df.reset_index(drop=True), BOI)
    else:
        dfiberbehav_clean = dfiberbehav_df.reset_index(drop=True)

    sr = pp.samplerate(dfiberbehav_clean)

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

        ## Plot joint PETHs

        ## Compute cross-correlation
        lags_s, mean_xcorr, sem_xcorr, peak_lag_s, _ = corr.compute_peth_crosscorr(
            PETH_list, PETH_list_560, sr
        )

        fig_corr = corr.plot_peth_crosscorr(lags_s, mean_xcorr, sem_xcorr, peak_lag_s,
                         BOI, exp, group, MAXBOUTSNUMBER,
                         color='cornflowerblue', fill_alpha=0.25)

        ## Compute deconvolved correlation (matches R-GECO signal to GRAB-ACh dynamics)
        peth_560_deconv_list = [
            np.array([corr.deconvolve_rgeco(row, sr) for row in peth])
            for peth in PETH_list_560
        ]
        lags_s, mean_xcorr, sem_xcorr, peak_lag_s, _ = corr.compute_peth_crosscorr(
            PETH_list, peth_560_deconv_list, sr
        )

        ## Compute Granger causality
        granger_df = corr.test_granger_causality(dfiberbehav_clean, BOI, max_lag_s=3,
                            sr=None, alpha=0.05)
        granger_df.to_excel(corr_path / f'{BOI}_granger.xlsx')
