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

from scripts.loader import analysis_path, data_path, proto_df, subjects_df, batches

#%%

dual_color = True

#filter characteristics
ORDER = 4
CUT_FREQ = None #in Hz
#threshold to fuse behaviour if bouts are too close, in secs
THRESH_S = 5
#threshold for PETH : if events are too short do not plot them and do not include them in PETH, in seconds
EVENT_TIME_THRESHOLD = 0

#%% Plot PETH for each mouse

# PETH parameters 
baseline = False # parameter to know how the z-score in calculated (mean and sd on short timewindow before event or wholetrace)
MAXBOUTSNUMBER = None
if baseline:
    tag = "windowedbaseline"
else:
    tag = "wholetrace"

# Plot parameters
EVENT_LIST = ['onset','withdrawal']
TIME_WINDOWS = [[3, 8],[3, 8]]  # Time window for PETH calculation (pre, post), for each event
Y_LIM = [-2,10]
Y_LIM_DUAL = [-2,5]
behaviors_of_interest = ['Licks_filtered']

#['Licks_filtered','Airpuffs']
#['Licks_filtered']
#['Open arm','Closed arm','Head dipping','Center']

for exp in ['RewardHab']: #[f.name for f in analysis_path.iterdir() if f.is_dir()]:
    exp_path = analysis_path / exp
    datapath_exp_dict = nom.get_experiment_data_path(batches, proto_df, data_path, exp)

    # Loop over each session folder in the experiment path
    print('##########################################')
    print(f'EXPERIMENT : {exp}')
    print('##########################################')

    # Create the repository and PETH paths
    repo_path = exp_path / f'length{EVENT_TIME_THRESHOLD}_interbout{THRESH_S}_o{ORDER}f{CUT_FREQ}'
    peth_path = repo_path / f'PETH_{tag}'
    peth_path.mkdir(parents=True, exist_ok=True)  # Create directory if it doesn't exist

    # Loop over each mouse in the subjects DataFrame
    for mouse, batch, group in zip(subjects_df['Subject'], subjects_df['Batch'], subjects_df['Group']):
        fiberbehav_path = repo_path / f'{batch}_{mouse}_fiberbehav.csv'

        if fiberbehav_path.exists():  # Check if fiber behavior file exists for this mouse
            print("--------------")
            print(f'MOUSE : {mouse} {batch}')
            print("--------------")

            try:
                # Read the fiber behavior file
                dfiberbehav_df = pd.read_csv(fiberbehav_path, index_col=0)
            except Exception as e:
                warnings.warn(f"Failed to read file {fiberbehav_path}: {e}")
                continue
            
            for behavior in behaviors_of_interest:

                if behavior == 'Airpuffs':
                    dfiberbehav_clean = bp.remove_first_bout(dfiberbehav_df.reset_index(drop=True), behavior)
                else:
                    dfiberbehav_clean = dfiberbehav_df.reset_index(drop=True)

                for event, time_window in zip(EVENT_LIST, TIME_WINDOWS):  
                    # Generate the PETH data for the current behavior, event, and time window
                    print(f"Getting PETH data for {behavior} {event} 465nm")
                    peth_data = bp.PETH(dfiberbehav_clean, behavior, event, time_window, 
                                        EVENT_TIME_THRESHOLD, baselinewindow = baseline, 
                                        maxboutsnumber=MAXBOUTSNUMBER)
                    
                    # Create a DataFrame from the PETH data
                    sr = round(pp.samplerate(dfiberbehav_df))
                    PRE_TIME, POST_TIME = time_window
                    n_timepoints = (PRE_TIME + POST_TIME) * sr + 1
                    time_index = np.linspace(-PRE_TIME, POST_TIME, n_timepoints)

                    peth_df = pd.DataFrame(np.transpose(peth_data), index=time_index)
                    
                    # Plot the PETH and save the figure 
                    print(f"Plotting PETH 465nm")
                    peth_plot = bp.plot_PETH(peth_data, behavior, event, time_window, exp, batch, mouse, group, ylim=Y_LIM)
                    peth_plot.savefig(peth_path / f'{batch}_{mouse}_{behavior}_465_{event[0]}{time_window[0] - time_window[1]}_PETH.png')
                    peth_plot.savefig(peth_path / f'{batch}_{mouse}_{behavior}_465_{event[0]}{time_window[0] - time_window[1]}_PETH.pdf')
                    plt.close(peth_plot)
            
                    if dual_color:
                        # Generate the PETH data for the current behavior, event, and time window
                        print(f"Getting PETH data for {behavior} {event} 560nm")
                        peth_data = bp.PETH(dfiberbehav_clean, behavior, event, time_window, 
                                            EVENT_TIME_THRESHOLD, baselinewindow = baseline, 
                                            maxboutsnumber=MAXBOUTSNUMBER, dFF_column = '560 dFF')
                        
                        # Create a DataFrame from the PETH data
                        sr = round(pp.samplerate(dfiberbehav_df))
                        PRE_TIME, POST_TIME = time_window
                        n_timepoints = (PRE_TIME + POST_TIME) * sr + 1
                        time_index = np.linspace(-PRE_TIME, POST_TIME, n_timepoints)

                        peth_df = pd.DataFrame(np.transpose(peth_data), index=time_index)
                        
                        # Plot the PETH and save the figure
                        print(f"Plotting PETH 560nm")
                        peth_plot = bp.plot_PETH(peth_data, behavior, event, time_window, 
                                                 exp, batch, mouse, group, ylim=Y_LIM_DUAL,
                                                 dff_column = '560')
                        peth_plot.savefig(peth_path / f'{batch}_{mouse}_{behavior}_560_{event[0]}{time_window[0] - time_window[1]}_PETH.png')
                        peth_plot.savefig(peth_path / f'{batch}_{mouse}_{behavior}_560_{event[0]}{time_window[0] - time_window[1]}_PETH.pdf')
                        plt.close(peth_plot)

print(f"All plots saved to {peth_path}")
        
#%% Plot PETH for each group and extract mean and max Z-scored data

# ----------------------------- #
# PETH parameters
exp = 'RewardHab'
BOI = 'Licks_filtered'
baseline = False
MAXBOUTSNUMBER = None
event = 'onset'

# Plot parameters
TIME_WINDOW = [3, 8]
Y_LIM = [-2,10]
Y_LIM_DUAL = [-2,5]

if baseline:
    tag = "windowedbaseline"
else:
    tag = "wholetrace"

# Set groups
subjects_df['Group'] = subjects_df['Group'].fillna('')
included_groups = set(subjects_df['Group'])
# ----------------------------- #

print('##########################################')
print(f'EXPERIMENT: {exp}')
print('##########################################')

repo_path = exp_path / f'length{EVENT_TIME_THRESHOLD}_interbout{THRESH_S}_o{ORDER}f{CUT_FREQ}'
peth_path = repo_path / f'PETH_grouped_{tag}'
peth_path.mkdir(parents=True, exist_ok=True)

# Initialize data storage lists
subject_list = []
group_list = []
PETH_array = None
PETH_array_560 = None
PETH_mean_list = []
PETH_max_list = []
PETH_mean_list_560 = []
PETH_max_list_560 = []

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

        # --- 465 channel ---
        PETH_mouse = bp.PETH(
            dfiberbehav_clean, BOI, event, TIME_WINDOW, EVENT_TIME_THRESHOLD,
            baselinewindow=baseline, maxboutsnumber=MAXBOUTSNUMBER
        )
        PETH_mouse_mean = np.mean(PETH_mouse, axis=0, keepdims=True)

        if PETH_array is None:
            PETH_array = PETH_mouse_mean
            print('Initialized PETH_array successfully')
        else:
            PETH_array = np.concatenate((PETH_array, PETH_mouse_mean))

        # --- 560 channel ---
        if dual_color:
            PETH_mouse_560 = bp.PETH(
                dfiberbehav_clean, BOI, event, TIME_WINDOW, EVENT_TIME_THRESHOLD,
                baselinewindow=baseline, maxboutsnumber=MAXBOUTSNUMBER, dFF_column='560 dFF'
            )
            PETH_mouse_mean_560 = np.mean(PETH_mouse_560, axis=0, keepdims=True)

            if PETH_array_560 is None:
                PETH_array_560 = PETH_mouse_mean_560
                print('Initialized PETH_array_560 successfully')
            else:
                PETH_array_560 = np.concatenate((PETH_array_560, PETH_mouse_mean_560))

        # --- 465 channel: mean and max dFF before and after the event ---
        mean_before = np.mean(PETH_mouse[:TIME_WINDOW[0]])
        mean_after  = np.mean(PETH_mouse[TIME_WINDOW[0]:])
        max_before  = np.max(PETH_mouse[:TIME_WINDOW[0]])
        max_after   = np.max(PETH_mouse[TIME_WINDOW[0]:])

        PETH_mean_list.append((mean_before, mean_after))
        PETH_max_list.append((max_before, max_after))

        # --- 560 channel: mean and max dFF before and after the event ---
        if dual_color:
            mean_before_560 = np.mean(PETH_mouse_560[:TIME_WINDOW[0]])
            mean_after_560  = np.mean(PETH_mouse_560[TIME_WINDOW[0]:])
            max_before_560  = np.max(PETH_mouse_560[:TIME_WINDOW[0]])
            max_after_560   = np.max(PETH_mouse_560[TIME_WINDOW[0]:])

            PETH_mean_list_560.append((mean_before_560, mean_after_560))
            PETH_max_list_560.append((max_before_560, max_after_560))

        # Export mean/max PETH data to Excel
        export_dict = {
            'Subject': subject_list,
            'Group': group_list,
            f'465 Mean dFF before {BOI}': [x[0] for x in PETH_mean_list],
            f'465 Mean dFF after {BOI}':  [x[1] for x in PETH_mean_list],
            f'465 Max dFF before {BOI}':  [x[0] for x in PETH_max_list],
            f'465 Max dFF after {BOI}':   [x[1] for x in PETH_max_list],
        }
        if dual_color:
            export_dict.update({
                f'560 Mean dFF before {BOI}': [x[0] for x in PETH_mean_list_560],
                f'560 Mean dFF after {BOI}':  [x[1] for x in PETH_mean_list_560],
                f'560 Max dFF before {BOI}':  [x[0] for x in PETH_max_list_560],
                f'560 Max dFF after {BOI}':   [x[1] for x in PETH_max_list_560],
            })

        meanmaxPETH_df = pd.DataFrame(export_dict)
        meanmaxPETH_df.to_excel(peth_path / f'{BOI}_{TIME_WINDOW[0]}_{TIME_WINDOW[1]}_PETHmeanmax.xlsx')

# Plot PETH for each group

for group in included_groups:
    group_indices = [i for i, g in enumerate(group_list) if g == group]

    # --- 465 channel ---
    PETH_array_group = PETH_array[group_indices]
    print(f"Group {group} 465 PETH data size: {PETH_array_group.shape}")

    fig_PETHpooled = bp.plot_PETH_pooled(PETH_array_group, BOI, event, TIME_WINDOW, exp, group, ylim=Y_LIM)
    fig_PETHpooled.savefig(peth_path / f'{group}_{BOI}_{TIME_WINDOW[1]}_PETH.pdf')
    fig_PETHpooled.savefig(peth_path / f'{group}_{BOI}_{TIME_WINDOW[1]}_PETH.png')
    plt.close(fig_PETHpooled)

    # --- 560 channel ---
    if dual_color and PETH_array_560 is not None:
        PETH_array_group_560 = PETH_array_560[group_indices]
        print(f"Group {group} 560 PETH data size: {PETH_array_group_560.shape}")

        fig_PETHpooled_560 = bp.plot_PETH_pooled(PETH_array_group_560, BOI, event, TIME_WINDOW, exp, group,
                                                  ylim=Y_LIM_DUAL, dff_column='560')
        fig_PETHpooled_560.savefig(peth_path / f'{group}_{BOI}_{TIME_WINDOW[1]}_560_PETH.pdf')
        fig_PETHpooled_560.savefig(peth_path / f'{group}_{BOI}_{TIME_WINDOW[1]}_560_PETH.png')
        plt.close(fig_PETHpooled_560)
# %%
