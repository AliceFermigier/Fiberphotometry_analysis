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

#%%
from scripts.loader import analysis_path, data_path, proto_df, subjects_df, batches

#filter characteristics
ORDER = 4
CUT_FREQ = None #in Hz
#threshold to fuse behaviour if bouts are too close, in secs
THRESH_S = 10
#threshold for PETH : if events are too short do not plot them and do not include them in PETH, in seconds
EVENT_TIME_THRESHOLD = 0

#%% Plot PETH for each mouse

# PETH parameters 
baseline = True
MAXBOUTSNUMBER = 10
if baseline:
    tag = "windowedbaseline"
else:
    tag = "wholetrace"
for exp in ['Reward_Hab']: #[f.name for f in analysis_path.iterdir() if f.is_dir()]:
    exp_path = analysis_path / exp
    datapath_exp_dict = nom.get_experiment_data_path(batches, proto_df, data_path, exp)

    EVENT_LIST = ['onset']  # Event triggers, e.g., onset, withdrawal
    TIME_WINDOWS = [[1, 3]]  # Time window for PETH calculation (pre, post)

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

            # List all behaviors of interest (BOI) by excluding specific behaviors
            behaviors_of_interest = ['Licks_filtered','Nose_in_any_airport']
            
            for behavior in behaviors_of_interest:
                for event, time_window in zip(EVENT_LIST, TIME_WINDOWS):  
                    try:
                        # Generate the PETH data for the current behavior, event, and time window
                        peth_data = bp.PETH(dfiberbehav_df, behavior, event, time_window, EVENT_TIME_THRESHOLD, baselinewindow = baseline, maxboutsnumber=MAXBOUTSNUMBER)
                        
                        # Create a DataFrame from the PETH data
                        sr = round(pp.samplerate(dfiberbehav_df))
                        PRE_TIME, POST_TIME = time_window
                        n_timepoints = (PRE_TIME + POST_TIME) * sr + 1
                        time_index = np.linspace(-PRE_TIME, POST_TIME, n_timepoints)

                        peth_df = pd.DataFrame(np.transpose(peth_data), index=time_index)
                        
                        # Plot the PETH and save the figure 
                        peth_plot = bp.plot_PETH(peth_data, behavior, event, time_window, exp, mouse, group)
                        plot_filename = f'{mouse}_{behavior}_{event[0]}{time_window[0] - time_window[1]}_PETH.png'
                        peth_plot_path = peth_path / plot_filename
                        peth_plot.savefig(peth_plot_path)
                        plt.close(peth_plot)
                    
                    except Exception as e:
                        print(f'Error computing PETH for {behavior}, {mouse} : {e}')
                    
                        
#%% Plot PETH for each group and extract mean and max Z-scored data

# ----------------------------- #
# Parameters
BOI = 'Licks_filtered'
#'Licks_filtered' 'Nose_in_any_airport'
TIME_WINDOW = [1, 3]  # In seconds
MAXBOUTSNUMBER = 10
# ----------------------------- #

print('##########################################')
print(f'EXPERIMENT: {exp}')
print('##########################################')

repo_path = exp_path / f'length{EVENT_TIME_THRESHOLD}_interbout{THRESH_S}_o{ORDER}f{CUT_FREQ}'
peth_path = repo_path / 'PETH'
peth_path.mkdir(parents=True, exist_ok=True)  # Create PETH directory if it doesn't exist

# Initialize data storage lists
subject_list = []
group_list = []
PETH_array = None
PETH_mean_list = []
PETH_max_list = []


# Loop over each subject (mouse)
for mouse, batch, group in zip(subjects_df['Subject'], subjects_df['Batch'], subjects_df['Group']):
    print("--------------")
    print(f'MOUSE: {mouse} {batch}')
    print("--------------")
    
    fiberbehav_file = repo_path / f'{batch}_{mouse}_fiberbehav.csv'
    
    if not fiberbehav_file.exists():
        print(f"File not found: {fiberbehav_file}")
        continue
    
    fiberbehav_df = pd.read_csv(fiberbehav_file, index_col=0)
    sr = pp.samplerate(fiberbehav_df)
    
    if BOI in fiberbehav_df.columns[2:].tolist():
        subject_list.append(mouse)
        group_list.append(group)
        print(f'PETH {BOI} for {mouse}')
        
        # Calculate mean PETH for the current mouse
        PETH_mouse = bp.PETH(
        fiberbehav_df, BOI, 'onset', TIME_WINDOW, EVENT_TIME_THRESHOLD,
        maxboutsnumber=MAXBOUTSNUMBER
        )
        #print(f'PETH mouse : {PETH_mouse}, lenght = {len(PETH_mouse)}')
        PETH_mouse_mean = np.mean(PETH_mouse, axis=0, keepdims=True)
        #PETH_mouse_mean = gaussian_filter1d(PETH_mouse_mean, sigma=0.7)


        if PETH_array is None:
            PETH_array = PETH_mouse_mean
            print('Initialized PETH_array successfully')
        else:
            PETH_array = np.concatenate((PETH_array, PETH_mouse_mean))  # Stack new data
        
        # Calculate mean and max dFF before and after the event (PETH)
        mean_before = np.mean(PETH_mouse[:TIME_WINDOW[0]])  # Mean before event
        mean_after = np.mean(PETH_mouse[TIME_WINDOW[0]:])   # Mean after event
        max_before = np.max(PETH_mouse[:TIME_WINDOW[0]])    # Max before event
        max_after = np.max(PETH_mouse[TIME_WINDOW[0]:])     # Max after event
        
        PETH_mean_list.append((mean_before, mean_after))
        PETH_max_list.append((max_before, max_after))

# Export mean/max PETH data to Excel
meanmaxPETH_df = pd.DataFrame({
    'Subject': subject_list,
    'Group': group_list,
    f'Mean dFF before {BOI}': [x[0] for x in PETH_mean_list],
    f'Mean dFF after {BOI}': [x[1] for x in PETH_mean_list],
    f'Max dFF before {BOI}': [x[0] for x in PETH_max_list],
    f'Max dFF after {BOI}': [x[1] for x in PETH_max_list]
})
meanmaxPETH_df.to_excel(peth_path / f'{BOI}_{TIME_WINDOW[0]}_{TIME_WINDOW[1]}_PETHmeanmax.xlsx')

# Plot PETH for each group
included_groups = ['Saline', 'MEC 20uM']
for group in included_groups:
    # Filter PETH data for the current group
    group_indices = [i for i, g in enumerate(group_list) if g == group]
    PETH_array_group = PETH_array[group_indices]
    
    print(f"Group {group} PETH data size: {PETH_array_group.shape}")

    # Plot pooled PETH for the group
    fig_PETHpooled = bp.plot_PETH_pooled(PETH_array_group, BOI, 'onset', TIME_WINDOW, exp, group)
    fig_PETHpooled.savefig(peth_path / f'{group}_{BOI}_{TIME_WINDOW[1]}_PETH.pdf')
    fig_PETHpooled.savefig(peth_path / f'{group}_{BOI}_{TIME_WINDOW[1]}_PETH.png')
# %%
