#%%IMPORTED
###########

import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import os
import importlib

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

from scripts.loader import analysis_path, data_path, exp, proto_df, subjects_df, batches

#filter characteristics
ORDER = 4
CUT_FREQ = None #in Hz
#threshold to fuse behaviour if bouts are too close, in secs
THRESH_S = 0
#threshold for PETH : if events are too short do not plot them and do not include them in PETH, in seconds
EVENT_TIME_THRESHOLD = 0

#%% 2.2 - Calculate mean, max, and delta dFF within behavioural states (for state behaviours)
for session_path in [Path(f.path) for f in os.scandir(exp_path) if f.is_dir()]:
    session = session_path.name  # Extract session name
    print('###################')
    print(f'EXPERIMENT: {exp}')
    print('###################')
    
    # Create necessary paths
    repo_path = session_path / f'length{EVENT_TIME_THRESHOLD}_interbout{THRESH_S}_o{ORDER}f{CUT_FREQ}'
    groupanalysis_path = repo_path / 'Group_analysis'
    groupanalysis_path.mkdir(exist_ok=True)  # Create directory if it doesn't exist
    
    # Create lists to store data for export 
    mean_dFFs_list = []
    diffmeanmaxdFF_list = []
    diffmeanmaxdFF_perbout_list = []
    
    for mouse, batch, group in zip(subjects_df['Subject'], subjects_df['Batch'], subjects_df['Group']):
        
        # Paths to the relevant fiberbehav files
        fiberbehav_path = repo_path / f'{mouse}_fiberbehav.csv'
        fiberbehav_notderived_path = repo_path / f'{mouse}_fiberbehavnotderived.csv'
        
        # Check if required data exists for this mouse
        if fiberbehav_path.exists():
            print("-------------------")
            print(f'PROCESSING MOUSE: {mouse} batch {batch}')
            print("-------------------")
            
            # 1 Load data from CSV files
            dfiberbehav_df = pd.read_csv(fiberbehav_path, index_col=0)
            fiberbehav_df = pd.read_csv(fiberbehav_notderived_path, index_col=0)
            
            # 2 Calculate metrics (mean, max, and delta dFF)
            mean_dFF_result = sc.meandFF_behav(list_BOI, dfiberbehav_df, exp, session, mouse, group)
            diffmeanmax_dFF_result = sc.diffmeanmaxdFF_behav(fiberbehav_df, list_BOI, mouse, group)
            diffmeanmax_dFF_perbout_result = sc.diffmeanmaxdFF_behav_perbout(fiberbehav_df, list_BOI, mouse, group)
            
            # 3 Add the results to the list for later export
            mean_dFFs_list.append(mean_dFF_result)
            diffmeanmaxdFF_list.append(diffmeanmax_dFF_result)
    
    # 4 Concatenate results and export to Excel
    if mean_dFFs_list:
        meandFFs_allmice = pd.concat(mean_dFFs_list, ignore_index=True)
        meandFFs_allmice_path = groupanalysis_path / f'{exp}_{session}_globmeandFFs.xlsx'
        meandFFs_allmice.to_excel(meandFFs_allmice_path, index=False)
        print(f'Saved global mean dFFs to {meandFFs_allmice_path}')
    else:
        print('No mean dFF data to export.')

    if diffmeanmaxdFF_list:
        diffdFFs_allmice = pd.concat(diffmeanmaxdFF_list, ignore_index=True)
        diffdFFs_allmice_path = groupanalysis_path / f'{exp}_{session}_diffmaxmeandFFs.xlsx'
        diffdFFs_allmice.to_excel(diffdFFs_allmice_path, index=False)
        print(f'Saved delta max mean dFFs to {diffdFFs_allmice_path}')
    else:
        print('No delta max mean dFF data to export.')
        
#%% 2.3 - Calculate mean and max before and after behaviour onset for whole group

# ----------------------------- #
BOI = 'Shock'  # Behavior of Interest
TIME_MEANMAX = 5  # Time window for mean and max calculation, before and after behaviour (seconds)
MAX_BOUTS_NUMBER = None  # Limit the number of bouts to analyze, if specified
# ----------------------------- #

for session_path in [Path(f.path) for f in os.scandir(exp_path) if f.is_dir()]:
    session = session_path.name  # Get the session name from the path
    print('##########################################')
    print(f'EXPERIMENT: {exp} - SESSION: {session}')
    print('##########################################')
    
    code = gp.session_code(session, exp)
    repo_path = session_path / f'length{EVENT_TIME_THRESHOLD}_interbout{THRESH_S}_o{ORDER}f{CUT_FREQ}'
    groupanalysis_path = repo_path / 'Group_analysis'
    
    if len(os.listdir(repo_path)) > 1:
        if not groupanalysis_path.exists():
            groupanalysis_path.mkdir()
        
        # Initialize result Dataframe
        results_df = pd.DataFrame(columns=['Subject',  
                                           'Group', 
                                           f'Mean dFF before {BOI}', 
                                           f'Mean dFF after {BOI}', 
                                           f'Max dFF before {BOI}', 
                                           f'Max dFF after {BOI}'])
        
        for mouse, batch, group in zip(subjects_df['Subject'], subjects_df['Batch'], subjects_df['Group']):
            file_prefix = f'{batch}_{mouse}_{code}'
            print("------------------")
            print(f'MOUSE: {mouse} {batch}')
            print("------------------")
            
            fiberbehav_file = repo_path / f'{file_prefix}_fiberbehav.csv'
            if not fiberbehav_file.exists():
                print(f"File not found: {fiberbehav_file}")
                continue
            
            try:
                # Read in the CSV file, get mouse group
                fiberbehav_df = pd.read_csv(fiberbehav_file)
                group = subjects_df.loc[subjects_df['Subject'] == mouse, 'Group'].values[0]
                # Calculate results for mouse
                results_mouse_df = sc.process_mouse_meanmax_beforeafter(fiberbehav_df, mouse, group, BOI, MAX_BOUTS_NUMBER)
                # Put results in general Dataframe
                results_df = pd.concat([results_df, results_mouse_df], ignore_index=True)
                
            except Exception as e:
                print(f"Error processing file {fiberbehav_file}: {e}")
                continue
        
        # Export data to Excel
        output_file = groupanalysis_path / f'{exp}_{session}_{TIME_MEANMAX}sbeforeafter{BOI}_maxmeandFF.xlsx'
        
        try:
            results_df.to_excel(output_file, index=False)
            print(f"Data successfully exported to {output_file}")
        except Exception as e:
            print(f"Error exporting data to {output_file}: {e}")



#%% 2.6 - Compute variance and transients on whole trace and pre/post baseline
# ----------------------------- #
# Parameters
exp = 'EPM_1'
if 'EPM' in exp:
    list_BOI = ['Open arm', 'Closed arm', 'Center']
else :
    list_BOI = []
lowcut = None  # Lowcut frequency for bandpass filter (Hz)
highcut = None    # Highcut frequency for bandpass filter (Hz)
threshold = 'two_MAD'
# ----------------------------- #

exp_path = analysis_path / exp
repo_path = exp_path / f'length{EVENT_TIME_THRESHOLD}_interbout{THRESH_S}_o{ORDER}f{CUT_FREQ}'
groupanalysis_path = repo_path / 'Group_analysis'
groupanalysis_path.mkdir(parents=True, exist_ok=True)

var_transients_list = []

# Loop over each subject (mouse)
for mouse, batch in zip(subjects_df['Subject'], subjects_df['Batch']):
    try :
        print("--------------")
        print(f'MOUSE: {mouse}, BATCH: {batch}')
        print("--------------")
        
        # Set file path for the fiber behavior CSV
        fiber_file = repo_path / f'{batch}_{mouse}_fiberbehavnotderived.csv'
        
        # Check if the file exists
        if not fiber_file.exists():
            print(f"File not found: {fiber_file}")
            continue

        # Load data
        dfiber_df = pd.read_csv(fiber_file)
        
        # Apply bandpass filter to the dFF signal
        dfiber_df['Filtered dFF'] = tr.bandpass_filter(dfiber_df, lowcut, highcut)
        
        # Plot signal and spectrum
        tr.plot_signal_and_spectrum(dfiber_df)
        
        # Retrieve group information
        group = subjects_df.loc[subjects_df['Subject'] == mouse, 'Group'].values[0]
        
        #Calculate variance and transients characteristics
        mouse_df, transients_fig = sc.variance_transients(dfiber_df, list_BOI, mouse, group, exp, batch, threshold)
        var_transients_list.append(mouse_df)
        transients_fig.savefig(groupanalysis_path / f'{mouse}_1o{ORDER}f{lowcut}_{highcut}_{threshold}.png')
    except Exception as e:
        print(f'Error while processing mouse {mouse} : {e}')
    
# Concatenate results and export to Excel
variability_df = pd.concat(var_transients_list, ignore_index=True)
output_file = groupanalysis_path / f'Variability_1o{ORDER}f{lowcut}_{highcut}_{threshold}.xlsx'
variability_df.to_excel(output_file, index=False)

print(f"Variability data saved to {output_file}")
 # %%
