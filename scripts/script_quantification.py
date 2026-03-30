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

#%% 1 - Compute variance and transients on whole trace
# ----------------------------- #
# Parameters
exp = 'EPM'
list_BOI = ['Open arm', 'Closed arm', 'Center']

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
for mouse, batch, group in zip(subjects_df['Subject'], subjects_df['Batch'], subjects_df['Group']):
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
