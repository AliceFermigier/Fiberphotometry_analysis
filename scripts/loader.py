# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 18:45:55 2025

Loader

@author: alice
"""

import os
import pandas as pd
from pathlib import Path
import sys

#path to other scripts in sys.path
path_to_gitrepo=r'C:\Users\afermigier\Documents\GitHub\Fiberphotometry_analysis'
if path_to_gitrepo not in sys.path:
    sys.path.append(path_to_gitrepo)
    
import modules.nomenclature as nom

########
#LOADER#
########

experiment_path = Path(r'C:\Users\alice\Desktop\Papier Bedding\Data\Data_Fiber')
analysis_path = experiment_path / 'Analysis' 
data_path = experiment_path / 'Data'
os.chdir(experiment_path)
os.getcwd()

#import ID and groups of all mice
subjects_df = pd.read_excel(experiment_path / 'subjects.xlsx', sheet_name='Included')
#import tasks in protocol
proto_df = pd.read_excel(experiment_path / 'protocol.xlsx')
#create artifacts file if not existent
artifact_file = experiment_path / 'artifacts.xlsx' # File to store artifact timestamps
nom.create_or_load_artifacts_file(artifact_file, option='create_only')

############
#PARAMETERS#
############

#time to crop at the beginning of the trial for, in seconds
TIME_BEGIN = 60
#filter characteristics
ORDER = 4
CUT_FREQ = 1 #in Hz

#------------------#
exp = 'Fear'
batches = list(set(subjects_df['Batch']))
#------------------#

#threshold to fuse behaviour if bouts are too close, in secs
THRESH_S = 0.5
#threshold for PETH : if events are too short do not plot them and do not include them in PETH, in seconds
EVENT_TIME_THRESHOLD = 1

