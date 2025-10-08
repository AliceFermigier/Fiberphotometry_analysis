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
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if project_root not in sys.path:
    sys.path.insert(0, project_root)
    
import modules.common.nomenclature as nom

########
#LOADER#
########

experiment_path = Path(r'E:\202510_FiberMEC')
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
TIME_BEGIN = 0
#------------------#
batches = list(set(subjects_df['Batch']))
#------------------#



