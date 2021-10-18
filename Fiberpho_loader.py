 # -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 09:30:55 2021

Analysis of Doric fiber photometry data 

DeltaF/F data for 470nm and 405nm are calculated using Doric Neuroscience Studio
Behavioral scoring is done using Boris software

@author: Alice Fermigier
"""
#%%
##########
#IMPORTED#
##########

import os
import pandas as pd
from pathlib import Path
#with pathlib.Path : file_to_open = data_folder / "raw_data.txt"

#%%
########
#LOADER#
########
experiment_path = Path('/Users/alice/Desktop/PhD/Manips/Fiber/202110_CA2dbatch2/')
analysis_path = experiment_path / '202110_Analysis'
os.chdir(experiment_path)
os.getcwd()

#import ID and groups of all mice
subjects_df = pd.read_excel('subjects.xlsx')

#import tasks in protocol
proto_df = pd.read_excel('protocol.xlsx')

#all analysis files in Batch/Experiment/Analysis/Subject_Trial
#each trial of each subject has a specific directory called CD1_0, CD1_2, etc...
#for example, for ORM : CD1_0 = training of subject CD1 ; CD1_2 = test of subject CD1

#PARAMETERS

#all analysis files in experimentpath/date_Analysis/date_Experiment/Session/Mouse
#example : /Volumes/My Passport/Alice/Fiber/202103_CDHFDdvHPC/202103_Analysis/20210302_ORM/Test 1h/CD1
    
SAMPLERATE = 10 #in samples per second

#list of behaviours on which to do peri event time histograms (PETH)
list_EVENT = ['onset', 'withdrawal']
list_TIMEWINDOW = [[4,5],[4,7]] 

#time before behaviour for calculation of PETH baseline, in seconds
PRE_EVENT_TIME = 1

#time to crop at the beginning of the trial to remove artifact, in seconds
TIME_BEGIN = 60

#CAUTION : with the code as it is, can process up to 4 different behaviours
#if more, add elif to align_behav function

#%%
########
#SCRIPT#
########

#create analysis folder nomenclature
for task in proto_df['Task']:
    if not os.path.exists(analysis_path / task):
        os.mkdir(analysis_path / task)
    for session in proto_df.loc[proto_df['Task']==task,'Sessions'].values[0].split(','):
        if not os.path.exists(analysis_path / task / session):
            os.mkdir(analysis_path / task / session)
        for subject in subjects_df['Subject']:
            if not os.path.exists(analysis_path / task / session / subject):
                os.mkdir(analysis_path / task / session / subject)