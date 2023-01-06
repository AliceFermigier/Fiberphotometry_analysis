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
import numpy as np
#with pathlib.Path : file_to_open = data_folder / "raw_data.txt"

#%%
########
#LOADER#
########
experiment_path = Path('K:\\Alice\\Fiber\\202209_CA2db4')
analysis_path = experiment_path / 'Analysis'
data_path = experiment_path / 'Data'
os.chdir(experiment_path)
os.getcwd()

#import ID and groups of all mice
subjects_df = pd.read_excel(experiment_path / 'subjects.xlsx', sheet_name='Included')

#import tasks in protocol
proto_df = pd.read_excel(experiment_path / 'protocol.xlsx')

#all analysis files in Batch/Experiment/Analysis/Subject_Trial
#each trial of each subject has a specific directory called CD1_0, CD1_2, etc...
#for example, for ORM : CD1_0 = training of subject CD1 ; CD1_2 = test of subject CD1

#PARAMETERS

#all analysis files in experimentpath/date_Analysis/date_Experiment/Session/Mouse
#example : /Volumes/My Passport/Alice/Fiber/202103_CDHFDdvHPC/202103_Analysis/20210302_ORM/Test 1h/CD1
    
SAMPLERATE = 10 #in Hz

#list of behaviours on which to do peri event time histograms (PETH)
list_EVENT = ['onset', 'withdrawal']
list_TIMEWINDOW = [[4,5],[4,7]] 

#FIBER_BEHAV################################################################
#time before behaviour for calculation of PETH baseline, in seconds
PRE_EVENT_TIME = 1

#time to crop at the beginning of the trial to remove artifact, in seconds
TIME_BEGIN = 60

#threshold to fuse behaviour if bouts are too close, in secs
THRESH_S = 2

#threshold for PETH : if events are too short, not in PETH
EVENT_TIME_THRESHOLD = 0

#filter caractéristics
ORDER = 3
CUT_FREQ = 1 #in Hz

#FIBER TRACKING#######################
#number of pixels in fiberpho heatmap
RES = 40 
#behaviour used as 'start'
BEHAV_START = 'Gate opens'
#time baseline if behaviour hasn't been scored
TIME_BASELINE = 120.0 #in seconds

#CAUTION : with the code as it is, can process up to 4 different behaviours
#if more, add elif to align_behav function

#%%
########
#SCRIPT#
########

# #create analysis folder nomenclature
# for task in ['OdDis_1']:
#     if not os.path.exists(analysis_path / task):
#         os.mkdir(analysis_path / task)
#     for session in proto_df.loc[proto_df['Task']==task,'Sessions'].values[0].split(','):
#         if not os.path.exists(analysis_path / task / session):
#             os.mkdir(analysis_path / task / session)
#         for subject in subjects_df['Subject']:
#             if not os.path.exists(analysis_path / task / session / subject):
#                 os.mkdir(analysis_path / task / session / subject)
                
                
# #%%                
# #extract camera df from rawdata file
# df_camera = pd.read_csv('/Volumes/My Passport/Alice/Fiber/202110_CA2db2/20211004_AliceF_CA2b2bedding/HFDm1_0_camera.csv')
# df_rawdata = pd.read_csv('/Volumes/My Passport/Alice/Fiber/202110_CA2db2/20211004_AliceF_CA2b2bedding/HFDm1_0.csv')
# ind_list_cam = np.where(df_camera['Digital I/O | Ch.3 DI/O-3'] == 1)
# ind_list_raw = np.where(df_rawdata['Digital I/O | Ch.3'] == 1)
# start_cam = df_camera.at[ind_list_cam[0][0], 'Time(s)']
# start = df_rawdata.at[ind_list_raw[0][0], '---']


# df_camera_2 = pd.read_csv('/Volumes/My Passport/Alice/Fiber/202110_CA2db2/20211004_AliceF_CA2b2bedding/HFDm1_1_camera.csv')
# df_rawdata_2 = pd.read_csv('/Volumes/My Passport/Alice/Fiber/202110_CA2db2/20211004_AliceF_CA2b2bedding/HFDm1_1.csv')
# ind_list_cam2 = np.where(df_camera_2['Digital I/O | Ch.3 DI/O-3'] == 1)
# ind_list_raw2 = np.where(df_rawdata_2['Digital I/O | Ch.3'] == 1)
# start_cam_2 = df_camera_2.at[ind_list_cam2[0][0], 'Time(s)']
# start_2 = df_rawdata_2.at[ind_list_raw2[0][0], '---']