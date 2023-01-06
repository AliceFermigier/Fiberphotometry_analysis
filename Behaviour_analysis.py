# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 10:21:22 2021

Score behaviour of fiberpho behav data

@author: alice fermigier
"""

##########
#IMPORTED#
##########

import pandas as pd
import os
import numpy as np
from pathlib import Path
import glob

########
#LOADER#
########

#%%
#############
# A REMPLIR #
######################################################################################
analysis_dir = Path('K:\\Alice\\Fiber\\202110_CA2db2')
data_path = analysis_dir / 'Data' / 'To_do'
TIMEBIN = 30 #in seconds
SAMPLERATE = 0,1 #in seconds
######################################################################################

os.chdir(analysis_dir)
os.getcwd()

#%%
#############
# A REMPLIR #
######################################################################################
subjects_path = 'subjects.xlsx'
######################################################################################

df_subjects = pd.read_excel(subjects_path)

###################
#DEFINED FUNCTIONS#
###################
    
def get_data(subject, behav10Sps):
    """
    Plots data all time for eachgroup
    Need Boris binary file
    """
    
    list_behav = behav10Sps.columns[1:].tolist()
    list_behav.remove('Entry in arena')
    list_behav.remove('Gate opens')
    list_behav.remove('Climbing')
    
    ind_start = np.where(behav10Sps['Gate opens'] == 1)[0].tolist()[0]
    time_start = behav10Sps.at[ind_start,'time']
    time_stop = behav10Sps.at[len(behav10Sps)-1,'time']
    
    list_time = []
    list_time_behavs = [[]]*len(list_behav)
    for i,time in enumerate(range(time_start+TIMEBIN,time_stop+1,TIMEBIN)):
        df_timebin = behav10Sps[behav10Sps['time']>=time_start & behav10Sps['time']<time]
        list_time.append((i+1)*TIMEBIN)
        for i,behav in enumerate(list_behav):
            list_time_behavs[i].append(df_timebin[behav].sum())
    
    return([subject]*len(list_time),
           [df_subjects.loc[df_subjects['Subject']==subject,['Group']]]*len(list_time),
           list_time,
           list_time_behavs)
    
def plot_data_indiv(df_subjects, df_data):
    """
    Plots data all time for each subject
    """
    
########
#SCRIPT#
########

behav_test = 'D:\\Alice\\Fiber\\202110_CA2db2\\Data\\To_do\\20211006_AliceF_CA2b2SRM\\behav_2_HFDm3.csv'
behav_df = pd.read_csv(behav_test)

for exp_path in os.listdir(data_path):
    exp = exp_path.split('_')[-1]
    for file in glob.glob(f'{data_path}\\{exp_path}\\behav_*'):
        code 
        
        

subject = 'HFDm3'
get_data(subject, behav_df, TIMEBIN, SAMPLERATE)
