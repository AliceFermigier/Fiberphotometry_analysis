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

########
#LOADER#
########

#%%
#############
# A REMPLIR #
######################################################################################
analysis_dir = Path('K:\\Alice\\Fiber\\202110_CA2db2')
data_path = analysis_dir / 'Data' / 'To_do'
TIMEBIN = 60 #in seconds
SAMPLERATE = 0,1 #in seconds
######################################################################################

os.chdir(analysis_dir)
os.getcwd()

#%%
#############
# A REMPLIR #
######################################################################################
subjects_path = 'subjects.xlsx'

#if Boris binary file:
SAMPLERATE = 0.1 #(in seconds)
######################################################################################

df_subjects = pd.read_excel(subjects_path)

###################
#DEFINED FUNCTIONS#
###################
    
def gen_data(subject, behav10Sps, TIMEBIN, SAMPLERATE):
    """
    Plots data all time for eachgroup
    Need Boris binary file
    """
    
    list_behav = behav10Sps.columns[1:].tolist()
    list_behav.remove('Entry in arena')
    list_behav.remove('Gate opens')
    
    ind_start = np.where(behav10Sps['Gate opens'] == 1)[0].tolist()[0]
    time_start = behav10Sps.at[ind_start,'time']
    time_stop = behav10Sps.at[len(behav10Sps)-1,'time']
    
    list_sum = []
    for time in range(time_start,time_stop,TIMEBIN):
        sum_behav = behav10Sps[behav10Sps['time']>=time_start & behav10Sps['time']<time].sum()
        list_sum.append(sum_behav)
        
        
    
def plot_data_indiv(df_subjects, df_data):
    """
    Plots data all time for each subject
    """
    
########
#SCRIPT#
########

behav_test = 'D:\\Alice\\Fiber\\202110_CA2db2\\Data\\To_do\\20211006_AliceF_CA2b2SRM\\behav_2_HFDm3.csv'
behav_df = pd.read_csv(behav_test)

subject = 'HFDm3'
gen_data(subject, behav_df, TIMEBIN, SAMPLERATE)
