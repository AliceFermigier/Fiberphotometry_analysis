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
analysis_dir = Path('D:/Alice/Fiber/202110_CA2db2')
data_path = analysis_dir / 'Data' / 'Behav_todo'
######################################################################################

os.chdir(analysis_dir)
os.getcwd()

#%%
#############
# A REMPLIR #
######################################################################################
session = 'Test'

subjects_path = 'subjects.xlsx'

#if Boris binary file:
SAMPLERATE = 0.1 #(in seconds)
######################################################################################

df_subjects = pd.read_excel(subjects_path)

df_data = pd.read_excel(behav_path)

set_behav = set(df_data['Behavior'])
list_behav = list(set_behav)
list_behav.sort()
print(list_behav)

###################
#DEFINED FUNCTIONS#
###################

def tri_data(df_subjects, df_data, list_behav):
    """
    Sorts data at 5 and 10 minutes
    Takes Boris file with duration of each behaviour
    """
    data_array=np.zeros((4,len(df_subjects)))
    list_subjects = df_subjects['Subject'].tolist()
    for (subject, group) in zip(df_subjects['Subject'], df_subjects['Group']):
        df_subj = df_data.loc[df_data['Subject'] == subject]
        first_explo = float(df_subj['Start (s)'].tolist()[0])
        for behav in list_behav:
            df_behav = df_subj.loc[df_subj['Behavior'] == behav]
            df_5min = df_behav[df_behav['Stop (s)']<300+first_explo]
            df_10min = df_behav.loc[df_behav['Stop (s)']<600+first_explo]
            total5min = df_5min['Duration (s)'].sum()
            total10min = df_10min['Duration (s)'].sum()
            data_array[list_behav.index(behav)][list_subjects.index(subject)] = total5min
            data_array[list_behav.index(behav)+2][list_subjects.index(subject)] = total10min
    df_analysis_5min = pd.DataFrame({'Subject':list_subjects,
                                'Group':df_subjects['Group'].tolist(),
                                list_behav[0]+' 5min (s)':data_array[0],
                                list_behav[1]+' 5min (s)':data_array[1],
                                list_behav[0]+' 5min (%)':data_array[0]/(data_array[0]+data_array[1])*100,
                                list_behav[1]+' 5min (%)':data_array[1]/(data_array[0]+data_array[1])*100,
                                'Total exploration 5min (s)':data_array[0]+data_array[1],
                                'Index 5min':(data_array[1]-data_array[0])/(data_array[0]+data_array[1])})
    df_analysis_10min = pd.DataFrame({'Subject':list_subjects,
                                'Group':df_subjects['Group'].tolist(),
                                list_behav[0]+' 10min (s)':data_array[2],
                                list_behav[1]+' 10min (s)':data_array[3],
                                list_behav[0]+' 10min (%)':data_array[2]/(data_array[2]+data_array[3])*100,
                                list_behav[1]+' 10min (%)':data_array[3]/(data_array[2]+data_array[3])*100,
                                'Total exploration 10min (s)':data_array[2]+data_array[3],
                                'Index 10min':(data_array[3]-data_array[2])/(data_array[2]+data_array[3])})
    data_array_510_1 = data_array[2]-data_array[0]
    data_array_510_2 = data_array[3]-data_array[1]
    df_analysis_510min = pd.DataFrame({'Subject':list_subjects,
                                'Group':df_subjects['Group'].tolist(),
                                list_behav[0]+' 5-10min (s)':data_array_510_1,
                                list_behav[1]+' 5-10min (s)':data_array_510_2,
                                list_behav[0]+' 5-10min (%)':data_array_510_1/(data_array_510_1+data_array_510_2)*100,
                                list_behav[1]+' 5-10min (%)':data_array_510_2/(data_array_510_1+data_array_510_2)*100,
                                'Total exploration 5-10min (s)':data_array_510_1+data_array_510_2,
                                'Index 5-10min':(data_array_510_2-data_array_510_1)/(data_array_510_1+data_array_510_2)})
    
    
    with pd.ExcelWriter(data_file[:-5]+'_analysis.xlsx') as writer:
        df_analysis_5min.to_excel(writer, sheet_name='5 min')
        df_analysis_10min.to_excel(writer, sheet_name='10 min')
        df_analysis_510min.to_excel(writer, sheet_name='5-10 min')
    
    return()
    
def plot_data(df_subjects, df_data):
    """
    Plots data all time for eachgroup
    Need Boris binary file
    """
    
    
    
def plot_data_indiv(df_subjects, df_data):
    """
    Plots data all time for each subject
    """
    
########
#SCRIPT#
########

tri_data(df_subjects, df_data, list_behav)