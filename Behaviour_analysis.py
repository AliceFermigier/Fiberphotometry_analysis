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
import matplotlib.pyplot as plt

########
#LOADER#
########

#%%
#############
# A REMPLIR #
######################################################################################
analysis_dir = Path('K:\\Alice\\Fiber\\202209_CA2db4')
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

df_subjects = pd.read_excel(subjects_path, sheet_name='Included')

###################
#DEFINED FUNCTIONS#
###################
    
def get_data_timedep(subject, behav10Sps, list_behav):
    """
    Gets data all time for eachgroup
    Need Boris binary file
    """
    
    ind_start = np.where(behav10Sps['Gate opens'] == 1)[0].tolist()[0]
    time_start = behav10Sps.at[ind_start,'time']
    time_stop = behav10Sps.at[len(behav10Sps)-1,'time']
    
    list_time = []
    list_time_behavs = [[]]*len(list_behav)
    for i,time in enumerate(range(time_start+TIMEBIN,time_stop+1,TIMEBIN)):
        df_timebin = behav10Sps[behav10Sps['time']>=time_start & behav10Sps['time']<time]
        list_time.append(time)
        for i,behav in enumerate(list_behav):
            list_time_behavs[i].append(df_timebin[behav].sum()/10)
    
    return(list_time_behavs)

def get_data(subject, behav10Sps, list_behav):
    """
    Gets data exploration on whole trial
    Need Boris binary file
    """
    
    list_time_behavs = []
    for behav in list_behav:
        list_time_behavs.append(behav10Sps[behav].sum()/10)
        
    return(list_time_behavs)
    
    
########
#%%SCRIPT#
########

exp_path = data_path / '20220928_AliceF_OdDispostFear'
exp = str(exp_path).split('_')[-1]

session = 'test'
code = 1

behav10Sps = pd.read_csv(glob.glob(f'{exp_path}\\behav_{code}_*')[0])
list_behav = behav10Sps.columns[1:].tolist()
list_behav.remove('Entry in arena')
list_behav.remove('Gate opens')
if 'Climbing' in list_behav:
    list_behav.remove('Climbing')
print(list_behav)

if code == 0 :
    cmap = 'tab20'
    
if code == 1 :
    cmap = 'viridis' 
    #cmap = 'tab20'
    
if code == 2 :
    cmap = 'summer' 
    #cmap = 'tab20b'
     
#generate dataframe with all subjects
list_subj = []
list_group = []
list_times = []*len(list_behav)

for file in glob.glob(f'{exp_path}\\behav_{code}_*'):
    behav10Sps = pd.read_csv(file)
    subject = file.split('_')[-1][:-4]
    list_subj.append(subject)
    list_group.append(df_subjects.loc[df_subjects['Subject']==subject,'Group'].values[0])
    list_time_behavs = get_data(subject, behav10Sps, list_behav)
    list_times.append(list_time_behavs)
    print(list_times)
        
df_exp=pd.DataFrame({'Subject':list_subj,
                     'Group':list_group,
                     list_behav[0]:[time[0] for time in list_times],
                     list_behav[1]:[time[1] for time in list_times]})

df_exp.to_excel(exp_path / f'{exp}_{code}_count.xlsx')
                        
df_plot=df_exp.groupby('Group')
means=df_plot.mean()
errors=df_plot.std()

#plot figure
fig_exp, axmean = plt.subplots(1, 1, figsize=(7, 6))
labels = df_exp.columns[2:]

means.plot.bar(y=labels, yerr=errors[labels], capsize=2, rot=0, ax=axmean, linewidth=.1, colormap=cmap)

axmean.set_ylabel('Exploration (s)')
axmean.set_title(f'Time exploring - {exp} {session}')
axmean.legend(loc='upper right')

fig_exp.savefig(exp_path / f'{exp}_{code}_histogram.png')


#%%generate dataframe with all timepoints and all subjects
list_subj = []
list_timebins = []

for file in glob.glob(f'{data_path}\\{exp_path}\\behav_{code}_*'):
    subject = file.split('_')[-1]
    group = df_subjects.loc[df_subjects['Subject']==subject,'Group']
    list_time_behavs = get_data_timedep(subject, file, list_behav)
    list_subj.append([subject]*len(list_time_behavs[0]))
    [df_subjects.loc[df_subjects['Subject']==subject,['Group']]]*len(list_time_behavs[0])


