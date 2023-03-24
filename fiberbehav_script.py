# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 13:39:29 2023

To run fiberphotometry analysis with behaviour

@author: alice fermigier
"""

#%%
##########
#IMPORTED#
##########

import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
plt.rcParams.update({'figure.max_open_warning' : 0})
from pathlib import Path
import os
import sys

#Custom
#put path to directory where python files are stored

if 'D:\Profil\Documents\GitHub\Fiberphotometry_analysis' not in sys.path:
    sys.path.append('D:\Profil\Documents\GitHub\Fiberphotometry_analysis')
    #sys.path.append('/Users/alice/Documents/GitHub/Fiberphotometry_analysis')

import preprocess as pp
import genplot as gp
import behavplot as bp
import visualize as vis
import statcalc as sc

#%%
########
#LOADER#
########

experiment_path = Path('K:\\Alice\\Fiber\\202301_CA2b5')
analysis_path = experiment_path / 'Analysis'
data_path = experiment_path / 'Data'
preprocess_path = experiment_path / 'Preprocessing'
os.chdir(experiment_path)
os.getcwd()

#import ID and groups of all mice
subjects_df = pd.read_excel(experiment_path / 'subjects.xlsx', sheet_name='Included')
#import tasks in protocol
proto_df = pd.read_excel(experiment_path / 'protocol.xlsx')

#PETH parameters
list_EVENT = ['onset', 'withdrawal']
list_TIMEWINDOW = [[4,5],[4,7]] 
#time before behaviour for calculation of PETH baseline, in seconds
PRE_EVENT_TIME = 1
#time to crop at the beginning of the trial for , in seconds
TIME_BEGIN = 60
#threshold to fuse behaviour if bouts are too close, in secs
THRESH_S = 0
#threshold for PETH : if events are too short do not plot them and do not include them in PETH, in seconds
EVENT_TIME_THRESHOLD = 0
#filter characteristics
ORDER = 3
CUT_FREQ = 2 #in Hz
# ideal samplerate : deinterleaved data samplerate should be 10 but is 9.99 instead ==> change index to align with behaviour
SAMPLERATE = 10 #in Hz

########
#SCRIPT#
########

#%% 1 - PREPROCESSING
#####################

#------------------#
exp =  'Essai1'

exp_path = analysis_path / exp
pp_path = exp_path / 'Preprocessing' # put all raw plots and deinterleaved_df in separate 
if not os.path.exists(pp_path):
    os.mkdir(pp_path)

#get data path related to the task in protocol excel file
data_path_exp = data_path / proto_df.loc[proto_df['Task']==exp, 'Data_path'].values[0]
#------------------#

#%% 1.1 - Deinterleave data and save in separate file.

for session_path in [Path(f.path) for f in os.scandir(exp_path) if f.is_dir()]:
    session = str(session_path).split('\\')[-1]
    print('##########################################')
    print(f'EXPERIMENT : {exp} - SESSION : {session}')
    print('##########################################')
    code = gp.session_code(session,exp)
    for mouse in subjects_df['Subject']:
        print("--------------")
        print(f'MOUSE : {mouse}')
        print("--------------")
        rawdata_df = pd.read_csv(data_path_exp/f'{mouse}_{code}.csv',skiprows=1,usecols=['Time(s)','AIn-1','DI/O-1','DI/O-2'])
        deinterleaved_df = pp.deinterleave(rawdata_df)
        deinterleaved_df.to_csv(pp_path/f'{mouse}_{code}_deinterleaved.csv')
        
# 1.2 - Look at rawdata and check for artifacts

        fig_raw = gp.plot_rawdata(rawdata_df,exp,session,mouse)
        fig_raw.savefig(pp_path/f'{mouse}_{code}_rawdata.pdf')
        
#%% 1.3 - Open artifacted data in Dash using plotly and enter artifact timestamps in excel file

#------------------#
session = 'Test'
mouse = 'A1f'
#------------------#
# in excel 'Filecode', put '{exp}_{session}_{mouse}'

deinterleaved_df = pd.read_csv(pp_path/f'{mouse}_{code}_deinterleaved.csv')
vis.visualize(deinterleaved_df,'fiberpho',exp,session,mouse)

#%% 1.4 - Artifact correction and dFF calculation

#import artifacts boundaries
artifacts_df = pd.read_excel(experiment_path / 'artifacts.xlsx')

for session_path in [Path(f.path) for f in os.scandir(exp_path) if f.is_dir()]:
    session = str(session_path).split('\\')[-1]
    print('##########################################')
    print(f'EXPERIMENT : {exp} - SESSION : {session}')
    print('##########################################')
    code = gp.session_code(session,exp)
    for mouse in subjects_df['Subject']:
        print("--------------")
        print(f'MOUSE : {mouse}')
        print("--------------")
        deinterleaved_df = pd.read_csv(pp_path/f'{mouse}_{code}_deinterleaved.csv')
        filecode = f'{exp}_{session}_{mouse}'
        sr = pp.samplerate(deinterleaved_df)
        
        # calculate dFF with artifacts removal, then interpolate missing data
        dFFdata_df = pp.dFF(deinterleaved_df,artifacts_df,filecode,sr,method='fit')
        interpdFFdata_df = pp.interpolate_dFFdata(dFFdata_df, method='linear')
        
        # smooth data with butterworth filter or simple moving average (SMA)
        smoothdFF_df=pp.smoothing_SMA(interpdFFdata_df,win_size=7)
        smoothdFF_df=pp.butterfilt(interpdFFdata_df, ORDER, CUT_FREQ)
        smoothdFF_df.to_csv(pp_path/f'{mouse}_{code}_dFFfilt.csv')

#%% 1.4 - Align with behaviour, create corresponding excel, plot fiberpho data with behaviour

for session_path in [Path(f.path) for f in os.scandir(exp_path) if f.is_dir()]:
    session = str(session_path).split('\\')[-1]
    print('##########################################')
    print(f'EXPERIMENT : {exp} - SESSION : {session}')
    print('##########################################')
    code = gp.session_code(session,exp)
    repo_path = session_path /  f'length{EVENT_TIME_THRESHOLD}_interbout{THRESH_S}_o{ORDER}f{CUT_FREQ}'
    if not os.path.exists(repo_path):
            os.mkdir(repo_path)
    for mouse in subjects_df['Subject']:
        print("--------------")
        print(f'MOUSE : {mouse}')
        print("--------------")
        rawdata_path = data_path_exp / f'{mouse}_{code}.csv'
        behav_path = data_path_exp / f'behav_{code}_{mouse}.csv'
        fiberpho_path = pp_path / f'{mouse}_{code}_dFFfilt.csv'
        
        #begin analysis only if behaviour has been scored and plots haven't been done
        ready = False
        if os.path.exists(behav_path):
            ready = True
        print(f'ready? {ready}')
        done = True
        if not os.path.exists(repo_path / f'{mouse}_{code}_fiberbehav.csv'):
            done = False
        print(f'done? {done}')
        
        if ready == True and done == False:
            rawdata_cam_df = pd.read_csv(rawdata_path, skiprows=1, usecols=['Time(s)','DI/O-3'])
            behav10Sps = pd.read_csv(behav_path)
            fiberpho = pd.read_csv(fiberpho_path)
            print(exp, session, mouse)
            
            #list of behaviours to analyze
            list_BOI = behav10Sps.columns[1:].tolist()
            print(f'Behaviours = {list_BOI}')
                
            #align behaviour and fiberpho data, create fiberbehav.csv
            timevector = gp.time_vector(fiberpho, SAMPLERATE)
            timestart_camera = gp.timestamp_camera(rawdata_cam_df)[0]
            print(f'start camera : {timestart_camera}')
            fiberbehav_df = bp.align_behav(behav10Sps, fiberpho, timevector, timestart_camera)
            dfiberbehav_df = bp.behav_process(fiberbehav_df, list_BOI, THRESH_S, EVENT_TIME_THRESHOLD)
            dfiberbehav_df.to_csv(repo_path / f'{mouse}_{code}_fiberbehav.csv')
            
            #plot fiberpho data with behaviour
            bp.plot_fiberpho_behav(dfiberbehav_df)
            
#%% 1.5 - Calculate mean and diff within behaviours

for session_path in [Path(f.path) for f in os.scandir(exp_path) if f.is_dir()]:
    session = str(session_path).split('\\')[-1]
    print('##########################################')
    print(f'EXPERIMENT : {exp} - SESSION : {session}')
    print('##########################################')
    code = gp.session_code(session,exp)
    repo_path = session_path /  f'length{EVENT_TIME_THRESHOLD}_interbout{THRESH_S}_o{ORDER}f{CUT_FREQ}'
    #create list of mean_dFFs and diff_dFFs
    mean_dFFs_list = []
    subject_list = []
    diffdFF_list = []
    for mouse in subjects_df['Subject']:
        print("--------------")
        print(f'MOUSE : {mouse}')
        print("--------------")
        subject_list.append(mouse)
        group = subjects_df.loc[subjects_df['Subject']==mouse, 'Group'].values[0]
        dfiberbehav_df = pd.read_csv(repo_path / f'{mouse}_{code}_fiberbehav.csv')
        
        mean_dFF_df = sc.meandFF_behav(list_BOI, fiberbehav_df, exp, session, mouse, group)
        mean_dFFs_list.append(mean_dFF_df)

#%% 1.6 - Plot PETH

#%%Plots for entry in zone
for exp_path in [Path(f.path) for f in os.scandir(analysis_path) if f.is_dir()]:
    exp = str(exp_path).split('\\')[-1]
    for session_path in [Path(f.path) for f in os.scandir(exp_path) if f.is_dir()]:
        session = str(session_path).split('\\')[-1]
        print('##########################################')
        print(f'EXPERIMENT : {exp} - SESSION : {session}')
        print('##########################################')
        code = session_code(session,exp)
        repo_path = session_path /  f'length{EVENT_TIME_THRESHOLD/10}_interbout{THRESH_S}_o{ORDER}f{CUT_FREQ}'
        
        #create list of mean_dFFs, max_dFFs
        meandFF_list = []
        subject_list = []
        group_list = []
        maxdFF_list = []
        PETH_array = None
        for mouse_path in [Path(f.path) for f in os.scandir(repo_path) if f.is_dir()]:
            # '/' on mac, '\\' on windows
            mouse = str(mouse_path).split('\\')[-1]
            print("################")
            print(f'MOUSE : {mouse}')
            if os.path.exists(mouse_path / f'{mouse}_fbprocess.xlsx'):
                fbprocess_df = pd.read_excel(mouse_path / f'{mouse}_fbprocess.xlsx')
                group = subjects_df.loc[subjects_df['Subject']==mouse, 'Group'].values[0]
                
                #1 - subjects
                subject_list.append(mouse)
                group_list.append(group)
                #2 - mean and max dFF
                ind_event = np.where(fbprocess_df['Entry in arena'] == 1)[0][0]
                ind_event = fbprocess_df.loc[ind_event-3*SAMPLERATE:ind_event, 'Denoised dFF'].idxmin()
                meandFF_list.append(fbprocess_df.loc[ind_event:ind_event+30*SAMPLERATE, 'Denoised dFF'].mean())
                maxdFF_list.append(max(fbprocess_df.loc[ind_event:ind_event+30*SAMPLERATE, 'Denoised dFF']))
                #3 - PETH data
                if PETH_array is None:
                    PETH_array = PETH(fbprocess_df, 'Entry in arena', 'onset', [6,10])
                else:
                    PETH_array = np.concatenate((PETH_array,PETH(fbprocess_df, 'Entry in arena', 'onset', [6,10])))
                    
        #plot PETH
        included_groups = ['CD','HFD']
        PETHarray_list=[]
        for group in included_groups:
            PETH_array_group = PETH_array
            list_todelete = []
            for (i,group_mouse) in enumerate(group_list): 
                if group not in group_mouse:
                    list_todelete.append(i)
            np.delete(PETH_array_group,(list_todelete),axis=0)
            PETHarray_list.append(np.delete(PETH_array_group,(list_todelete),axis=0))
        plot_PETH_pooled(included_groups, PETHarray_list, 'Entry in arena', 'onset', [6,10])
        
        #export data to excel
        if not os.path.exists(repo_path / f'{exp}_{session}_maxmeandFFentry.xlsx'):
            meanmaxdFFs_df = pd.DataFrame(data={'Subject' : subject_list, 'Group' : group_list, 
                                                'Mean dFF entry' : meandFF_list, 'Max dFF entry' : maxdFF_list})
            meanmaxdFFs_df.to_excel(repo_path / f'{exp}_{session}_maxmeandFFentry.xlsx')
            

#%%Run for all 

# for EVENT_TIME_THRESHOLD in [0, 1, 2]:
#     for THRESH_S in [0, 0.5, 1, 2, 3]:
#         for CUT_FREQ in [1,2,3,4]:

for exp_path in [Path(f.path) for f in os.scandir(analysis_path) if f.is_dir()]:
    exp = str(exp_path).split('\\')[-1]
    for session_path in [Path(f.path) for f in os.scandir(exp_path) if f.is_dir()]:
        session = str(session_path).split('\\')[-1]
        print('##########################################')
        print(f'EXPERIMENT : {exp} - SESSION : {session}')
        print('##########################################')
        code = session_code(session,exp)
        #get data path related to the task in protocol excel file
        data_path_exp = data_path / proto_df.loc[proto_df['Task']==exp, 'Data_path'].values[0]
        #create repository for values of thresholds : length and interbout
        repo_path = session_path / f'length{EVENT_TIME_THRESHOLD/10}_interbout{THRESH_S}_o{ORDER}f{CUT_FREQ}'
        if not os.path.exists(repo_path):
            os.mkdir(repo_path)
            for subject in subjects_df['Subject']:
                os.mkdir(repo_path / subject)
        #create list of mean_dFFs and diff_dFFs
        mean_dFFs_list = []
        subject_list = []
        diffdFF_list = []
        for mouse_path in [Path(f.path) for f in os.scandir(repo_path) if f.is_dir()]:
            # '/' on mac, '\\' on windows
            mouse = str(mouse_path).split('\\')[-1]
            subject_list.append(mouse)
            print("################")
            print(f'MOUSE : {mouse}')
            group = subjects_df.loc[subjects_df['Subject']==mouse, 'Group'].values[0]
            
            #get data
            behav_path = data_path_exp / f'behav_{code}_{mouse}.csv'
            fiberpho_path = data_path_exp / f'{mouse}_{code}_dFFfilt.csv'
            camera_path = data_path_exp / f'{mouse}_{code}_camera.csv'
            rawdata_path = data_path_exp / f'{mouse}_{code}.csv'
            
            #begin analysis only if behaviour has been scored and plots haven't been done
            ready = False
            if os.path.exists(behav_path):
                ready = True
            print(f'ready? {ready}')
            
            done = True
            if not os.path.exists(mouse_path / f'{mouse}_{code}_fbprocess.xlsx'):
                done = False
            print(f'done? {done}')
            
            if os.path.exists(camera_path) and ready == True and done == False:
                camera = pd.read_csv(camera_path)
                
            elif os.path.exists(rawdata_path) and ready == True and done == False:
                rawdata_cam_df = pd.read_csv(rawdata_path, skiprows=1, usecols=['Time(s)','DI/O-3'])
            
            if os.path.exists(behav_path) and os.path.exists(fiberpho_path) and ready == True and done == False:
                behav10Sps = pd.read_csv(behav_path)
                fiberpho = pd.read_csv(fiberpho_path)
                print(exp, session, mouse)
                
                #list of behaviours to analyze
                list_BOI = behav10Sps.columns[1:].tolist()
                
                #align behaviour and fiberpho data, create fbprocess.xslx
                print('timevector')
                timevector = time_vector(fiberpho, SAMPLERATE)
                print('timestamp')
                if os.path.exists(camera_path):
                    print('---------> from camera')
                    timestart_camera = timestamp_camera(camera)[0]
                else:
                    print('---------> from rawdata')
                    timestart_camera = timestamp_camera_fromraw(rawdata_cam_df)[0]
                print('start camera : ', timestart_camera)
                print('aligning')
                fiberbehav_df = align_behav(behav10Sps, fiberpho, timevector, timestart_camera)
                fiberbehav_df = filter_dFF(fiberbehav_df, ORDER, CUT_FREQ)
                
                mean_dFF_df = meandFF_behav(list_BOI, fiberbehav_df)
                mean_dFFs_list.append(mean_dFF_df)
                
                if EVENT_TIME_THRESHOLD==0 and THRESH_S==0:
                    #plot isosbestic and gcamp data
                    plot_fiberpho(fiberbehav_df)
                print('processing')
                (fiberbehav2_df, fbprocess_df) = behav_process(fiberbehav_df, list_BOI)
                
                #add diffdFFs to the list
                diffdFF_df = diff_dFF(fiberbehav_df, fbprocess_df, list_BOI)
                diffdFF_list.append(diffdFF_df)
                
                # #plot raw data
                #plot_rawdata(rawdata_df)
                
                #plot fiberpho data aligned with behav
                plot_fiberpho_behav(fbprocess_df)
                #plot_fiberpho_behav_snip(fiberbehav2_df, timestart_camera)
                
                for BOI in list_BOI:
                    if BOI in ['Entry in arena','Gate opens']:
                        PETH_data = PETH(fbprocess_df, BOI, 'onset', [6,10])
                        plot_PETH_average(PETH_data, BOI, 'onset', [6,10])
                    else:
                        for event, timewindow in zip(list_EVENT,list_TIMEWINDOW): #CAUTION : adapt function!!
                            PETH_data = PETH(fbprocess_df, BOI, event, timewindow)
                            plot_PETH(PETH_data, BOI, event, timewindow)
       
        # if mean_dFFs_list != []: 
        #     meandFFs_allmice = pd.concat(mean_dFFs_list)
        #     meandFFs_allmice.to_excel(repo_path / f'{exp}_{session}_length{EVENT_TIME_THRESHOLD/10}_interbout{THRESH_S}_filtero{ORDER}f{CUT_FREQ}_meandFFs.xlsx')
            
        #     df_plotmean=meandFFs_allmice.groupby('Group')
        #     means=df_plotmean.mean()
        #     errors=df_plotmean.std()
            
        #     #plot figure
        #     fig_meandFF, axmean = plt.subplots(1, 1, figsize=(7, 6))
        #     labels = meandFFs_allmice.columns[2:]
            
        #     means.plot.bar(y=labels, yerr=errors[labels], capsize=2, rot=0, 
        #                    ax=axmean, linewidth=.1, colormap='viridis')
                
        #     axmean.set_ylabel('Mean dFF')
        #     axmean.set_title(f'Mean dFF - {exp} {session} - length{EVENT_TIME_THRESHOLD/10} interbout{THRESH_S} - filtero{ORDER}f{CUT_FREQ}')
            
        #     fig_meandFF.savefig(repo_path / f'{exp}_{session}_length{EVENT_TIME_THRESHOLD/10}_interbout{THRESH_S}_filtero{ORDER}f{CUT_FREQ}_meandFFs.png')
            
        # if diffdFF_list != []:
        #     diffdFFs_allmice = pd.concat(diffdFF_list)
        #     diffdFFsmean = diffdFFs_allmice.groupby(['Subject','Behaviour'], as_index=False).mean()
        #     diffdFFsmean = diffdFFsmean[diffdFFsmean['Bout']>2]
            
        #     #group diffdFFsmean dataframe to plot it!
        #     cols = ['Group']
        #     list_behav = set(diffdFFsmean['Behaviour'])
        #     for behav in list_behav:
        #         cols.append(behav)
        #     list_subjects = set(diffdFFsmean['Subject'])
            
        #     diffdFFsmeanplot = pd.DataFrame(len(list_subjects)*[[0]*3],
        #                                     columns=cols, index=list_subjects)
        #     for subject in list_subjects:
        #         diffdFFsmeanplot.loc[subject, 'Group']=subject[:-1]
        #         df_subj = diffdFFsmean[diffdFFsmean['Subject']==subject]
        #         for behav in list_behav:
        #             if behav in df_subj['Behaviour'].values:
        #                 value = df_subj.loc[df_subj['Behaviour']==behav,'Delta dFF'].values[0]
        #                 diffdFFsmeanplot.loc[subject, behav]=value
             
        #     diffdFFs_allmice.to_excel(repo_path / f'{exp}_{session}_length{EVENT_TIME_THRESHOLD/10}_interbout{THRESH_S}_filtero{ORDER}f{CUT_FREQ}_diffdFFsall.xlsx')
        #     diffdFFsmeanplot.to_excel(repo_path / f'{exp}_{session}_length{EVENT_TIME_THRESHOLD/10}_interbout{THRESH_S}_filtero{ORDER}f{CUT_FREQ}_diffdFFsmean.xlsx')
            
        #     df_plotdiff=diffdFFsmeanplot.groupby('Group')
        #     means_diff=df_plotdiff.mean()
        #     errors_diff=df_plotdiff.std()
            
        #     #plot figure
        #     fig_diffdFF, axdiff = plt.subplots(1, 1, figsize=(7, 6))
            
        #     means_diff.plot.bar(y=list_behav, yerr=errors_diff[list_behav], rot=0,
        #                         ax=axdiff, capsize=2, colormap='summer')
            
        #     axdiff.set_ylabel('Diff dFF')
        #     axdiff.set_title(f'Diff dFF - {exp} {session} - length{EVENT_TIME_THRESHOLD/10} interbout{THRESH_S} - filtero{ORDER}f{CUT_FREQ}')
        #     #save figure
        #     fig_diffdFF.savefig(repo_path / f'{exp}_{session}_length{EVENT_TIME_THRESHOLD/10}_interbout{THRESH_S}_filtero{ORDER}f{CUT_FREQ}_diffdFFs.png')

