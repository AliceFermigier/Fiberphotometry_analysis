# -*- coding: utf-8 -*-
"""
Created on Fri Jun 18 11:41:03 2021

Plot tracking

@author: Alice Fermigier
"""
#%%
##########
#IMPORTED#
##########

import matplotlib.pyplot as plt
import pandas as pd
import os


#%%
########
#LOADER#
########

from Fiberpho_loader import experiment_path, analysis_path, subjects_df, SAMPLERATE
from Fiberpho_loader import list_EVENT, list_TIMEWINDOW, PRE_EVENT_TIME, TIME_BEGIN

os.chdir(experiment_path)


#%%
########
#SCRIPT#
########

def tracking(df_tracking):

    fig1 = plt.figure(figsize=(20,20))
    ax1 = fig1.add_subplot(111)
    
    p1, = ax1.plot('Item2.X', 'Item2.Y', linewidth=.5, color='black', data=df_tracking)
    
def align_dFFtrack(df_tracking, fiberpho, camera):
    

def track_fiberpho(df_process):
    
    fig2 = plt.figure(figsize=(20,20))
    ax2 = fig2.add_subplot(111)
    
    p2, = ax2.pcolormesh('Item2.X', 'Item2.Y', 'dFFfilt', cmap=cmap, norm=norm, data=df_process)
    
#for plotly animations : https://plotly.com/python/sliders/


#%%Run test

mouse_path = analysis_path / '20211004_OdDis/Test/HFDm3'
mouse = str(mouse_path).split('/')[-1]

fiberpho_path = str(mouse_path) + '/' + mouse + '_dFFfilt.csv'
camera_path = str(mouse_path) + '/' + mouse + '_camera.csv'
rawdata_path = str(mouse_path) + '/' + mouse +  '_rawdata.csv'
    
fiberpho = pd.read_csv(fiberpho_path)
camera = pd.read_csv(camera_path)
#rawdata_df = pd.read_csv(rawdata_path)

print('timevector')
timevector = time_vector(fiberpho, SAMPLERATE)
print('timestamp')
timestart_camera = timestamp_camera(camera)[0]
print('start camera : ', timestart_camera)
print('aligning')
fiberbehav_df = align_behav(behav10Sps, fiberpho, timevector, timestart_camera)
print('processing')
(fiberbehav2_df, fbprocess_df) = behav_process(fiberbehav_df, list_BOI)

#%%Run for all 

#list experiments, sessions and mice:
for exp_path in Path(analysis_path).iterdir():
    if exp_path.is_dir():
        for session_path in Path(exp_path).iterdir():
            if session_path.is_dir():
                for mouse_path in Path(session_path).iterdir():
                    #get the name of the mouse, session and experiment
                    # '/' on mac, '\\' on windows
                    exp = str(mouse_path).split('/')[-3]
                    session = str(mouse_path).split('/')[-2]
                    mouse = str(mouse_path).split('/')[-1]
                    print(exp, session, mouse)
                    
                    #get data
                    tracking_path = str(mouse_path) + '/tracking_' + mouse + '.csv'
                    fiberpho_path = str(mouse_path) + '/' + mouse + '_dFFfilt.csv'
                    camera_path = str(mouse_path) + '/' + mouse + '_camera.csv'
                    
                    if os.path.exists(behav_path):
                        fiberpho = pd.read_csv(fiberpho_path)
                        camera = pd.read_csv(camera_path)
                        tracking = pd.read_csv(tracking_path)
                
                    

