# -*- coding: utf-8 -*-
"""
Created on Fri Jun 18 11:41:03 2021

Plot tracking from Bonsai csv file
Bonsai file used to generate file : tracking_fiberOF.bonsai

@author: Alice Fermigier
"""
#%%
##########
#IMPORTED#
##########

import matplotlib.pyplot as plt
import pandas as pd
import os
import numpy as np


#%%
########
#LOADER#
########

from Fiberpho_loader import experiment_path, analysis_path, data_path, subjects_df, SAMPLERATE, proto_df
from Fiberpho_loader import list_EVENT, list_TIMEWINDOW, PRE_EVENT_TIME, TIME_BEGIN, THRESH_S, EVENT_TIME_THRESHOLD

from Fiberpho_plots import session_code, truncate, time_vector, timestamp_camera, timestamp_camera_fromraw

os.chdir(experiment_path)


#%%
########
#SCRIPT#
########

def tracking(df_tracking):

    fig1 = plt.figure(figsize=(20,20))
    ax1 = fig1.add_subplot(111)
    
    p1, = ax1.plot('Item2.X', 'Item2.Y', linewidth=.5, color='black', data=df_tracking)
    
    plt.savefig(mouse_path / f'{mouse}_tracking.pdf')
    
    return()
    
def align_dFFtrack(df_tracking, fiberpho, timevector, camera_start, camera_stop):
    
    # #trouver samplerate df_tracking en temps avec camera_start - camera_stop
    # samplerate_tracking = len(df_tracking)/(camera_stop-camera_start)
    
    #resample df_tracking
    df_tracking_resampled = df_tracking.asfreq(SAMPLERATE)
    
    #index of df_tracking where camera starts
    list_indstart = np.where(round(timevector,1) == timestart_camera)
    indstart = list_indstart[0].tolist()[0]
    
    #align tracking data
    tracking_x = [0]*indstart
    tracking_x.extend(df_tracking_resampled['Item2.X'])
    
    tracking_y = [0]*indstart
    tracking_y.extend(df_tracking_resampled['Item2.Y'])
    
    # creates list of denoised fiberpho data    
    denoised_fiberpho = fiberpho['Analog In. | Ch.1 470 nm (Deinterleaved)_dF/F0' +
                                 '-Analog In. | Ch.1 405 nm (Deinterleaved)_dF/F0_LowPass'].dropna()
    denoised_fiberpho_list = denoised_fiberpho.tolist()
    
    # creates list of isosbestic (405) and Ca dependent (470) data
    dff_405nm = fiberpho['Analog In. | Ch.1 405 nm (Deinterleaved)_dF/F0_LowPass'].dropna()
    dff_405nm_list = dff_405nm.tolist()
    dff_470nm = fiberpho['Analog In. | Ch.1 470 nm (Deinterleaved)_dF/F0_LowPass'].dropna()
    dff_470nm_list = dff_470nm.tolist()
    
    # makes timevector into a list
    timelist = timevector.tolist()
    
    # crops lists so that all lengths match
    min_length = min([len(timelist), len(denoised_fiberpho_list), len(tracking_x),
                      len(dff_405nm_list), len(dff_470nm_list)])
    timelist = timelist[:min_length]
    denoised_fiberpho_list = denoised_fiberpho_list[:min_length]
    dff_405nm_list = dff_405nm_list[:min_length]
    dff_470nm_list = dff_470nm_list[:min_length]
    tracking_x = tracking_x[:min_length]
    tracking_y = tracking_y[:min_length]

    # create dataframe
    fibertracking_df = pd.DataFrame(data = {'Time(s)':timelist, 'Denoised dFF' : denoised_fiberpho_list,
                                             '405nm deltaF/F' : dff_405nm_list, '470nm deltaF/F' : dff_470nm_list, 
                                             'Track_x' : tracking_x, 'Track_y' : tracking_y})
    
    return(fibertracking_df)

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
                
                    

