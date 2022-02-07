#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  4 11:05:52 2022

@author: alice
"""
#%%
##########
#IMPORTED#
##########

import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy import signal
plt.rcParams.update({'figure.max_open_warning' : 0})
from pathlib import Path
import os
import sys


#Custom
#put path to directory where python files are stored
#sys.path.append('/Users/alice/Documents/GitHub/Fiberphotometry_analysis')
sys.path.append('C:\\Users\\afermigier\\Documents\\GitHub\\Fiberphotometry_analysis')

#%%
########
#LOADER#
########

from Fiberpho_loader import experiment_path, analysis_path, data_path, subjects_df, SAMPLERATE, proto_df
from Fiberpho_loader import list_EVENT, list_TIMEWINDOW, PRE_EVENT_TIME, TIME_BEGIN, ORDER #THRESH_S, EVENT_TIME_THRESHOLD, CUT_FREQ

from Func_fiberplots import session_code, truncate, time_vector, timestamp_camera, timestamp_camera_fromraw

os.chdir(experiment_path)

#%%
###################
#DEFINED FUNCTIONS#
###################
def align_behav(behav10Sps, fiberpho, timestart_camera):
    """
    Aligns fiberpho data with behavioural data from Boris on timevector
    --> Parameters 
        behav10Sps : pd df, binary behavioural data from Boris, with sample rate = 10Sps
        fiberpho : pd df, pre-processed fiber photometry data (deltaF/F 470nm - 405nm)
            to be aligned with behav data
        timevector : pd df, timevector to plot data
        timestart_camera : int, index of when camera starts
    --> Returns
        fiberbehav_df : pd df of aligned data

    """
    #change names in behaviour file
    behav10Sps.rename(columns={"time": "Time(s)", "e": 'Animal entry',
                               "f" : 'Running', 's': 'Feeding', 'z': 'Door opening'}, 
                          inplace= True)
    #scored behaviours in Boris
    list_behav = behav10Sps.columns[1:].tolist()
    behav_comp = [0]*len(list_behav)
    print(list_behav)
    
    #index where camera starts
    indstart = int(timestart_camera*10)

    # create lists of behaviour data for each scored behaviour
    # aligned with start of the camera
    for (ind, behav) in enumerate(list_behav) :
        behav_comp[ind] = [0]*(indstart-1)
        behav_comp[ind].extend(behav10Sps[behav].tolist()[:len(fiberpho)-indstart+1])
    
        
    if len(list_behav) == 1 :
        fiberbehav_df = pd.DataFrame(data = {'Time(s)':fiberpho['Time(s)'], 'Denoised dFF' : fiberpho['Denoised dFF'],
                                             list_behav[0] : behav_comp[0]})
        
    elif len(list_behav) == 2 :
        fiberbehav_df = pd.DataFrame(data = {'Time(s)':fiberpho['Time(s)'], 'Denoised dFF' : fiberpho['Denoised dFF'],
                                             list_behav[0] : behav_comp[0], list_behav[1] : behav_comp[1]})
        
    elif len(list_behav) == 3 :
        fiberbehav_df = pd.DataFrame(data = {'Time(s)':fiberpho['Time(s)'], 'Denoised dFF' : fiberpho['Denoised dFF'],
                                             list_behav[0] : behav_comp[0], list_behav[1] : behav_comp[1],
                                             list_behav[2] : behav_comp[2]})
    
    elif len(list_behav) == 4 :
        fiberbehav_df = pd.DataFrame(data = {'Time(s)':fiberpho['Time(s)'], 'Denoised dFF' : fiberpho['Denoised dFF'],
                                             list_behav[0] : behav_comp[0], list_behav[1] : behav_comp[1],
                                             list_behav[2] : behav_comp[2], list_behav[3] : behav_comp[3] })
        
    else :
        print('Error : too many behaviours in Boris binary file. Score up to 4 or add line to function')
    
    return(fiberbehav_df)

def filter_dFF(fiberbehav_df, ORDER, CUT_FREQ):
    """
    Apply additional filter to dFF data
    """
    fiberpho = fiberbehav_df['Denoised dFF']
    samplingrate = 1000/(fiberbehav_df.loc[1000,'Time(s)']-fiberbehav_df.loc[0,'Time(s)'])
    sos = signal.butter(ORDER, CUT_FREQ, btype='low', analog=False, output='sos', fs=samplingrate)
    filtered_data = signal.sosfilt(sos, fiberpho)
    
    filtered_df = fiberbehav_df
    filtered_df['Denoised dFF'] = filtered_data
    
    return(filtered_df)

def behav_process(fiberbehav_df, list_BOI):

    #1 fuse explorations that are too close
    
    THRESHOLD = THRESH_S*10 #set threshold in 0.1seconds (because samplerate = 10Sps)
    
    for BOI in list_BOI:
        i_1 = 0
        count = 0
        for (ind,i) in zip(fiberbehav_df.index, fiberbehav_df[BOI]):
            if i==i_1 and i==0:
                count += 1
            elif i!=i_1 and i==0:
                count = 1
            elif i!=i_1 and i==1:
                if count <= THRESHOLD:
                    #print(fiberbehav_df.loc[ind-count:ind-1, BOI])
                    #print([1]*count)
                    fiberbehav_df.loc[ind-count:ind-1, BOI] = [1]*count
            i_1 = i
                
    #2 calculate the derivative of behav of interest and put in new df
    # that way, it will show 1 when behaviour starts and -1 when it stops
    
    if len(list_BOI) == 1:
        behavprocess_df = pd.DataFrame(data = {'Time(s)': fiberbehav_df['Time(s)'], 
                                               'Denoised dFF' : fiberbehav_df['Denoised dFF'],
                                               list_BOI[0] : fiberbehav_df[list_BOI[0]].diff()})
    elif len(list_BOI) == 2:
        behavprocess_df = pd.DataFrame(data = {'Time(s)': fiberbehav_df['Time(s)'], 
                                               'Denoised dFF' : fiberbehav_df['Denoised dFF'],
                                               list_BOI[0] : fiberbehav_df[list_BOI[0]].diff(),
                                               list_BOI[1] : fiberbehav_df[list_BOI[1]].diff()})
    elif len(list_BOI) == 3:
        behavprocess_df = pd.DataFrame(data = {'Time(s)': fiberbehav_df['Time(s)'], 
                                               'Denoised dFF' : fiberbehav_df['Denoised dFF'],
                                               list_BOI[0] : fiberbehav_df[list_BOI[0]].diff(),
                                               list_BOI[1] : fiberbehav_df[list_BOI[1]].diff(),
                                               list_BOI[2] : fiberbehav_df[list_BOI[2]].diff()})
    elif len(list_BOI) == 4:
        behavprocess_df = pd.DataFrame(data = {'Time(s)': fiberbehav_df['Time(s)'], 
                                               'Denoised dFF' : fiberbehav_df['Denoised dFF'],
                                               list_BOI[0] : fiberbehav_df[list_BOI[0]].diff(),
                                               list_BOI[1] : fiberbehav_df[list_BOI[1]].diff(),
                                               list_BOI[2] : fiberbehav_df[list_BOI[2]].diff(),
                                               list_BOI[3] : fiberbehav_df[list_BOI[3]].diff()})
        
    behavprocess_df.to_excel(mouse_path / f'{mouse}_fbprocess.xlsx')

    
    return(fiberbehav_df, behavprocess_df)

def plot_fiberpho_behav(behavprocess_df, fiberbehav_df):
    """
    Plots denoised deltaF/F aligned with behaviour (includes baseline)
    """
    
    behavprocesssnip_df = behavprocess_df[behavprocess_df['Time(s)'] > 40]
    
    fig2 = plt.figure(figsize=(20,4))
    ax1 = fig2.add_subplot(111)
    
    #plots fiberpho trace and behaviour
    p1, = ax1.plot('Time(s)', 'Denoised dFF', linewidth=.6, color='black', label='_GCaMP', data = behavprocesssnip_df)
    # p2, = ax1.plot(behavprocesssnip_df['Time(s)'], behavprocesssnip_df['Exploration fam'], linewidth=1, 
    #                color='moccasin', label='Exploration fam')
    
    #makes areas corresponding to behaviours
    i = 0
    j = 0
    for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Running'].tolist()):
        if y == 1:
            x_start = x
        if y == -1 and x_start!=0:
            ax1.axvspan(x_start, x, facecolor='gold', alpha=0.5, label = '_'*i + 'Running')
            x_start=0
            i+=1
    for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Feeding'].tolist()):
        if y == 1:
            x_start = x
        if y == -1 and x_start!=0:
            ax1.axvspan(x_start, x, facecolor='purple', alpha=0.5, label = '_'*j + 'Feeding')
            x_start=0
            j+=1

    #makes vertical line for entry opening of gate
    k = 0
    for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Door opening'].tolist()):
        if y == 1:
            x_door = x
            ax1.axvline(x_door, color='lightsteelblue', ls = '--', label = '_'*k + 'Door opening' )
            k += 1
            
    #makes vertical line for entry in open field
    l = 0
    for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Animal entry'].tolist()):
        if y == 1:
            x_entry = x
            ax1.axvline(x_entry, color='slategrey', ls = '--', label = '_'*l + 'Animal entry' )
            l += 1
    
    ax1.set_ylabel(r'$\Delta$F/F')
    ax1.set_xlabel('Seconds')
    ax1.set_ylim(-2,3)
    ax1.legend(loc = 'upper right')
    ax1.margins(0.01,0.03)
    ax1.set_title(f'dFF with Behavioural Scoring - {mouse} - interbout {THRESH_S}')
    
    #save figure
    fig2.savefig(mouse_path / f'{mouse}_threshold{THRESH_S}s_fiberbehav_scaled.pdf')
    
    return

def PETH(behavprocess_df, BOI, event, timewindow):
    """
    Creates dataframe of fiberpho data centered on bout event for BOI
    --> Parameters
        behavprocess_df : pd dataframe, aligned fiberpho and behav data for 1 mouse
        BOI : str, behaviour of interest (must have the same name as behavprocess column name)
        event : str, onset or withdrawal
        timewindow : list, time pre and post event
    --> Returns
        PETHo_array : np array, normalized fiberpho data centered on event
    """
    print(f'PETH {BOI}')
    #set time window relative to event
    PRE_TIME = timewindow[0]
    POST_TIME = timewindow[1]
    
    #if behaviour not recognized
    if BOI not in behavprocess_df.columns[2:].tolist() :
        print('BOI not recognized')
        
    list_ind_event_o = np.where(behavprocess_df[BOI] == 1)[0].tolist()
    list_ind_event_w = np.where(behavprocess_df[BOI] == -1)[0].tolist()
    
    #if event too short, don't process it
    for (start,end) in zip(list_ind_event_o, list_ind_event_w):
        #print(start,end)
        if end-start<EVENT_TIME_THRESHOLD and end-start>1:
            #print(end-start)
            list_ind_event_o.remove(start)
            list_ind_event_w.remove(end)
    
    #creates array of fiberpho data centered on event
    if event == 'onset':
        list_ind_event = list_ind_event_o
    elif event == 'withdrawal':
        list_ind_event = list_ind_event_w
        
    #if event happens too late in dataframe, don't process it
    if list_ind_event[-1]+POST_TIME*SAMPLERATE >= len(behavprocess_df):
        list_ind_event.pop(-1)
        
    list_i = np.arange(len(list_ind_event), dtype = int).tolist()
    PETH_array = np.zeros((len(list_ind_event),(POST_TIME+PRE_TIME)*SAMPLERATE+1))
    for (ind_event, i) in zip(list_ind_event, list_i) :
        #calculates baseline F0 on time window before event (from PRE_TIME to PRE_EVENT_TIME)
        F0 = np.mean(behavprocess_df.loc[ind_event-PRE_TIME*SAMPLERATE:ind_event-PRE_EVENT_TIME*SAMPLERATE, 'Denoised dFF'])
        std0 = np.std(behavprocess_df.loc[ind_event-PRE_TIME*SAMPLERATE:ind_event-PRE_EVENT_TIME*SAMPLERATE, 'Denoised dFF'])
        #creates array of z-scored dFF : z = (dFF-meandFF_baseline)/stddFF_baseline
        if len(behavprocess_df.loc[ind_event-PRE_TIME*SAMPLERATE:ind_event+(POST_TIME*SAMPLERATE), 'Denoised dFF']) == len(PETH_array[0]):
            PETH_array[i] = (behavprocess_df.loc[ind_event-PRE_TIME*SAMPLERATE:ind_event+(POST_TIME*SAMPLERATE), 'Denoised dFF']-F0)/std0
    
    return(PETH_array)

#%%

for BOI in list_BOI:
    # if BOI == 'Entry in arena' or BOI == 'Gate opens':
    #     PETH_data = PETH(behavprocess_df, BOI, 'onset', [6,10])
    #     plot_PETH_average(PETH_data, BOI, 'onset', [6,10])
    # else:
    for event, timewindow in zip(list_EVENT,list_TIMEWINDOW): #CAUTION : adapt function!!
        PETH_data = PETH(behavprocess_df, BOI, event, timewindow)
        plot_PETH(PETH_data, BOI, event, timewindow)
    
            