#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 07:58:17 2022

Code fiber workshop

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


#%%
###################
#DEFINED FUNCTIONS#
###################

def filter_dFF(fiberbehav_df, ORDER, CUT_FREQ):
    """
    Apply additional filter to dFF data
    """
    fiberpho = fiberbehav_df['Denoised dFF']
    fiber470nm = fiberbehav_df['470nm deltaF/F']
    fiber405nm = fiberbehav_df['405nm deltaF/F']
    samplingrate = 1000/(fiberbehav_df.loc[1000,'Time(s)']-fiberbehav_df.loc[0,'Time(s)'])
    sos = signal.butter(ORDER, CUT_FREQ, btype='high', analog=False, output='sos', fs=samplingrate)
    filtered_data = signal.sosfilt(sos, fiberpho)
    filtered_470nm = signal.sosfilt(sos, fiber470nm)
    filtered_405nm = signal.sosfilt(sos, fiber405nm)
    
    filtered_df = fiberbehav_df
    filtered_df['Denoised dFF'] = filtered_data
    filtered_df['470nm deltaF/F'] = filtered_470nm
    filtered_df['405nm deltaF/F'] = filtered_405nm
    
    return(filtered_df)

def session_code(session):
    """
    Generate session code in file name
    """
    if session in ['Habituation','Training','S1']:
        code = '0'
    if session in ['S2','Test 1h','Test']:
        code = '1'
    elif session in ['S3','Test 24h']:
        code = '2'
        
    return(code)

def truncate(n, decimals=0):
    multiplier = 10 ** decimals
    return int(n * multiplier) / multiplier

def time_vector(fiberpho, SAMPLERATE) :
    """
    Creates timevector on which to plot the data, in pd format
    --> Parameters :
        fiberpho = float, duration of trial in secs
        samplerate = int, in Sps (for processed fiberpho data in Doric Neuroscience Studio, samplerate = 10Sps)
    --> Returns :
        timevector = pd series
    """
    #denoised_fiberpho = fiberpho['Analog In. | Ch.1 470 nm (Deinterleaved)_dF/F0_LowPass' + 
    #                             '-Analog In. | Ch.1 405 nm (Deinterleaved)_dF/F0_LowPass'].dropna()
    #--> if better timevector of the exact same lenght as fiberpho data
    
    duration =  math.ceil(fiberpho.at[len(fiberpho)-2,'Time(s)'])
    #duration = len(fiberpho) #for Roman data
    return pd.Series(np.linspace(0.0, duration, num = int(duration*SAMPLERATE)+1))

def timestamp_camera_fromraw(rawdata_df) :
    """
    Function to extract the timestamps where the camera starts and stops
    --> Parameters
        camera : pd dataframe, camera I/O with sample rate = 12kSps
    --> Returns
        (camera_start, camera_stop) = timestamp when camera starts and stops in seconds (truncated to 0,1s) #camera_stop à enlever si pas besoin
    """
    ind_list = np.where(rawdata_df['DI/O-3'] == 1)
    ind_list = ind_list[0].tolist()
    (ind_start, ind_stop) = (ind_list[0],ind_list[len(ind_list)-1])
    return (truncate(rawdata_df.at[ind_start, 'Time(s)'], 1),
            truncate(rawdata_df.at[ind_stop, 'Time(s)'], 1))

def align_behav(behav10Sps, fiberpho, timevector, timestart_camera):
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
    #scored behaviours in Boris
    list_behav = behav10Sps.columns[1:].tolist()
    list_ind = np.arange(len(list_behav), dtype = int).tolist()
    behav_comp = [0]*len(list_behav)
    
    #index where camera starts
    list_indstart = np.where(round(timevector,1) == timestart_camera)
    indstart = list_indstart[0].tolist()[0]

    # create lists of behaviour data for each scored behaviour
    # aligned with start of the camera
    for (behav, ind) in zip(list_behav, list_ind) :
        behav_comp[ind] = [0]*indstart
        behav_comp[ind].extend(behav10Sps[behav].tolist())
       
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
    min_length = min([len(timelist), len(denoised_fiberpho_list), len(behav_comp[0]),
                      len(dff_405nm_list), len(dff_470nm_list)])
    timelist = timelist[:min_length]
    denoised_fiberpho_list = denoised_fiberpho_list[:min_length]
    dff_405nm_list = dff_405nm_list[:min_length]
    dff_470nm_list = dff_470nm_list[:min_length]
    behav_crop = []
    for behav in behav_comp:
        behav_crop.append(behav[:min_length])
        
    if len(list_behav) == 1 :
        fiberbehav_df = pd.DataFrame(data = {'Time(s)':timelist, 'Denoised dFF' : denoised_fiberpho_list,
                                             '405nm deltaF/F' : dff_405nm_list, '470nm deltaF/F' : dff_470nm_list, 
                                             list_behav[0] : behav_crop[0]})
        
    elif len(list_behav) == 2 :
        fiberbehav_df = pd.DataFrame(data = {'Time(s)':timelist, 'Denoised dFF' : denoised_fiberpho_list,
                                             '405nm deltaF/F' : dff_405nm_list, '470nm deltaF/F' : dff_470nm_list,
                                             list_behav[0] : behav_crop[0], list_behav[1] : behav_crop[1]})
        
    elif len(list_behav) == 3 :
        fiberbehav_df = pd.DataFrame(data = {'Time(s)':timelist, 'Denoised dFF' : denoised_fiberpho_list,
                                             '405nm deltaF/F' : dff_405nm_list, '470nm deltaF/F' : dff_470nm_list,
                                             list_behav[0] : behav_crop[0], list_behav[1] : behav_crop[1],
                                             list_behav[2] : behav_crop[2]})
    
    elif len(list_behav) == 4 :
        fiberbehav_df = pd.DataFrame(data = {'Time(s)':timelist, 'Denoised dFF' : denoised_fiberpho_list,
                                             '405nm deltaF/F' : dff_405nm_list, '470nm deltaF/F' : dff_470nm_list,
                                             list_behav[0] : behav_crop[0], list_behav[1] : behav_crop[1],
                                             list_behav[2] : behav_crop[2], list_behav[3] : behav_crop[3] })
        
    else :
        print('Error : too many behaviours in Boris binary file. Score up to 4 or add line to function')
        
    return(fiberbehav_df)
        
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

def plot_rawdata(rawdata_path):
    """
    Plots raw isosbestic and GCaMP traces
    """
    
    rawdata_df_470 = pd.read_csv(rawdata_path, skiprows=1, usecols=['Time(s)', '470 nm (Deinterleaved)'])
    rawdata_df_405 = pd.read_csv(rawdata_path, skiprows=1, usecols=['Time(s)', '405 nm (Deinterleaved)'])
    
    rawdata_df_470.dropna(inplace=True)
    rawdata_df_405.dropna(inplace=True)
    
    
    fig5 = plt.figure(figsize=(10,6))
    ax7 = fig5.add_subplot(211)

    p1, = ax7.plot('Time(s)', '470 nm (Deinterleaved)', 
                   linewidth=1, color='deepskyblue', label='GCaMP', data = rawdata_df_470) 
    p2, = ax7.plot('Time(s)', '405 nm (Deinterleaved)', 
                   linewidth=1, color='blueviolet', label='ISOS', data = rawdata_df_405)
    
    ax7.set_ylabel('V')
    ax7.set_xlabel('Seconds')
    ax7.legend(handles=[p1,p2], loc='upper right')
    ax7.margins(0.01,0.3)
    ax7.set_title(f'GCaMP and Isosbestic raw traces - {exp} {session} {mouse}')
    
    #save figure
    fig5.savefig(mouse_path / f'{mouse}_rawdata.pdf')
    
    return

def plot_fiberpho(fiberbehav_df):
    """
    Plots isosbestic and Ca dependent deltaF/F
    """
    
    #crops starting artifacts
    fiberbehavsnip_df = fiberbehav_df[fiberbehav_df['Time(s)'] > 40]
       
    fig1 = plt.figure(figsize=(20,6))
    ax0 = fig1.add_subplot(111)
    
    p1, = ax0.plot('Time(s)', '470nm deltaF/F', linewidth=1, color='deepskyblue', label='GCaMP', data = fiberbehavsnip_df) 
    p2, = ax0.plot('Time(s)', '405nm deltaF/F', linewidth=1, color='blueviolet', label='ISOS', data = fiberbehavsnip_df)
    
    ax0.set_ylabel(r'$\Delta$F/F')
    ax0.set_xlabel('Seconds')
    ax0.legend(handles=[p1,p2], loc='upper right')
    ax0.margins(0.01,0.2)
    ax0.set_title(f'GCaMP and Isosbestic dFF - {exp} {session} {mouse}')
    
    #save figure
    fig1.savefig(mouse_path / f'{mouse}_GCaMP_ISOS.pdf')
    
    return

def plot_fiberpho_behav(behavprocess_df):
    """
    Plots denoised deltaF/F aligned with behaviour (includes baseline)
    """
    
    behavprocesssnip_df = behavprocess_df[behavprocess_df['Time(s)'] > 10]
    
    fig2 = plt.figure(figsize=(20,4))
    ax1 = fig2.add_subplot(111)
    
    if session in ['Test 1h','Test 24h','Test','S3']:
        #plots fiberpho trace and behaviour
        p1, = ax1.plot('Time(s)', 'Denoised dFF', linewidth=1, color='black', label='_GCaMP', data = behavprocesssnip_df)
        
        #makes areas corresponding to behaviours
        i = 0
        j = 0
        k = 0
        for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Exploration familiar'].tolist()):
            if y == 1:
                x_start = x
            if y == -1 and x_start!=0:
                ax1.axvspan(x_start, x, facecolor='gold', alpha=0.3, label = '_'*i + 'Exploration familiar')
                x_start=0
                i+=1
        for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Exploration new'].tolist()):
            if y == 1:
                x_start = x
            if y == -1 and x_start!=0:
                ax1.axvspan(x_start, x, facecolor='purple', alpha=0.3, label = '_'*j + 'Exploration new')
                x_start=0
                j+=1
        for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Climbing'].tolist()):
            if y == 1:
                x_start = x
            if y == -1 and x_start!=0:
                ax1.axvspan(x_start, x, facecolor='cornflowerblue', alpha=0.3, label = '_'*k + 'Climbing')
                x_start=0
                k+=1
            
        #makes vertical line for entry in arena
        x_entry = fiberbehav_df.at[int(np.where(fiberbehav_df['Entry'] == 1)[0][0]), 'Time(s)']
        ax1.axvline(x_entry, color='slategrey', ls = '--', label = 'Entry' )
        
    ###########################################################
    ##EXERCISE : write a code to plot homecage / new arena data
    ###########################################################
    
    if session in ['S1']:
        #plots fiberpho trace and behaviour
        p1, = ax1.plot('Time(s)', 'Denoised dFF', linewidth=1, color='black', label='_GCaMP', data = behavprocesssnip_df)
        
        #makes areas corresponding to behaviours
        i = 0
        j = 0
        x_start = 0.1
        for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Homecage'].tolist()):
            if y == 1:
                x_start = x
            if y == -1 and x_start!=0:
                ax1.axvspan(x_start, x, facecolor='gold', alpha=0.3, label = '_'*i + 'Homecage')
                x_start=0
                i+=1
        for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['New arena'].tolist()):
            if y == 1:
                x_start = x
            if y == -1 and x_start!=0:
                ax1.axvspan(x_start, x, facecolor='purple', alpha=0.3, label = '_'*j + 'New arena')
                x_start=0
                j+=1
           
    ###########################################################
    
    ax1.set_ylabel(r'$\Delta$F/F')
    ax1.set_xlabel('Seconds')
    ax1.set_ylim(-2,3)
    ax1.legend(loc = 'upper right')
    ax1.margins(0.01,0.03)
    ax1.set_title(f'dFF with Behavioural Scoring - {exp} {session} {mouse} - interbout {THRESH_S} - filtero{ORDER}f{CUT_FREQ}')
    
    #save figure
    fig2.savefig(mouse_path / f'{mouse}_threshold{THRESH_S}s_filtero{ORDER}f{CUT_FREQ}_fiberbehav_scaled.pdf')
    
    return

def meandFF_behav(list_BOI, fiberbehav_df):
    """
    Calculates mean dFF during each behaviour
    """
    list_behav_analyzed = []
    list_meandFF = []
    list_stddFF = []
    
    #create list of behaviours and list of correspondind mean dFFs
    for behav in list_BOI:
        if fiberbehav_df[behav].sum() > 2:
            list_behav_analyzed.append(behav)
            meandFF_behav_df = fiberbehav_df.groupby([behav], as_index=False).mean()
            stddFF_behav_df = fiberbehav_df.groupby([behav], as_index=False).std()
            print(meandFF_behav_df)
            list_meandFF.append(meandFF_behav_df.loc[meandFF_behav_df[behav]==1, 'Denoised dFF'].values[0])
            list_stddFF.append(stddFF_behav_df.loc[stddFF_behav_df[behav]==1, 'Denoised dFF'].values[0])
         
    if 'Entry' in list_BOI:
        #get index of when the animal enters
        ind_start_trial = fiberbehav_df.index[fiberbehav_df['Entry'] == 1].tolist()[0]
                                    
        #calculate mean dFF during baseline
        
        meandFF_baseline = fiberbehav_df.loc[:ind_start_trial, 'Denoised dFF'].mean()
        stddFF_baseline = fiberbehav_df.loc[:ind_start_trial, 'Denoised dFF'].std()
        
        #calculate mean dFF when no exploration and after gate opens
        
        meandFF_postbaseline_df = fiberbehav_df.loc[ind_start_trial:]
        stddFF_postbaseline_df = fiberbehav_df.loc[ind_start_trial:]
        for behav in list_behav_analyzed:
            meandFF_postbaseline_df = meandFF_postbaseline_df.loc[meandFF_postbaseline_df[behav]==0]
            stddFF_postbaseline_df = stddFF_postbaseline_df.loc[stddFF_postbaseline_df[behav]==0]
            
        meandFF_postbaseline = meandFF_postbaseline_df['Denoised dFF'].mean()
        stddFF_postbaseline = stddFF_postbaseline_df['Denoised dFF'].std()
        
        #calculate mean dFF when no exploration on total trial
        meandFF_df = fiberbehav_df
        stddFF_df = fiberbehav_df
        for behav in list_behav_analyzed:
            meandFF_df = meandFF_df.loc[meandFF_df[behav]==0]
            stddFF_df = stddFF_df.loc[stddFF_df[behav]==0]
            
        meandFF = meandFF_df['Denoised dFF'].mean()
        stddFF = stddFF_df['Denoised dFF'].std()
    
        #create dataframe with values
    
        list_meandFFs = [mouse, meandFF,meandFF_baseline,meandFF_postbaseline]
        list_stddFFs = [mouse, stddFF,stddFF_baseline,stddFF_postbaseline]
        list_columns = ['Subject','Mean dFF','Baseline', 'Post_baseline']
    else:
        list_meandFFs = [mouse]
        list_stddFFs = [mouse]
        list_columns = ['Subject']
    for (behav,meandFF) in zip(list_behav_analyzed,list_meandFF):
        list_meandFFs.append(meandFF)
        list_columns.append(behav)
    for stddFF in list_stddFF:
        list_stddFFs.append(stddFF)
        
    meandFFs_df = pd.DataFrame(data=[list_meandFFs,list_stddFFs], columns=list_columns)
    print(meandFFs_df)
        
    #Create figure of mean dFFs
    fig7 = plt.figure(figsize=(7,6))
    ax71 = fig7.add_subplot(111)
    print('x = ',meandFFs_df.columns.to_list()[1:])
    print('y = ',meandFFs_df.iloc[0, 1:].to_list())
    p71 = ax71.bar(meandFFs_df.columns.to_list()[1:], meandFFs_df.iloc[0, 1:], yerr=meandFFs_df.iloc[1, 1:])
    ax71.axhline(y=0, linewidth=.6, color='black')
    ax71.set_ylabel(r'$\Delta$F/F')
    ax71.set_title(f'Mean dFF - {exp} {session} {mouse}')
    plt.xticks(rotation=35)
    
    #save figure
    fig7.savefig(mouse_path / f'{mouse}_meandFFs.pdf')
    
    meandFFs_df.to_excel(mouse_path / f'{mouse}_meandFF.xlsx')
    
    return(meandFFs_df)

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
    EVENT_TIME_THRESHOLDms = EVENT_TIME_THRESHOLD*10
    for (start,end) in zip(list_ind_event_o, list_ind_event_w):
        #print(start,end)
        if end-start<EVENT_TIME_THRESHOLDms and end-start>1:
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

def plot_PETH(PETH_data, BOI, event, timewindow):
    """
    Plots PETH average and heatmap
    """
    PRE_TIME = timewindow[0]
    POST_TIME = timewindow[1]
    
    #create figure
    fig4 = plt.figure(figsize=(6,10))
    
    #create time vector
    peri_time = np.arange(-PRE_TIME, POST_TIME+0.1, 0.1)
    
    #calculate mean dFF and std error
    mean_dFF_snips = np.mean(PETH_data, axis=0)
    std_dFF_snips = np.std(PETH_data, axis=0)
        
    #plot individual traces and mean
    ax5 = fig4.add_subplot(212)
    for snip in PETH_data:
        p1, = ax5.plot(peri_time, snip, linewidth=.5,
                       color=[.7, .7, .7], label='Individual trials')
    p2, = ax5.plot(peri_time, mean_dFF_snips, linewidth=2,
                   color='green', label='Mean response')
    
    #plot standard error bars
    p3 = ax5.fill_between(peri_time, mean_dFF_snips+std_dFF_snips,
                      mean_dFF_snips-std_dFF_snips, facecolor='green', alpha=0.2)
    p4 = ax5.axvline(x=0, linewidth=2, color='slategray', ls = '--', label='Exploration '+event)
    
    #ax5.axis('tight')
    ax5.set_xlabel('Seconds')
    ax5.set_ylabel('z-scored $\Delta$F/F')
    ax5.legend(handles=[p1, p2, p4], loc='upper left', fontsize = 'small')
    ax5.margins(0,0.01)
    
    #add heatmap
    ax6 = fig4.add_subplot(211)
    cs = ax6.imshow(PETH_data, cmap='magma', aspect = 'auto',
                    interpolation='none', extent=[-PRE_TIME, POST_TIME, len(PETH_data), 0],
                    vmin = -6, vmax = 9)
    ax6.set_ylabel('Bout Number')
    ax6.set_yticks(np.arange(.5, len(PETH_data), 2))
    ax6.set_yticklabels(np.arange(0, len(PETH_data), 2))
    ax6.axvline(x=0, linewidth=2, color='black', ls = '--')
    ax6.set_title(BOI+' '+event+' - '+exp+' '+session+' '+mouse)
    
    fig4.subplots_adjust(right=0.8, hspace = 0.1)
    cbar_ax = fig4.add_axes([0.85, 0.54, 0.02, 0.34])
    fig4.colorbar(cs, cax=cbar_ax)
    
    #save figure
    fig4.savefig(mouse_path / f'{mouse}{BOI}_PETH{event[0]}.pdf')
    
    return

#%%
########
#LOADER#
########

SAMPLERATE = 10 #in Hz

#list of behaviours on which to do peri event time histograms (PETH)
list_EVENT = ['onset', 'withdrawal']
list_TIMEWINDOW = [[4,5],[4,7]] 

#time before behaviour for calculation of PETH baseline, in seconds
PRE_EVENT_TIME = 1

#time to crop at the beginning of the trial to remove artifact, in seconds
TIME_BEGIN = 60

#threshold to fuse behaviour if bouts are too close, in secs
THRESH_S = 0

#threshold for PETH : if events are too short, not in PETH
EVENT_TIME_THRESHOLD = 0

#filter caracteristics
ORDER = 3
CUT_FREQ = 1 #in Hz

EVENT_TIME_THRESHOLD = 0
THRESH_S = 0

#%%
########
#SCRIPT#
########

#set the general analysis path
experiment_path = Path('/Users/alice/Desktop/FiberWorkshop') #Path where you stored all your fiber data and analysis files
analysis_path = experiment_path / 'Analysis'
data_path = experiment_path / 'Data'
os.chdir(experiment_path)
os.getcwd()

#set the nomenclature
exp = 'SRM'
session = 'S3'
mouse = 'HFDf2'

#get data files
data_path_exp = data_path / exp
code = session_code(session)

#get analysis paths (where you are going to store your analysed files)
exp_path = analysis_path / exp
session_path = exp_path / session
repo_path = session_path / f'length{EVENT_TIME_THRESHOLD/10}_interbout{THRESH_S}_o{ORDER}f{CUT_FREQ}'
if not os.path.exists(repo_path):
    os.mkdir(repo_path)
    os.mkdir(repo_path / mouse)
mouse_path = repo_path / mouse
    
#get data
rawdata_path = data_path_exp / f'{mouse}_{code}.csv'
behav_path = data_path_exp / f'{exp}_{code}_{mouse}.csv'
fiberpho_path = data_path_exp / f'{mouse}_{code}_dFFfilt.csv'

#get pandas dataframes
rawdata_df = pd.read_csv(rawdata_path, skiprows=1)
behav10Sps = pd.read_csv(behav_path)
fiberpho = pd.read_csv(fiberpho_path)

#list of behaviours to analyze
list_BOI = behav10Sps.columns[1:].tolist()
print(list_BOI)

#PROCESSING
###########

print('timevector')
timevector = time_vector(fiberpho, SAMPLERATE)

print('timestamp')
timestart_camera = timestamp_camera_fromraw(rawdata_df)[0]
print('start camera : ', timestart_camera)

print('aligning')
fiberbehav_df = align_behav(behav10Sps, fiberpho, timevector, timestart_camera)

print('filtering')
fiberbehav_df = filter_dFF(fiberbehav_df, ORDER, 0.02)

print('processing')
(fiberbehav2_df, fbprocess_df) = behav_process(fiberbehav_df, list_BOI)


#PLOTTING
#########

#plot raw data
plot_rawdata(rawdata_path)

#plot dFFs
plot_fiberpho(fiberbehav_df)

#plot dFF with behaviour
plot_fiberpho_behav(fbprocess_df)

#MEAN
#####
mean_dFF_df = meandFF_behav(list_BOI, fiberbehav_df)

#PETH (PERI-EVENT TIME HISTOGRAM)
#################################
if exp == 'SRM':
    for BOI in list_BOI:
        for event, timewindow in zip(list_EVENT,list_TIMEWINDOW):
            PETH_data = PETH(fbprocess_df, BOI, event, timewindow)
            plot_PETH(PETH_data, BOI, event, timewindow)




