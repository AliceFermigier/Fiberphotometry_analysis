 # -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 09:30:55 2021

Analysis of Doric fiber photometry data 

DeltaF F data for 470nm and 405nm are calculated using Doric Neuroscience Studio
Behavioral scoring is done using Boris software

@author: Alice Fermigier
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
#sys.path.append('/Users/alice/Documents/GitHub/Fiberphotometry_analysis')
sys.path.append('D:\Profil\Documents\GitHub\Fiberphotometry_analysis')

#%%
########
#LOADER#
########

from Fiberpho_loader import experiment_path, analysis_path, data_path, subjects_df, SAMPLERATE, proto_df
from Fiberpho_loader import list_EVENT, list_TIMEWINDOW, PRE_EVENT_TIME, TIME_BEGIN, ORDER, THRESH_S, EVENT_TIME_THRESHOLD, CUT_FREQ

from Func_fiberplots import session_code, truncate, time_vector, timestamp_camera, timestamp_camera_fromraw, filter_dFF
from Func_fiberplots import meandFF_behav, diff_dFF

os.chdir(experiment_path)

#%%
###################
#DEFINED FUNCTIONS#
###################
    
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
    
    elif len(list_behav) == 5 :
        fiberbehav_df = pd.DataFrame(data = {'Time(s)':timelist, 'Denoised dFF' : denoised_fiberpho_list,
                                             '405nm deltaF/F' : dff_405nm_list, '470nm deltaF/F' : dff_470nm_list,
                                             list_behav[0] : behav_crop[0], list_behav[1] : behav_crop[1],
                                             list_behav[2] : behav_crop[2], list_behav[3] : behav_crop[3],
                                             list_behav[4] : behav_crop[4]})
    else :
        print('Error : too many behaviours in Boris binary file. Score up to 5 or add line to function')
        
    if exp=='NewContext':
        fiberbehav_df.loc[:405,'Homecage']=0
    
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
        
    elif len(list_BOI) == 5:
        behavprocess_df = pd.DataFrame(data = {'Time(s)': fiberbehav_df['Time(s)'], 
                                               'Denoised dFF' : fiberbehav_df['Denoised dFF'],
                                               list_BOI[0] : fiberbehav_df[list_BOI[0]].diff(),
                                               list_BOI[1] : fiberbehav_df[list_BOI[1]].diff(),
                                               list_BOI[2] : fiberbehav_df[list_BOI[2]].diff(),
                                               list_BOI[3] : fiberbehav_df[list_BOI[3]].diff(),
                                               list_BOI[4] : fiberbehav_df[list_BOI[4]].diff()})
        
    behavprocess_df.to_excel(mouse_path / f'{mouse}_fbprocess.xlsx')

    return(fiberbehav_df, behavprocess_df)
    
def plot_rawdata(rawdata_df):
    """
    Plots raw isosbestic and GCaMP traces
    """
    
    fig5 = plt.figure(figsize=(10,6))
    ax7 = fig5.add_subplot(211)

    p1, = ax7.plot('Time(s)', 'Analog In. | Ch.1 470 nm (Deinterleaved)', 
                   linewidth=1, color='deepskyblue', label='GCaMP', data = rawdata_df) 
    p2, = ax7.plot('Time(s)', 'Analog In. | Ch.1 405 nm (Deinterleaved)', 
                   linewidth=1, color='blueviolet', label='ISOS', data = rawdata_df)
    
    ax7.set_ylabel('V')
    ax7.set_xlabel('Time(s)')
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
    
    behavprocesssnip_df = behavprocess_df[behavprocess_df['Time(s)'] > 40]
    
    fig2 = plt.figure(figsize=(20,5))
    ax1 = fig2.add_subplot(111)
    
    if exp=='Fear':
        if session =='Test':
            #plots fiberpho trace and behaviourb 
            p1, = ax1.plot('Time(s)', 'Denoised dFF', linewidth=.6, color='black', label='_GCaMP', data = behavprocesssnip_df)
    
            #makes areas corresponding to behaviours
            i = 0
            j = 0
            for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Homecage'].tolist()):
                if y == 1:
                    x_start = x
                if y == -1 and x_start!=0:
                    ax1.axvspan(x_start, x, facecolor='gold', alpha=0.3, label = '_'*i + 'Homecage')
                    x_start=0
                    i+=1
            for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Fear cage'].tolist()):
                if y == 1:
                    x_start = x
                if y == -1 and x_start!=0:
                    ax1.axvspan(x_start, x, facecolor='blue', alpha=0.3, label = '_'*j + 'Fear cage')
                    x_start=0
                    j+=1
            
        elif session =='Conditioning':
            #plots fiberpho trace and behaviourb 
            p1, = ax1.plot('Time(s)', 'Denoised dFF', linewidth=.6, color='black', label='_GCaMP', data = behavprocesssnip_df)
            x_entry = fiberbehav_df.at[int(np.where(fiberbehav_df['Shock'] == 1)[0][0]), 'Time(s)']
            ax1.axvline(x_entry, color='yellow', ls = '-', lw=2, label = 'Shock' )
            
    elif exp=='NewContext':
        #plots fiberpho trace and behaviourb 
        p1, = ax1.plot('Time(s)', 'Denoised dFF', linewidth=.6, color='black', label='_GCaMP', data = behavprocesssnip_df)

        #makes areas corresponding to behaviours
        i = 0
        j = 0
        for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Homecage'].tolist()):
            if y == 1:
                x_start = x
            if y == -1 and x_start!=0:
                ax1.axvspan(x_start, x, facecolor='gold', alpha=0.3, label = '_'*i + 'Homecage')
                x_start=0
                i+=1
        for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['New context'].tolist()):
            if y == 1:
                x_start = x
            if y == -1 and x_start!=0:
                ax1.axvspan(x_start, x, facecolor='darkturquoise', alpha=0.3, label = '_'*j + 'New context')
                x_start=0
                j+=1

    else:
        if 'Essai' in exp:
            p1, = ax1.plot('Time(s)', 'Denoised dFF', linewidth=1, color='black', label='_GCaMP', data = behavprocesssnip_df)
            
        elif session in ['Test 1h','Test 24h','Test','S3']:
            #plots fiberpho trace and behaviour
            p1, = ax1.plot('Time(s)', 'Denoised dFF', linewidth=1, color='black', label='_GCaMP', data = behavprocesssnip_df)
            # p2, = ax1.plot(behavprocesssnip_df['Time(s)'], behavprocesssnip_df['Exploration fam'], linewidth=1, 
            #                color='moccasin', label='Exploration fam')
            
            #makes areas corresponding to behaviours
            i = 0
            j = 0
            m = 0
            n = 0
            if 'Exploration fam' in list_BOI:
                for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Exploration fam'].tolist()):
                    if y == 1:
                        x_start = x
                    if y == -1 and x_start!=0:
                        ax1.axvspan(x_start, x, facecolor='gold', alpha=0.3, label = '_'*i + 'Exploration fam')
                        x_start=0
                        i+=1
            for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Exploration new'].tolist()):
                if y == 1:
                    x_start = x
                if y == -1 and x_start!=0:
                    ax1.axvspan(x_start, x, facecolor='purple', alpha=0.3, label = '_'*j + 'Exploration new')
                    x_start=0
                    j+=1
            if 'Climbing' in list_BOI:
                for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Climbing'].tolist()):
                    if y == 1:
                        x_start = x
                    if y == -1 and x_start!=0:
                        ax1.axvspan(x_start, x, facecolor='cornflowerblue', alpha=0.3, label = '_'*m + 'Climbing')
                        x_start=0
                        m+=1
            if 'Rearing' in list_BOI:
                for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Rearing'].tolist()):
                    if y == 1:
                        x_start = x
                    if y == -1 and x_start!=0:
                        ax1.axvspan(x_start, x, facecolor='mediumpurple', alpha=0.3, label = '_'*n + 'Rearing')
                        x_start=0
                        n+=1
            
        elif session in ['Habituation','Training','S1']:
            #plots fiberpho trace and behaviourb 
            p1, = ax1.plot('Time(s)', 'Denoised dFF', linewidth=.6, color='black', label='_GCaMP', data = behavprocesssnip_df)
    
            #makes areas corresponding to behaviours
            i = 0
            j = 0
            m = 0
            n = 0
            for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Exploration left'].tolist()):
                if y == 1:
                    x_start = x
                if y == -1 and x_start!=0:
                    ax1.axvspan(x_start, x, facecolor='orange', alpha=0.3, label = '_'*i + 'Exploration left')
                    x_start=0
                    i+=1
            for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Exploration right'].tolist()):
                if y == 1:
                    x_start = x
                if y == -1 and x_start!=0:
                    ax1.axvspan(x_start, x, facecolor='darkturquoise', alpha=0.3, label = '_'*j + 'Exploration right')
                    x_start=0
                    j+=1
            if 'Climbing' in list_BOI:
                for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Climbing'].tolist()):
                    if y == 1:
                        x_start = x
                    if y == -1 and x_start!=0:
                        ax1.axvspan(x_start, x, facecolor='cornflowerblue', alpha=0.3, label = '_'*m + 'Climbing')
                        x_start=0
                        m+=1
            if 'Rearing' in list_BOI:
                for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Rearing'].tolist()):
                    if y == 1:
                        x_start = x
                    if y == -1 and x_start!=0:
                        ax1.axvspan(x_start, x, facecolor='mediumpurple', alpha=0.3, label = '_'*n + 'Rearing')
                        x_start=0
                        n+=1
        
        elif session in ['S2']:
            #plots fiberpho trace and behaviour
            p1, = ax1.plot('Time(s)', 'Denoised dFF', linewidth=.6, color='black', label='_GCaMP', data = behavprocesssnip_df)
        
            #makes areas corresponding to behaviours
            i = 0
            j = 0
            m = 0
            n = 0
            for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Exploration non social'].tolist()):
                if y == 1:
                    x_start = x
                if y == -1 and x_start!=0:   
                    ax1.axvspan(x_start, x, facecolor='grey', alpha=0.3, label = '_'*i + 'Exploration non social')
                    x_start=0
                    i += 1
            for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Exploration social'].tolist()):
                if y == 1:
                    x_start = x
                if y == -1 and x_start!=0:   
                    ax1.axvspan(x_start, x, facecolor='mediumvioletred', alpha=0.3, label = '_'*j + 'Exploration social')
                    x_start=0
                    j += 1
            if 'Climbing' in list_BOI:
                for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Climbing'].tolist()):
                    if y == 1:
                        x_start = x
                    if y == -1 and x_start!=0:
                        ax1.axvspan(x_start, x, facecolor='cornflowerblue', alpha=0.3, label = '_'*m + 'Climbing')
                        x_start=0
                        m+=1
                        
            if 'Rearing' in list_BOI:
                for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Rearing'].tolist()):
                    if y == 1:
                        x_start = x
                    if y == -1 and x_start!=0:
                        ax1.axvspan(x_start, x, facecolor='mediumpurple', alpha=0.3, label = '_'*n + 'Rearing')
                        x_start=0
                        n+=1
        
        #makes vertical line for entry opening of gate
        if 'Gate opens' in list_BOI:
            x_entry = fiberbehav_df.at[int(np.where(fiberbehav_df['Gate opens'] == 1)[0][0]), 'Time(s)']
            ax1.axvline(x_entry, color='lightsteelblue', ls = '--', label = 'Gate opens' )
                
        #makes vertical line for entry in open field
        x_entry = fiberbehav_df.at[int(np.where(fiberbehav_df['Entry in arena'] == 1)[0][0]), 'Time(s)']
        ax1.axvline(x_entry, color='slategrey', ls = '--', label = 'Entry in arena' )
    
    ax1.set_ylabel(r'$\Delta$F/F')
    ax1.set_xlabel('Seconds')
    #ax1.set_ylim(-1.5,2.5)
    ax1.legend(loc = 'upper right')
    ax1.margins(0.01,0.03)
    ax1.set_title(f'dFF with Behavioural Scoring - {exp} {session} {mouse} - interbout {THRESH_S} - filtero{ORDER}f{CUT_FREQ}')
    
    #save figure
    fig2.savefig(mouse_path / f'{mouse}_threshold{THRESH_S}s_filtero{ORDER}f{CUT_FREQ}_fiberbehav_scaled.pdf')
    fig2.savefig(mouse_path / f'{mouse}_threshold{THRESH_S}s_filtero{ORDER}f{CUT_FREQ}_fiberbehav_scaled.png')
    
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
    
    #take the beginning of the behaviour as the beginning of the slope : min(dFF) on 3 secs before entry
    if BOI in ['Gate opens','Entry in arena']:
        for i,ind_onset in enumerate(list_ind_event):
            list_ind_event[i] = behavprocess_df.loc[ind_onset-3*SAMPLERATE:ind_onset, 'Denoised dFF'].idxmin()
        
    PETH_array = np.zeros((len(list_ind_event),(POST_TIME+PRE_TIME)*SAMPLERATE+1))
    for i, ind_event in enumerate(list_ind_event) :
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

def plot_PETH_average(PETH_data, BOI, event, timewindow):
    """
    Plots only average of PETH, no heatmap
    For BOIs only happening once in the trial
    """
    
    PRE_TIME = timewindow[0]
    POST_TIME = timewindow[1]
    
    #create figure
    fig4 = plt.figure(figsize=(6,4))
    ax5 = fig4.add_subplot(111)
    
    #create time vector
    peri_time = np.arange(-PRE_TIME, POST_TIME+0.1, 0.1) 
    
    #plot trace
    p1, = ax5.plot(peri_time, PETH_data[0], linewidth=2, color='green')
    p2 = ax5.axvline(x=0, linewidth=2, color='slategray', ls = '--', label=BOI)
    
    ax5.set_xlabel('Seconds')
    ax5.set_ylabel('z-scored $\Delta$F/F')
    ax5.legend(handles=[p2], loc='upper left', fontsize = 'small')
    ax5.margins(0, 0.1)
    ax5.set_title(f'{BOI} - {exp} {session} {mouse}')
    
    #save figure
    fig4.savefig(mouse_path / f'{mouse}{BOI}_PETH.pdf')
    
def plot_PETH_pooled(included_groups, PETHarray_list, BOI, event, timewindow):
    
    PRE_TIME = timewindow[0]
    POST_TIME = timewindow[1]
    
    #create figure
    fig4 = plt.figure(figsize=(6,4))
    ax5 = fig4.add_subplot(111)
    
    #create time vector
    peri_time = np.arange(-PRE_TIME, POST_TIME+0.1, 0.1)       
    
    listmean_dFF_snips = []
    listsem_dFF_snips = []
    #calculate mean dFF and std error
    for (i,PETH_data) in enumerate(PETHarray_list):
        listmean_dFF_snips.append(np.mean(PETH_data, axis=0))
        listsem_dFF_snips.append(np.std(PETH_data, axis=0))
        
    #plot individual traces and mean CD
    # for snip in PETHarray_list[0]:
    #     p1, = ax5.plot(peri_time, snip, linewidth=.5,
    #                    color='cornflowerblue', alpha=.3)
    p2, = ax5.plot(peri_time, listmean_dFF_snips[0], linewidth=2,
                   color='cornflowerblue', label=included_groups[0])   
    #plot standard error bars
    p3 = ax5.fill_between(peri_time, listmean_dFF_snips[0]+(listsem_dFF_snips[0]/np.sqrt(len(PETHarray_list[0]))),
                      listmean_dFF_snips[0]-(listsem_dFF_snips[0]/np.sqrt(len(PETHarray_list[0]))), facecolor='cornflowerblue', alpha=0.2)
    
    #plot individual traces and mean HFD
    # for snip in PETHarray_list[1]:
    #     p4, = ax5.plot(peri_time, snip, linewidth=.5,
    #                    color='coral', alpha=.3)
    p5, = ax5.plot(peri_time, listmean_dFF_snips[1], linewidth=2,
                   color='orangered', label=included_groups[1])   
    #plot standard error bars
    p6 = ax5.fill_between(peri_time, listmean_dFF_snips[1]+(listsem_dFF_snips[1]/np.sqrt(len(PETHarray_list[1]))),
                      listmean_dFF_snips[1]-(listsem_dFF_snips[1]/np.sqrt(len(PETHarray_list[1]))), facecolor='orangered', alpha=0.2)
    
    p8 = ax5.axvline(x=0, linewidth=2, color='slategray', ls = '--', label=BOI)
    
    ax5.set_xlabel('Seconds')
    ax5.set_ylabel('z-scored $\Delta$F/F')
    ax5.legend(handles=[p2,p5,p8], loc='upper left', fontsize = 'small')
    ax5.margins(0, 0.1)
    ax5.set_title(f'{BOI} - {exp} {session} {included_groups[0]} {included_groups[1]}')
    
    #save figure
    fig4.savefig(repo_path / f'{included_groups[0]}{included_groups[1]}{BOI}_PETHred.pdf')
    fig4.savefig(repo_path / f'{included_groups[0]}{included_groups[1]}{BOI}_PETHred.png')
    
########
#SCRIPT#
########

#%%Plots for entry in zone
#for exp_path in [Path(f.path) for f in os.scandir(analysis_path) if f.is_dir()]:
for exp_path in ['K:\\Alice\\Fiber\\202301_CA2b5\\Analysis\\Essai1', 'K:\\Alice\\Fiber\\202301_CA2b5\\Analysis\\Essai2',
                 'K:\\Alice\\Fiber\\202301_CA2b5\\Analysis\\OdDis1', 'K:\\Alice\\Fiber\\202301_CA2b5\\Analysis\\OdDis_postshock',
                 'K:\\Alice\\Fiber\\202301_CA2b5\\Analysis\\SRM']:
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
            if not os.path.exists(mouse_path / f'{mouse}_threshold{THRESH_S}s_filtero{ORDER}f{CUT_FREQ}_fiberbehav_scaled.pdf'):
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
   

# #%%Run test
           
# EVENT_TIME_THRESHOLD = 0
# THRESH_S = 0
# CUT_FREQ = 1
            
# exp_path = Path('K://Alice//Fiber//202207_CA2db3//Analysis//OdDis_3')
# exp = str(exp_path).split('\\')[-1]
# for session_path in [Path(f.path) for f in os.scandir(exp_path) if f.is_dir()]:
#     session = str(session_path).split('\\')[-1]
#     print('##########################################')
#     print(f'EXPERIMENT : {exp} - SESSION : {session}')
#     print('##########################################')
#     code = session_code(session)
#     #get data path related to the task in protocol excel file
#     data_path_exp = data_path / proto_df.loc[proto_df['Task']==exp, 'Data_path'].values[0]
#     #create repository for values of thresholds : length and interbout
#     repo_path = session_path / f'length{EVENT_TIME_THRESHOLD/10}_interbout{THRESH_S}_o{ORDER}f{CUT_FREQ}'
#     if not os.path.exists(repo_path):
#         os.mkdir(repo_path)
#         for subject in subjects_df['Subject']:
#             os.mkdir(repo_path / subject)
#     #create list of mean_dFFs and diff_dFFs
#     mean_dFFs_list = []
#     subject_list = []
#     diffdFF_list = []
#     for mouse_path in [Path(f.path) for f in os.scandir(repo_path) if f.is_dir()]:
#         # '/' on mac, '\\' on windows
#         mouse = str(mouse_path).split('\\')[-1]
#         print("################")
#         print(f'MOUSE : {mouse}')
#         group = mouse[:-1]
        
#         #get data
#         behav_path = data_path_exp / f'behav_{code}_{mouse}.csv'
#         fiberpho_path = data_path_exp / f'{mouse}_{code}_dFFfilt.csv'
#         camera_path = data_path_exp / f'{mouse}_{code}_camera.csv'
#         rawdata_path = data_path_exp / f'{mouse}_{code}.csv'
        
#         #begin analysis only if behaviour has been scored and plots haven't been done
#         ready = False
#         if os.path.exists(behav_path):
#             ready = True
#         print(f'ready? {ready}')
        
#         done = True
#         if not os.path.exists(mouse_path / f'{mouse}_threshold{THRESH_S}s_filtero{ORDER}f{CUT_FREQ}_fiberbehav_scaled.pdf'):
#             done = False
#         print(f'done? {done}')
        
#         if os.path.exists(camera_path) and ready == True and done == False:
#             camera = pd.read_csv(camera_path)
            
#         elif os.path.exists(rawdata_path) and ready == True and done == False:
#             rawdata_cam_df = pd.read_csv(rawdata_path, skiprows=1, usecols=['Time(s)','DI/O-3'])
        
#         if os.path.exists(behav_path) and os.path.exists(fiberpho_path) and ready == True and done == False:
#             behav10Sps = pd.read_csv(behav_path)
#             fiberpho = pd.read_csv(fiberpho_path)
#             print(exp, session, mouse)
            
#               #list of behaviours to analyze
#             list_BOI = behav10Sps.columns[1:].tolist()
            
#             #align behaviour and fiberpho data, create fbprocess.xslx
#             print('timevector')
#             timevector = time_vector(fiberpho, SAMPLERATE)
#             print('timestamp')
#             if os.path.exists(camera_path):
#                 print('---------> from camera')
#                 timestart_camera = timestamp_camera(camera)[0]
#             else:
#                 print('---------> from rawdata')
#                 timestart_camera = timestamp_camera_fromraw(rawdata_cam_df)[0]
#             print('start camera : ', timestart_camera)
#             print('aligning')
#             fiberbehav_df = align_behav(behav10Sps, fiberpho, timevector, timestart_camera)
#             fiberbehav_df = filter_dFF(fiberbehav_df, ORDER, CUT_FREQ)
#             if EVENT_TIME_THRESHOLD==0 and THRESH_S==0:
#                 mean_dFF_df = meandFF_behav(list_BOI, fiberbehav_df)
#                 mean_dFFs_list.append(mean_dFF_df)
#                 #plot isosbestic and gcamp data
#                 plot_fiberpho(fiberbehav_df)
#             print('processing')
#             (fiberbehav2_df, fbprocess_df) = behav_process(fiberbehav_df, list_BOI)
            
#             #add diffdFFs to the list
#             diffdFF_df = diff_dFF(fiberbehav_df, fbprocess_df, list_BOI)
#             diffdFF_list.append(diffdFF_df)
            
#             # #plot raw data
#             #plot_rawdata(rawdata_df)
            
#             #plot fiberpho data aligned with behav
#             plot_fiberpho_behav(fbprocess_df)
#             #plot_fiberpho_behav_snip(fiberbehav2_df, timestart_camera)
            
#             # for BOI in list_BOI:
#             #     if BOI == 'Entry in arena' or BOI == 'Gate opens':
#             #         PETH_data = PETH(fbprocess_df, BOI, 'onset', [6,10])
#             #         plot_PETH_average(PETH_data, BOI, 'onset', [6,10])
#             #     else:
#             #         for event, timewindow in zip(list_EVENT,list_TIMEWINDOW): #CAUTION : adapt function!!
#             #             PETH_data = PETH(fbprocess_df, BOI, event, timewindow)
#             #             plot_PETH(PETH_data, BOI, event, timewindow)
   
#     if EVENT_TIME_THRESHOLD==0 and THRESH_S==0 and mean_dFFs_list != []: 
#         meandFFs_allmice = pd.concat(mean_dFFs_list)
#         meandFFs_allmice.to_excel(repo_path / f'{exp}_{session}_length{EVENT_TIME_THRESHOLD/10}_interbout{THRESH_S}_filtero{ORDER}f{CUT_FREQ}_meandFFs.xlsx')
        
#         df_plotmean=meandFFs_allmice.groupby('Group')
#         means=df_plotmean.mean()
#         errors=df_plotmean.std()
        
#         #plot figure
#         fig_meandFF, axmean = plt.subplots(1, 1, figsize=(7, 6))
#         labels = meandFFs_allmice.columns[2:]
        
#         means.plot.bar(y=labels, yerr=errors[labels], capsize=2, rot=0, 
#                        ax=axmean, linewidth=.1, colormap='viridis')
            
#         axmean.set_ylabel('Mean dFF')
#         axmean.set_title(f'Mean dFF - {exp} {session} - length{EVENT_TIME_THRESHOLD/10} interbout{THRESH_S} - filtero{ORDER}f{CUT_FREQ}')
        
#         fig_meandFF.savefig(repo_path / f'{exp}_{session}_length{EVENT_TIME_THRESHOLD/10}_interbout{THRESH_S}_filtero{ORDER}f{CUT_FREQ}_meandFFs.png')
        
#     if diffdFF_list != []:
#         diffdFFs_allmice = pd.concat(diffdFF_list)
#         diffdFFsmean = diffdFFs_allmice.groupby(['Subject','Behaviour'], as_index=False).mean()
#         diffdFFsmean = diffdFFsmean[diffdFFsmean['Bout']>2]
        
#         #group diffdFFsmean dataframe to plot it!
#         cols = ['Group']
#         list_behav = set(diffdFFsmean['Behaviour'])
#         for behav in list_behav:
#             cols.append(behav)
#         list_subjects = set(diffdFFsmean['Subject'])
        
#         diffdFFsmeanplot = pd.DataFrame(len(list_subjects)*[[0]*3],
#                                         columns=cols, index=list_subjects)
#         for subject in list_subjects:
#             diffdFFsmeanplot.loc[subject, 'Group']=subject[:-1]
#             df_subj = diffdFFsmean[diffdFFsmean['Subject']==subject]
#             for behav in list_behav:
#                 if behav in df_subj['Behaviour'].values:
#                     value = df_subj.loc[df_subj['Behaviour']==behav,'Delta dFF'].values[0]
#                     diffdFFsmeanplot.loc[subject, behav]=value
         
#         diffdFFs_allmice.to_excel(repo_path / f'{exp}_{session}_length{EVENT_TIME_THRESHOLD/10}_interbout{THRESH_S}_filtero{ORDER}f{CUT_FREQ}_diffdFFsall.xlsx')
#         diffdFFsmeanplot.to_excel(repo_path / f'{exp}_{session}_length{EVENT_TIME_THRESHOLD/10}_interbout{THRESH_S}_filtero{ORDER}f{CUT_FREQ}_diffdFFsmean.xlsx')
        
#         df_plotdiff=diffdFFsmeanplot.groupby('Group')
#         means_diff=df_plotdiff.mean()
#         errors_diff=df_plotdiff.std()
        
#         #plot figure
#         fig_diffdFF, axdiff = plt.subplots(1, 1, figsize=(7, 6))
        
#         means_diff.plot.bar(y=list_behav, yerr=errors_diff[list_behav], rot=0,
#                             ax=axdiff, capsize=2, colormap='summer')
        
#         axdiff.set_ylabel('Diff dFF')
#         axdiff.set_title(f'Diff dFF - {exp} {session} - length{EVENT_TIME_THRESHOLD/10} interbout{THRESH_S} - filtero{ORDER}f{CUT_FREQ}')
#         #save figure
#         fig_diffdFF.savefig(repo_path / f'{exp}_{session}_length{EVENT_TIME_THRESHOLD/10}_interbout{THRESH_S}_filtero{ORDER}f{CUT_FREQ}_diffdFFs.png')
                                         
# #%%delete files
# repo_path = Path('D:\\Alice\\Fiber\\202110_CA2db2\\Analysis\\SRM\\S2\\length2_interbout2')
# for subject in Path(repo_path).iterdir():
#     files_in_directory = os.listdir(subject)
#     filtered_files = [file for file in files_in_directory if file.endswith(".csv")]
#     for file in filtered_files:
#         path_to_file = os.path.join(subject, file)
#         os.remove(path_to_file)
        
# #%%copy files
# new_path = Path('D:\\Alice\\Fiber\\202110_CA2db2\\Data\\20211006_AliceF_CA2b2SRMmales1')
# session_path = Path('D:\\Alice\\Fiber\\202110_CA2db2\\Analysis\\SRM\\S2\\length2_interbout2')
# for subject in Path(session_path).iterdir():
#     files_in_directory = os.listdir(subject)
#     filtered_files = [file for file in files_in_directory if file.endswith(".csv")]
#     for file in filtered_files:
#         path_to_file = os.path.join(subject, file)
#         newpath_to_file = os.path.join(new_path, file)
#         os.rename(path_to_file, newpath_to_file)

# #%%test filter
# fiberdata = pd.read_csv('D:\\Alice\\Fiber\\202110_CA2db2\\Data\\To_do\\20211110_AliceF_CA2b2beddingRimo\\CDf1_1_dFFfilt.csv')
# fiberpho = fiberdata['Analog In. | Ch.1 470 nm (Deinterleaved)_dF/F0-Analog In. | Ch.1 405 nm (Deinterleaved)_dF/F0_LowPass']
# samplingrate = 1000/(fiberdata.loc[1000,'Time(s)']-fiberdata.loc[0,'Time(s)'])
# sos = signal.butter(5, 0.1, btype='low', analog=False, output='sos', fs=samplingrate)
# filtered = signal.sosfilt(sos, fiberpho)

# fiberdata['Filtered_dFF'] = filtered

# fig8 = plt.figure(figsize=(20,4))
# ax81 = fig8.add_subplot(111)
# p81, = ax81.plot('Time(s)', 'Filtered_dFF',
#                 linewidth=.6, color='black', label='_GCaMP', data = fiberdata)