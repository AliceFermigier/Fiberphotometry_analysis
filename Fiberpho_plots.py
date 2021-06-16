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
sys.path.append('C:\\Users\\Alice\\Documents\\GitHub\\Fiberphotometry_analysis')

#%%
########
#LOADER#
########

from Fiberpho_loader import experiment_path, analysis_path, subjects_df, SAMPLERATE
from Fiberpho_loader import list_EVENT, list_TIMEWINDOW, PRE_EVENT_TIME, TIME_BEGIN

os.chdir(experiment_path)

#%%
###################
#DEFINED FUNCTIONS#
###################

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
    return pd.Series(np.linspace(0.0, duration, num = int(duration*SAMPLERATE)+1))

def timestamp_camera(camera) :
    """
    Function to extract the timestamps where the camera starts and stops
    --> Parameters
        camera : pd dataframe, camera I/O with sample rate = 12kSps
    --> Returns
        (camera_start, camera_stop) = timestamp when camera starts and stops in seconds (truncated to 0,1s) #camera_stop à enlever si pas besoin
    """
    ind_list = np.where(camera['Digital I/O | Ch.3 DI/O-3'] == 1)
    ind_list = ind_list[0].tolist()
    (ind_start, ind_stop) = (ind_list[0],ind_list[len(ind_list)-1])
    return (truncate(camera.at[ind_start, 'Time(s)'], 1),
            truncate(camera.at[ind_stop, 'Time(s)'], 1))
    
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
    
    return fiberbehav_df

def behav_process(fiberbehav_df, list_BOI):

    #1 fuse explorations that are too close
    
    THRESHOLD = 10 #set threshold in 0.1seconds (because samplerate = 10Sps)
    
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
        
    behavprocess_df.to_excel(mouse_path / 'fbprocess.xlsx')

    
    return(fiberbehav_df, behavprocess_df)

def plot_rawdata(rawdata_df):
    """
    Plots raw isosbestic and GCaMP traces
    """
    
    # if windows : \\
    # if mac : /
    mouse = str(mouse_path).split('\\')[-1]
    session = str(mouse_path).split('\\')[-2]
    exp = str(mouse_path).split('\\')[-3]
    
    time_stop = rawdata_df.loc[len(rawdata_df)-1,'Time(s)']
    time_stop = np.ceil(time_stop)
    
    fig5 = plt.figure(figsize=(10,6))
    ax7 = fig5.add_subplot(211)

    p1, = ax7.plot('Time(s)', 'Analog In. | Ch.1 470 nm (Deinterleaved)', 
                   linewidth=1, color='deepskyblue', label='GCaMP', data = rawdata_df) 
    p2, = ax7.plot('Time(s)', 'Analog In. | Ch.1 405 nm (Deinterleaved)', 
                   linewidth=1, color='blueviolet', label='ISOS', data = rawdata_df)
    
    ax7.set_ylabel('V')
    ax7.set_xlabel('Seconds')
    ax7.legend(handles=[p1,p2], loc='upper right')
    ax7.margins(0.01,0.3)
    ax7.set_title('GCaMP and Isosbestic raw traces - '+exp+' '+session+' '+mouse)
    
    #save figure
    fig5.savefig(str(mouse_path)+'/'+mouse+'_rawdata.png')
    
    return

def plot_fiberpho(fiberbehav_df):
    """
    Plots isosbestic and Ca dependent deltaF/F
    """
    
    # if windows : \\
    # if mac : /
    mouse = str(mouse_path).split('\\')[-1]
    session = str(mouse_path).split('\\')[-2]
    exp = str(mouse_path).split('\\')[-3]
    
    #crops starting artifacts
    fiberbehavsnip_df = fiberbehav_df[fiberbehav_df['Time(s)'] > 40]
       
    fig1 = plt.figure(figsize=(10,6))
    ax0 = fig1.add_subplot(211)
    
    p1, = ax0.plot('Time(s)', '470nm deltaF/F', linewidth=1, color='deepskyblue', label='GCaMP', data = fiberbehavsnip_df) 
    p2, = ax0.plot('Time(s)', '405nm deltaF/F', linewidth=1, color='blueviolet', label='ISOS', data = fiberbehavsnip_df)
    
    ax0.set_ylabel(r'$\Delta$F/F')
    ax0.set_xlabel('Seconds')
    ax0.legend(handles=[p1,p2], loc='upper right')
    ax0.margins(0.01,0.2)
    ax0.set_title('GCaMP and Isosbestic dFF - '+exp+' '+session+' '+mouse)
    
    #save figure
    fig1.savefig(str(mouse_path)+'/'+mouse+'_GCaMP_ISOS.png')
    
    return

def plot_fiberpho_behav(fiberbehav_df):
    """
    Plots denoised deltaF/F aligned with behaviour (includes baseline)
    """
    
    fiberbehavsnip_df = fiberbehav_df[fiberbehav_df['Time(s)'] > 40]
    
    # if windows : \\
    # if mac : /
    mouse = str(mouse_path).split('\\')[-1]
    session = str(mouse_path).split('\\')[-2]
    exp = str(mouse_path).split('\\')[-3]
    
    fig2 = plt.figure(figsize=(20,12))
    ax1 = fig2.add_subplot(311)
    
    if session in ['Test 1h','Test 24h','Test','S3']:
        #plots fiberpho trace and behaviour
        p1, = ax1.plot('Time(s)', 'Denoised dFF', linewidth=1, color='black', label='_GCaMP', data = fiberbehavsnip_df)
        # p2, = ax1.plot(fiberbehavsnip_df['Time(s)'], fiberbehavsnip_df['Exploration fam'], linewidth=1, 
        #                color='moccasin', label='Exploration fam')
        
        #makes areas corresponding to behaviours
        i = 0
        j = 0
        for (x,y) in zip(fiberbehavsnip_df['Time(s)'].tolist(), fiberbehavsnip_df['Exploration fam'].tolist()):
            if y == 1:
                ax1.axvspan(x, x+0.1, facecolor='moccasin', alpha=0.5, label = '_'*i + 'Exploration fam')
                i+=1
        for (x,y) in zip(fiberbehavsnip_df['Time(s)'].tolist(), fiberbehavsnip_df['Exploration new'].tolist()):
            if y == 1:
                ax1.axvspan(x, x+0.1, facecolor='cornflowerblue', alpha=0.5, label = '_'*j + 'Exploration new')
                j+=1
        # for (x,y) in zip(fiberbehavsnip_df['Time(s)'].tolist(), fiberbehavsnip_df['Climbing'].tolist()):
        #     if y == 1:
        #         ax1.axvspan(x, x+0.1, facecolor='cornflowerblue', alpha=0.5)
        
    if session in ['Habituation','Training','S1','S2']:
        #plots fiberpho trace and behaviourb 
        p1, = ax1.plot('Time(s)', 'Denoised dFF', linewidth=1, color='black', label='_GCaMP', data = fiberbehavsnip_df)

        #makes areas corresponding to behaviours
        i = 0
        j = 0
        for (x,y) in zip(fiberbehavsnip_df['Time(s)'].tolist(), fiberbehavsnip_df['Exploration left'].tolist()):
            if y == 1:
                ax1.axvspan(x, x+0.1, facecolor='orange', alpha=0.5, label = '_'*i + 'Exploration left')
                i+=1
        for (x,y) in zip(fiberbehavsnip_df['Time(s)'].tolist(), fiberbehavsnip_df['Exploration right'].tolist()):
            if y == 1:
                ax1.axvspan(x, x+0.1, facecolor='darkturquoise', alpha=0.5, label = '_'*j + 'Exploration right')
                j+=1
        # for (x,y) in zip(fiberbehavsnip_df['Time(s)'].tolist(), fiberbehavsnip_df['Climbing'].tolist()):
        #     if y == 1:
        #         ax1.axvspan(x, x+0.1, facecolor='cornflowerblue', alpha=0.5)
            
    #makes vertical line for entry in open field
    x_entry = fiberbehav_df.at[int(np.where(fiberbehav_df['Entry in open field'] == 1)[0][0]), 'Time(s)']
    ax1.axvline(x_entry, color='slategrey', ls = '--', label = 'Entry OF' )
    
    ax1.set_ylabel(r'$\Delta$F/F')
    ax1.set_xlabel('Seconds')
    ax1.set_ylim(-15,35)
    ax1.legend(loc = 'upper right')
    ax1.margins(0.01,0.03)
    ax1.set_title('dFF with Behavioural Scoring - '+exp+' '+session+' '+mouse)
    
    #save figure
    fig2.savefig(str(mouse_path)+'/'+mouse+'_fiberbehav_scaled.png')
    
    return

def plot_fiberpho_behav_snip(fiberbehav_df, timestart_camera):
    """
    Plots denoised deltaF/F aligned with behaviour (starts when camera starts)
    """
    
    # if windows : \\
    # if mac : /
    mouse = str(mouse_path).split('\\')[-1]
    session = str(mouse_path).split('\\')[-2]
    exp = str(mouse_path).split('\\')[-3]
    
    fig3 = plt.figure(figsize=(20,12))
    ax2 = fig3.add_subplot(311)
    
    #y_scale = 7 #adjust according to data needs
    #y_shift = 30 #scale and shift are just for aesthetics
    
    #crops data when camera starts
    fiberbehavsnip_df = fiberbehav_df[fiberbehav_df['Time(s)'] > timestart_camera]
    
    if session in ['Test 1h','Test 24h','Test','S3']:
        #plots fiberpho trace and behaviour
        p1, = ax2.plot('Time(s)', 'Denoised dFF', linewidth=1, color='black', label='_GCaMP', data = fiberbehavsnip_df)
        
        #makes areas corresponding to behaviours
        i = 0
        j = 0
        for (x,y) in zip(fiberbehavsnip_df['Time(s)'].tolist(), fiberbehavsnip_df['Exploration fam'].tolist()):
            if y == 1:
                ax2.axvspan(x, x+0.1, facecolor='moccasin', alpha=0.5, label = '_'*i + 'Exploration fam')
                i += 1
        for (x,y) in zip(fiberbehavsnip_df['Time(s)'].tolist(), fiberbehavsnip_df['Exploration new'].tolist()):
            if y == 1:
                ax2.axvspan(x, x+0.1, facecolor='cornflowerblue', alpha=0.5, label = '_'*j + 'Exploration new')
                j += 1
        # for (x,y) in zip(fiberbehavsnip_df['Time(s)'].tolist(), fiberbehavsnip_df['Climbing'].tolist()):
        #     if y == 1:
        #         ax2.axvspan(x, x+0.1, facecolor='cornflowerblue', alpha=0.5)
        
    if session in ['Habituation','Training','S1','S2']:
        #plots fiberpho trace and behaviour
        p1, = ax2.plot('Time(s)', 'Denoised dFF', linewidth=1, color='black', label='_GCaMP', data = fiberbehavsnip_df)
    
        #makes areas corresponding to behaviours
        i = 0
        j = 0
        for (x,y) in zip(fiberbehavsnip_df['Time(s)'].tolist(), fiberbehavsnip_df['Exploration left'].tolist()):
            if y == 1:
                ax2.axvspan(x, x+0.1, facecolor='orange', alpha=0.5, label = '_'*i + 'Exploration left')
                i += 1
        for (x,y) in zip(fiberbehavsnip_df['Time(s)'].tolist(), fiberbehavsnip_df['Exploration right'].tolist()):
            if y == 1:
                ax2.axvspan(x, x+0.1, facecolor='darkturquoise', alpha=0.5, label = '_'*j + 'Exploration right')
                j += 1
        # for (x,y) in zip(fiberbehavsnip_df['Time(s)'].tolist(), fiberbehavsnip_df['Climbing'].tolist()):
        #     if y == 1:
        #         ax2.axvspan(x, x+0.1, facecolor='cornflowerblue', alpha=0.5)
            
    #makes vertical line for entry in open field
    x_entry = fiberbehav_df.at[int(np.where(fiberbehav_df['Entry in open field'] == 1)[0][0]), 'Time(s)']
    ax2.axvline(x_entry, color='slategrey', ls = '--', label = 'Entry OF' )
       
    ax2.set_ylabel(r'$\Delta$F/F')
    ax2.set_xlabel('Seconds')
    ax2.set_ylim(-15,35)
    ax2.legend(loc = 'upper right')
    ax2.margins(0.01,0.03)
    ax2.set_title('dFF with Behavioural Scoring - '+exp+' '+session+' '+mouse)
    
    #save figure
    fig3.savefig(str(mouse_path)+'/'+mouse+'_fiberbehavsnip_scaled.png')
    
    return
    
#%%
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
    print('PETH')
    #set time window relative to event
    PRE_TIME = timewindow[0]
    POST_TIME = timewindow[1]
    
    #if behaviour not recognized
    if BOI not in behavprocess_df.columns[2:].tolist() :
        print('BOI not recognized')
    
    #creates array of fiberpho data centered on event
    if event == 'onset':
        list_ind_event = np.where(behavprocess_df[BOI] == 1)[0].tolist()
    elif event == 'withdrawal':
        list_ind_event = np.where(behavprocess_df[BOI] == -1)[0].tolist()
        
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
    
    # if windows : \\
    # if mac : /
    mouse = str(mouse_path).split('\\')[-1]
    session = str(mouse_path).split('\\')[-2]
    exp = str(mouse_path).split('\\')[-3]
    
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
    fig4.savefig(str(mouse_path)+'/'+mouse+'_'+BOI+'_PETH'+event[0]+'.png')
    
    return

def plot_PETH_average(PETH_data, BOI, event, timewindow):
    """
    Plots only average of PETH, no heatmap
    For BOIs only happening once in the trial
    """
    # if windows : \\
    # if mac : /
    mouse = str(mouse_path).split('\\')[-1]
    session = str(mouse_path).split('\\')[-2]
    exp = str(mouse_path).split('\\')[-3]
    
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
    ax5.set_title(BOI+' - '+exp+' '+session+' '+mouse)
    
    #save figure
    fig4.savefig(str(mouse_path)+'/'+mouse+'_'+BOI+'_PETH.png')
    
########
#SCRIPT#
########

#%%Run test

mouse_path = analysis_path / '20210304_SRM/S3/CD3'
mouse = str(mouse_path).split('\\')[-1]

behav_path = str(mouse_path) + '/behav_' + mouse + '.csv'
fiberpho_path = str(mouse_path) + '/' + mouse + '_dFFfilt.csv'
camera_path = str(mouse_path) + '/' + mouse + '_camera.csv'
rawdata_path = str(mouse_path) + '/' + mouse +  '_rawdata.csv'
    
behav10Sps = pd.read_csv(behav_path)
fiberpho = pd.read_csv(fiberpho_path)
camera = pd.read_csv(camera_path)
rawdata_df = pd.read_csv(rawdata_path)


#list of behaviours to analyze
list_BOI = behav10Sps.columns[1:].tolist()

print('timevector')
timevector = time_vector(fiberpho, SAMPLERATE)
print('timestamp')
timestart_camera = timestamp_camera(camera)[0]
print('start camera : ', timestart_camera)
print('aligning')
fiberbehav_df = align_behav(behav10Sps, fiberpho, timevector, timestart_camera)
print('processing')
(fiberbehav2_df, fbprocess_df) = behav_process(fiberbehav_df, list_BOI)

#plot raw data
plot_rawdata(rawdata_df)

#plot isosbestic and gcamp data
plot_fiberpho(fiberbehav_df)

#plot fiberpho data aligned with behav
plot_fiberpho_behav(fiberbehav2_df)
plot_fiberpho_behav_snip(fiberbehav2_df, timestart_camera)

for BOI in list_BOI:
    if BOI == 'Entry in open field':
        PETH_data = PETH(fbprocess_df, BOI, 'onset', [10,10])
        plot_PETH_average(PETH_data, BOI, 'onset', [10,10])
    else:
        for event, timewindow in zip(list_EVENT,list_TIMEWINDOW): #CAUTION : adapt function!!
            PETH_data = PETH(fbprocess_df, BOI, event, timewindow)
            plot_PETH(PETH_data, BOI, event, timewindow)

#%%Run for all 

#list experiments, sessions and mice:
for exp_path in Path(analysis_path).iterdir():
    for session_path in Path(exp_path).iterdir():
        for mouse_path in Path(session_path).iterdir():
            #get the name of the mouse, session and experiment
            # '/' on mac, '\\' on windows
            exp = str(mouse_path).split('\\')[-3]
            session = str(mouse_path).split('\\')[-2]
            mouse = str(mouse_path).split('\\')[-1]
            print(exp, session, mouse)
            
            if Path(str(mouse_path)+'/'+mouse+'_fiberbehavsnip.png').is_file():
                #get data
                behav_path = str(mouse_path) + '/behav_' + mouse + '.csv'
                fiberpho_path = str(mouse_path) + '/' + mouse + '_dFFfilt.csv'
                camera_path = str(mouse_path) + '/' + mouse + '_camera.csv'
                rawdata_path = str(mouse_path) + '/' + mouse +  '_rawdata.csv'
                    
                behav10Sps = pd.read_csv(behav_path)
                fiberpho = pd.read_csv(fiberpho_path)
                camera = pd.read_csv(camera_path)
                rawdata_df = pd.read_csv(rawdata_path)
                
                 #list of behaviours to analyze
                list_BOI = behav10Sps.columns[1:].tolist()
                
                #align behaviour and fiberpho data, create fbprocess.xslx
                print('timevector')
                timevector = time_vector(fiberpho, SAMPLERATE)
                print('timestamp')
                timestart_camera = timestamp_camera(camera)[0]
                print('start camera : ', timestart_camera)
                print('aligning')
                fiberbehav_df = align_behav(behav10Sps, fiberpho, timevector, timestart_camera)
                print('processing')
                (fiberbehav2_df, fbprocess_df) = behav_process(fiberbehav_df, list_BOI)
                
                # #plot raw data
                # plot_rawdata(rawdata_df)
                
                # #plot isosbestic and gcamp data
                # plot_fiberpho(fiberbehav_df)
                
                #plot fiberpho data aligned with behav
                plot_fiberpho_behav(fiberbehav2_df)
                plot_fiberpho_behav_snip(fiberbehav2_df, timestart_camera)
                
                # for BOI in list_BOI:
                #     if BOI == 'Entry in open field':
                #         PETH_data = PETH(fbprocess_df, BOI, 'onset', [6,10])
                #         plot_PETH_average(PETH_data, BOI, 'onset', [6,10])
                #     else:
                #         for event, timewindow in zip(list_EVENT,list_TIMEWINDOW): #CAUTION : adapt function!!
                #             PETH_data = PETH(fbprocess_df, BOI, event, timewindow)
                #             plot_PETH(PETH_data, BOI, event, timewindow)