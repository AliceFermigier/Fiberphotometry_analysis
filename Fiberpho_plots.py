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
import matplotlib.cm as cm
import matplotlib.colors as mcolors
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
from Fiberpho_loader import list_EVENT, list_TIMEWINDOW, PRE_EVENT_TIME, TIME_BEGIN, THRESH_S, EVENT_TIME_THRESHOLD

from Func_fiberplots import session_code, truncate, time_vector, timestamp_camera, timestamp_camera_fromraw

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
        
    else :
        print('Error : too many behaviours in Boris binary file. Score up to 4 or add line to function')
    
    return(fiberbehav_df)

def meandFF_behav(list_BOI, fiberbehav_df):
    """
    Calculates mean dFF during each behaviour
    """
    list_behav_analyzed = []
    list_meandFF = []
    
    #get index of when the gate opens
    ind_start_trial = fiberbehav_df.index[fiberbehav_df['Gate opens'] == 1].tolist()[0]
    
    #create list of behaviours and list of correspondind mean dFFs
    for behav in list_BOI:
        if fiberbehav_df[behav].sum() > 1:
            list_behav_analyzed.append(behav)
            meandFF_behav_df = fiberbehav_df.groupby([behav], as_index=False).mean()
            print(meandFF_behav_df)
            list_meandFF.append(meandFF_behav_df.loc[meandFF_behav_df[behav]==1, 'Denoised dFF'].values[0])
                                
    #calculate mean dFF during baseline
    
    meandFF_baseline = fiberbehav_df.loc[:ind_start_trial, 'Denoised dFF'].mean()
    
    #calculate mean dFF when no exploration and after gate opens
    
    meandFF_postbaseline_df = fiberbehav_df.loc[ind_start_trial:]
    for behav in list_behav_analyzed:
        meandFF_postbaseline_df = meandFF_postbaseline_df.loc[meandFF_postbaseline_df[behav]==0]
        
    meandFF_postbaseline = meandFF_postbaseline_df['Denoised dFF'].mean()
    
    #calculate mean dFF when no exploration on total trial
    meandFF_df = fiberbehav_df
    for behav in list_behav_analyzed:
        meandFF_df = meandFF_df.loc[meandFF_df[behav]==0]
        
    meandFF = meandFF_df['Denoised dFF'].mean()
    
    #create dataframe with values
    
    list_dFFs = [mouse, meandFF,meandFF_baseline,meandFF_postbaseline]
    list_columns = ['Subject','Mean dFF','Baseline', 'Post_baseline']
    for (behav,meandFF) in zip(list_behav_analyzed,list_meandFF):
        list_dFFs.append(meandFF)
        list_columns.append(behav)
    meandFFs_df = pd.DataFrame(data=[list_dFFs], columns=list_columns)
    print(meandFFs_df)
    
    #Create figure of mean dFFs
    fig7 = plt.figure(figsize=(7,6))
    ax71 = fig7.add_subplot(111)
    print('x = ',meandFFs_df.columns.to_list()[1:])
    print('y = ',meandFFs_df.iloc[0, 1:].to_list())
    p71 = ax71.bar(meandFFs_df.columns.to_list()[1:], meandFFs_df.iloc[0, 1:])
    ax71.axhline(y=0, linewidth=.6, color='black')
    ax71.set_ylabel(r'$\Delta$F/F')
    ax71.set_title(f'Mean dFF - {exp} {session} {mouse}')
    #save figure
    fig7.savefig(mouse_path / f'{mouse}_meandFFs.pdf')
    
    #meandFFs_df.to_excel(mouse_path / f'{mouse}_meandFFglob.xlsx')
    
    return(meandFFs_df)

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

def filter_dFF(fiberbehav_df):
    """
    Apply additional filter to dFF data
    """
    fiberpho = fiberbehav_df['Denoised dFF']
    samplingrate = 1000/(fiberbehav_df.loc[1000,'Time(s)']-fiberbehav_df.loc[0,'Time(s)'])
    sos = signal.butter(5, 1, btype='low', analog=False, output='sos', fs=samplingrate)
    filtered_data = signal.sosfilt(sos, fiberpho)
    
    filtered_df = fiberbehav_df
    filtered_df['Denoised dFF'] = filtered_data
    
    return(filtered_df)
    

def diff_dFF(fiberbehav_df, behavprocess_df, list_BOI):
    """
    Calculates dFF difference between beginning and end of bouts
    Also calculates mean dFF during each bout
    """
    list_behav_analyzed = []
    for behav in list_BOI:
        if fiberbehav_df[behav].sum() > 1:
            list_behav_analyzed.append(behav)
    list_deltadFF = []
    list_meandFF = []
    listind_behav = []
    listind_bouts = []
    
    for i, behav in enumerate(list_BOI):
        list_starts = np.where(behavprocess_df[behav]==1)[0].tolist()
        list_stops = np.where(behavprocess_df[behav]==-1)[0].tolist()
        bout_n = 1
        for (start, stop) in zip(list_starts, list_stops):
            delta = behavprocess_df.loc[stop, 'Denoised dFF']-behavprocess_df.loc[start, 'Denoised dFF']
            mean = behavprocess_df.loc[start:stop, 'Denoised dFF'].mean()
            list_deltadFF.append(delta)
            list_meandFF.append(mean)
            listind_behav.append(behav)
            listind_bouts.append(bout_n)
            bout_n+=1
    
    diffdFF_df = pd.DataFrame(data = {'Subject':[mouse]*len(listind_behav), 'Behaviour':listind_behav, 
                                      'Bout':listind_bouts, 'Mean dFF':list_meandFF, 'Delta dFF':list_deltadFF})
    
    return(diffdFF_df)
    
def plot_meandFF(meandFFs_df):
    
    fig6 = plt.figure(figsize=(7,6))
    ax61 = fig6.add_subplot(111)

    # Choose a colormap
    colormap = cm.viridis
    #colormap = cm.jet 
    
    colorparams = [i+1 for i in range(len(meandFFs_df.columns))]
    normalize = mcolors.Normalize(vmin=1, vmax=colorparams[-1])
    list_colors = []
    
    #Create color list
    for param in colorparams:
        color = colormap(normalize(param))
        list_colors.append(color)
        
    #Plot histogram
    p61, = ax61.bar(meandFFs_df.columns, meandFFs_df[0], list_colors)
    ax61.set_ylabel(r'$\Delta$F/F')
    
    #save figure
    fig6.savefig(mouse_path / f'{mouse}_mean.pdf')
    
    return()
    
def plot_rawdata(rawdata_df):
    """
    Plots raw isosbestic and GCaMP traces
    """
    
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
    
    fig2 = plt.figure(figsize=(20,4))
    ax1 = fig2.add_subplot(111)
    
    if session in ['Test 1h','Test 24h','Test','S3']:
        #plots fiberpho trace and behaviour
        p1, = ax1.plot('Time(s)', 'Denoised dFF', linewidth=.6, color='black', label='_GCaMP', data = behavprocesssnip_df)
        # p2, = ax1.plot(behavprocesssnip_df['Time(s)'], behavprocesssnip_df['Exploration fam'], linewidth=1, 
        #                color='moccasin', label='Exploration fam')
        
        #makes areas corresponding to behaviours
        i = 0
        j = 0
        for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Exploration fam'].tolist()):
            if y == 1:
                x_start = x
            if y == -1 and x_start!=0:
                ax1.axvspan(x_start, x, facecolor='gold', alpha=0.5, label = '_'*i + 'Exploration fam')
                x_start=0
                i+=1
        for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Exploration new'].tolist()):
            if y == 1:
                x_start = x
            if y == -1 and x_start!=0:
                ax1.axvspan(x_start, x, facecolor='purple', alpha=0.5, label = '_'*j + 'Exploration new')
                x_start=0
                j+=1
        # for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Climbing'].tolist()):
        #     if y == 1:
        #         ax1.axvspan(x, x+0.1, facecolor='cornflowerblue', alpha=0.5)
        
    elif session in ['Habituation','Training','S1']:
        #plots fiberpho trace and behaviourb 
        p1, = ax1.plot('Time(s)', 'Denoised dFF', linewidth=.6, color='black', label='_GCaMP', data = behavprocesssnip_df)

        #makes areas corresponding to behaviours
        i = 0
        j = 0
        for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Exploration left'].tolist()):
            if y == 1:
                x_start = x
            if y == -1 and x_start!=0:
                ax1.axvspan(x_start, x, facecolor='orange', alpha=0.5, label = '_'*i + 'Exploration left')
                x_start=0
                i+=1
        for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Exploration right'].tolist()):
            if y == 1:
                x_start = x
            if y == -1 and x_start!=0:
                ax1.axvspan(x_start, x, facecolor='darkturquoise', alpha=0.5, label = '_'*j + 'Exploration right')
                x_start=0
                j+=1
        # for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Climbing'].tolist()):
        #     if y == 1:
        #         ax1.axvspan(x, x+0.1, facecolor='cornflowerblue', alpha=0.5)
    
    elif session in ['S2']:
        #plots fiberpho trace and behaviour
        p1, = ax1.plot('Time(s)', 'Denoised dFF', linewidth=.6, color='black', label='_GCaMP', data = behavprocesssnip_df)
    
        #makes areas corresponding to behaviours
        i = 0
        j = 0
        for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Exploration non social'].tolist()):
            if y == 1:
                x_start = x
            if y == -1 and x_start!=0:   
                ax1.axvspan(x_start, x, facecolor='grey', alpha=0.5, label = '_'*i + 'Exploration non social')
                x_start=0
                i += 1
        for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Exploration social'].tolist()):
            if y == 1:
                x_start = x
            if y == -1 and x_start!=0:   
                ax1.axvspan(x_start, x, facecolor='mediumvioletred', alpha=0.5, label = '_'*j + 'Exploration social')
                x_start=0
                j += 1
        # for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Climbing'].tolist()):
        #     if y == 1:
        #         ax2.axvspan(x, x+0.1, facecolor='cornflowerblue', alpha=0.5)
        
    #makes vertical line for entry opening of gate
    x_entry = fiberbehav_df.at[int(np.where(fiberbehav_df['Gate opens'] == 1)[0][0]), 'Time(s)']
    ax1.axvline(x_entry, color='lightsteelblue', ls = '--', label = 'Gate opens' )
            
    #makes vertical line for entry in open field
    x_entry = fiberbehav_df.at[int(np.where(fiberbehav_df['Entry in arena'] == 1)[0][0]), 'Time(s)']
    ax1.axvline(x_entry, color='slategrey', ls = '--', label = 'Entry OF' )
    
    ax1.set_ylabel(r'$\Delta$F/F')
    ax1.set_xlabel('Seconds')
    ax1.set_ylim(-1,2)
    ax1.legend(loc = 'upper right')
    ax1.margins(0.01,0.03)
    ax1.set_title(f'dFF with Behavioural Scoring - {exp} {session} {mouse} - interbout {THRESH_S}')
    
    #save figure
    fig2.savefig(mouse_path / f'{mouse}_threshold{THRESH_S}s_fiberbehav_scaled.pdf')
    
    return

def plot_fiberpho_behav_snip(behavprocess_df, timestart_camera):
    """
    Plots denoised deltaF/F aligned with behaviour (starts when camera starts)
    """
    
    fig3 = plt.figure(figsize=(20,12))
    ax2 = fig3.add_subplot(311)
    
    #y_scale = 7 #adjust according to data needs
    #y_shift = 30 #scale and shift are just for aesthetics
    
    #crops data when camera starts
    behavprocesssnip_df = behavprocess_df[behavprocess_df['Time(s)'] > timestart_camera]
    
    if session in ['Test 1h','Test 24h','Test','S3']:
        #plots fiberpho trace and behaviour
        p1, = ax2.plot('Time(s)', 'Denoised dFF', linewidth=.6, color='black', label='_GCaMP', data = behavprocesssnip_df)
        
        #makes areas corresponding to behaviours
        i = 0
        j = 0
        for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Exploration fam'].tolist()):
            if y == 1:
                x_start = x
            if y == -1 and x_start!=0:
                ax2.axvspan(x_start, x, facecolor='gold', alpha=0.5, label = '_'*i + 'Exploration fam')
                x_start=0
                i+=1
        for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Exploration new'].tolist()):
            if y == 1:
                x_start = x
            if y == -1 and x_start!=0:
                ax2.axvspan(x_start, x, facecolor='purple', alpha=0.5, label = '_'*j + 'Exploration new')
                x_start=0
                j+=1
        # for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Climbing'].tolist()):
        #     if y == 1:
        #         ax2.axvspan(x, x+0.1, facecolor='cornflowerblue', alpha=0.5)
        
    if session in ['Habituation','Training','S1']:
        #plots fiberpho trace and behaviour
        p1, = ax2.plot('Time(s)', 'Denoised dFF', linewidth=.6, color='black', label='_GCaMP', data = behavprocesssnip_df)
    
        #makes areas corresponding to behaviours
        i = 0
        j = 0
        for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Exploration left'].tolist()):
            if y == 1:
                x_start = x
            if y == -1 and x_start!=0:
                ax2.axvspan(x_start, x, facecolor='orange', alpha=0.5, label = '_'*i + 'Exploration left')
                x_start=0
                i+=1
        for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Exploration right'].tolist()):
            if y == 1:
                x_start = x
            if y == -1 and x_start!=0:
                ax2.axvspan(x_start, x, facecolor='darkturquoise', alpha=0.5, label = '_'*j + 'Exploration right')
                x_start=0
                j+=1
                j += 1
        # for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Climbing'].tolist()):
        #     if y == 1:
        #         ax2.axvspan(x, x+0.1, facecolor='cornflowerblue', alpha=0.5)
        
    if session in ['S2']:
        #plots fiberpho trace and behaviour
        p1, = ax2.plot('Time(s)', 'Denoised dFF', linewidth=.6, color='black', label='_GCaMP', data = behavprocesssnip_df)
    
        #makes areas corresponding to behaviours
        i = 0
        j = 0
        for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Exploration non social'].tolist()):
            if y == 1:
                x_start = x
            if y == -1 and x_start!=0:   
                ax2.axvspan(x_start, x, facecolor='grey', alpha=0.5, label = '_'*i + 'Exploration non social')
                x_start=0
                i += 1
        for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Exploration social'].tolist()):
            if y == 1:
                x_start = x
            if y == -1 and x_start!=0:   
                ax2.axvspan(x_start, x, facecolor='mediumvioletred', alpha=0.5, label = '_'*j + 'Exploration social')
                x_start=0
                j += 1
        # for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Climbing'].tolist()):
        #     if y == 1:
        #         ax2.axvspan(x, x+0.1, facecolor='cornflowerblue', alpha=0.5)
            
    #makes vertical line for opening of gate
    x_entry = fiberbehav_df.at[int(np.where(fiberbehav_df['Gate opens'] == 1)[0][0]), 'Time(s)']
    ax2.axvline(x_entry, color='lightsteelblue', ls = '--', label = 'Gate opens' )
            
    #makes vertical line for entry in open field
    x_entry = fiberbehav_df.at[int(np.where(fiberbehav_df['Entry in arena'] == 1)[0][0]), 'Time(s)']
    ax2.axvline(x_entry, color='slategrey', ls = '--', label = 'Entry OF' )
       
    ax2.set_ylabel(r'$\Delta$F/F')
    ax2.set_xlabel('Seconds')
    ax2.set_ylim(-1,2)
    ax2.legend(loc = 'upper right')
    ax2.margins(0.01,0.03)
    ax2.set_title(f'dFF with Behavioural Scoring - {exp} {session} {mouse} - interbout {THRESH_S}')
    
    #save figure
    fig3.savefig(mouse_path / f'{mouse}_threshold{THRESH_S}s_fiberbehavsnip_scaled.pdf')
    
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
    print('PETH')
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
        print(start,end)
        if end-start<EVENT_TIME_THRESHOLD and end-start>1:
            print(end-start)
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
    
########
#SCRIPT#
########

#%%Run test

mouse_path = analysis_path / '20211004_OdDis/Test/HFDm3'
mouse = str(mouse_path).split('/')[-1]

behav_path = f'{mouse_path}\\behav_{mouse}.csv'
fiberpho_path = f'{mouse_path}\\{mouse}_dFFfilt.csv'
camera_path = f'{mouse_path}\\{mouse}_camera.csv'
rawdata_path = f'{mouse_path}\\{mouse}.csv'
    
behav10Sps = pd.read_csv(behav_path)
fiberpho = pd.read_csv(fiberpho_path)
camera = pd.read_csv(camera_path)
#rawdata_df = pd.read_csv(rawdata_path)


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
#plot_rawdata(rawdata_df)

#plot isosbestic and gcamp data
plot_fiberpho(fiberbehav_df)

#plot fiberpho data aligned with behav
plot_fiberpho_behav(fiberbehav2_df)
#plot_fiberpho_behav_snip(fiberbehav2_df, timestart_camera)

for BOI in list_BOI:
    if BOI == 'Entry in arena' or BOI == 'Gate opens':
        PETH_data = PETH(fbprocess_df, BOI, 'onset', [6,10])
        plot_PETH_average(PETH_data, BOI, 'onset', [6,10])
    else:
        for event, timewindow in zip(list_EVENT,list_TIMEWINDOW): #CAUTION : adapt function!!
            PETH_data = PETH(fbprocess_df, BOI, event, timewindow)
            plot_PETH(PETH_data, BOI, event, timewindow)

#%%Run for all 

#list experiments, sessions and mice:
for exp_path in Path(analysis_path).iterdir():
    if exp_path.is_dir():
        exp = str(exp_path).split('\\')[-1]
        for session_path in Path(exp_path).iterdir():
            if session_path.is_dir():
                session = str(session_path).split('\\')[-1]
                print(exp, session)
                code = session_code(session)
                #get data path related to the task in protocol excel file
                data_path_exp = data_path / proto_df.loc[proto_df['Task']==exp, 'Data_path'].values[0]
                #create repository for values of thresholds : length and interbout
                repo_path = session_path / f'length{EVENT_TIME_THRESHOLD/10}_interbout{THRESH_S}'
                if not os.path.exists(repo_path):
                    os.mkdir(repo_path)
                    for subject in subjects_df['Subject']:
                        os.mkdir(repo_path / subject)
                #create list of mean_dFFs and diff_dFFs
                mean_dFFs_list = []
                subject_list = []
                diffdFF_list = []
                for mouse_path in Path(repo_path).iterdir():
                    # '/' on mac, '\\' on windows
                    mouse = str(mouse_path).split('\\')[-1]
                    print(mouse)
                    
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
                    if not os.path.exists(mouse_path / f'{mouse}_threshold{THRESH_S}s_fiberbehav_scaled.pdf'):
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
                        if EVENT_TIME_THRESHOLD==0 and THRESH_S==0:
                            mean_dFF_df = meandFF_behav(list_BOI, fiberbehav_df)
                            mean_dFFs_list.append(mean_dFF_df)
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
                            if BOI == 'Entry in arena' or BOI == 'Gate opens':
                                PETH_data = PETH(fbprocess_df, BOI, 'onset', [6,10])
                                plot_PETH_average(PETH_data, BOI, 'onset', [6,10])
                            else:
                                for event, timewindow in zip(list_EVENT,list_TIMEWINDOW): #CAUTION : adapt function!!
                                    PETH_data = PETH(fbprocess_df, BOI, event, timewindow)
                                    plot_PETH(PETH_data, BOI, event, timewindow)
               
                if EVENT_TIME_THRESHOLD==0 and THRESH_S==0 and mean_dFFs_list != []: 
                    meandFFs_allmice = pd.concat(mean_dFFs_list)
                    meandFFs_allmice.to_excel(repo_path / f'{exp}_{session}_meandFFs.xlsx')
                    
                if diffdFF_list != []:
                    diffdFFs_allmice = pd.concat(diffdFF_list)
                    diffdFFsmean = diffdFFs_allmice.groupby(['Subject','Behaviour'], as_index=False).mean()
                
                    diffdFFs_allmice.to_excel(repo_path / f'{exp}_{session}_diffdFFsall.xlsx')
                    diffdFFsmean.to_excel(repo_path / f'{exp}_{session}_diffdFFsmean.xlsx')
                                    
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

#%%test filter
fiberdata = pd.read_csv('D:\\Alice\\Fiber\\202110_CA2db2\\Data\\To_do\\20211110_AliceF_CA2b2beddingRimo\\CDf1_1_dFFfilt.csv')
fiberpho = fiberdata['Analog In. | Ch.1 470 nm (Deinterleaved)_dF/F0-Analog In. | Ch.1 405 nm (Deinterleaved)_dF/F0_LowPass']
samplingrate = 1000/(fiberdata.loc[1000,'Time(s)']-fiberdata.loc[0,'Time(s)'])
sos = signal.butter(5, 1, btype='low', analog=False, output='sos', fs=samplingrate)
filtered = signal.sosfilt(sos, fiberpho)

fiberdata['Filtered_dFF'] = filtered

fig8 = plt.figure(figsize=(20,4))
ax81 = fig8.add_subplot(111)
p81, = ax81.plot('Time(s)', 'Filtered_dFF',
                linewidth=.6, color='black', label='_GCaMP', data = fiberdata)