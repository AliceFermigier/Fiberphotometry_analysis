#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  4 11:53:48 2021

Useful functions for fiberpho analysis

@author: alice
"""

#%%
##########
#IMPORTED#
##########

import pandas as pd
import numpy as np
import math
from scipy import signal

#%%
###################
#DEFINED FUNCTIONS#
###################

def session_code(session,exp):
    """
    Generate session code in file name
    """
    if exp=='NewContext':
        code = '0'
    else:
        if session in ['Habituation','Training','S1','Conditioning']:
            code = '0'
        elif session in ['S2','Test 1h','Test']:
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

# def timestamp_camera(camera) :
#     """
#     Function to extract the timestamps where the camera starts and stops
#     --> Parameters
#         camera : pd dataframe, camera I/O with sample rate = 12kSps
#     --> Returns
#         (camera_start, camera_stop) = timestamp when camera starts and stops in seconds (truncated to 0,1s) #camera_stop à enlever si pas besoin
#     """
#     ind_list = np.where(camera['Digital I/O | Ch.3 DI/O-3'] == 1)
#     ind_list = ind_list[0].tolist()
#     (ind_start, ind_stop) = (ind_list[0],ind_list[len(ind_list)-1])
#     return (truncate(camera.at[ind_start, 'Time(s)'], 1),
#             truncate(camera.at[ind_stop, 'Time(s)'], 1))

def timestamp_camera(rawdata_df) :
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

# def filter_dFF(fiberbehav_df, ORDER, CUT_FREQ):
#     """
#     Apply additional filter to dFF data
#     """
#     fiberpho = fiberbehav_df['Denoised dFF']
#     samplingrate = 1000/(fiberbehav_df.loc[1000,'Time(s)']-fiberbehav_df.loc[0,'Time(s)'])
#     sos = signal.butter(ORDER, CUT_FREQ, btype='low', analog=False, output='sos', fs=samplingrate)
#     filtered_data = signal.sosfilt(sos, fiberpho)
    
#     filtered_df = fiberbehav_df
#     filtered_df['Denoised dFF'] = filtered_data
    
#     return(filtered_df)

def meandFF_behav(list_BOI, fiberbehav_df, exp, session, mouse, group):
    """
    Calculates mean dFF during each behaviour
    Output : dataframe with 'Subject','Group','Mean dFF','Baseline', 'Post_baseline', and mean dFF for each behaviour
    """
    list_behav_analyzed = []
    list_meandFF = []
    
    fiberbehavsnip_df = fiberbehav_df[fiberbehav_df['Time(s)'] > 15]
    
    #get index of when the gate opens
    if exp == 'NewContext':
        ind_start_trial = fiberbehavsnip_df.index[fiberbehavsnip_df['New context'] == 1].tolist()[0]
    elif exp == 'Fear' and session == 'Test':
        ind_start_trial = fiberbehavsnip_df.index[fiberbehavsnip_df['Fear cage'] == 1].tolist()[0]
    elif exp == 'Fear' and session == 'Conditioning':
        ind_start_trial = fiberbehavsnip_df.index[fiberbehavsnip_df['Shock'] == 1].tolist()[0]
    else:
        ind_start_trial = fiberbehavsnip_df.index[fiberbehavsnip_df['Entry in arena'] == 1].tolist()[0]
    
    #create list of behaviours and list of corresponding mean dFFs
    for behav in list_BOI:
        if fiberbehavsnip_df[behav].sum() > 2:
            list_behav_analyzed.append(behav)
            meandFF_behav_df = fiberbehavsnip_df.groupby([behav], as_index=False).mean()
            print(meandFF_behav_df)
            list_meandFF.append(meandFF_behav_df.loc[meandFF_behav_df[behav]==1, 'Denoised dFF'].values[0])
                                
    #calculate mean dFF during baseline
    
    meandFF_baseline = fiberbehavsnip_df.loc[:ind_start_trial, 'Denoised dFF'].mean()
    
    #calculate mean dFF when no exploration and after gate opens
    
    meandFF_postbaseline_df = fiberbehavsnip_df.loc[ind_start_trial:]
    for behav in list_behav_analyzed:
        meandFF_postbaseline_df = meandFF_postbaseline_df.loc[meandFF_postbaseline_df[behav]==0]
        
    meandFF_postbaseline = meandFF_postbaseline_df['Denoised dFF'].mean()
    
    #calculate mean dFF when no exploration on total trial
    meandFF_df = fiberbehavsnip_df
    for behav in list_behav_analyzed:
        meandFF_df = meandFF_df.loc[meandFF_df[behav]==0]
        
    meandFF = meandFF_df['Denoised dFF'].mean()
    
    #create dataframe with values
    
    list_dFFs = [mouse, group, meandFF,meandFF_baseline,meandFF_postbaseline]
    list_columns = ['Subject','Group','Mean dFF','Baseline', 'Post_baseline']
    for (behav,meandFF) in zip(list_behav_analyzed,list_meandFF):
        list_dFFs.append(meandFF)
        list_columns.append(behav)
    meandFFs_df = pd.DataFrame(data=[list_dFFs], columns=list_columns)
    print(meandFFs_df)
    
    # #Create figure of mean dFFs
    # fig7 = plt.figure(figsize=(7,6))
    # ax71 = fig7.add_subplot(111)
    # print('x = ',meandFFs_df.columns.to_list()[1:])
    # print('y = ',meandFFs_df.iloc[0, 1:].to_list())
    # p71 = ax71.bar(meandFFs_df.columns.to_list()[1:], meandFFs_df.iloc[0, 1:])
    # ax71.axhline(y=0, linewidth=.6, color='black')
    # ax71.set_ylabel(r'$\Delta$F/F')
    # ax71.set_title(f'Mean dFF - {exp} {session} {mouse}')
    # #save figure
    # fig7.savefig(mouse_path / f'{mouse}_meandFFs.pdf')
    
    meandFFs_df.to_excel(mouse_path / f'{mouse}_meandFFglob.xlsx')
    
    return(meandFFs_df)

def diff_dFF(fiberbehav_df, behavprocess_df, list_BOI):
    """
    Calculates dFF difference between beginning and end of bouts
    Also calculates mean dFF during each bout
    """
    list_behav_analyzed = []
    for behav in list_BOI:
        if fiberbehav_df[behav].sum() > 3:
            list_behav_analyzed.append(behav)
    list_deltadFF = []
    list_meandFF = []
    listind_behav = []
    listind_bouts = []
    
    for i, behav in enumerate(list_BOI):
        list_starts = np.where(behavprocess_df[behav]==1)[0].tolist()
        list_stops = np.where(behavprocess_df[behav]==-1)[0].tolist()
        for (start,end) in zip(list_starts, list_stops):
            if end-start<EVENT_TIME_THRESHOLD*10 and end-start>1:
                list_starts.remove(start)
                list_stops.remove(end)
        bout_n = 1
        for (start, stop) in zip(list_starts, list_stops):
            mean_start = behavprocess_df.loc[start-5:start+5, 'Denoised dFF'].mean()
            mean_stop = behavprocess_df.loc[stop-5:stop+5, 'Denoised dFF'].mean()
            delta = mean_stop-mean_start
            #delta = behavprocess_df.loc[stop, 'Denoised dFF']-behavprocess_df.loc[start, 'Denoised dFF']
            mean = behavprocess_df.loc[start:stop, 'Denoised dFF'].mean()
            list_deltadFF.append(delta)
            list_meandFF.append(mean)
            listind_behav.append(behav)
            listind_bouts.append(bout_n)
            bout_n+=1
    
    diffdFF_df = pd.DataFrame(data = {'Subject':[mouse]*len(listind_behav), 'Group':[group]*len(listind_behav),
                                      'Behaviour':listind_behav, 'Bout':listind_bouts, 'Mean dFF':list_meandFF,
                                      'Delta dFF':list_deltadFF})
    
    return(diffdFF_df)