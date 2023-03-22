# -*- coding: utf-8 -*-
"""
Created on Fri Mar 19 13:58:32 2021

Processing of fiberpho data
Full trace dFF according to Lerner et al. 2015
https://www.tdt.com/support/python-sdk/offline-analysis-examples/licking-bout-epoc-filtering/

@author: Alice Fermigier
"""

#%%
##########
#IMPORTED#
##########

import pandas as pd
import numpy as np
import math
import os
import sys
from scipy import signal
from ast import literal_eval
from pathlib import Path

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

def deinterleave(rawdata_df):
    """
    deinterleave signal
    """
    derivative405 = rawdata_df['DI/O-1'].diff()
    derivative470 = rawdata_df['DI/O-2'].diff()
    
    list_405 = []
    list_470 = []
    for i in np.where(derivative405==1)[0]:
        if i+250 > len(rawdata_df):
            list_405.append(rawdata_df.loc[len(rawdata_df)-1,'AIn-1'])
        else:
            list_405.append(rawdata_df.loc[i+250,'AIn-1'])
    for i in np.where(derivative470==1)[0]:
        if i+250 > len(rawdata_df):
            list_405.append(rawdata_df.loc[len(rawdata_df)-1,'AIn-1'])
        else:
            list_470.append(rawdata_df.loc[i+250,'AIn-1'])
        
    while len(list_470) > len(list_405):
        list_470.pop()
        
    while len(list_405) > len(list_470):
        list_405.pop()

    timevector = np.linspace(0,rawdata_df['Time(s)'].max(),len(list_405))
    
    deinterleaved_df = pd.DataFrame(data = {'Time(s)' : timevector, '405 Deinterleaved' : list_405, '470 Deinterleaved' : list_470})
    deinterleaved_df.replace(0,value=None,inplace=True)
        
    return(deinterleaved_df)

def filter_dFF(data_df, ORDER, CUT_FREQ, columns):
    """
    Apply additional filter to dFF data
    """
    filtered_df = data_df
    for col in columns:
        fiberpho = data_df[col]
        samplingrate = 1000/(data_df.loc[1000,'Time(s)']-data_df.loc[0,'Time(s)'])
        sos = signal.butter(ORDER, CUT_FREQ, btype='low', analog=False, output='sos', fs=samplingrate)
        filtered_data = signal.sosfilt(sos, fiberpho)
        
        filtered_df[col] = filtered_data

    return(filtered_df)

def samplerate(data_df):
    
    sr = len(data_df)/(data_df['Time(s)'].max()-data_df['Time(s)'].min())
    
    return(sr)

def rem_artifacts(data_df,artifacts_df,filecode,sr):
    """
    rawdata_df = pandas dataframe from csv file
    artifacts_df = df from excel file with artifacts to remove
    columns = list of columns with signals to process
    """
    for col in data_df.columns[1:] :
        #remove 5 first seconds of data
        data_df.loc[0:round(5*sr),col]=np.nan
        if filecode in artifacts_df['File']:
            list_artifacts = artifacts_df.loc[artifacts_df['File']==filecode,'Artifacts'].values
            for [x_start, x_stop] in literal_eval(list_artifacts[0]):
                data_df.loc[round(x_start*sr):round(x_stop*sr),col]=np.nan
    
    return(data_df)

def dFF(data_df,artifacts_df,filecode,sr,method='mean'):
    """
    columns = list of columns with signals to process, 405 and 470 in that order
    """
    dFFdata = np.full([3,len(data_df)], np.nan)
    for (i,col) in enumerate(data_df.columns[1:]):
        #compute dFF specifically for separate tranches if artifacts are present
        if filecode in artifacts_df['File']:
            list_artifacts = artifacts_df.loc[artifacts_df['File']==filecode,'Artifacts'].values
            for [x_start, x_stop] in literal_eval(list_artifacts[0]):
                meanF = np.nanmean(data_df.loc[x_start:x_stop,col])
                dFFdata[i][x_start:x_stop] = [(j-meanF)/meanF for j in data_df.loc[x_start:x_stop,col] if j!=np.nan]
        else:
            meanF = np.nanmean(data_df[col])
            dFFdata[i] = [(j-meanF)/meanF for j in data_df[col] if j!=np.nan]
    
    if method == 'mean':
        dFFdata[2]=dFFdata[1]-dFFdata[0]
    if method == 'mean fit':
        dFFdata[2]=dFFdata[1]-dFFdata[0]
    
    dFFdata_df = pd.DataFrame(data = {'Time(s)':data_df['Time(s)'],'405 dFF':dFFdata[0],
                                      '470 dFF':dFFdata[1],'Denoised dFF':dFFdata[2]})
    
    return(dFFdata_df)

def interpolate_dFFdata(data_df, method='linear'):
    """
    data_df = dFFdata_df with nan values after removing artifacts
    """
    for (i,col) in enumerate(data_df.columns[1:]):
        meandFF = np.nanmean(data_df[col])
        if method == 'linear':
            data_df.loc[0,col]=meandFF
            data_df.interpolate(method='linear',inplace=True) #interpolate data linearly
        elif method == 'mean pad':
            data_df.fillna(meandFF,inplace=True)
        
    return(data_df)
    
#%%
########
#SCRIPT#
########
    
# data_path = Path('K:\\Alice\\Fiber\\202301_CA2b5\\Data\\20230222_AliceF_CA2b5OdDispostshock')
# for mouse in subjects_df['Subject']:
#     print(mouse)
#     rawdata_df = pd.read_csv(data_path/f'{mouse}_1.csv',skiprows=1,usecols=['Time(s)','AIn-1','DI/O-1','DI/O-2'])
#     deinterleaved_df = deinterleave(rawdata_df)
    
#     fig5 = plt.figure(figsize=(10,6))
#     ax7 = fig5.add_subplot(211)

#     p1, = ax7.plot('Time(s)', '470 Deinterleaved', 
#                    linewidth=1, color='deepskyblue', label='GCaMP', data = deinterleaved_df[50:]) 
#     p2, = ax7.plot('Time(s)', '405 Deinterleaved', 
#                    linewidth=1, color='blueviolet', label='ISOS', data = deinterleaved_df[50:])
    
#     ax7.set_ylabel('V')
#     ax7.set_xlabel('Time(s)')
#     ax7.legend(handles=[p1,p2], loc='upper right')
#     ax7.margins(0.01,0.3)
#     ax7.set_title(f'GCaMP and Isosbestic raw traces - {mouse}')


