# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 15:51:56 2023

Functions for preprocessing fiberphotometry data

@author: alice fermigier
"""

#%%
##########
#IMPORTED#
##########

import pandas as pd
import numpy as np
from scipy import signal
from ast import literal_eval


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
        
    return deinterleaved_df


def samplerate(data_df):
    
    sr = len(data_df)/(data_df['Time(s)'].max()-data_df['Time(s)'].min())
    
    return sr

def controlFit(control, signal):
    
	p = np.polyfit(control, signal, 1)
	arr = (p[0]*control)+p[1]
	return arr

def dFF(data_df,artifacts_df,filecode,sr,method='mean'):
    """
    Calculates dFF and removes artifacts
    Output : pandas dataframe with 405 dFF, 470 dFF and Denoised dFF
    """   
    sr = round(sr)
    if method == 'mean':
        dFFdata = np.full([3,len(data_df)], np.nan)
        for (i,col) in enumerate(data_df.columns[1:]):
            #compute dFF specifically for separate tranches if artifacts are present
            if filecode in artifacts_df['Filecode'].values:
                print('Removing artifacts in dFF')
                list_artifacts = artifacts_df.loc[artifacts_df['Filecode']==filecode,'Artifacts'].values
                begin=round(5*sr)
                for k in range(len(literal_eval(list_artifacts[0]))):
                    [x_start, x_stop] = literal_eval(list_artifacts[0])[k]
                    end = data_df.loc[data_df['Time(s)']>x_start].index[0]
                    meanF = np.nanmean(data_df.loc[begin+1:end,col])
                    dFFdata[i][begin+1:end] = [(j-meanF)/meanF*100 for j in data_df.loc[begin+1:end-1,col]]
                    begin= data_df.loc[data_df['Time(s)']<x_stop].index[-1]
                end=len(data_df)
                meanF = np.nanmean(data_df.loc[begin+1:end,col])
                dFFdata[i][begin+1:end] = [(j-meanF)/meanF*100 for j in data_df.loc[begin+1:end-1,col]]
            else:
                meanF = np.nanmean(data_df[col])
                dFFdata[i] = [(j-meanF)/meanF*100 for j in data_df[col] if j!=np.nan]
        
        dFFdata[2]=dFFdata[1]-dFFdata[0]
    
    elif method == 'fit':
        dFFdata = np.full([3,len(data_df)], np.nan)
        #compute dFF specifically for separate tranches if artifacts are present
        if filecode in artifacts_df['Filecode'].values:
            print('Removing artifacts in dFF')
            list_artifacts = artifacts_df.loc[artifacts_df['Filecode']==filecode,'Artifacts'].values
            begin=round(5*sr)
            for k in range(len(literal_eval(list_artifacts[0]))):
                [x_start, x_stop] = literal_eval(list_artifacts[0])[k]
                end = data_df.loc[data_df['Time(s)']>x_start].index[0]
                print(data_df.loc[begin+1:end-1,'470 Deinterleaved'])
                dFFdata[0][begin+1:end] = controlFit(data_df.loc[begin+1:end-1,'405 Deinterleaved'], 
                                                     data_df.loc[begin+1:end-1,'470 Deinterleaved'])
                dFFdata[1][begin+1:end] = data_df.loc[begin+1:end-1,'470 Deinterleaved']
                print(begin,end,dFFdata[0][begin+1:end])
                begin= data_df.loc[data_df['Time(s)']<x_stop].index[-1]
            end=len(data_df)
            dFFdata[0][begin+1:end] = controlFit(data_df.loc[begin+1:end-1,'405 Deinterleaved'], 
                                                 data_df.loc[begin+1:end-1,'470 Deinterleaved'])
            dFFdata[1][begin+1:end] = data_df.loc[begin+1:end-1,'470 Deinterleaved']
        else:
            [begin,end] = [round(5*sr),len(data_df)]
            dFFdata[0][begin+1:end] = controlFit(data_df.loc[begin+1:end-1,'405 Deinterleaved'], 
                                                 data_df.loc[begin+1:end-1,'470 Deinterleaved'])
            dFFdata[1][begin+1:end] = data_df.loc[begin+1:end-1,'470 Deinterleaved']
        
        res = np.subtract(dFFdata[1], dFFdata[0])
        normData = np.divide(res, dFFdata[0])
        dFFdata[2] = normData*100
    
    dFFdata_df = pd.DataFrame(data = {'Time(s)':data_df['Time(s)'],'405 dFF':dFFdata[0],
                                      '470 dFF':dFFdata[1],'Denoised dFF':dFFdata[2]})
    
    return dFFdata_df

def interpolate_dFFdata(data_df, method='linear'):
    """
    data_df = dFFdata_df with nan values after removing artifacts
    method : linear will connect datapoints with segment ; mean pad will fill missing data with mean dFF
    """
    for (i,col) in enumerate(data_df.columns[1:]):
        meandFF = np.nanmean(data_df[col])
        if method == 'linear':
            data_df.loc[0,col]=meandFF
            data_df.interpolate(method='linear',inplace=True) #interpolate data linearly
        elif method == 'mean pad':
            data_df.fillna(meandFF,inplace=True) #pad missing data with mean dFF
        
    return data_df

def butterfilt(data_df, ORDER, CUT_FREQ):
    """
    Butterworth filter
    """
    filtered_df = data_df
    for col in data_df.columns[1:]:
        fiberpho = data_df[col]
        samplingrate = 1000/(data_df.loc[1000,'Time(s)']-data_df.loc[0,'Time(s)'])
        sos = signal.butter(ORDER, CUT_FREQ, btype='low', output='sos', fs=samplingrate)
        filtered_data = signal.sosfilt(sos, fiberpho)
        
        filtered_df[col] = filtered_data

    return filtered_df

def smoothing_SMA(data_df,win_size):
    """
    Simple moving average (SMA) smoothing
    win_size : size of the moving window
    """
    for col in data_df.columns[1:]:
        data_df[col] = data_df[col].rolling(win_size, min_periods=1).mean()
        
    return data_df