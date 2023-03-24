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
from scipy.sparse import csc_matrix, eye, diags
from scipy.sparse.linalg import spsolve
from ast import literal_eval
from pathlib import Path

#Custom
#put path to directory where python files are stored
sys.path.append('C:\\Users\\Alice\\Documents\\GitHub\\Fiberphotometry_analysis')

#%%
########
#LOADER#
########

from Fiberpho_loader import experiment_path, analysis_path, subjects_df, SAMPLERATE, artifacts_df
from Fiberpho_loader import list_EVENT, list_TIMEWINDOW, PRE_EVENT_TIME, TIME_BEGIN

os.chdir(experiment_path)

#%%
###################
#DEFINED FUNCTIONS#
###################

###Smoothing methods

#airPLS

def WhittakerSmooth(x,w,lambda_,differences=1):
    '''
    Penalized least squares algorithm for background fitting
    airPLS.py Copyright 2014 Renato Lombardo - renato.lombardo@unipa.it
    
    input
        x: input data (i.e. chromatogram of spectrum)
        w: binary masks (value of the mask is zero if a point belongs to peaks and one otherwise)
        lambda_: parameter that can be adjusted by user. The larger lambda is,  the smoother the resulting background
        differences: integer indicating the order of the difference of penalties
    
    output
        the fitted background vector
    '''
    X=np.matrix(x)
    m=X.size
    E=eye(m,format='csc')
    for i in range(differences):
        E=E[1:]-E[:-1] # numpy.diff() does not work with sparse matrix. This is a workaround.
    W=diags(w,0,shape=(m,m))
    A=csc_matrix(W+(lambda_*E.T*E))
    B=csc_matrix(W*X.T)
    background=spsolve(A,B)
    return np.array(background)

def airPLS(x, lambda_=100, porder=1, itermax=15):
    '''
    Adaptive iteratively reweighted penalized least squares for baseline fitting
    airPLS.py Copyright 2014 Renato Lombardo - renato.lombardo@unipa.it
    
    input
        x: input data (i.e. chromatogram of spectrum)
        lambda_: parameter that can be adjusted by user. The larger lambda is,  the smoother the resulting background, z
        porder: adaptive iteratively reweighted penalized least squares for baseline fitting
    
    output
        the fitted background vector
    '''
    m=x.shape[0]
    w=np.ones(m)
    for i in range(1,itermax+1):
        z=WhittakerSmooth(x,w,lambda_, porder)
        d=x-z
        dssn=np.abs(d[d<0].sum())
        if(dssn<0.001*(abs(x)).sum() or i==itermax):
            if(i==itermax): print('WARING max iteration reached!')
            break
        w[d>=0]=0 # d>0 means that this point is part of a peak, so its weight is set to 0 in order to ignore it
        w[d<0]=np.exp(i*np.abs(d[d<0])/dssn)
        w[0]=np.exp(i*(d[d<0]).max()/dssn) 
        w[-1]=w[0]
    return z

def smoothing_airPLS(data_df):
    for col in data_df.columns[1:]:
        data_df[col] = airPLS(data_df[col], lambda_=2, porder=1, itermax=15)
        
    return(data_df)

#butterworth filter

def filter_dFF(data_df, ORDER, CUT_FREQ):
    """
    Apply additional filter to dFF data
    """
    filtered_df = data_df
    for col in data_df.columns[1:]:
        fiberpho = data_df[col]
        samplingrate = 1000/(data_df.loc[1000,'Time(s)']-data_df.loc[0,'Time(s)'])
        sos = signal.butter(ORDER, CUT_FREQ, btype='low', analog=False, output='sos', fs=samplingrate)
        filtered_data = signal.sosfilt(sos, fiberpho)
        
        filtered_df[col] = filtered_data

    return(filtered_df)

#simple moving average (SMA) smoothing

def smoothing_SMA(data_df):
    for col in data_df.columns[1:]:
        data_df[col] = data_df[col].rolling(7, min_periods=1).mean()
        
    return(data_df)

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
        if filecode in artifacts_df['Filecode'].values:
            list_artifacts = artifacts_df.loc[artifacts_df['Filecode']==filecode,'Artifacts'].values
            for [x_start, x_stop] in literal_eval(list_artifacts[0]):
                data_df.loc[(data_df['Time(s)']>x_start)&(data_df['Time(s)']<x_stop),col]=np.nan
    
    return(data_df)

def detrend(data_df,artifacts_df,filecode,sr):
    """
    
    """   
    detrended_data = np.full([2,len(data_df)], np.nan)
    for (i,col) in enumerate(data_df.columns[1:]):
        #compute dFF specifically for separate tranches if artifacts are present
        if filecode in artifacts_df['Filecode'].values:
            print('Removing artifacts in detrending')
            list_artifacts = artifacts_df.loc[artifacts_df['Filecode']==filecode,'Artifacts'].values
            begin=round(5*sr)
            for k in range(len(literal_eval(list_artifacts[0]))):
                [x_start, x_stop] = literal_eval(list_artifacts[0])[k]
                end = data_df.loc[data_df['Time(s)']>x_start].index[0]
                detrended_data[i][begin+1:end] = airPLS(data_df.loc[begin+1:end-1,col], lambda_=100, porder=1, itermax=15)
                begin= data_df.loc[data_df['Time(s)']<x_stop].index[-1]
            end=len(data_df)
            detrended_data[i][begin+1:end] = airPLS(data_df.loc[begin+1:end-1,col], lambda_=100, porder=1, itermax=15)
        else:
            detrended_data[i] = airPLS(data_df.loc[begin+1:end-1,col], lambda_=100, porder=1, itermax=15)
    
    detrended_df = pd.DataFrame(data = {'Time(s)':data_df['Time(s)'],'405 Deinterleaved':detrended_data[0],
                                      '470 Deinterleaved':detrended_data[1]})
    
    return(detrended_df)

def controlFit(control, signal):
    
	p = np.polyfit(control, signal, 1)
	arr = (p[0]*control)+p[1]
	return arr

def dFF(data_df,artifacts_df,filecode,sr,method='mean'):
    """
    
    """   
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
                    dFFdata[i][begin+1:end] = [((j-meanF)/meanF)*100 for j in data_df.loc[begin+1:end-1,col]]
                    begin= data_df.loc[data_df['Time(s)']<x_stop].index[-1]
                end=len(data_df)
                meanF = np.nanmean(data_df.loc[begin+1:end,col])
                dFFdata[i][begin+1:end] = [((j-meanF)/meanF)*100 for j in data_df.loc[begin+1:end-1,col]]
            else:
                meanF = np.nanmean(data_df[col])
                dFFdata[i] = [((j-meanF)/meanF)*100 for j in data_df[col] if j!=np.nan]
        
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
                print(data_df.loc[begin+1:end-1,'405 Deinterleaved'])
                dFFdata[0][begin+1:end] = controlFit(data_df.loc[begin+1:end-1,'405 Deinterleaved'], data_df.loc[begin+1:end-1,'470 Deinterleaved'])
                dFFdata[1][begin+1:end] = data_df.loc[begin+1:end-1,'470 Deinterleaved']
                begin= data_df.loc[data_df['Time(s)']<x_stop].index[-1]
            end=len(data_df)
            dFFdata[0][begin+1:end] = controlFit(data_df.loc[begin+1:end-1,'405 Deinterleaved'], data_df.loc[begin+1:end-1,'470 Deinterleaved'])
            dFFdata[1][begin+1:end] = data_df.loc[begin+1:end-1,'470 Deinterleaved']
        else:
            [begin,end] = [round(5*sr),len(data_df)]
            dFFdata[0][begin+1:end] = controlFit(data_df.loc[begin+1:end-1,'405 Deinterleaved'], data_df.loc[begin+1:end-1,'470 Deinterleaved'])
            dFFdata[1][begin+1:end] = data_df.loc[begin+1:end-1,'470 Deinterleaved']
        
        res = np.subtract(dFFdata[1], dFFdata[0])
        normData = np.divide(res, dFFdata[0])
        dFFdata[2] = normData*100
    
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
            data_df.fillna(meandFF,inplace=True) #pad missing data with mean dFF
        
    return(data_df)
    
#%%
########
#SCRIPT#
########

import matplotlib.pyplot as plt
    
data_path = Path('K:\\Alice\\Fiber\\202301_CA2b5\\Data\\20230222_AliceF_CA2b5OdDispostshock')
exp = 'OdDispostshock'
session = 'Test'
for mouse in ['A1f']: #subjects_df['Subject']:
    print(mouse)
    rawdata_df = pd.read_csv(data_path/f'{mouse}_1.csv',skiprows=1,usecols=['Time(s)','AIn-1','DI/O-1','DI/O-2'])
    print('Deinterleaving')
    deinterleaved_df = deinterleave(rawdata_df)
    
    filecode = f'{exp}_{session}_{mouse}'
    sr = samplerate(deinterleaved_df)

    print('dFF')
    dFFdata_df = dFF(deinterleaved_df,artifacts_df,filecode,sr,method='fit')
    print('interpolate')
    interpdFFdata_df = interpolate_dFFdata(dFFdata_df, method='linear')
    print('smoothing')
    smootheddFFdata_df =smoothing_SMA(interpdFFdata_df)
    
fig3 = plt.figure(figsize=(10,6))
ax4 = fig3.add_subplot(211)
p1, = ax4.plot('Time(s)', '470 Deinterleaved', 
               linewidth=1, color='deepskyblue', label='GCaMP', data = deinterleaved_df) 
p2, = ax4.plot('Time(s)', '405 Deinterleaved', 
               linewidth=1, color='blueviolet', label='ISOS', data = deinterleaved_df)

ax4.set_ylabel('V')
ax4.set_xlabel('Time(s)')
ax4.legend(handles=[p1,p2], loc='upper right')
ax4.margins(0.01,0.3)
ax4.set_title(f'GCaMP and Isosbestic deinterleaved - {exp} {session} {mouse}')


fig5 = plt.figure(figsize=(10,6))
ax7 = fig5.add_subplot(211)

p1, = ax7.plot('Time(s)', '470 dFF', 
               linewidth=1, color='deepskyblue', label='GCaMP', data = smootheddFFdata_df) 
p2, = ax7.plot('Time(s)', '405 dFF', 
               linewidth=1, color='blueviolet', label='ISOS', data = smootheddFFdata_df)

ax7.set_ylabel('V')
ax7.set_xlabel('Time(s)')
ax7.legend(handles=[p1,p2], loc='upper right')
ax7.margins(0.01,0.3)
ax7.set_title(f'GCaMP and Isosbestic dFF traces - {exp} {session} {mouse}')

fig2 = plt.figure(figsize=(20,5))
ax1 = fig2.add_subplot(111)
p1, = ax1.plot('Time(s)', 'Denoised dFF', linewidth=.6, color='black', label='_GCaMP', data = smootheddFFdata_df)


