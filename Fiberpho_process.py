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
    
    old_samplerate = len(rawdata_df)/rawdata_df.loc[len(rawdata_df)-1,'Time(s)']
    new_samplerate = len(np.where(derivativeAIn1==1))/rawdata_df.loc[len(rawdata_df)-1,'Time(s)']
    
    
    
    for i in range (0,len(rawdata_df),new_samplerate):
        
        
    return(deinterleaved_df)

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

def rem_artifacts(rawdata_df):
    """
    rawdata_df = pandas dataframe from csv file
    """
    
    
    return(data)

def rem_photobleach(data):
    """
    data = np array
    """
    
    return(filt_data)

def dFF(filt_data):
    """
    filt_data = np array
    """
    
    return(filt_dFF)


#%%
########
#SCRIPT#
########
import matplotlib.pyplot as plt

data_df = pd.read_csv('K:\\Alice\\Fiber\\202301_CA2b5\\Data\\20230112_AliceF_CA2b5Essais2\\A1f_0.csv', skiprows=1)

data_df['AIn-1 x DI/O-1'].replace(0.0,inplace=True)
data_df['AIn-1 x DI/O-2'].replace(0.0,inplace=True)

plt.plot('Time(s)','AIn-1 x DI/O-2',data = data_df[300:1800])

data_df = data_df.loc[[i for i in range(0,len(data_df),1207)]]

derivative405 = data_df['DI/O-1'].diff()
old_samplerate = len(data_df)/data_df['Time(s)'].max()
new_samplerate = len(np.where(derivative405==1)[0])/data_df['Time(s)'].max()

#remove deinterleaving artifacts
for i in data_df['AIn-1 x DI/O-1']:
    




