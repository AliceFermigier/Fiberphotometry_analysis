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
    
    list_405 = []
    list_470 = []
    for i in np.where(derivative405==1)[0]:
        list_405.append(rawdata_df.loc[i+250,'AIn-1'])
    for i in np.where(derivative470==1)[0]:
        list_470.append(rawdata_df.loc[i+250,'AIn-1'])
        
    while len(list_470) > len(list_405):
        list_470.pop()
        
    while len(list_405) > len(list_470):
        list_405.pop()

    timevector = np.linspace(0,rawdata_df['Time(s)'].max(),len(list_405))
    
    deinterleaved_df = pd.DataFrame(data = {'Time(s)' : timevector, '405 Deinterleaved' : list_405, '470 Deinterleaved' : list_470})
    deinterleaved_df.replace(0,value=None,inplace=True)
        
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

deinterleaved_df = deinterleave(data_df)
plt.plot('Time(s)','405 Deinterleaved',data=deinterleaved_df)
    




