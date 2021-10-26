# -*- coding: utf-8 -*-
"""
Created on Tue Oct 26 13:52:47 2021

Process and plot plethysmograph and fiberphotometry data

@author: alice fermigier
"""

#%%
##########
#IMPORTED#
##########

#Libraries
import numpy as np
import statistics as stat
from pathlib import Path
import glob
from scipy import signal
import os
import sys
# citing pyabf : 
# Harden, SW (2020). pyABF 2.2.3. [Online]. Available: https://pypi.org/project/pyabf/

#Custom
#put path to directory where python files are stored
sys.path.append('C:\\Users\\afermigier\\Documents\\GitHub\\Plethysmograph_analysis')
sys.path.append('C:\\Users\\afermigier\\Documents\\GitHub\\Fiberphotometry_analysis')

from PlethyL_process import r_slowvar_crop, find_stim, find_inspirations, resp_freq, resp_amp_AUC


#%%
########
#LOADER#
########

#Imported from param file
from Plethyfiber_params import Shift, SamplingRate, LOW_CUTOFF_FREQ, HIGH_CUTOFF_FREQ, StimStart_list, StimDurTh
from Plethyfiber_params import PostStimDur, minpeakdist, thresh, ArTfact, CROP

#Imported from fiberpho file
from Fiberpho_loader import experiment_path, analysis_path, data_path, subjects_df, SAMPLERATE, proto_df
from Fiberpho_loader import list_EVENT, list_TIMEWINDOW, PRE_EVENT_TIME, TIME_BEGIN, THRESH_S, EVENT_TIME_THRESHOLD


#%%
###################
#DEFINED FUNCTIONS#
###################

def r_highartifacts_fiberplethy(plethys_df):
    """
    Removes high artifacts
    --> Parameters : raw plethysmo abf data
    --> Return : corrected data
    """
    
    list_artefacts = []
    signal_median = stat.median(plethys_df['AIn-4'][0:StimStart_list[0]*SamplingRate])
    signal_std = np.std(plethys_df['AIn-4'][0:StimStart_list[0]*SamplingRate])
    print(signal_median)
    print(signal_std)
    
    for point in plethys_df['AIn-4']:
        if point > signal_median + 3*signal_std or point < signal_median - 3*signal_std:
            list_artefacts.append(False)
        else:
            list_artefacts.append(True)
            
    plethys_df_clean = pd.DataFrame({'Time(s)':plethys_df['Time(s)'], 'WBP(V)':plethys_df['AIn-4'].where(list_artefacts, signal_median)})
    
    return(plethys_df_clean)

def r_slowvar_crop(plethys_df):
    """
    Removes slow variations (~0.1Hz)
    Crops signal to remove artifact at beginning of timecourse
    --> Parameters : plethysmo abf data ; CROP : int
    --> Return : filtered data
    """
    ORDER = 6
    n=round(ORDER/2) #why divide the order?
    Wn=[LOW_CUTOFF_FREQ/(SamplingRate/2), HIGH_CUTOFF_FREQ/(SamplingRate/2)]
    
    sos = signal.butter(n, Wn, btype='bandpass', output='sos')
    plethys_ddp_fi = signal.sosfilt(sos, plethys_df['WBP(V)'])
    plethys_ddp_f = plethys_ddp_fi - plethys_df['WBP(V)'] #out of bandpass data
    
    return(plethys_ddp_fi, plethys_ddp_f)


def stim_list(plethy_df):
    """
    Find stims on traces and give timepoints (truncated to 0,1ms to match fiber data)
    """
    for stim_start in StimStart_list:
        

def plethyfiber_plot(fiberpho_df, plethy_df):
    """
    Take fiberpho data, plethysmo data and stims
    Plots data
    """


#%%
########
#SCRIPT#
########

#test_data
rawdata_path = data_path / '20211022_AliceF_CA2b2plethysmoNovel\\CDf1_0.csv'
fiberpho_path = data_path / '20211022_AliceF_CA2b2plethysmoNovel\\CDf1_0_dFFfilt.csv'

plethys_df = pd.read_csv(rawdata_path, skiprows=1, usecols=['Time(s)','AIn-4'])

#remove artifacts
plethys_df_clean = r_highartifacts_fiberplethy(plethys_df)
#filter and crop
(plethys_ddp_fi, plethys_ddp_f) = r_slowvar_crop(plethys_ddp)
#find start and stop of stims
(stim_start, stim_stop, stim_dur) = find_stim(plethys_ddp_f)
#find inspirations
(Ilocs, Ipks) = find_inspirations(plethys_ddp_fi)
#calculate frequency, amplitude and area under curve
(IFreqi, IFreqim_preStim, IFreqim_Stim, IFreqm_preStim, IFreqm_Stim) = resp_freq(plethys_ddp_fi, stim_start, stim_stop, Ilocs)
(Ipks_preStim, mAmp_preStim, IArea_preStim, Ipks_Stim, mAmp_Stim, IArea_Stim) = resp_amp_AUC(plethys_ddp_fi, stim_start, stim_stop, Ilocs, Ipks)
#Plot data with computed stim period
if not os.path.exists(f'{analysis_path}\\{subject}\\{session}_WBPstim.pdf'):
    fig3 = plot_procdatastim(plethys_ddp, stim_start, stim_stop, subject, session)
    #save figure
    fig3.savefig(f'{analysis_path}\\{subject}\\{session}_WBPstim.pdf')
if not os.path.exists(f'{analysis_path}\\{subject}\\{session}_WBPpeaksfreq.pdf'):
        #Plot with peaks
    fig4 = plot_peaksfreq(plethys_ddp_fi, Ipks, Ilocs, IFreqi, subject, session, stim_start, stim_stop)
    #save figure
    fig4.savefig(f'{analysis_path}\\{subject}\\{session}_WBPpeaksfreq.pdf')
    
