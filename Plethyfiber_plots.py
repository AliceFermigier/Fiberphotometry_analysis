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
import matplotlib.pyplot as plt
import pandas as pd
# citing pyabf : 
# Harden, SW (2020). pyABF 2.2.3. [Online]. Available: https://pypi.org/project/pyabf/

#put path to directory where python files are stored
if 'D:\\Profil\\Documents\\GitHub\\Fiberphotometry_analysis' not in sys.path:
    sys.path.append('D:\\Profil\\Documents\\GitHub\\Fiberphotometry_analysis')
if 'D:\\Profil\\Documents\\GitHub\\Plethysmograph_analysis' not in sys.path:
    sys.path.append('D:\\Profil\\Documents\\GitHub\\Plethysmograph_analysis')

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
        if point > signal_median + 5*signal_std or point < signal_median - 5*signal_std:
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


# def stim_list(plethy_df):
#     """
#     Find stims on traces and give timepoints (truncated to 0,1ms to match fiber data)
#     """
#     for stim_start in StimStart_list:
        

def plethyfiber_plot_raw(fiberpho_df, plethys_df):
    """
    Take fiberpho data, plethysmo data and stims
    Plots data
    """
    fig4 = plt.figure(figsize=(20,9))
    ax3 = fig4.add_subplot(311)
    
    #plot plethysmo data
    p1, = ax3.plot('Time(s)', 'AIn-4', linewidth=.5, color='black', data=plethys_df)    
    ax3.set_ylabel('WBP (A.U.)')
    ax3.set_title(f'Whole body plethysmography and fiberphotometry - {session} {mouse}')
    #ax3.legend(loc='lower left')
    ax3.margins(0,0.1)
    
    #plot denoised fiberphotometry data
    ax4 = fig4.add_subplot(312, sharex=ax3)
    p2, = ax4.plot('Time(s)', 'Analog In. | Ch.1 470 nm (Deinterleaved)_dF/F0-Analog In. | Ch.1 405 nm (Deinterleaved)_dF/F0_LowPass',
                   linewidth=.6, color='black', label='GCaMP-ISOS', data = fiberpho_df)
    ax4.set_ylabel(r'$\Delta$F/F')
    ax4.legend(handles=[p2], loc='upper right')
    #ax4.set_ylim(-1,2)
    ax4.margins(0,0.1)
    
    #plot fiberphotometry data GCaMP and Isos
    ax7 = fig4.add_subplot(313, sharex=ax3)
    p3, = ax7.plot('Time(s)', 'Analog In. | Ch.1 470 nm (Deinterleaved)_dF/F0_LowPass', 
                   linewidth=.5, color='deepskyblue', label='GCaMP', data = fiberpho_df) 
    p4, = ax7.plot('Time(s)', 'Analog In. | Ch.1 405 nm (Deinterleaved)_dF/F0_LowPass', 
                   linewidth=.5, color='blueviolet', label='ISOS', data = fiberpho_df)
    
    ax7.set_ylabel(r'$\Delta$F/F')
    ax7.legend(handles=[p3,p4], loc='upper right')
    ax7.margins(0,0.1)
    #ax7.set_ylim(-1,2)
    ax7.set_xlabel('Time(s)')
    
    return(fig4)


#%%
########
#SCRIPT#
########

# for exp_path in Path(analysis_path).iterdir():
#     if exp_path.is_dir():
#         exp = str(exp_path).split('\\')[-1]
#         for session_path in Path(exp_path).iterdir():
#             if session_path.is_dir():
#                 session = str(session_path).split('\\')[-1]
#                 #get data path related to the task in protocol excel file
#                 data_path_exp = data_path / proto_df.loc[proto_df['Task']==exp, 'Data_path'].values[0]

exp_path = analysis_path / 'Plethysmo'
session_path = exp_path / 'Novel'
session = str(session_path).split('\\')[-1]
data_path_exp = data_path / '20211022_AliceF_CA2b2plethysmoNovel'

for mouse_path in Path(session_path).iterdir():
    # '/' on mac, '\\' on windows
    mouse = str(mouse_path).split('\\')[-1]
    print(mouse)
    
    #get data
    rawdata_path = data_path_exp / f'{mouse}_0.csv'
    fiberpho_path = data_path_exp / f'{mouse}_0_dFFfilt.csv'
    plethys_df = pd.read_csv(rawdata_path, skiprows=1, usecols=['Time(s)','AIn-4'])
    fiberpho_df = pd.read_csv(fiberpho_path)
    if len(fiberpho_df.columns) == 5:
        fiberpho_df.drop(columns='Unnamed: 4', inplace = True)
        fiberpho_df.interpolate(methode = 'nearest', inplace = True)
    
    #plot figure
    #if not (mouse_path / f'{mouse}_WBPfiberpho_raw.pdf').is_file():
    fig4 = plethyfiber_plot_raw(fiberpho_df, plethys_df)
    fig4.savefig(mouse_path / f'{mouse}_WBPfiberpho_raw.png') 
    fig4.savefig(mouse_path / f'{mouse}_WBPfiberpho_raw.pdf')

# #remove artifacts
# plethys_df_clean = r_highartifacts_fiberplethy(plethys_df)
# #filter and crop
# (plethys_ddp_fi, plethys_ddp_f) = r_slowvar_crop(plethys_df_clean)
# plt.plot(plethys_ddp_fi)
# plt.plot(plethys_ddp_f)
# #find start and stop of stims
# (stim_start, stim_stop, stim_dur) = find_stim(plethys_ddp_f)
# #find inspirations
# (Ilocs, Ipks) = find_inspirations(plethys_ddp_fi)
# #calculate frequency, amplitude and area under curve
# (IFreqi, IFreqim_preStim, IFreqim_Stim, IFreqm_preStim, IFreqm_Stim) = resp_freq(plethys_ddp_fi, stim_start, stim_stop, Ilocs)
# (Ipks_preStim, mAmp_preStim, IArea_preStim, Ipks_Stim, mAmp_Stim, IArea_Stim) = resp_amp_AUC(plethys_ddp_fi, stim_start, stim_stop, Ilocs, Ipks)
# #Plot data with computed stim period
# if not os.path.exists(f'{analysis_path}\\{subject}\\{session}_WBPstim.pdf'):
#     fig3 = plot_procdatastim(plethys_ddp, stim_start, stim_stop, subject, session)
#     #save figure
#     fig3.savefig(f'{analysis_path}\\{subject}\\{session}_WBPstim.pdf')
# if not os.path.exists(f'{analysis_path}\\{subject}\\{session}_WBPpeaksfreq.pdf'):
#         #Plot with peaks
#     fig4 = plot_peaksfreq(plethys_ddp_fi, Ipks, Ilocs, IFreqi, subject, session, stim_start, stim_stop)
#     #save figure
#     fig4.savefig(f'{analysis_path}\\{subject}\\{session}_WBPpeaksfreq.pdf')
    
