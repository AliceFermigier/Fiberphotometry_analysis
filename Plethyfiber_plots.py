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
from ast import literal_eval
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
from Fiberpho_loader import list_EVENT, list_TIMEWINDOW, PRE_EVENT_TIME, TIME_BEGIN, ORDER, THRESH_S, EVENT_TIME_THRESHOLD, CUT_FREQ

from Func_fiberplots import session_code, truncate, time_vector, timestamp_camera, timestamp_camera_fromraw, filter_dFF


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
    ax3.set_title(f'Whole body plethysmography and fiberphotometry - {mouse}')
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

def PETH(fiberpho_df, odor, sniffs_df, event, timewindow):
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
    print(f'PETH {odor}')
    #set time window relative to event
    PRE_TIME = timewindow[0]
    POST_TIME = timewindow[1]
    SAMPLERATE = round(len(fiberpho_df)/fiberpho_df['Time(s)'].values[-1])
    print(SAMPLERATE)
    
    #if behaviour not recognized
    if odor not in set(sniffs_df['Odor']) :
        print('Odor not recognized')
    
    list_t_event_o = []
    list_t_event_w = []
    for list_sniffs in sniffs_df.loc[(sniffs_df['Subject']==mouse) & (sniffs_df['Odor']==odor),['Start_Stop']].values:
        for [x_start, x_stop] in literal_eval(list_sniffs[0]): 
            if [x_start, x_stop] != [0,0]: #if length behav < threshold
                list_t_event_o.append(x_start)
                list_t_event_w.append(x_stop)

    #if event too short, don't process it
    for (start,end) in zip(list_t_event_o, list_t_event_w):
        if end-start<EVENT_TIME_THRESHOLD:
            list_t_event_o.remove(start)
            list_t_event_w.remove(end)
    
    #creates array of fiberpho data centered on event
    if event == 'onset':
        list_t_event = list_t_event_o
    elif event == 'withdrawal':
        list_t_event = list_t_event_w
        
    print(list_t_event)
    PETH_array = np.zeros((len(list_t_event),(POST_TIME+PRE_TIME)*SAMPLERATE+1))
    for (i, ind_event) in enumerate(list_t_event) :
        #calculates baseline F0 on time window before event (from PRE_TIME to PRE_EVENT_TIME)
        F0 = np.mean(fiberpho_df.loc[(ind_event-PRE_TIME)*SAMPLERATE:(ind_event-PRE_EVENT_TIME)*SAMPLERATE, 'Analog In. | Ch.1 470 nm (Deinterleaved)_dF/F0-Analog In. | Ch.1 405 nm (Deinterleaved)_dF/F0_LowPass'])
        std0 = np.std(fiberpho_df.loc[(ind_event-PRE_TIME)*SAMPLERATE:(ind_event-PRE_EVENT_TIME)*SAMPLERATE, 'Analog In. | Ch.1 470 nm (Deinterleaved)_dF/F0-Analog In. | Ch.1 405 nm (Deinterleaved)_dF/F0_LowPass'])
        #creates array of z-scored dFF : z = (dFF-meandFF_baseline)/stddFF_baseline
        PETH_array[i] = (fiberpho_df.loc[(ind_event-PRE_TIME)*SAMPLERATE:(ind_event+POST_TIME)*SAMPLERATE, 'Analog In. | Ch.1 470 nm (Deinterleaved)_dF/F0-Analog In. | Ch.1 405 nm (Deinterleaved)_dF/F0_LowPass']-F0)/std0
    
    return(PETH_array)

def plethyfiber_plot_sniffs(fiberpho_df, plethys_df, sniffs_df):
    """
    Take fiberpho data, plethysmo data and stims
    Plots data
    """
    fig5 = plt.figure(figsize=(20,7))
    ax8 = fig5.add_subplot(211)
    
    #plot plethysmo data
    p1, = ax8.plot('Time(s)', 'AIn-4', linewidth=.5, color='black', data=plethys_df.loc[plethys_df['Time(s)']>300])    
    ax8.set_ylabel('WBP (A.U.)')
    ax8.set_title(f'Whole body plethysmography and fiberphotometry - {mouse}')
    #ax8.legend(loc='lower left')
    ax8.margins(0,0.1)
    
    #plot denoised fiberphotometry data
    ax9 = fig5.add_subplot(212, sharex=ax8)
    p2, = ax9.plot('Time(s)', 'Analog In. | Ch.1 470 nm (Deinterleaved)_dF/F0-Analog In. | Ch.1 405 nm (Deinterleaved)_dF/F0_LowPass',
                   linewidth=.6, color='black', label='GCaMP-ISOS', data = fiberpho_df.loc[fiberpho_df['Time(s)']>300])
    ax9.set_ylabel(r'$\Delta$F/F')
    ax9.legend(handles=[p2], loc='upper right')
    #ax9.set_ylim(-1,2)
    ax9.margins(0,0.1)
    ax9.set_xlabel('Time(s)')

    #makes areas corresponding to behaviours
    i=0
    j=0
    k=0
    for list_sniffs in sniffs_df.loc[(sniffs_df['Subject']==mouse) & (sniffs_df['Odor']=='Clean'),['Start_Stop']].values:
        for [x_start, x_stop] in literal_eval(list_sniffs[0]):
            if [x_start, x_stop] != [0,0]:
                ax8.axvspan(x_start, x_stop, facecolor='grey', alpha=0.5, label = '_'*i + 'Sniff Clean')
                ax9.axvspan(x_start, x_stop, facecolor='grey', alpha=0.5)
                i+=1
    for list_sniffs in sniffs_df.loc[(sniffs_df['Subject']==mouse) & (sniffs_df['Odor']=='HC'),['Start_Stop']].values:
        for [x_start, x_stop] in literal_eval(list_sniffs[0]):
            if [x_start, x_stop] != [0,0]:
                ax8.axvspan(x_start, x_stop, facecolor='gold', alpha=0.5, label = '_'*j + f'Sniff HC')
                ax9.axvspan(x_start, x_stop, facecolor='gold', alpha=0.5)
                j+=1
    for list_sniffs in sniffs_df.loc[(sniffs_df['Subject']==mouse) & (sniffs_df['Odor']=='Novel'),['Start_Stop']].values:
        for [x_start, x_stop] in literal_eval(list_sniffs[0]):
            if [x_start, x_stop] != [0,0]:
                ax8.axvspan(x_start, x_stop, facecolor='purple', alpha=0.5, label = '_'*k + f'Sniff Novel')
                ax9.axvspan(x_start, x_stop, facecolor='purple', alpha=0.5)
                k+=1
        
    #makes vertical lines for stims
    for stim in sniffs_df.loc[(sniffs_df['Subject']==mouse),['Stim']].values:
        [x_start,x_stop] = literal_eval(stim[0])
        ax8.axvline(x_start, color='lightsteelblue', ls = '--')
        ax8.axvline(x_stop, color='slategrey', ls = '--')
        ax9.axvline(x_start, color='lightsteelblue', ls = '--')
        ax9.axvline(x_stop, color='slategrey', ls = '--')
    
    return(fig5)

def plot_PETH(PETH_data, odor, event, timewindow):
    """
    Plots PETH average and heatmap
    """
    PRE_TIME = timewindow[0]
    POST_TIME = timewindow[1]
    
    SAMPLERATE = round(len(fiberpho_df)/fiberpho_df['Time(s)'].values[-1])
    
    #create figure
    fig7 = plt.figure(figsize=(6,10))
    
    #create time vector
    peri_time = np.arange(-PRE_TIME, POST_TIME+1/SAMPLERATE, 1/SAMPLERATE)
    print(len(peri_time))
    
    #calculate mean dFF and std error
    mean_dFF_snips = np.mean(PETH_data, axis=0)
    std_dFF_snips = np.std(PETH_data, axis=0)
        
    #plot individual traces and mean
    ax5 = fig7.add_subplot(212)
    for snip in PETH_data:
        print(len(snip))
        p1, = ax5.plot(peri_time, snip, linewidth=.5,
                       color=[.7, .7, .7], label='Individual trials')
    p2, = ax5.plot(peri_time, mean_dFF_snips, linewidth=2,
                   color='green', label='Mean response')
    
    #plot standard error bars
    p3 = ax5.fill_between(peri_time, mean_dFF_snips+std_dFF_snips,
                      mean_dFF_snips-std_dFF_snips, facecolor='green', alpha=0.2)
    p4 = ax5.axvline(x=0, linewidth=2, color='slategray', ls = '--', label=f'Sniff {event}')
    
    #ax5.axis('tight')
    ax5.set_xlabel('Time (s)')
    ax5.set_ylabel('z-scored $\Delta$F/F')
    ax5.legend(handles=[p1, p2, p4], loc='upper left', fontsize = 'small')
    ax5.margins(0,0.01)
    
    #add heatmap
    ax6 = fig7.add_subplot(211)
    cs = ax6.imshow(PETH_data, cmap='magma', aspect = 'auto',
                    interpolation='none', extent=[-PRE_TIME, POST_TIME, len(PETH_data), 0],
                    vmin = -6, vmax = 9)
    ax6.set_ylabel('Bout Number')
    ax6.set_yticks(np.arange(.5, len(PETH_data), 2))
    ax6.set_yticklabels(np.arange(0, len(PETH_data), 2))
    ax6.axvline(x=0, linewidth=2, color='black', ls = '--')
    ax6.set_title(f'{odor} {event} - Plethysmo {mouse}')
    
    fig7.subplots_adjust(right=0.8, hspace = 0.1)
    cbar_ax = fig7.add_axes([0.85, 0.54, 0.02, 0.34])
    fig7.colorbar(cs, cax=cbar_ax)
    
    return(fig7)

def meanmax_dFF_sniffs(fiberpho_df, sniffs_df):
    """
    calculates mean and max dFF during sniffs
    """
    
    sniffs_mouse = sniffs_df[sniffs_df['Subject']==mouse]
    
    list_mean_sniffs=[]
    list_max_sniffs=[]
    list_odors = []
    list_counts = []
    list_mouse = []
    list_group = []
    for odor in set(sniffs_mouse['Odor']):
        for count in sniffs_mouse.loc[sniffs_mouse['Odor']==odor,'Count']:
            for list_sniffs in sniffs_mouse.loc[(sniffs_mouse['Count']==count)|(sniffs_mouse['Odor']==odor),['Start_Stop']].values:
                for [x_start, x_stop] in literal_eval(list_sniffs[0]):
                    if [x_start, x_stop] != [0,0]:
                        list_mean_sniffs.append(fiberpho_df.loc[(fiberpho_df['Time(s)']<x_stop)&(fiberpho_df['Time(s)']>x_start),
                                                               'Analog In. | Ch.1 470 nm (Deinterleaved)_dF/F0-Analog In. | Ch.1 405 nm (Deinterleaved)_dF/F0_LowPass'].mean())
                        list_max_sniffs.append(fiberpho_df.loc[(fiberpho_df['Time(s)']<x_stop)&(fiberpho_df['Time(s)']>x_start), 
                                                              'Analog In. | Ch.1 470 nm (Deinterleaved)_dF/F0-Analog In. | Ch.1 405 nm (Deinterleaved)_dF/F0_LowPass'].max())
                        list_odors.append(odor)
                        list_counts.append(count)
                        list_mouse.append(mouse)
                        list_group.append(sniffs_mouse['Group'].values[0])
                        
    list_data = [list_mouse, list_group, list_odors, list_counts, list_mean_sniffs, list_max_sniffs]
    list_columns = ['Subject','Group','Odor','Count', 'Mean dFF sniff', 'Max dFF sniff']
    meandFFsniffs_df = pd.DataFrame(data=[list_data], columns=list_columns)
                        
    return(meandFFsniffs_df)
                        
def meanmax_dFF_stims(fiberpho_df, sniffs_df):
    """
    calculates mean and max dFF during stims
    """
    
    sniffs_mouse = sniffs_df[sniffs_df['Subject']==mouse]
    
    list_mean_stim=[]
    list_max_stim=[]
    list_odors = []
    list_counts = []
    list_mouse = []
    list_group = []
    for odor in set(sniffs_mouse['Odor']):
        for count in sniffs_mouse.loc[sniffs_mouse['Odor']==odor,'Count']:
            [x_start, x_stop]=literal_eval(sniffs_mouse.loc[(sniffs_mouse['Count']==count)|(sniffs_mouse['Odor']==odor),['Stim']].values[0][0])
            list_mean_stim.append(fiberpho_df.loc[(fiberpho_df['Time(s)']<x_stop)&(fiberpho_df['Time(s)']>x_start),
                                                   'Analog In. | Ch.1 470 nm (Deinterleaved)_dF/F0-Analog In. | Ch.1 405 nm (Deinterleaved)_dF/F0_LowPass'].mean())
            list_max_stim.append(fiberpho_df.loc[(fiberpho_df['Time(s)']<x_stop)&(fiberpho_df['Time(s)']>x_start), 
                                                  'Analog In. | Ch.1 470 nm (Deinterleaved)_dF/F0-Analog In. | Ch.1 405 nm (Deinterleaved)_dF/F0_LowPass'].max())
            list_odors.append(odor)
            list_counts.append(count)
            list_mouse.append(mouse)
            list_group.append(sniffs_mouse['Group'].values[0])
            
    list_data = [list_mouse, list_group, list_odors, list_counts, list_mean_stim, list_max_stim]
    list_columns = ['Subject','Group','Odor','Count', 'Mean dFF stim', 'Max dFF stim']
    meandFFstim_df = pd.DataFrame(data=[list_data], columns=list_columns)
                        
    return(meandFFstim_df)


#%%
########
#SCRIPT#
########

exp_path = analysis_path / 'Plethysmo'
exp = str(exp_path).split('\\')[-1]
data_path_exp = data_path / proto_df.loc[proto_df['Task']==exp, 'Data_path'].values[0]
sniffs_df = pd.read_excel(data_path_exp / 'Sniffs.xlsx') #converters={'Stim': pd.eval, 'Start_Stop': pd.eval})

#create repository for values of thresholds : length and interbout
repo_path = exp_path / f'length{EVENT_TIME_THRESHOLD/10}_interbout{THRESH_S}_o{ORDER}f{CUT_FREQ}'
if not os.path.exists(repo_path):
    os.mkdir(repo_path)
    for subject in subjects_df['Subject']:
        os.mkdir(repo_path / subject)
        
mean_dFFsniff_list = []    
mean_dFFstim_list = []    
for mouse_path in Path(repo_path).iterdir():
    # '/' on mac, '\\' on windows
    mouse = str(mouse_path).split('\\')[-1]
    print(mouse)
    if mouse in set(sniffs_df['Subject']) and mouse != 'A3f':
        #get data
        rawdata_path = data_path_exp / f'{mouse}_1.csv'
        fiberpho_path = data_path_exp / f'{mouse}_1_dFFfilt.csv'
        if mouse == 'A3f':
            plethys_df = pd.read_csv(rawdata_path, usecols=['Time(s)','AIn-4'])
        else :
            plethys_df = pd.read_csv(rawdata_path, skiprows=1, usecols=['Time(s)','AIn-4'])
        fiberpho_df = pd.read_csv(fiberpho_path)
        if len(fiberpho_df.columns) == 5:
            fiberpho_df.drop(columns='Unnamed: 4', inplace = True)
            fiberpho_df.interpolate(methode = 'nearest', inplace = True)
        #plot figure
        if not (mouse_path / f'{mouse}_WBPfiberpho_raw.pdf').is_file():
            fig4 = plethyfiber_plot_raw(fiberpho_df, plethys_df)
            fig4.savefig(mouse_path / f'{mouse}_WBPfiberpho_raw.png') 
            fig4.savefig(mouse_path / f'{mouse}_WBPfiberpho_raw.pdf')
            plt.close('fig4') 
        if not (mouse_path / f'{mouse}_WBPfiberpho_sniffs.pdf').is_file():
            fig5 = plethyfiber_plot_sniffs(fiberpho_df, plethys_df, sniffs_df)
            fig5.savefig(mouse_path / f'{mouse}_WBPfiberpho_sniffs.png') 
            fig5.savefig(mouse_path / f'{mouse}_WBPfiberpho_sniffs.pdf')
        for odor in set(sniffs_df['Odor']):
            print(odor)
            for (event, timewindow) in zip(list_EVENT, list_TIMEWINDOW):
                print(event)
                if not (mouse_path / f'{mouse}{odor}_PETH{event[0]}.pdf').is_file():
                    PETH_data = PETH(fiberpho_df, odor, sniffs_df, event, timewindow)
                    fig7 = plot_PETH(PETH_data, odor, event, timewindow)
                    fig7.savefig(mouse_path / f'{mouse}{odor}_PETH{event[0]}.png')
                    fig7.savefig(mouse_path / f'{mouse}{odor}_PETH{event[0]}.pdf')
                    plt.close('fig7') 
        
        mean_dFFsniffs_df = meanmax_dFF_sniffs(fiberpho_df, sniffs_df)
        mean_dFFsniff_list.append(mean_dFFsniffs_df)
       
        mean_dFFstim_df = meanmax_dFF_stims(fiberpho_df, sniffs_df)
        mean_dFFstim_list.append(mean_dFFstim_df)
        
meandFFsniff_allmice = pd.concat(mean_dFFsniff_list)
meandFFsniff_allmice.to_excel(repo_path / f'{exp}_meandFFsniffs.xlsx')
meandFFstim_allmice = pd.concat(mean_dFFstim_list)
meandFFstim_allmice.to_excel(repo_path / f'{exp}_meandFFstims.xlsx')
                    
                

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
    
