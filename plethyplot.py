# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 11:30:12 2023

Functions for plotting with plethysmography data

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
from scipy import signal
import os
import sys
import matplotlib.pyplot as plt
import pandas as pd
from ast import literal_eval
# citing pyabf : 
# Harden, SW (2020). pyABF 2.2.3. [Online]. Available: https://pypi.org/project/pyabf/

#%%
###################
#DEFINED FUNCTIONS#
###################

def plethyfiber_plot_raw(fiberpho_df, plethys_df, mouse):
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
    p2, = ax4.plot('Time(s)', 'Denoised dFF',
                   linewidth=.6, color='black', label='GCaMP-ISOS', data = fiberpho_df)
    ax4.set_ylabel(r'$\Delta$F/F')
    ax4.legend(handles=[p2], loc='upper right')
    #ax4.set_ylim(-1,2)
    ax4.margins(0,0.1)
    
    #plot fiberphotometry data GCaMP and Isos
    ax7 = fig4.add_subplot(313, sharex=ax3)
    p3, = ax7.plot('Time(s)', '470 dFF', 
                   linewidth=.5, color='deepskyblue', label='GCaMP', data = fiberpho_df) 
    p4, = ax7.plot('Time(s)', '405 dFF', 
                   linewidth=.5, color='blueviolet', label='ISOS', data = fiberpho_df)
    
    ax7.set_ylabel(r'$\Delta$F/F')
    ax7.legend(handles=[p3,p4], loc='upper right')
    ax7.margins(0,0.1)
    #ax7.set_ylim(-1,2)
    ax7.set_xlabel('Time(s)')
    
    return fig4

def plethyfiber_plot_sniffs(fiberpho_df, plethys_df, sniffs_df,mouse):
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
    p2, = ax9.plot('Time(s)', 'Denoised dFF',
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
                ax8.axvspan(x_start, x_stop, facecolor='gold', alpha=0.5, label = '_'*j + 'Sniff HC')
                ax9.axvspan(x_start, x_stop, facecolor='gold', alpha=0.5)
                j+=1
    for list_sniffs in sniffs_df.loc[(sniffs_df['Subject']==mouse) & (sniffs_df['Odor']=='Novel'),['Start_Stop']].values:
        for [x_start, x_stop] in literal_eval(list_sniffs[0]):
            if [x_start, x_stop] != [0,0]:
                ax8.axvspan(x_start, x_stop, facecolor='purple', alpha=0.5, label = '_'*k + 'Sniff Novel')
                ax9.axvspan(x_start, x_stop, facecolor='purple', alpha=0.5)
                k+=1
        
    #makes vertical lines for stims
    for stim in sniffs_df.loc[(sniffs_df['Subject']==mouse),['Stim']].values:
        [x_start,x_stop] = literal_eval(stim[0])
        ax8.axvline(x_start, color='lightsteelblue', ls = '--')
        ax8.axvline(x_stop, color='slategrey', ls = '--')
        ax9.axvline(x_start, color='lightsteelblue', ls = '--')
        ax9.axvline(x_stop, color='slategrey', ls = '--')
    
    return fig5

def PETH_sniff(fiberpho_df, odor, sniffs_df, event, timewindow, mouse, sr, PRE_EVENT_TIME):
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
    #round samplerate because has to be an int
    sr = round(sr)
    
    #set time window relative to event
    PRE_TIME = timewindow[0]
    POST_TIME = timewindow[1]
    
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
    
    #creates array of fiberpho data centered on event
    if event == 'onset':
        list_t_event = list_t_event_o
    elif event == 'withdrawal':
        list_t_event = list_t_event_w
        
        #take the beginning of the behaviour as the beginning of the slope : min(dFF) on 3 secs before entry
    if event == 'onset':
        for i,ind_onset in enumerate(list_t_event):
            list_t_event[i] = fiberpho_df.loc[ind_onset-3*sr:ind_onset, 'Denoised dFF'].idxmin()
        
    print(list_t_event)
    PETH_array = np.zeros((len(list_t_event),(POST_TIME+PRE_TIME)*sr+1))
    for (i, ind_event) in enumerate(list_t_event) :
        #calculates baseline F0 on time window before event (from PRE_TIME to PRE_EVENT_TIME)
        F0 = np.mean(fiberpho_df.loc[(ind_event-PRE_TIME)*sr:(ind_event-PRE_EVENT_TIME)*sr, 'Denoised dFF'])
        std0 = np.std(fiberpho_df.loc[(ind_event-PRE_TIME)*sr:(ind_event-PRE_EVENT_TIME)*sr, 'Denoised dFF'])
        #creates array of z-scored dFF : z = (dFF-meandFF_baseline)/stddFF_baseline
        PETH_array[i] = (fiberpho_df.loc[(ind_event-PRE_TIME)*sr:(ind_event+POST_TIME)*sr, 'Denoised dFF']-F0)/std0
    
    return(PETH_array)

def PETH_stim(fiberpho_df, odor, sniffs_df, event, timewindow, mouse, sr, PRE_EVENT_TIME):
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
    #round samplerate because has to be an int
    sr = round(sr)
    
    #set time window relative to event
    PRE_TIME = timewindow[0]
    POST_TIME = timewindow[1]
    
    #if behaviour not recognized
    if odor not in set(sniffs_df['Odor']) :
        print('Odor not recognized')
    
    list_t_event_o = []
    list_t_event_w = []
    for [x_start, x_stop] in literal_eval(sniffs_df.loc[(sniffs_df['Subject']==mouse) & (sniffs_df['Odor']==odor),['Stim']].values[0]):
        list_t_event_o.append(x_start)
        list_t_event_w.append(x_stop)
    
    #creates array of fiberpho data centered on event
    if event == 'onset':
        list_t_event = list_t_event_o
    elif event == 'withdrawal':
        list_t_event = list_t_event_w
        
    PETH_array = np.zeros((len(list_t_event),(POST_TIME+PRE_TIME)*sr+1))
    for (i, ind_event) in enumerate(list_t_event) :
        #calculates baseline F0 on time window before event (from PRE_TIME to PRE_EVENT_TIME)
        F0 = np.mean(fiberpho_df.loc[(ind_event-PRE_TIME)*sr:(ind_event-PRE_EVENT_TIME)*sr, 'Denoised dFF'])
        std0 = np.std(fiberpho_df.loc[(ind_event-PRE_TIME)*sr:(ind_event-PRE_EVENT_TIME)*sr, 'Denoised dFF'])
        #creates array of z-scored dFF : z = (dFF-meandFF_baseline)/stddFF_baseline
        PETH_array[i] = (fiberpho_df.loc[(ind_event-PRE_TIME)*sr:(ind_event+POST_TIME)*sr, 'Denoised dFF']-F0)/std0
    
    return(PETH_array)

def plot_PETH(PETH_data, odor, event, timewindow, BOI, sr, mouse):
    """
    Plots PETH average and heatmap
    BOI = 'Sniff' or 'Stim'
    """
    PRE_TIME = timewindow[0]
    POST_TIME = timewindow[1]
    
    sr = round(sr)
    
    #create figure
    fig7 = plt.figure(figsize=(6,10))
    
    #create time vector
    peri_time = np.arange(-PRE_TIME, POST_TIME+1/sr, 1/sr)
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
    p4 = ax5.axvline(x=0, linewidth=2, color='slategray', ls = '--', label=f'{BOI} {event}')
    
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