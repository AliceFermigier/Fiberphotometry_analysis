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
import matplotlib.pyplot as plt
import pandas as pd
from ast import literal_eval
# citing pyabf : 
# Harden, SW (2020). pyABF 2.2.3. [Online]. Available: https://pypi.org/project/pyabf/

#%%
###################
#DEFINED FUNCTIONS#
###################

def align_sniffs(fiberpho_df, plethys_df, sniffs_df, sr, mouse):
    """
    Aligns fiberpho data with plethysmo and sniffs
    """
    sniffmouse_df = sniffs_df.loc[(sniffs_df['Subject']==mouse)]
    fibersniff_df = pd.DataFrame({'Time(s)':fiberpho_df['Time(s)'], 'Denoised dFF':fiberpho_df['Denoised dFF']})
    plethy_list = []
    for time in fiberpho_df['Time(s)']:
        plethy_list.append(plethys_df.loc[plethys_df['Time(s)']==time,['AIn-4']])
    fibersniff_df.insert(len(fibersniff_df.columns),'Plethysmograph',plethy_list,allow_duplicates = False)
    for odor in set(sniffs_df['Odor']):
        for (count,stim) in enumerate(sniffmouse_df.loc[sniffmouse_df['Odor']==odor,'Stim'].values):
            fibersniff_df.insert(len(fibersniff_df.columns),f'Stim {odor} {count}',0,allow_duplicates = False)
            [x_start,x_stop] = literal_eval(stim[0])
            fibersniff_df.loc[round(x_start*sr):round(x_stop*sr), f'Stim {odor} {count}'] = 1
        for (count,list_sniffs) in enumerate(sniffmouse_df.loc[sniffmouse_df['Odor']==odor,'Start_Stop'].values):
            fibersniff_df.insert(len(fibersniff_df.columns),f'Sniff {odor} {count}',0,allow_duplicates = False)
            for [x_start, x_stop] in literal_eval(list_sniffs[0]):
                if [x_start, x_stop] != [0,0]:
                    fibersniff_df.loc[round(x_start*sr):round(x_stop*sr), f'Sniff {odor} {count}'] = 1
    
    return fibersniff_df

def process_fibersniff(fibersniff_df, EVENT_TIME_THRESHOLD, THRESH_S, sr):
    """
    Fuse sniffs that are too close and remove sniffs that are too short
    """
    sr = round(sr)
    
    #1 fuse sniffs that are too close
    for col in fibersniff_df.columns[3:]:
        if col.split()[0] == 'Sniff':
            i_1 = 0
            count = 0
            for (ind,i) in zip(fibersniff_df.index, fibersniff_df[col]):
                if i==i_1 and i==0:
                    count += 1
                elif i!=i_1 and i==0:
                    count = 1
                elif i!=i_1 and i==1:
                    if count <= THRESH_S*sr:
                        fibersniff_df.loc[ind-count:ind-1, col] = [1]*count
                i_1 = i
            
    #2 remove explorations that are too short
    for col in fibersniff_df.columns[2:]:
        if col.split()[0] == 'Sniff':
            i_1 = 0
            count = 0
            for (ind,i) in zip(fibersniff_df.index, fibersniff_df[col]):
                if i==i_1 and i==1:
                    count += 1
                elif i!=i_1 and i==1:
                    count = 1
                elif i!=i_1 and i==0:
                    if count <= EVENT_TIME_THRESHOLD*sr:
                        fibersniff_df.loc[ind-count:ind-1, col] = [0]*count
                i_1 = i
    
    return fibersniff_df

def derive(fibersniff_df):
    """
    calculate the derivative of behav of interest and put in new df
    that way, it will show 1 when behaviour starts and -1 when it stops
    """
    for col in fibersniff_df.columns[3:]:
        fibersniff_df[col] = fibersniff_df[col].diff()

    return fibersniff_df
    
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

def plethyfiber_plot_sniffs(dfibersniff_df,sniffs_df,mouse):
    """
    Take fiberpho data, aligned with plethysmo data and stims
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
    i = 0
    j = 0
    m = 0
    for count in set(sniffs_df.loc[sniffs_df['Odor']=='Clean','Count']):
        for (x,y) in zip(dfibersniff_df['Time(s)'].tolist(), dfibersniff_df['Sniff Clean {count}'].tolist()):
            if y == 1:
                x_start = x
            if y == -1 and x_start!=0:
                ax9.axvspan(x_start, x, facecolor='grey', alpha=.5, label = '_'*i + 'Sniff Clean')
                x_start=0
                i+=1
    for count in set(sniffs_df.loc[sniffs_df['Odor']=='HC','Count']):
        for (x,y) in zip(dfibersniff_df['Time(s)'].tolist(), dfibersniff_df['Sniff Clean {count}'].tolist()):
            if y == 1:
                x_start = x
            if y == -1 and x_start!=0:
                ax9.axvspan(x_start, x, facecolor='gold', alpha=.5, label = '_'*j + 'Sniff HC')
                x_start=0
                j+=1
    for count in set(sniffs_df.loc[sniffs_df['Odor']=='HC','Count']):
        for (x,y) in zip(dfibersniff_df['Time(s)'].tolist(), dfibersniff_df['Sniff Clean {count}'].tolist()):
            if y == 1:
                x_start = x
            if y == -1 and x_start!=0:
                ax9.axvspan(x_start, x, facecolor='gold', alpha=.5, label = '_'*j + 'Sniff HC')
                x_start=0
                j+=1


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
    
    return PETH_array

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
    
    return PETH_array

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
    ax5.set_xlabel('Time(s)')
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
    
    return fig7

def plot_PETH_pooled(included_groups, colorscheme, PETHarray_list, BOI, event, timewindow, exp, session):
    """
    Plots PETH averaged over 2 groups

    included_groups : list of included groups (['CD','HFD'])
    PETHarray_list : list of PETH arrays
    BOI : behaviour of interest
    event : onset or withdrawal
    timewindow : time before and after behaviour
    """
    
    PRE_TIME = timewindow[0]
    POST_TIME = timewindow[1]
    
    #create figure
    fig4 = plt.figure(figsize=(6,4))
    ax5 = fig4.add_subplot(111)
    
    #create time vector
    peri_time = np.arange(-PRE_TIME, POST_TIME+0.1, 0.1)       
    
    listmean_dFF_snips = []
    listsem_dFF_snips = []
    #calculate mean dFF and std error
    for (i,PETH_data) in enumerate(PETHarray_list):
        listmean_dFF_snips.append(np.mean(PETH_data, axis=0))
        listsem_dFF_snips.append(np.std(PETH_data, axis=0))
        
    #plot individual traces and mean CD
    # for snip in PETHarray_list[0]:
    #     p1, = ax5.plot(peri_time, snip, linewidth=.5,
    #                    color='cornflowerblue', alpha=.3)
    p2, = ax5.plot(peri_time, listmean_dFF_snips[0], linewidth=2,
                   color=colorscheme[0], label=included_groups[0])   
    #plot standard error bars
    p3 = ax5.fill_between(peri_time, listmean_dFF_snips[0]+(listsem_dFF_snips[0]/np.sqrt(len(PETHarray_list[0]))),
                      listmean_dFF_snips[0]-(listsem_dFF_snips[0]/np.sqrt(len(PETHarray_list[0]))), facecolor=colorscheme[0], alpha=0.2)
    
    #plot individual traces and mean HFD
    # for snip in PETHarray_list[1]:
    #     p4, = ax5.plot(peri_time, snip, linewidth=.5,
    #                    color='coral', alpha=.3)
    p5, = ax5.plot(peri_time, listmean_dFF_snips[1], linewidth=2,
                   color=colorscheme[1], label=included_groups[1])   
    #plot standard error bars
    p6 = ax5.fill_between(peri_time, listmean_dFF_snips[1]+(listsem_dFF_snips[1]/np.sqrt(len(PETHarray_list[1]))),
                      listmean_dFF_snips[1]-(listsem_dFF_snips[1]/np.sqrt(len(PETHarray_list[1]))), colorscheme[1], alpha=0.2)
    
    p8 = ax5.axvline(x=0, linewidth=2, color='slategray', ls = '--', label=BOI)
    
    ax5.set_xlabel('Time(s)')
    ax5.set_ylabel('z-scored $\Delta$F/F')
    ax5.legend(handles=[p2,p5,p8], loc='upper left', fontsize = 'small')
    ax5.margins(0, 0.1)
    ax5.set_title(f'{BOI} - {exp} {session} {included_groups[0]} {included_groups[1]}')
    
    return fig4