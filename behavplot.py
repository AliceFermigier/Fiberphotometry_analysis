# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 16:44:33 2023

Functions for plotting with behavioural data

@author: alice fermigier
"""

#%%
##########
#IMPORTED#
##########

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#%%
###################
#DEFINED FUNCTIONS#
###################

def align_behav(behav10Sps, fiberpho, timevector, timestart_camera, exp):
    """
    Aligns fiberpho data with behavioural data from Boris on timevector
    --> Parameters 
        behav10Sps : pd df, binary behavioural data from Boris, with sample rate = 10Sps
        fiberpho : pd df, pre-processed fiber photometry data (deltaF/F 470nm - 405nm)
            to be aligned with behav data
        timevector : pd df, timevector to plot data
        timestart_camera : int, index of when camera starts
    --> Returns
        fiberbehav_df : pd df of aligned data

    """
    #scored behaviours in Boris
    list_behav = behav10Sps.columns[1:].tolist()
    behav_comp = [0]*len(list_behav)
    
    #index where camera starts
    indstart = np.where(round(timevector,1) == timestart_camera)[0].tolist()[0]
    print(indstart)

    # create lists of behaviour data for each scored behaviour
    # aligned with start of the camera
    for (ind, behav) in enumerate(list_behav):
        behav_comp[ind] = [0]*indstart
        print(len(behav_comp[ind]))
        behav_comp[ind].extend(behav10Sps[behav].tolist())
        print(len(behav10Sps[behav]),len(behav_comp[ind]))
       
    # creates list of denoised fiberpho data    
    denoised_fiberpho = fiberpho['Denoised dFF'].dropna()
    denoised_fiberpho_list = denoised_fiberpho.tolist()
    
    # creates list of isosbestic (405) and Ca dependent (470) data
    dff_405nm = fiberpho['405 dFF'].dropna()
    dff_405nm_list = dff_405nm.tolist()
    dff_470nm = fiberpho['470 dFF'].dropna()
    dff_470nm_list = dff_470nm.tolist()
    
    # makes timevector into a list
    timelist = timevector.tolist()
    
    # crops lists so that all lengths match
    min_length = min([len(timelist), len(denoised_fiberpho_list), len(behav_comp[0]),
                      len(dff_405nm_list), len(dff_470nm_list)])
    print(min_length)
    timelist = timelist[:min_length]
    denoised_fiberpho_list = denoised_fiberpho_list[:min_length]
    dff_405nm_list = dff_405nm_list[:min_length]
    dff_470nm_list = dff_470nm_list[:min_length]
    behav_crop = []
    for behav in behav_comp:
        behav_crop.append(behav[:min_length])
        
    print(len(behav_crop[0])) 
    fiberbehav_df = pd.DataFrame(data = {'Time(s)':timelist, 'Denoised dFF' : denoised_fiberpho_list,
                                         '405nm dFF' : dff_405nm_list, '470nm dFF' : dff_470nm_list})
        
    for (behav,data) in zip(list_behav,behav_crop):
        fiberbehav_df.insert(len(fiberbehav_df.columns),behav,data,allow_duplicates=False)
        
    if exp=='NewContext':
        fiberbehav_df.loc[:405,'Homecage']=0
    
    return(fiberbehav_df)

def behav_process(fiberbehav_df, list_BOI, THRESH_S, EVENT_TIME_THRESHOLD, sr):
    
    sr = round(sr)

    #1 fuse explorations that are too close
    for BOI in list_BOI:
        i_1 = 0
        count = 0
        for (ind,i) in zip(fiberbehav_df.index, fiberbehav_df[BOI]):
            if i==i_1 and i==0:
                count += 1
            elif i!=i_1 and i==0:
                count = 1
            elif i!=i_1 and i==1:
                if count <= THRESH_S*sr and count!=1:
                    fiberbehav_df.loc[ind-count:ind-1, BOI] = [1]*count
            i_1 = i
            
    #2 remove explorations that are too short
    for BOI in list_BOI:
        i_1 = 0
        count = 0
        for (ind,i) in zip(fiberbehav_df.index, fiberbehav_df[BOI]):
            if i==i_1 and i==1:
                count += 1
            elif i!=i_1 and i==1:
                count = 1
            elif i!=i_1 and i==0:
                if count <= EVENT_TIME_THRESHOLD*sr and count!=1:
                    fiberbehav_df.loc[ind-count:ind-1, BOI] = [0]*count
            i_1 = i
            
    return fiberbehav_df 
                    
def derive(fiberbehav_df):
    """
    calculate the derivative of behav of interest and put in new df
    that way, it will show 1 when behaviour starts and -1 when it stops
    """
    for col in fiberbehav_df.columns[5:]:
        fiberbehav_df[col] = fiberbehav_df[col].diff()

    return fiberbehav_df 

def plot_fiberpho_behav(behavprocess_df,list_BOI,exp,session,mouse,THRESH_S,EVENT_TIME_THRESHOLD):
    """
    Plots denoised deltaF/F aligned with behaviour (includes baseline)
    """
    
    behavprocesssnip_df = behavprocess_df[behavprocess_df['Time(s)'] > 40]
    
    fig2 = plt.figure(figsize=(20,5))
    ax1 = fig2.add_subplot(111)
    
    if exp=='Fear':
        if session =='Test':
            #plots fiberpho trace and behaviourb 
            p1, = ax1.plot('Time(s)', 'Denoised dFF', linewidth=.6, color='black', label='_GCaMP', data = behavprocesssnip_df)
    
            #makes areas corresponding to behaviours
            i = 0
            j = 0
            for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Homecage'].tolist()):
                if y == 1:
                    x_start = x
                if y == -1 and x_start!=0:
                    ax1.axvspan(x_start, x, facecolor='gold', alpha=0.3, label = '_'*i + 'Homecage')
                    x_start=0
                    i+=1
            for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Fear cage'].tolist()):
                if y == 1:
                    x_start = x
                if y == -1 and x_start!=0:
                    ax1.axvspan(x_start, x, facecolor='blue', alpha=0.3, label = '_'*j + 'Fear cage')
                    x_start=0
                    j+=1
            
        elif session =='Conditioning':
            #plots fiberpho trace and behaviourb 
            p1, = ax1.plot('Time(s)', 'Denoised dFF', linewidth=.6, color='black', label='_GCaMP', data = behavprocesssnip_df)
            x_entry = fiberbehav_df.at[int(np.where(behavprocesssnip_df['Shock'] == 1)[0][0]), 'Time(s)']
            ax1.axvline(x_entry, color='yellow', ls = '-', lw=2, label = 'Shock' )
            
    elif exp=='NewContext':
        #plots fiberpho trace and behaviourb 
        p1, = ax1.plot('Time(s)', 'Denoised dFF', linewidth=.6, color='black', label='_GCaMP', data = behavprocesssnip_df)

        #makes areas corresponding to behaviours
        i = 0
        j = 0
        for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Homecage'].tolist()):
            if y == 1:
                x_start = x
            if y == -1 and x_start!=0:
                ax1.axvspan(x_start, x, facecolor='gold', alpha=0.3, label = '_'*i + 'Homecage')
                x_start=0
                i+=1
        for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['New context'].tolist()):
            if y == 1:
                x_start = x
            if y == -1 and x_start!=0:
                ax1.axvspan(x_start, x, facecolor='darkturquoise', alpha=0.3, label = '_'*j + 'New context')
                x_start=0
                j+=1

    else:
        if 'Essai' in exp:
            p1, = ax1.plot('Time(s)', 'Denoised dFF', linewidth=1, color='black', label='_GCaMP', data = behavprocesssnip_df)
            
        elif session in ['Test 1h','Test 24h','Test','S3']:
            #plots fiberpho trace and behaviour
            p1, = ax1.plot('Time(s)', 'Denoised dFF', linewidth=1, color='black', label='_GCaMP', data = behavprocesssnip_df)
            # p2, = ax1.plot(behavprocesssnip_df['Time(s)'], behavprocesssnip_df['Exploration fam'], linewidth=1, 
            #                color='moccasin', label='Exploration fam')
            
            #makes areas corresponding to behaviours
            i = 0
            j = 0
            m = 0
            n = 0
            if 'Exploration fam' in list_BOI:
                for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Exploration fam'].tolist()):
                    if y == 1:
                        x_start = x
                    if y == -1 and x_start!=0:
                        ax1.axvspan(x_start, x, facecolor='gold', alpha=0.3, label = '_'*i + 'Exploration fam')
                        x_start=0
                        i+=1
            for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Exploration new'].tolist()):
                if y == 1:
                    x_start = x
                if y == -1 and x_start!=0:
                    ax1.axvspan(x_start, x, facecolor='purple', alpha=0.3, label = '_'*j + 'Exploration new')
                    x_start=0
                    j+=1
            if 'Climbing' in list_BOI:
                for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Climbing'].tolist()):
                    if y == 1:
                        x_start = x
                    if y == -1 and x_start!=0:
                        ax1.axvspan(x_start, x, facecolor='cornflowerblue', alpha=0.3, label = '_'*m + 'Climbing')
                        x_start=0
                        m+=1
            if 'Rearing' in list_BOI:
                for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Rearing'].tolist()):
                    if y == 1:
                        x_start = x
                    if y == -1 and x_start!=0:
                        ax1.axvspan(x_start, x, facecolor='mediumpurple', alpha=0.3, label = '_'*n + 'Rearing')
                        x_start=0
                        n+=1
            
        elif session in ['Habituation','Training','S1']:
            #plots fiberpho trace and behaviourb 
            p1, = ax1.plot('Time(s)', 'Denoised dFF', linewidth=.6, color='black', label='_GCaMP', data = behavprocesssnip_df)
    
            #makes areas corresponding to behaviours
            i = 0
            j = 0
            m = 0
            n = 0
            for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Exploration left'].tolist()):
                if y == 1:
                    x_start = x
                if y == -1 and x_start!=0:
                    ax1.axvspan(x_start, x, facecolor='orange', alpha=0.3, label = '_'*i + 'Exploration left')
                    x_start=0
                    i+=1
            for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Exploration right'].tolist()):
                if y == 1:
                    x_start = x
                if y == -1 and x_start!=0:
                    ax1.axvspan(x_start, x, facecolor='darkturquoise', alpha=0.3, label = '_'*j + 'Exploration right')
                    x_start=0
                    j+=1
            if 'Climbing' in list_BOI:
                for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Climbing'].tolist()):
                    if y == 1:
                        x_start = x
                    if y == -1 and x_start!=0:
                        ax1.axvspan(x_start, x, facecolor='cornflowerblue', alpha=0.3, label = '_'*m + 'Climbing')
                        x_start=0
                        m+=1
            if 'Rearing' in list_BOI:
                for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Rearing'].tolist()):
                    if y == 1:
                        x_start = x
                    if y == -1 and x_start!=0:
                        ax1.axvspan(x_start, x, facecolor='mediumpurple', alpha=0.3, label = '_'*n + 'Rearing')
                        x_start=0
                        n+=1
        
        elif session in ['S2']:
            #plots fiberpho trace and behaviour
            p1, = ax1.plot('Time(s)', 'Denoised dFF', linewidth=.6, color='black', label='_GCaMP', data = behavprocesssnip_df)
        
            #makes areas corresponding to behaviours
            i = 0
            j = 0
            m = 0
            n = 0
            for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Exploration non social'].tolist()):
                if y == 1:
                    x_start = x
                if y == -1 and x_start!=0:   
                    ax1.axvspan(x_start, x, facecolor='grey', alpha=0.3, label = '_'*i + 'Exploration non social')
                    x_start=0
                    i += 1
            for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Exploration social'].tolist()):
                if y == 1:
                    x_start = x
                if y == -1 and x_start!=0:   
                    ax1.axvspan(x_start, x, facecolor='mediumvioletred', alpha=0.3, label = '_'*j + 'Exploration social')
                    x_start=0
                    j += 1
            if 'Climbing' in list_BOI:
                for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Climbing'].tolist()):
                    if y == 1:
                        x_start = x
                    if y == -1 and x_start!=0:
                        ax1.axvspan(x_start, x, facecolor='cornflowerblue', alpha=0.3, label = '_'*m + 'Climbing')
                        x_start=0
                        m+=1
                        
            if 'Rearing' in list_BOI:
                for (x,y) in zip(behavprocesssnip_df['Time(s)'].tolist(), behavprocesssnip_df['Rearing'].tolist()):
                    if y == 1:
                        x_start = x
                    if y == -1 and x_start!=0:
                        ax1.axvspan(x_start, x, facecolor='mediumpurple', alpha=0.3, label = '_'*n + 'Rearing')
                        x_start=0
                        n+=1
        
        #makes vertical line for entry opening of gate
        if 'Gate opens' in list_BOI:
            x_entry = behavprocess_df.at[int(np.where(behavprocess_df['Gate opens'] == 1)[0][0]), 'Time(s)']
            ax1.axvline(x_entry, color='lightsteelblue', ls = '--', label = 'Gate opens' )
                
        #makes vertical line for entry in open field
        x_entry = behavprocess_df.at[int(np.where(behavprocess_df['Entry in arena'] == 1)[0][0]), 'Time(s)']
        ax1.axvline(x_entry, color='slategrey', ls = '--', label = 'Entry in arena' )
    
    ax1.set_ylabel(r'$\Delta$F/F')
    ax1.set_xlabel('Seconds')
    #ax1.set_ylim(-1.5,2.5)
    ax1.legend(loc = 'upper right')
    ax1.margins(0.01,0.03)
    ax1.set_title(f'dFF with Behavioural Scoring - {exp} {session} {mouse} - interbout {THRESH_S} - cut {EVENT_TIME_THRESHOLD}')
    
    return fig2

def PETH(behavprocess_df, BOI, event, timewindow, EVENT_TIME_THRESHOLD, sr, PRE_EVENT_TIME):
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
    print(f'PETH {BOI}')
    #set time window relative to event
    PRE_TIME = timewindow[0]
    POST_TIME = timewindow[1]
    
    #round samplerate
    sr = round(sr)
    
    #if behaviour not recognized
    if BOI not in behavprocess_df.columns[2:].tolist() :
        print('BOI not recognized')
        
    list_ind_event_o = np.where(behavprocess_df[BOI] == 1)[0].tolist()
    list_ind_event_w = np.where(behavprocess_df[BOI] == -1)[0].tolist()
    
    #if event too short, don't process it
    EVENT_TIME_THRESHOLD = EVENT_TIME_THRESHOLD*sr
    for (start,end) in zip(list_ind_event_o, list_ind_event_w):
        #print(start,end)
        if end-start<EVENT_TIME_THRESHOLD and end-start>1:
            #print(end-start)
            list_ind_event_o.remove(start)
            list_ind_event_w.remove(end)
    
    #creates array of fiberpho data centered on event
    if event == 'onset':
        list_ind_event = list_ind_event_o
    elif event == 'withdrawal':
        list_ind_event = list_ind_event_w
        
    #if event happens too late in dataframe, don't process it
    if list_ind_event[-1]+POST_TIME*sr >= len(behavprocess_df):
        list_ind_event.pop(-1)
    
    #take the beginning of the behaviour as the beginning of the slope : min(dFF) on 3 secs before entry
    if event == 'onset':
        for i,ind_onset in enumerate(list_ind_event):
            list_ind_event[i] = behavprocess_df.loc[ind_onset-3*sr:ind_onset, 'Denoised dFF'].idxmin()
        
    PETH_array = np.zeros((len(list_ind_event),(POST_TIME+PRE_TIME)*sr+1))
    for i, ind_event in enumerate(list_ind_event) :
        #calculates baseline F0 on time window before event (from PRE_TIME to PRE_EVENT_TIME)
        F0 = np.mean(behavprocess_df.loc[ind_event-PRE_TIME*sr:ind_event-PRE_EVENT_TIME*sr, 'Denoised dFF'])
        std0 = np.std(behavprocess_df.loc[ind_event-PRE_TIME*sr:ind_event-PRE_EVENT_TIME*sr, 'Denoised dFF'])
        #creates array of z-scored dFF : z = (dFF-meandFF_baseline)/stddFF_baseline
        if len(behavprocess_df.loc[ind_event-PRE_TIME*sr:ind_event+(POST_TIME*sr), 'Denoised dFF']) == len(PETH_array[0]):
            PETH_array[i] = (behavprocess_df.loc[ind_event-PRE_TIME*sr:ind_event+(POST_TIME*sr), 'Denoised dFF']-F0)/std0
    
    return PETH_array

def plot_PETH(PETH_data, BOI, event, timewindow, exp, session, mouse):
    """
    Plots PETH average and heatmap
    """
    PRE_TIME = timewindow[0]
    POST_TIME = timewindow[1]
    
    #create figure
    fig4 = plt.figure(figsize=(6,10))
    
    #create time vector
    peri_time = np.arange(-PRE_TIME, POST_TIME+0.1, 0.1)
    
    #calculate mean dFF and std error
    mean_dFF_snips = np.mean(PETH_data, axis=0)
    std_dFF_snips = np.std(PETH_data, axis=0)
        
    #plot individual traces and mean
    ax5 = fig4.add_subplot(212)
    for snip in PETH_data:
        p1, = ax5.plot(peri_time, snip, linewidth=.5,
                       color=[.7, .7, .7], label='Individual trials')
    p2, = ax5.plot(peri_time, mean_dFF_snips, linewidth=2,
                   color='green', label='Mean response')
    
    #plot standard error bars
    p3 = ax5.fill_between(peri_time, mean_dFF_snips+std_dFF_snips,
                      mean_dFF_snips-std_dFF_snips, facecolor='green', alpha=0.2)
    p4 = ax5.axvline(x=0, linewidth=2, color='slategray', ls = '--', label='Exploration '+event)
    
    #ax5.axis('tight')
    ax5.set_xlabel('Time(s)')
    ax5.set_ylabel('z-scored $\Delta$F/F')
    ax5.legend(handles=[p1, p2, p4], loc='upper left', fontsize = 'small')
    ax5.margins(0,0.01)
    
    #add heatmap
    ax6 = fig4.add_subplot(211)
    cs = ax6.imshow(PETH_data, cmap='magma', aspect = 'auto',
                    interpolation='none', extent=[-PRE_TIME, POST_TIME, len(PETH_data), 0],
                    vmin = -6, vmax = 9)
    ax6.set_ylabel('Bout Number')
    ax6.set_yticks(np.arange(.5, len(PETH_data), 2))
    ax6.set_yticklabels(np.arange(0, len(PETH_data), 2))
    ax6.axvline(x=0, linewidth=2, color='black', ls = '--')
    ax6.set_title(f'{BOI} {event} {exp} {session} {mouse}')
    
    fig4.subplots_adjust(right=0.8, hspace = 0.1)
    cbar_ax = fig4.add_axes([0.85, 0.54, 0.02, 0.34])
    fig4.colorbar(cs, cax=cbar_ax)
    
    return fig4

def plot_PETH_average(PETH_data, BOI, event, timewindow, exp, session, mouse):
    """
    Plots only average of PETH, no heatmap
    For BOIs only happening once in the trial
    """
    
    PRE_TIME = timewindow[0]
    POST_TIME = timewindow[1]
    
    #create figure
    fig4 = plt.figure(figsize=(6,4))
    ax5 = fig4.add_subplot(111)
    
    #create time vector
    peri_time = np.arange(-PRE_TIME, POST_TIME+0.1, 0.1) 
    
    #plot trace
    p1, = ax5.plot(peri_time, PETH_data[0], linewidth=2, color='green')
    p2 = ax5.axvline(x=0, linewidth=2, color='slategray', ls = '--', label=BOI)
    
    ax5.set_xlabel('Seconds')
    ax5.set_ylabel('z-scored $\Delta$F/F')
    ax5.legend(handles=[p2], loc='upper left', fontsize = 'small')
    ax5.margins(0, 0.1)
    ax5.set_title(f'{BOI} - {exp} {session} {mouse}')
    
    return fig4
    
def plot_PETH_pooled(included_groups, PETHarray_list, BOI, event, timewindow, exp, session):
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
                   color='grey', label=included_groups[0])   
    #plot standard error bars
    p3 = ax5.fill_between(peri_time, listmean_dFF_snips[0]+(listsem_dFF_snips[0]/np.sqrt(len(PETHarray_list[0]))),
                      listmean_dFF_snips[0]-(listsem_dFF_snips[0]/np.sqrt(len(PETHarray_list[0]))), facecolor='grey', alpha=0.2)
    
    #plot individual traces and mean HFD
    # for snip in PETHarray_list[1]:
    #     p4, = ax5.plot(peri_time, snip, linewidth=.5,
    #                    color='coral', alpha=.3)
    p5, = ax5.plot(peri_time, listmean_dFF_snips[1], linewidth=2,
                   color='purple', label=included_groups[1])   
    #plot standard error bars
    p6 = ax5.fill_between(peri_time, listmean_dFF_snips[1]+(listsem_dFF_snips[1]/np.sqrt(len(PETHarray_list[1]))),
                      listmean_dFF_snips[1]-(listsem_dFF_snips[1]/np.sqrt(len(PETHarray_list[1]))), facecolor='purple', alpha=0.2)
    
    p8 = ax5.axvline(x=0, linewidth=2, color='slategray', ls = '--', label=BOI)
    
    ax5.set_xlabel('Seconds')
    ax5.set_ylabel('z-scored $\Delta$F/F')
    ax5.legend(handles=[p2,p5,p8], loc='upper left', fontsize = 'small')
    ax5.margins(0, 0.1)
    ax5.set_title(f'{BOI} - {exp} {session} {included_groups[0]} {included_groups[1]}')
    
    return fig4