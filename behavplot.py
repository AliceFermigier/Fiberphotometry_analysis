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
from pathlib import Path
import os

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
    list_ind = np.arange(len(list_behav), dtype = int).tolist()
    behav_comp = [0]*len(list_behav)
    
    #index where camera starts
    list_indstart = np.where(round(timevector,1) == timestart_camera)
    indstart = list_indstart[0].tolist()[0]

    # create lists of behaviour data for each scored behaviour
    # aligned with start of the camera
    for (behav, ind) in zip(list_behav, list_ind) :
        behav_comp[ind] = [0]*indstart
        behav_comp[ind].extend(behav10Sps[behav].tolist())
       
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
    timelist = timelist[:min_length]
    denoised_fiberpho_list = denoised_fiberpho_list[:min_length]
    dff_405nm_list = dff_405nm_list[:min_length]
    dff_470nm_list = dff_470nm_list[:min_length]
    behav_crop = []
    for behav in behav_comp:
        behav_crop.append(behav[:min_length])
        
    if len(list_behav) == 1 :
        fiberbehav_df = pd.DataFrame(data = {'Time(s)':timelist, 'Denoised dFF' : denoised_fiberpho_list,
                                             '405nm deltaF/F' : dff_405nm_list, '470nm deltaF/F' : dff_470nm_list, 
                                             list_behav[0] : behav_crop[0]})
        
    elif len(list_behav) == 2 :
        fiberbehav_df = pd.DataFrame(data = {'Time(s)':timelist, 'Denoised dFF' : denoised_fiberpho_list,
                                             '405nm deltaF/F' : dff_405nm_list, '470nm deltaF/F' : dff_470nm_list,
                                             list_behav[0] : behav_crop[0], list_behav[1] : behav_crop[1]})
        
    elif len(list_behav) == 3 :
        fiberbehav_df = pd.DataFrame(data = {'Time(s)':timelist, 'Denoised dFF' : denoised_fiberpho_list,
                                             '405nm deltaF/F' : dff_405nm_list, '470nm deltaF/F' : dff_470nm_list,
                                             list_behav[0] : behav_crop[0], list_behav[1] : behav_crop[1],
                                             list_behav[2] : behav_crop[2]})
    
    elif len(list_behav) == 4 :
        fiberbehav_df = pd.DataFrame(data = {'Time(s)':timelist, 'Denoised dFF' : denoised_fiberpho_list,
                                             '405nm deltaF/F' : dff_405nm_list, '470nm deltaF/F' : dff_470nm_list,
                                             list_behav[0] : behav_crop[0], list_behav[1] : behav_crop[1],
                                             list_behav[2] : behav_crop[2], list_behav[3] : behav_crop[3] })
    
    elif len(list_behav) == 5 :
        fiberbehav_df = pd.DataFrame(data = {'Time(s)':timelist, 'Denoised dFF' : denoised_fiberpho_list,
                                             '405nm deltaF/F' : dff_405nm_list, '470nm deltaF/F' : dff_470nm_list,
                                             list_behav[0] : behav_crop[0], list_behav[1] : behav_crop[1],
                                             list_behav[2] : behav_crop[2], list_behav[3] : behav_crop[3],
                                             list_behav[4] : behav_crop[4]})
    else :
        print('Error : too many behaviours in Boris binary file. Score up to 5 or add line to function')
        
    if exp=='NewContext':
        fiberbehav_df.loc[:405,'Homecage']=0
    
    return(fiberbehav_df)

def behav_process(fiberbehav_df, list_BOI, THRESH_S, EVENT_TIME_THRESHOLD):

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
                if count <= THRESH_S*10:
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
                if count <= EVENT_TIME_THRESHOLD*10:
                    fiberbehav_df.loc[ind-count:ind-1, BOI] = [0]*count
            i_1 = i
            
    return fiberbehav_df 
                    
def derive(fiberbehav_df, list_BOI):
    """
    calculate the derivative of behav of interest and put in new df
    that way, it will show 1 when behaviour starts and -1 when it stops
    """
    for col in fiberbehav_df.columns[2:]:
        fiberbehav_df[col] = fiberbehav_df[col].diff()

    return fiberbehav_df 