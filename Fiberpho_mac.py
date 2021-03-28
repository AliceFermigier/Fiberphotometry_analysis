# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 09:30:55 2021

Analysis of Doric fiber photometry data 

DeltaF/F data for 470nm and 405nm are calculated using Doric Neuroscience Studio
Behavioral scoring is done using Boris software

@author: Alice Fermigier
"""
#%%
##########
#IMPORTED#
##########

import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
from pathlib import Path
#import os
#with pathlib.Path : file_to_open = data_folder / "raw_data.txt"

#%%
########
#LOADER#
########

analysis_path = Path('/Volumes/My Passport/Alice/Fiber/20210302_ORM/Analysis')



behav_path = analysis_path / 'ORM1_CD1_1.csv'
fiberpho_path = analysis_path / 'CD1_1_procbis_tracesfilt.csv'
camera_path = analysis_path / 'CD1_1_procbis_camera.csv'

behav10Sps = pd.read_csv(behav_path)
fiberpho = pd.read_csv(fiberpho_path)
camera = pd.read_csv(camera_path)

#all analysis files in Experiment/Analysis/Subject_Trial
#each trial of each subject has a specific directory called CD1_0, CD1_2, etc...
#for example, for ORM : CD1_0 = training of subject CD1 ; CD1_2 = test of subject CD1

# for trial in os.path.listdir(analysis_path):
#     behav_path = analysis_path / trial / str(trial) + '_behav.csv'
#     fiberpho_path = analysis_path / trial / str(trial) + '_proc_tracesfilt.csv'
#     camera_path = analysis_path / trial / str(trial) + '_proc_camera.csv'
    
#     behav10Sps = pd.read_csv(behav_path)
#     fiberpho = pd.read_csv(fiberpho_path)
#     camera = pd.read_csv(camera_path)
    
#     # to be created 
    
#     plot_fiberpho_path = analysis_path / trial / str(trial) + '_fiberpho_isosgcamp.png'
#     plot_fiberbehav_path = analysis_path / trial / str(trial) + '_fiberbehav.png'
#     plot_fiberbehavsnip_path = analysis_path / trial / str(trial) + '_fiberbehavsnip.png'
    

samplerate = 10 #in samples per second

#CAUTION : with the code as it is, can process up to 4 different behaviours
#if more, add elif to align_behav function



#%%
###################
#DEFINED FUNCTIONS#
###################

def truncate(n, decimals=0):
    multiplier = 10 ** decimals
    return int(n * multiplier) / multiplier

def time_vector(fiberpho, samplerate) :
    """
    Creates timevector on which to plot the data, in pd format
    --> Parameters :
        fiberpho = float, duration of trial in secs
        samplerate = int, in Sps (for processed fiberpho data in Doric Neuroscience Studio, samplerate = 10Sps)
    --> Returns :
        timevector = pd series
    """
    #denoised_fiberpho = fiberpho['Analog In. | Ch.1 470 nm (Deinterleaved)_dF/F0_LowPass' + 
    #                             '-Analog In. | Ch.1 405 nm (Deinterleaved)_dF/F0_LowPass'].dropna()
    #--> if better timevector of the exact same lenght as fiberpho data
    
    duration =  math.ceil(fiberpho.at[len(fiberpho)-1,'Time(s)'])
    return pd.Series(np.linspace(0.0, duration, num = int(duration*samplerate)+1))

def timestamp_camera(camera) :
    """
    Function to extract the timestamps where the camera starts and stops
    --> Parameters
        camera : pd dataframe, camera I/O with sample rate = 12kSps
    --> Returns
        (camera_start, camera_stop) = timestamp when camera starts and stops in seconds (truncated to 0,1s) #camera_stop à enlever si pas besoin
    """
    ind_list = np.where(camera['Digital I/O | Ch.3 DI/O-3'] == 1)
    ind_list = ind_list[0].tolist()
    (ind_start, ind_stop) = (ind_list[0],ind_list[len(ind_list)-1])
    return (truncate(camera.at[ind_start, 'Time(s)'], 1),
            truncate(camera.at[ind_stop, 'Time(s)'], 1))

# def indstart_camera(camera) : #on s'en fout, mauvais sample rate de toutes facons donc mauvais index
#     """
#     Function to extract the timestamps where the camera starts and stops
#     --> Parameters
#         camera : pd dataframe, camera I/O with sample rate = 12kSps
#     --> Returns
#         indstart_camera : int, index of df where camera starts
#     """
#     ind_list = np.where(camera['Digital I/O | Ch.3 DI/O-3'] == 1)
#     return (ind_list[0])
    
def align_behav(behav10Sps, fiberpho, timevector, timestart_camera):
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
    
    #index where camrera starts
    list_indstart = np.where(timevector == timestart_camera)
    indstart = list_indstart[0].tolist()[0]

    # create lists of behaviour data for each scored behaviour
    # aligned with start of the camera
    for (behav, ind) in zip(list_behav, list_ind) :
        behav_comp[ind] = [0]*indstart
        behav_comp[ind].extend(behav10Sps[behav].tolist())
       
    # create list of denoised fiberpho data    
    denoised_fiberpho = fiberpho['Analog In. | Ch.1 470 nm (Deinterleaved)_dF/F0_LowPass' + 
                                 '-Analog In. | Ch.1 405 nm (Deinterleaved)_dF/F0_LowPass'].dropna()
    denoised_fiberpho_list = denoised_fiberpho.tolist()
    
    # create list of isosbestic (405) and Ca dependent (470) data
    dff_405nm = fiberpho['Analog In. | Ch.1 405 nm (Deinterleaved)_dF/F0_LowPass'].dropna()
    dff_470nm = fiberpho['Analog In. | Ch.1 470 nm (Deinterleaved)_dF/F0_LowPass'].dropna()
    
    # makes timevector into a list
    timelist = timevector.tolist()
    
    # crops lists so that all lengths match
    min_length = min([len(timelist), len(denoised_fiberpho_list), len(behav_comp[0])])
    timelist = timelist[:min_length]
    denoised_fiberpho_list = denoised_fiberpho_list[:min_length]
    behav_crop = []
    for behav in behav_comp:
        behav_crop.append(behav[:min_length])
        
    print(len(timelist), len(denoised_fiberpho_list), len(behav_crop[0]))
        
    if len(list_behav) == 1 :
        fiberbehav_df = pd.DataFrame(data = {'Time(s)':timelist, 'Denoised deltaF/F' : denoised_fiberpho_list,
                                             '405nm deltaF/F' : dff_405nm, '470nm deltaF/F' : dff_470nm, 
                                             list_behav[0] : behav_crop[0]})
        
    elif len(list_behav) == 2 :
        fiberbehav_df = pd.DataFrame(data = {'Time(s)':timelist, 'Denoised deltaF/F' : denoised_fiberpho_list,
                                             '405nm deltaF/F' : dff_405nm, '470nm deltaF/F' : dff_470nm,
                                             list_behav[0] : behav_crop[0], list_behav[1] : behav_crop[1]})
        
    elif len(list_behav) == 3 :
        fiberbehav_df = pd.DataFrame(data = {'Time(s)':timelist, 'Denoised deltaF/F' : denoised_fiberpho_list,
                                             list_behav[0] : behav_crop[0], list_behav[1] : behav_crop[1],
                                             list_behav[2] : behav_crop[2]})
    
    elif len(list_behav) == 4 :
        fiberbehav_df = pd.DataFrame(data = {'Time(s)':timelist, 'Denoised deltaF/F' : denoised_fiberpho_list,
                                             '405nm deltaF/F' : dff_405nm, '470nm deltaF/F' : dff_470nm,
                                             list_behav[0] : behav_crop[0], list_behav[1] : behav_crop[1],
                                             list_behav[2] : behav_crop[2], list_behav[3] : behav_crop[3] })
        
    else :
        print('Error : too many behaviours in Boris binary file. Score up to 4 or add line to function')
        
    
    return fiberbehav_df

#%%

def plot_fiberpho(fiberbehav_df):
    """
    Plots isosbestic and Ca dependent deltaF/F
    """
    
    # Plot both unprocessed demodulated stream            
    fig1 = plt.figure(figsize=(10,6))
    ax0 = fig1.add_subplot(111)
    
    # Plotting the traces
    p1, = ax0.plot('Time(s)', '470nm deltaF/F', linewidth=2, color='dodgerblue', label='GCaMP', data = fiberbehav_df) 
    p2, = ax0.plot('Time(s)', '405nm deltaF/F', linewidth=2, color='blueviolet', label='ISOS', data = fiberbehav_df)
    
    ax0.set_ylabel(r'$\Delta$F/F')
    ax0.set_xlabel('Seconds')
    ax0.set_title('GCaMP and Isosbestic dFF')
    ax0.legend(handles=[p1,p2], loc='upper right')
    fig1.tight_layout()
    
    return

#%%
def plot_fiberpho_behav(fiberbehav_df):
    """
    Plots denoised deltaF/F aligned with behaviour (includes baseline)
    """
    
    fig3 = plt.figure(figsize=(20,12))
    ax2 = fig3.add_subplot(311)
    
    #y_scale = 10 #adjust according to data needs
    #y_shift = -30 #scale and shift are just for asthetics
    
    p1, = ax2.plot('Time(s)', 'Denoised deltaF/F', linewidth=2, color='black', label='GCaMP', data = fiberbehav_df)
    p2, = ax2.plot('Time(s)', 'Exploration fam', linewidth=0.01, color='moccasin', label='Exploration fam', data = fiberbehav_df)
    p3, = ax2.plot('Time(s)', 'Exploration new', linewidth=0.01, color='cornflowerblue', label='Exploration new', data = fiberbehav_df)
    for (x,y) in zip(fiberbehav_df['Time(s)'].tolist(), fiberbehav_df['Exploration fam'].tolist()):
        if y == 1:
            ax2.axvspan(x, x+0.1, facecolor='moccasin', alpha=0.5)
    for (x,y) in zip(fiberbehav_df['Time(s)'].tolist(), fiberbehav_df['Exploration new'].tolist()):
        if y == 1:
            ax2.axvspan(x, x+0.1, facecolor='cornflowerblue', alpha=0.5)
    entry = np.where()
    for (x,y) in zip(fiberbehav_df['Time(s)'].tolist(), fiberbehav_df['Exploration new'].tolist()):
        ax2.vlines()
    ax2.set_ylabel(r'$\Delta$F/F')
    ax2.set_xlabel('Seconds')
    ax2.set_title('dFF with Behavioural Scoring')
    ax2.legend(handles=[p1,p2,p3], loc='upper right')
    fig3.tight_layout()
    
    return

#%%
def plot_fiberpho_behav_snip(fiberbehav_df, timestart_camera):
    """
    Plots denoised deltaF/F aligned with behaviour (starts when camera starts)
    """
    
    fig3 = plt.figure(figsize=(20,12))
    ax2 = fig3.add_subplot(311)
    
    #y_scale = 10 #adjust according to data needs
    #y_shift = -30 #scale and shift are just for asthetics
    
    fiberbehavsnip_df = fiberbehav_df[fiberbehav_df['Time(s)'] > timestart_camera]
    
    p1, = ax2.plot('Time(s)', 'Denoised deltaF/F', linewidth=2, color='black', label='GCaMP', data = fiberbehavsnip_df)
    p2, = ax2.plot('Time(s)', 'Exploration fam', linewidth=0.01, color='moccasin', label='Exploration fam', data = fiberbehavsnip_df)
    p3, = ax2.plot('Time(s)', 'Exploration new', linewidth=0.01, color='cornflowerblue', label='Exploration new', data = fiberbehavsnip_df)
    for (x,y) in zip(fiberbehav_df['Time(s)'].tolist(), fiberbehav_df['Exploration fam'].tolist()):
        if y == 1:
            ax2.axvspan(x, x+0.1, facecolor='moccasin', alpha=0.5)
    for (x,y) in zip(fiberbehav_df['Time(s)'].tolist(), fiberbehav_df['Exploration new'].tolist()):
        if y == 1:
            ax2.axvspan(x, x+0.1, facecolor='cornflowerblue', alpha=0.5)
    ax2.set_ylabel(r'$\Delta$F/F')
    ax2.set_xlabel('Seconds')
    ax2.set_title('dFF with Behavioural Scoring')
    ax2.legend(handles=[p1,p2,p3], loc='upper right')
    fig3.tight_layout()
    
    return

#%%
def behav_averaging(fiberbehav_df):
    
    #1 turn separate sniffings into bouts of at least 1 sec
    
    #2 average dFF after bout onset and after bout ending
    
    #3 average dFF peak for each animal and then average? on calculate peak of average trace?
    
    #4 calculate latency to peak
    
    #5 calculate peak for each animal during bout and then average
    
    return

def plot_heatmaps(fiberbehav_df):

    #plot heatmaps of fluorescence on bout onset and bout ending

    return
    
#%%
########
#SCRIPT#
########


#Tests
######

timevector = time_vector(fiberpho, samplerate)

timestart_camera = timestamp_camera(camera)[0]

fiberbehav_df = align_behav(behav10Sps, fiberpho, timevector, timestart_camera)
fiberbehav_df.to_excel('/Users/alice/Desktop/Test.xlsx')
print(fiberbehav_df)
#plot_fiberpho_behav(fiberbehav_df)
