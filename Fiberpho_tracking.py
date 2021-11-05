# -*- coding: utf-8 -*-
"""
Created on Fri Jun 18 11:41:03 2021

Plot tracking from Bonsai csv file
Bonsai file used to generate file : tracking_fiberOF.bonsai

@author: Alice Fermigier
"""
#%%
##########
#IMPORTED#
##########

import matplotlib.pyplot as plt
plt.rcParams.update({'figure.max_open_warning' : 0})
from pathlib import Path
import pandas as pd
import os
import numpy as np
import plotly.express as px
import sys

#Custom
#put path to directory where python files are stored
sys.path.append('/Users/alice/Documents/GitHub/Fiberphotometry_analysis')
#sys.path.append('C:\\Users\\afermigier\\Documents\\GitHub\\Fiberphotometry_analysis')


#%%
########
#LOADER#
########

#from Fiberpho_loader import experiment_path, analysis_path, data_path, subjects_df
#from Fiberpho_loader import SAMPLERATE

from Func_fiberplots import session_code, truncate, time_vector, timestamp_camera, timestamp_camera_fromraw

#os.chdir(experiment_path)


#%%
########
#SCRIPT#
########

def plot_tracking(df_tracking):
    
    fig1 = plt.figure(figsize=(20,20))
    ax1 = fig1.add_subplot(111)
    
    p1, = ax1.plot('Item2.X', 'Item2.Y', linewidth=.5, color='black', data=df_tracking)
    
    #fig1.savefig(mouse_path / f'{mouse}_tracking.pdf')

#aborted
# def plot_trackingdensity(df_tracking):
#     """
#     Plots density of tracking for each animal
#     """
#     #lat +90 à -90
#     #lon -180 à +180
    
#     (min_x, max_x) = (min(df_tracking['Item2.X'].dropna()), max(df_tracking['Item2.X'].dropna()))
#     (min_y, max_y) = (min(df_tracking['Item2.Y'].dropna()), max(df_tracking['Item2.Y'].dropna()))
    
#     x_mid = (max_x+min_x)/2
#     y_mid = (max_y+min_y)/2
    
#     range_x = max_x-min_x
#     range_y = max_y-min_y
    
#     #normalized coordinates
#     lon_coords = [(x-x_mid)*90/range_x for x in df_tracking['Item2.X']]
#     lat_coords = [(y-y_mid)*90/range_y for y in df_tracking['Item2.Y']]
    
#     #put all the values to 1
#     values = np.ones(len(df_tracking['Item2.X']))
    
#     #plot figure
#     fig2 = px.density_mapbox(lat=lat_coords, lon=lon_coords, z=values, radius=10,
#                              zoom=0, opacity=1, title = f'Position density - {exp} {session} {mouse}',
#                              mapbox_style = 'white-bg')
#     fig2.show()
#     fig2.write_image(mouse_path / f'{mouse}_trackingdensity.pdf')
    
#     return()

def plot_densityheatmap(df_tracking):
    
    df_tracking_clean = df_tracking.dropna()
    
    #find space edges
    (min_x, max_x) = (min(df_tracking_clean['Item2.X']), max(df_tracking_clean['Item2.X']))
    (min_y, max_y) = (min(df_tracking_clean['Item2.Y']), max(df_tracking_clean['Item2.Y']))
    
    #define matrix
    space = np.zeros((RES, RES))-0.1
    
    #discretize position in dataframe
    dis_x = [int(((RES-1)*(x-min_x))/(max_x-min_x)) for x in df_tracking_clean['Item2.X']]
    dis_y = [int(((RES-1)*(y-min_y))/(max_y-min_y)) for y in df_tracking_clean['Item2.Y']]
    
    density_df = pd.DataFrame(data = {'Track_x' : dis_x, 'Track_y' : dis_y,
                                      'Count' : [1]*len(dis_x)})
    
    density = density_df.groupby(['Track_x','Track_y'], as_index=False).sum()
    density['Density'] = [c/max(density['Count']) for c in density['Count']]
    
    for (i,j) in zip(density['Track_x'], density['Track_y']):
        space[j][i] = density.loc[(density['Track_x']==i) & (density['Track_y']==j),'Density'].values[0]
    
    #plot space array
    fig2 = plt.figure(figsize=(30,30))
    ax2 = fig2.add_subplot(111)
    
    p2 = ax2.contourf(space, cmap='magma', vmin = -0.1, vmax = 1, levels=15)
    #fig2.colorbar(p2, ax=ax2)
    ax2.set_title(f'Position density - {exp} {session} {mouse}', fontsize=28)
    ax2.set_xticks([])
    ax2.set_yticks([])
    
    #for discrete visualization
    # ax2b = fig2.add_subplot(212,sharex=ax2)
    # p2b = ax2b.pcolormesh(space, cmap='inferno', vmin = -0.1, vmax = 1)
    
    fig2.savefig(mouse_path / f'{mouse}_{code}_densityheatmap.png', transparent=True)
    
    return()
    
def align_dFFtrack(df_tracking, fiberpho, timevector, camera_start, camera_stop):
    
    #add timeframe to tracking data
    timevector_tracking = np.linspace(camera_start, camera_stop, len(df_tracking))
    timevector_tracking_trunc = [truncate(i,1) for i in timevector_tracking]
    df_tracking['Time(s)']=timevector_tracking_trunc
    df_tracking_resampled = df_tracking.drop_duplicates('Time(s)')
    
    #index of df_tracking where camera starts
    list_indstart = np.where(round(timevector,1) == camera_start)
    indstart = list_indstart[0].tolist()[0]
    
    #align tracking data by adding zeros before camera starts
    tracking_x = [None]*indstart
    tracking_x.extend(df_tracking_resampled['Item2.X'])
    
    tracking_y = [None]*indstart
    tracking_y.extend(df_tracking_resampled['Item2.Y'])
    
    # creates list of denoised fiberpho data    
    denoised_fiberpho = fiberpho['Analog In. | Ch.1 470 nm (Deinterleaved)_dF/F0' +
                                 '-Analog In. | Ch.1 405 nm (Deinterleaved)_dF/F0_LowPass'].dropna()
    denoised_fiberpho_list = denoised_fiberpho.tolist()
    
    # creates list of isosbestic (405) and Ca dependent (470) data
    dff_405nm = fiberpho['Analog In. | Ch.1 405 nm (Deinterleaved)_dF/F0_LowPass'].dropna()
    dff_405nm_list = dff_405nm.tolist()
    dff_470nm = fiberpho['Analog In. | Ch.1 470 nm (Deinterleaved)_dF/F0_LowPass'].dropna()
    dff_470nm_list = dff_470nm.tolist()
    
    # makes timevector into a list
    timelist = timevector.tolist()
    
    # crops lists so that all lengths match
    min_length = min([len(timelist), len(denoised_fiberpho_list), len(tracking_x),
                      len(dff_405nm_list), len(dff_470nm_list)])
    timelist = timelist[:min_length]
    denoised_fiberpho_list = denoised_fiberpho_list[:min_length]
    dff_405nm_list = dff_405nm_list[:min_length]
    dff_470nm_list = dff_470nm_list[:min_length]
    tracking_x = tracking_x[:min_length]
    tracking_y = tracking_y[:min_length]

    # create dataframe
    fibertracking_df = pd.DataFrame(data = {'Time(s)':timelist, 'Denoised dFF' : denoised_fiberpho_list,
                                             'Track_x' : tracking_x, 'Track_y' : tracking_y})
    
    return(fibertracking_df)

def plot_fibertrack_heatmap(fibertracking_df):
    
    fibertracking_df_clean = fibertracking_df.dropna()
    
    #find space edges
    (min_x, max_x) = (min(fibertracking_df_clean['Track_x']), max(fibertracking_df_clean['Track_x']))
    (min_y, max_y) = (min(fibertracking_df_clean['Track_y']), max(fibertracking_df_clean['Track_y']))
    
    #discretize position in dataframe
    dis_x = [int(((RES-1)*(x-min_x))/(max_x-min_x)) for x in fibertracking_df_clean['Track_x']]
    dis_y = [int(((RES-1)*(y-min_y))/(max_y-min_y)) for y in fibertracking_df_clean['Track_y']]
    
    fibertracking_dis_df = pd.DataFrame(data = {'Denoised dFF' : fibertracking_df_clean['Denoised dFF'],
                                                'Track_x' : dis_x, 'Track_y' : dis_y})
    
    mean_dFF = fibertracking_dis_df.groupby(['Track_x','Track_y'], as_index=False).mean()
    (min_dFF, max_dFF) = (min(mean_dFF['Denoised dFF']), max(mean_dFF['Denoised dFF']))
    
    #define matrix
    space = np.zeros((RES, RES))-(abs(min_dFF)+0.05)
    
    for (i,j) in zip(mean_dFF['Track_x'], mean_dFF['Track_y']):
        space[j][i] = mean_dFF.loc[(mean_dFF['Track_x']==i) & (mean_dFF['Track_y']==j),'Denoised dFF']
    
    #plot space array contourf
    fig3 = plt.figure(figsize=(30,36))
    ax3 = fig3.add_subplot(111)
    
    p3 = ax3.contourf(space, cmap='turbo', vmin = min_dFF, vmax = max_dFF, levels=20)
    #p3 = ax3.pcolormesh(space, cmap='turbo', vmin = min_dFF, vmax = max_dFF)
    fig3.colorbar(p3, ax=ax3, aspect=60, orientation = 'horizontal', pad=0.01)
    ax3.set_title(f'dFF tracking - {exp} {session} {mouse}', fontsize=28)
    #fig3.suptitle(f'dFF tracking - {exp} {session} {mouse}', fontsize=20)
    ax3.set_xticks([])
    ax3.set_yticks([])
    
    #plot colormesh
    fig3b = plt.figure(figsize=(30,36))
    ax3b = fig3b.add_subplot(111)
    
    p3b = ax3b.pcolormesh(space, cmap='turbo', vmin = min_dFF, vmax = max_dFF)
    #p3 = ax3.pcolormesh(space, cmap='turbo', vmin = min_dFF, vmax = max_dFF)
    fig3b.colorbar(p3b, ax=ax3b, aspect=60, orientation = 'horizontal', pad=0.01)
    ax3b.set_title(f'dFF tracking - {exp} {session} {mouse}', fontsize=28)
    #fig3.suptitle(f'dFF tracking - {exp} {session} {mouse}', fontsize=20)
    ax3b.set_xticks([])
    ax3b.set_yticks([])
    
    fig3.savefig(mouse_path / f'{mouse}_{code}_dFFheatmap.png', transparent=True)
    fig3b.savefig(mouse_path / f'{mouse}_{code}_dFFheatmap_cmesh.png', transparent=True)
    
    return()

def plot_interactive_heatmap(fibertracking_df, arena):
    
    return()
    
#for plotly animations : https://plotly.com/python/sliders/


#%%Run test

SAMPLERATE = 10 #in Hz

mouse_path = Path('/Users/alice/Desktop/Data_test')
exp = 'OdDis'
session = 'Habituation'
mouse = 'CDf1'
code = 0
RES = 150 #number of pixels in fiberpho heatmap

fiberpho_path = mouse_path / f'{mouse}_{code}_dFFfilt.csv'
camera_path = mouse_path / f'{mouse}_{code}_camera.csv'
#rawdata_path = mouse_path / f'{mouse}_{code}.csv'
tracking_path = mouse_path / f'{mouse}_{code}_tracking.csv'
    
fiberpho = pd.read_csv(fiberpho_path)
camera = pd.read_csv(camera_path)
df_tracking = pd.read_csv(tracking_path)
#rawdata_df = pd.read_csv(rawdata_path)

print('timevector')
timevector = time_vector(fiberpho, SAMPLERATE)
print('timestamp')
(camera_start, camera_stop) = timestamp_camera(camera)
print('start camera : ', camera_start)
print('aligning')
fibertracking_df = align_dFFtrack(df_tracking, fiberpho, timevector, camera_start, camera_stop)

#%%plot tracking
plot_tracking(df_tracking)

#%%plot tracking density
plot_densityheatmap(df_tracking)

#%%plot 
plot_fibertrack_heatmap(fibertracking_df)

#%%Run for all 

#list experiments, sessions and mice:

                
                    

