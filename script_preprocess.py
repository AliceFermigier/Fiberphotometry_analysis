# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 13:39:29 2023

To run fiberphotometry analysis with behaviour or plethysmography data
1 - PREPROCESSING 

@author: alice fermigier
"""

#%%
##########
#IMPORTED# 
##########

import pandas as pd
from pathlib import Path
import os
from dash import Dash, dcc, html 
import plotly.express as px
import sys  
import matplotlib.pyplot as plt

#path to other scripts in sys.path
path_to_gitrepo=r'C:\Users\alice\Documents\GitHub\Fiberphotometry_analysis'
if path_to_gitrepo not in sys.path:
    sys.path.append(path_to_gitrepo)

#import functions
import preprocess as pp
import genplot as gp

#%%
########
#LOADER#
########

experiment_path = Path(r'C:\Users\alice\Desktop\Data_Fiber')
analysis_path = experiment_path / 'Analysis' 
data_path = experiment_path / 'Data'
os.chdir(experiment_path)
os.getcwd()

#import ID and groups of all mice
subjects_df = pd.read_excel(experiment_path / 'subjects.xlsx', sheet_name='Included')
#import tasks in protocol
proto_df = pd.read_excel(experiment_path / 'protocol.xlsx')

############
#PARAMETERS#
############

#time to crop at the beginning of the trial for , in seconds
TIME_BEGIN = 60
#filter characteristics
ORDER = 4
CUT_FREQ = 1 #in Hz

#%% 1 - PREPROCESSING
#####################

#------------------#
exp = 'OdDis1'
#------------------#

exp_path = analysis_path / exp
if not os.path.exists(exp_path):
    os.mkdir(exp_path)
for session_list in proto_df.loc[proto_df['Task']==exp,'Sessions'].values[0].split():
    if not os.path.exists(exp_path / session_list):
        os.mkdir(exp_path / session_list)
            
#get data path related to the task in protocol excel file
data_path_exp = data_path / proto_df.loc[proto_df['Task']==exp, 'Data_path'].values[0]

# put all raw plots and deinterleaved_df in separate folder
pp_path = data_path_exp / 'Preprocessing'
if not os.path.exists(pp_path):
    os.mkdir(pp_path)

#%% 1.1 - Deinterleave data and save in separate file

# Loop through each session directory
for session_path in [Path(f.path) for f in os.scandir(exp_path) if f.is_dir()]:
    session = session_path.name  # Extract session name
    print('##########################################')
    print(f'EXPERIMENT : {exp} - SESSION : {session}')
    print('##########################################')
    
    # Generate session code for the current session
    code = gp.session_code(session, exp)
    
    # Loop through each mouse in the subject DataFrame
    for mouse in subjects_df['Subject']:
        print("--------------")
        print(f'MOUSE : {mouse}')
        print("--------------")
        
        # Create a common file prefix
        file_prefix = f'{mouse}_{code}'
        
        # Paths for input raw data and output deinterleaved and plot files
        raw_data_path = data_path_exp / f'{file_prefix}.csv'
        deinterleaved_path = pp_path / f'{file_prefix}_deinterleaved.csv'
        raw_plot_path = pp_path / f'{file_prefix}_rawdata.png'
        
        # Check if raw data exists and deinterleaved data does not exist
        if raw_data_path.exists() and not deinterleaved_path.exists():
            
            #1 Load raw data (only necessary columns) and deinterleave
            rawdata_df = pd.read_csv(
                raw_data_path, 
                skiprows=1, 
                usecols=['Time(s)', 'AIn-1', 'DI/O-1', 'DI/O-2'], 
                encoding="ISO-8859-1"
            )
            
            #2 Deinterleave the data and save to CSV
            deinterleaved_df = pp.deinterleave(rawdata_df)
            deinterleaved_df.to_csv(deinterleaved_path, index=False)
            
            #3 Plot raw data and save as PNG
            fig_raw = gp.plot_rawdata(deinterleaved_df, exp, session, mouse)
            fig_raw.savefig(raw_plot_path)
            plt.close(fig_raw)

        
#%% 1.3 - Open artifacted data in Dash using plotly and enter artifact timestamps in excel file

#------------------#
session = 'Habituation'
mouse = 'A1m'
#------------------#
# in excel 'Filecode', put '{exp}_{session}_{mouse}'
code = gp.session_code(session,exp)
deinterleaved_df = pd.read_csv(pp_path/f'{mouse}_{code}_deinterleaved.csv')

app = Dash(__name__)
fig = px.line(deinterleaved_df[TIME_BEGIN:], x='Time(s)', y='405 Deinterleaved')
app.layout = html.Div([
    html.H4(f'{exp} {session} {mouse}'),
    dcc.Graph(
        id=f'{exp}',
        figure=fig
    )
])
if __name__ == '__main__':
    app.run_server(debug=False, use_reloader=False)


#%% 1.4 - Artifact correction and dFF calculation

#import artifacts boundaries
artifacts_df = pd.read_excel(experiment_path / 'artifacts.xlsx')
method = 'mean'

for session_path in [Path(f.path) for f in os.scandir(exp_path) if f.is_dir()]:
    session = str(session_path).split('\\')[-1]
    print('##########################################')
    print(f'EXPERIMENT : {exp} - SESSION : {session}')
    print('##########################################')
    code = gp.session_code(session,exp)
    for mouse in ['A2m']:#subjects_df['Subject']:
        if os.path.exists(pp_path/f'{mouse}_{code}_deinterleaved.csv'):
  #          try:
            print("--------------")
            print(f'MOUSE : {mouse}')
            print("--------------")
            deinterleaved_df = pd.read_csv(pp_path/f'{mouse}_{code}_deinterleaved.csv')
            filecode = f'{exp}_{session}_{mouse}'
            
            # calculate dFF with artifacts removal, then interpolate missing data
            dFFdata_df = pp.dFF(deinterleaved_df,artifacts_df,filecode,method)
            interpdFFdata_df = pp.interpolate_dFFdata(dFFdata_df, method='linear')
            #sometimes 1st timestamps=Nan instead of 0, raises an error
            interpdFFdata_df['Time(s)'] = interpdFFdata_df['Time(s)'].fillna(0) 
            
            # smooth data with butterworth filter or simple moving average (SMA)
            #smoothdFF_df=pp.smoothing_SMA(interpdFFdata_df,win_size=7)
            smoothdFF_df=pp.butterfilt(interpdFFdata_df, ORDER, CUT_FREQ)
            smoothdFF_df.to_csv(pp_path/f'{mouse}_{code}_dFFfilt.csv')
            
            #plotted GCaMP and isosbestic curves after dFF or fitting
            fig_dFF = gp.plot_fiberpho(smoothdFF_df,exp,session,mouse,method)
            fig_dFF.savefig(pp_path/f'{mouse}_{code}_{method}dFF.png')
            plt.close(fig_dFF) 
        #    except Exception as e:
         #       print(f'Problem in processing mouse {mouse} : {e}')
            