# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 13:39:29 2023

To run fiberphotometry analysis with behaviour or plethysmography data
1 - PREPROCESSING
2 - ANALYSIS WITH BEHAVIOUR BORIS FILE
3 - ANALYSIS WITH PLETHYSMOGRAPHY

@author: alice fermigier
"""

#%%
##########
#IMPORTED#
##########

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'figure.max_open_warning' : 0})
from pathlib import Path
import os
import sys
from dash import Dash, dcc, html
import plotly.express as px

#Custom
#put path to directory where python files are stored
if 'D:\Profil\Documents\GitHub\Fiberphotometry_analysis' not in sys.path:
    sys.path.append('D:\Profil\Documents\GitHub\Fiberphotometry_analysis')
if '/Users/alice/Documents/GitHub/Fiberphotometry_analysis' not in sys.path:
    sys.path.append('/Users/alice/Documents/GitHub/Fiberphotometry_analysis')

import preprocess as pp
import genplot as gp
import behavplot as bp
import plethyplot as plp
import statcalc as sc

#%%
########
#LOADER#
########

experiment_path = Path('D:\\Alice\\Fiber\\202301_CA2b5') #/Users/alice/Desktop/Fiberb5
analysis_path = experiment_path / 'Analysis'
data_path = experiment_path / 'Data'
os.chdir(experiment_path)
os.getcwd()

#import ID and groups of all mice
subjects_df = pd.read_excel(experiment_path / 'subjects.xlsx', sheet_name='Included')
#import tasks in protocol
proto_df = pd.read_excel(experiment_path / 'protocol.xlsx')

#time before behaviour for calculation of PETH baseline, in seconds
PRE_EVENT_TIME = 0
#time to crop at the beginning of the trial for , in seconds
TIME_BEGIN = 20
#threshold to fuse behaviour if bouts are too close, in secs
THRESH_S = 5
#threshold for PETH : if events are too short do not plot them and do not include them in PETH, in seconds
EVENT_TIME_THRESHOLD = 0
#filter characteristics
ORDER = 4
CUT_FREQ = 1 #in Hz
# ideal samplerate : deinterleaved data samplerate should be 10 but is 9.99 instead ==> change index to align with behaviour
SAMPLERATE = 10 #in Hz

########
#SCRIPT#
########

#%% 1 - PREPROCESSING
#####################

#------------------#
exp = 'Fear'
#------------------#

exp_path = analysis_path / exp
if not os.path.exists(exp_path):
    os.mkdir(exp_path)
for session_list in proto_df.loc[proto_df['Task']==exp,'Sessions'].values[0].split():
    if not os.path.exists(exp_path / session_list):
        os.mkdir(exp_path / session_list)
            
#get data path related to the task in protocol excel file
data_path_exp = data_path / proto_df.loc[proto_df['Task']==exp, 'Data_path'].values[0]

pp_path = data_path_exp / 'Preprocessing' # put all raw plots and deinterleaved_df in separate folder
if not os.path.exists(pp_path):
    os.mkdir(pp_path)

#%% 1.1 - Deinterleave data and save in separate file.

for session_path in [Path(f.path) for f in os.scandir(exp_path) if f.is_dir()]:#[exp_path / 'Test']:#
    session = str(session_path).split('\\')[-1]
    print('##########################################')
    print(f'EXPERIMENT : {exp} - SESSION : {session}')
    print('##########################################')
    code = gp.session_code(session,exp)
    for mouse in subjects_df['Subject']:
        print("--------------")
        print(f'MOUSE : {mouse}')
        print("--------------")
        if os.path.exists(data_path_exp/f'{mouse}_{code}.csv'):
            if not os.path.exists(pp_path/f'{mouse}_{code}_deinterleaved.csv'):
                rawdata_df = pd.read_csv(data_path_exp/f'{mouse}_{code}.csv',skiprows=1,usecols=['Time(s)','AIn-1','DI/O-1','DI/O-2'], encoding = "ISO-8859-1")
                deinterleaved_df = pp.deinterleave(rawdata_df)
                deinterleaved_df.to_csv(pp_path/f'{mouse}_{code}_deinterleaved.csv')
                
        # 1.2 - Look at rawdata and check for artifacts
        
                fig_raw = gp.plot_rawdata(deinterleaved_df,exp,session,mouse)
                fig_raw.savefig(pp_path/f'{mouse}_{code}_rawdata.png')
        
#%% 1.3 - Open artifacted data in Dash using plotly and enter artifact timestamps in excel file

#------------------#
session = 'Habituation'
mouse = 'A1m'
#------------------#
# in excel 'Filecode', put '{exp}_{session}_{mouse}'
code = gp.session_code(session,exp)
deinterleaved_df = pd.read_csv(pp_path/f'{mouse}_{code}_deinterleaved.csv')

app = Dash(__name__)
fig = px.line(deinterleaved_df[100:], x='Time(s)', y='405 Deinterleaved')
#fig = px.line(deinterleaved_df, x='Time(s)', y='405 Deinterleaved')
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
    for mouse in subjects_df['Subject']:
        if os.path.exists(pp_path/f'{mouse}_{code}_deinterleaved.csv'):
            print("--------------")
            print(f'MOUSE : {mouse}')
            print("--------------")
            deinterleaved_df = pd.read_csv(pp_path/f'{mouse}_{code}_deinterleaved.csv',index_col=0)
            filecode = f'{exp}_{session}_{mouse}'
            sr = pp.samplerate(deinterleaved_df)
            
            # calculate dFF with artifacts removal, then interpolate missing data
            dFFdata_df = pp.dFF(deinterleaved_df,artifacts_df,filecode,sr,method)
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

#%% 2 - ANALYSIS - BEHAVIOUR
############################

# 2.1 - Align with behaviour, create corresponding excel, plot fiberpho data with behaviour

for session_path in [Path(f.path) for f in os.scandir(exp_path) if f.is_dir()]:
    session = str(session_path).split('\\')[-1]
    print('##########################################')
    print(f'EXPERIMENT : {exp} - SESSION : {session}')
    print('##########################################')
    code = gp.session_code(session,exp)
    repo_path = session_path /  f'length{EVENT_TIME_THRESHOLD}_interbout{THRESH_S}_o{ORDER}f{CUT_FREQ}'
    if not os.path.exists(repo_path):
            os.mkdir(repo_path)
    for mouse in subjects_df['Subject']:
        print("--------------")
        print(f'MOUSE : {mouse}')
        print("--------------")
        rawdata_path = data_path_exp / f'{mouse}_{code}.csv'
        behav_path = data_path_exp / f'behav_{code}_{mouse}.csv'
        fiberpho_path = pp_path / f'{mouse}_{code}_dFFfilt.csv'
        
        #begin analysis only if behaviour has been scored and plots haven't been done
        ready = False
        if os.path.exists(behav_path):
            ready = True
        print(f'ready? {ready}')
        done = True
        if not os.path.exists(repo_path / f'{mouse}_{code}_fiberbehav.csv'):
            done = False
        print(f'done? {done}')
        
        if ready == True and done == False:
            rawdata_cam_df = pd.read_csv(rawdata_path, skiprows=1, usecols=['Time(s)','DI/O-3'], encoding = "ISO-8859-1")
            behav10Sps = pd.read_csv(behav_path)
            fiberpho = pd.read_csv(fiberpho_path,index_col=0)
            sr=pp.samplerate(fiberpho)
            print(exp, session, mouse)
            
            #list of behaviours to analyze
            list_BOI = behav10Sps.columns[1:].tolist()
            print(f'Behaviours = {list_BOI}')
                
            #align behaviour and fiberpho data, create fiberbehav.csv
            timevector = gp.time_vector(fiberpho, sr)
            timestart_camera = gp.timestamp_camera(rawdata_cam_df)[0]
            print(f'start camera : {timestart_camera}')
            print(f'len timevector : {len(timevector)}')
            fiberbehav_df = bp.align_behav(behav10Sps, fiberpho, timevector, timestart_camera, exp)
            fiberbehav_df = bp.behav_process(fiberbehav_df, list_BOI, THRESH_S, EVENT_TIME_THRESHOLD, sr)
            fiberbehav_df.to_csv(repo_path / f'{mouse}_{code}_fiberbehavnotderived.csv')
            dfiberbehav_df = bp.derive(fiberbehav_df)
            dfiberbehav_df.to_csv(repo_path / f'{mouse}_{code}_fiberbehav.csv')
            
            #plot fiberpho data with behaviour
            fig_fiberbehav = bp.plot_fiberpho_behav(dfiberbehav_df,list_BOI,exp,session,mouse,THRESH_S,EVENT_TIME_THRESHOLD)
            fig_fiberbehav.savefig(repo_path/f'{mouse}_{code}_fiberbehav.pdf')
            fig_fiberbehav.savefig(repo_path/f'{mouse}_{code}_fiberbehav.png')

#%% 2.2 - Plot PETH for each mouse and calculate mean and max z-scored dFF before and after behaviour
#the beginning of the PETH will be plotted at the beginning of an ascending curve

#PETH parameters 
list_EVENT = ['onset'] #, 'withdrawal']
list_TIMEWINDOW = [[5,30]] #,[4,15]] 

for session_path in [Path(f.path) for f in os.scandir(exp_path) if f.is_dir()]:
    session = str(session_path).split('\\')[-1]
    print('##########################################')
    print(f'EXPERIMENT : {exp} - SESSION : {session}')
    print('##########################################')
    code = gp.session_code(session,exp)
    repo_path = session_path /  f'length{EVENT_TIME_THRESHOLD}_interbout{THRESH_S}_o{ORDER}f{CUT_FREQ}'
    PETH_path = repo_path / 'PETH'
    if not os.path.exists(PETH_path):
            os.mkdir(PETH_path)
    for mouse in subjects_df['Subject']:
        group = subjects_df.loc[subjects_df['Subject']==mouse, 'Group'].values[0]
        if os.path.exists(repo_path /  f'{mouse}_{code}_fiberbehav.csv'):
            print("--------------")
            print(f'MOUSE : {mouse}')
            print("--------------")
            dfiberbehav_df = pd.read_csv(repo_path / f'{mouse}_{code}_fiberbehav.csv',index_col=0)
            list_BOI = dfiberbehav_df.columns[4:].tolist()
            for BOI in list_BOI:
                if BOI not in ['Entry in arena','Gate opens','Tail suspension']:
                    for event, timewindow in zip(list_EVENT,list_TIMEWINDOW): #CAUTION : adapt function!!
                        PETH_data = bp.PETH(dfiberbehav_df, BOI, event, timewindow, EVENT_TIME_THRESHOLD, PRE_EVENT_TIME)
                        PETH_df = pd.DataFrame(np.transpose(PETH_data),index=np.arange(-timewindow[0], timewindow[1]+0.1, 0.1))
                        PETH_df.to_excel(PETH_path / f'{mouse}_{code}_{BOI}{event[0]}{timewindow}_PETHdata.xlsx')
                        fig_PETH = bp.plot_PETH(PETH_data, BOI, event, timewindow, exp, session, mouse, group)
                        fig_PETH.savefig(PETH_path /  f'{mouse}_{code}_{BOI}{event[0]}{timewindow[0]-timewindow[1]}_PETH.png')
                        
#%% 2.3 - Calculate mean, max and diff within behaviours (for state behaviours)

for session_path in [Path(f.path) for f in os.scandir(exp_path) if f.is_dir()]:
    session = str(session_path).split('\\')[-1]
    print('##########################################')
    print(f'EXPERIMENT : {exp} - SESSION : {session}')
    print('##########################################')
    code = gp.session_code(session,exp)
    repo_path = session_path /  f'length{EVENT_TIME_THRESHOLD}_interbout{THRESH_S}_o{ORDER}f{CUT_FREQ}'
    groupanalysis_path = repo_path / 'Group_analysis'
    if not os.path.exists(groupanalysis_path):
        os.mkdir(groupanalysis_path)
    #create list of mean_dFFs and diff_dFFs
    mean_dFFs_list = []
    subject_list = []
    diffmeanmaxdFF_list = []
    group_list = []
    PETH_array = None
    for mouse in subjects_df['Subject']:
        if os.path.exists(repo_path /  f'{mouse}_{code}_fiberbehav.csv'):
            print("--------------")
            print(f'MOUSE : {mouse}')
            print("--------------")
            subject_list.append(mouse)
            group = subjects_df.loc[subjects_df['Subject']==mouse, 'Group'].values[0]
            dfiberbehav_df = pd.read_csv(repo_path / f'{mouse}_{code}_fiberbehav.csv',index_col=0)
            fiberbehav_df = pd.read_csv(repo_path / f'{mouse}_{code}_fiberbehavnotderived.csv',index_col=0)
            #1 - subjects
            subject_list.append(mouse)
            group_list.append(group)
            #2 - mean, max and delta dFF
            mean_dFFs_list.append(sc.meandFF_behav(list_BOI, dfiberbehav_df, exp, session, mouse, group))
            diffmeanmaxdFF_list.append(sc.diffmeanmax_dFF(fiberbehav_df, list_BOI, mouse, group))

    #export to excel
    meandFFs_allmice = pd.concat(mean_dFFs_list)
    meandFFs_allmice.to_excel(groupanalysis_path / f'{exp}_{session}_globmeandFFs.xlsx')
    
    diffdFFs_allmice = pd.concat(diffmeanmaxdFF_list)
    diffdFFs_allmice.to_excel(groupanalysis_path / f'{exp}_{session}_diffmaxmeandFFs.xlsx')
    
#%% 2.5 - Calculate mean and max before and after behaviour onset for whole group

#------------------------------#
BOI = 'Shock'
TIME_MEANMAX = 5 #in seconds
TIME_BEFORE = 1 #in seconds
maxboutsnumber = None
#------------------------------#

for session_path in [Path(f.path) for f in os.scandir(exp_path) if f.is_dir()]:
    session = str(session_path).split('\\')[-1]
    print('##########################################')
    print(f'EXPERIMENT : {exp} - SESSION : {session}')
    print('##########################################')
    code = gp.session_code(session,exp)
    repo_path = session_path /  f'length{EVENT_TIME_THRESHOLD}_interbout{THRESH_S}_o{ORDER}f{CUT_FREQ}'
    groupanalysis_path = repo_path / 'Group_analysis'
    if len(os.listdir(repo_path)) > 1:
        if not os.path.exists(groupanalysis_path): 
                os.mkdir(groupanalysis_path)
        #create list of mean_dFFs, max_dFFs
        meandFF_list = []
        meandFF_before_list = []
        subject_list = []
        group_list = []
        maxdFF_list = []
        maxdFF_before_list = []
        for mouse in subjects_df['Subject']:
            print("--------------")
            print(f'MOUSE : {mouse}')
            print("--------------")
            if os.path.exists(repo_path /  f'{mouse}_{code}_fiberbehav.csv'):
                fiberbehav_df = pd.read_csv(repo_path /  f'{mouse}_{code}_fiberbehav.csv',index_col=0)
                group = subjects_df.loc[subjects_df['Subject']==mouse, 'Group'].values[0]
                sr = pp.samplerate(fiberbehav_df)
                
                if BOI in fiberbehav_df.columns[2:].tolist():
                    ind_event_list = np.where(fiberbehav_df[BOI] == 1)[0].tolist()
                    if maxboutsnumber != None:
                        ind_event_list = ind_event_list[:maxboutsnumber] #ne mesurer que pour le nombre de bouts désiré
                    for ind_event in ind_event_list:
                        #1 - subjects 
                        subject_list.append(mouse)
                        group_list.append(group)
                        #2 - mean and max dFF
                        ind_event = fiberbehav_df.loc[ind_event-TIME_BEFORE*sr:ind_event, 'Denoised dFF'].idxmin()
                        meandFF_list.append(fiberbehav_df.loc[ind_event:ind_event+TIME_MEANMAX*sr, 'Denoised dFF'].mean()) #mean after entry
                        meandFF_before_list.append(fiberbehav_df.loc[ind_event-TIME_MEANMAX*sr:ind_event, 'Denoised dFF'].mean()) #mean before entry
                        maxdFF_list.append(fiberbehav_df.loc[ind_event:ind_event+TIME_MEANMAX*sr, 'Denoised dFF'].max()) #max after entry
                        maxdFF_before_list.append(fiberbehav_df.loc[ind_event-TIME_MEANMAX*sr:ind_event, 'Denoised dFF'].max()) #max before entry
            
        #export data to excel
        if not os.path.exists(groupanalysis_path / f'{exp}_{session}_{TIME_MEANMAX}s_maxmeandFFentry.xlsx'):
            meanmaxdFFs_df = pd.DataFrame(data={'Subject' : subject_list, 'Group' : group_list, 
                                                'Mean dFF before entry' : meandFF_before_list, 'Mean dFF entry' : meandFF_list, 
                                                'Max dFF before entry' : maxdFF_before_list, 'Max dFF entry' : maxdFF_list})
            meanmaxdFFs_df.to_excel(groupanalysis_path / f'{exp}_{session}_{TIME_MEANMAX}s_maxmeandFFentry.xlsx')
        
        
#%%Plot PETH for each group
#------------------------------#
BOI = 'Shock'
timewindow = [5,30]
maxboutsnumber = 1
#------------------------------#

for session_path in [Path(f.path) for f in os.scandir(exp_path) if f.is_dir()]:
    session = str(session_path).split('\\')[-1]
    print('##########################################')
    print(f'EXPERIMENT : {exp} - SESSION : {session}')
    print('##########################################')
    code = gp.session_code(session,exp)
    repo_path = session_path /  f'length{EVENT_TIME_THRESHOLD}_interbout{THRESH_S}_o{ORDER}f{CUT_FREQ}'
    if len(os.listdir(repo_path)) > 1:
        if not os.path.exists(PETH_path): 
                os.mkdir(PETH_path)
        #create list of mean_dFFs, max_dFFs
        subject_list = []
        group_list = []
        PETH_array = None
        PETH_mean_list = []
        PETH_max_list = []
        for mouse in subjects_df['Subject']:
            print("--------------")
            print(f'MOUSE : {mouse}')
            print("--------------")
            if os.path.exists(repo_path /  f'{mouse}_{code}_fiberbehav.csv'):
                fiberbehav_df = pd.read_csv(repo_path /  f'{mouse}_{code}_fiberbehav.csv',index_col=0)
                group = subjects_df.loc[subjects_df['Subject']==mouse, 'Group'].values[0]
                sr = pp.samplerate(fiberbehav_df)
                
                if BOI in fiberbehav_df.columns[2:].tolist():
                    subject_list.append(mouse)
                    group_list.append(group)
                    print(f'PETH {BOI}')
                    PETH_mouse = bp.PETH(fiberbehav_df, BOI, 'onset', timewindow, EVENT_TIME_THRESHOLD, sr, PRE_EVENT_TIME, maxboutsnumber)
                    if PETH_array is None:
                        PETH_array = PETH_mouse
                        print('initialize successful')
                    else:
                        PETH_array = np.concatenate((PETH_array,PETH_mouse))
                        print('stacking successful')
                    #calculate mean and max a z-scored dFF
                    PETH_mean_list.append(np.mean(PETH_mouse[:timewindow[0]]),np.mean(PETH_mouse[timewindow[0]:]))
                    PETH_max_list.append(np.max(PETH_mouse[:timewindow[0]]),np.max(PETH_mouse[timewindow[0]:]))
                    
        meanmaxPETH_df = pd.DataFrame(data={'Subject' : subject_list, 'Group' : group_list, 
                                                        f'Mean dFF before {BOI}' : list(list(zip(*PETH_mean_list))[0]), f'Mean dFF {BOI}' : list(list(zip(*PETH_mean_list))[1]), 
                                                        f'Max dFF before {BOI}' : list(list(zip(*PETH_max_list))[0]), f'Max dFF {BOI}' : list(list(zip(*PETH_max_list))[1])})
        meanmaxPETH_df.to_excel(PETH_path / f'{BOI}_{timewindow[0]}{timewindow[1]}_PETHmeanmax.xlsx')
                    
        included_groups = ['CD','HFD']
        for group in included_groups:
            PETH_array_group = PETH_array
            list_todelete = []
            for (i,group_mouse) in enumerate(group_list): 
                if group not in group_mouse:
                    list_todelete.append(i)
            PETH_array_group = np.delete(PETH_array_group,(list_todelete),axis=0)
            print(len(PETH_array_group))
            fig_PETHpooled = bp.plot_PETH_pooled(PETH_array_group, BOI, 'onset', timewindow, exp, session, group)
            
            fig_PETHpooled.savefig(PETH_path / f'{group}{BOI}{timewindow[1]}_PETH.pdf')
            fig_PETHpooled.savefig(PETH_path / f'{group}{BOI}{timewindow[1]}_PETH.png')
#%% 3 - ANALYSIS - PLETHYSMOGRAPH
#################################

# 3.1 - Visualize and score manually sniffs and stims in excel file

#------------------#
session = 'Test'
mouse = 'B2f'
#------------------#

code = gp.session_code(session,exp)
rawdata_path = data_path_exp / f'{mouse}_{code}.csv'
plethys_df = pd.read_csv(rawdata_path, skiprows=1, usecols=['Time(s)','AIn-4'])

app = Dash(__name__)
data_df = plethys_df.loc[[i for i in range(0,len(plethys_df),600)]] #downsample plethysmo data
fig = px.line(data_df, x='Time(s)', y='AIn-4')
app.layout = html.Div([
    html.H4(f'{exp} {session} {mouse}'),
    dcc.Graph(
        id=f'{exp}',
        figure=fig
    )
])
if __name__ == '__main__':
    app.run_server(debug=False, use_reloader=False)

#%% 3.2 - Align with sniffs, create corresponding csv, plot fiberpho data with sniffs and stims

sniffs_df = pd.read_excel(data_path_exp / 'Sniffs.xlsx')
session_path = exp_path / 'Test'
session = str(session_path).split('\\')[-1]

print('##########################################')
print(f'EXPERIMENT : {exp} - SESSION : {session}')
print('##########################################')
code = gp.session_code(session,exp)
repo_path = session_path /  f'length{EVENT_TIME_THRESHOLD}_interbout{THRESH_S}_o{ORDER}f{CUT_FREQ}'
raw_path = repo_path / 'Raw'
if not os.path.exists(repo_path):
    os.mkdir(repo_path)
if not os.path.exists(raw_path):
    os.mkdir(raw_path)
for mouse in subjects_df['Subject']:
    print("--------------")
    print(f'MOUSE : {mouse}')
    print("--------------")
    if mouse in set(sniffs_df['Subject']):
        fiberpho_df = pd.read_csv(pp_path / f'{mouse}_{code}_dFFfilt.csv',index_col=0)
        if not (repo_path / f'{mouse}_{code}_fibersniff.csv').is_file():
            srows=1  
            rawdata_path = data_path_exp / f'{mouse}_{code}.csv'
            sr = pp.samplerate(fiberpho_df)
            plethys_df = pd.read_csv(rawdata_path, skiprows=srows, usecols=['Time(s)','AIn-4'])
            fibersniff_df = plp.align_sniffs(fiberpho_df, plethys_df, sniffs_df, sr, mouse)
            fibersniff_df = plp.process_fibersniff(fibersniff_df, EVENT_TIME_THRESHOLD, THRESH_S, sr)
            dfibersniff_df = plp.derive(fibersniff_df)
            fibersniff_df.to_csv(repo_path / f'{mouse}_{code}_fibersniffnotderived.csv')
            dfibersniff_df.to_csv(repo_path / f'{mouse}_{code}_fibersniff.csv')
        #plot figures
        if not (raw_path / f'{mouse}_WBPfiberpho_raw.png').is_file():
            fig_raw = plp.plethyfiber_plot_raw(fiberpho_df, plethys_df,mouse)
            fig_raw.savefig(raw_path / f'{mouse}_WBPfiberpho_raw.png') 
        if not (repo_path / f'{mouse}_WBPfiberpho_sniffs.pdf').is_file():
            fig_sniffs = plp.plethyfiber_plot_sniffs(dfibersniff_df,sniffs_df,mouse)
            fig_sniffs.savefig(repo_path / f'{mouse}_WBPfiberpho_sniffs.png') 
            fig_sniffs.savefig(repo_path / f'{mouse}_WBPfiberpho_sniffs.pdf')
        else: print('Done')
            
                
                
#%% 3.3 - Plot PETH for each mouse, sniffs and stim   

sniffs_df = pd.read_excel(data_path_exp / 'Sniffs.xlsx')
session_path = exp_path / 'Test'
session = str(session_path).split('\\')[-1]           
                
print('##########################################')
print(f'EXPERIMENT : {exp} - SESSION : {session}')
print('##########################################')
code = gp.session_code(session,exp)
repo_path = session_path /  f'length{EVENT_TIME_THRESHOLD}_interbout{THRESH_S}_o{ORDER}f{CUT_FREQ}'
PETH_path = repo_path / 'PETH'
if not os.path.exists(repo_path):
        os.mkdir(repo_path)
if not os.path.exists(PETH_path):
        os.mkdir(PETH_path)
for mouse in subjects_df['Subject']:
    print("--------------")
    print(f'MOUSE : {mouse}')
    print("--------------")
    if mouse in set(sniffs_df['Subject']):
        dfibersniff_df = pd.read_csv(repo_path / f'{mouse}_{code}_fibersniff.csv',index_col=0)
        sr = pp.samplerate(dfibersniff_df)
        for odor in set(sniffs_df['Odor']):
            print(odor)
            for (event, timewindow) in zip(list_EVENT, list_TIMEWINDOW):
                print(event)
                if not (PETH_path / f'{mouse}{odor}_PETH{event[0]}.pdf').is_file():
                    #plot PETH for sniff
                    PETHsniff_data = plp.PETH_sniff(dfibersniff_df, odor, 0, event, timewindow, mouse, sr, PRE_EVENT_TIME)
                    fig_PETHsniff = plp.plot_PETH(PETHsniff_data, odor, event, timewindow, 'Sniff', sr, mouse)
                    fig_PETHsniff.savefig(PETH_path / f'{mouse}{odor}_PETHsniff{event[0]}.png')
                    fig_PETHsniff.savefig(PETH_path / f'{mouse}{odor}_PETHsniff{event[0]}.pdf')
                    #plot PETH for stim
                    PETHstim_data = plp.PETH_stim(dfibersniff_df, odor, 0, event, timewindow, mouse, sr, PRE_EVENT_TIME)
                    fig_PETHstim = plp.plot_PETH(PETHstim_data, odor, event, timewindow, 'Stim', sr, mouse)
                    fig_PETHstim.savefig(PETH_path / f'{mouse}{odor}_PETHstim{event[0]}.png')
                    fig_PETHstim.savefig(PETH_path / f'{mouse}{odor}_PETHstim{event[0]}.pdf')
        plt.close('all')
        
#%% 3.4.1 - Calculate mean and max before and during stims for whole group

#------------------------------#
TIME_MEANMAX = 15 #in seconds
#------------------------------#

print('##########################################')
print(f'EXPERIMENT : {exp} - SESSION : {session}')
print('##########################################')
code = gp.session_code(session,exp)
repo_path = session_path /  f'length{EVENT_TIME_THRESHOLD}_interbout{THRESH_S}_o{ORDER}f{CUT_FREQ}'
groupanalysis_path = repo_path / 'Group analysis'
if not os.path.exists(groupanalysis_path):
        os.mkdir(groupanalysis_path)
        
meanmaxdFF_list = []

for mouse in subjects_df['Subject']:
    print("--------------")
    print(f'MOUSE : {mouse}')
    print("--------------")
    if mouse in set(sniffs_df['Subject']):
        group = subjects_df.loc[subjects_df['Subject']==mouse, 'Group'].values[0]
        dfibersniff_df = pd.read_csv(repo_path / f'{mouse}_{code}_fibersniff.csv',index_col=0)
        sr = pp.samplerate(dfibersniff_df)
        list_BOI = [i for i in dfibersniff_df.columns[3:] if 'Stim' in i]

        meanmaxdFF_list.append(sc.meanmax_dFF_stims(dfibersniff_df, list_BOI, mouse, group, TIME_MEANMAX, sr))
        
#export to excel
meanmaxdFF_allmice = pd.concat(meanmaxdFF_list)
meanmaxdFF_allmice.to_excel(groupanalysis_path / f'{exp}_{session}_{TIME_MEANMAX}s_maxmeandFFstims.xlsx')

        
#%% 3.4.2 - Calculate mean, max and diff dFF for sniffs

#note : ne donne rien avec sniffs pour batch 5

print('##########################################')
print(f'EXPERIMENT : {exp} - SESSION : {session}')
print('##########################################')
code = gp.session_code(session,exp)
repo_path = session_path /  f'length{EVENT_TIME_THRESHOLD}_interbout{THRESH_S}_o{ORDER}f{CUT_FREQ}'
groupanalysis_path = repo_path / 'Group analysis'
if not os.path.exists(groupanalysis_path):
        os.mkdir(groupanalysis_path)

#create list of mean_dFFs and diff_dFFs
mean_dFFs_list = []
diffmeanmaxdFF_list = []

for mouse in subjects_df['Subject']:
    print("--------------")
    print(f'MOUSE : {mouse}')
    print("--------------")
    if mouse in set(sniffs_df['Subject']):
        group = subjects_df.loc[subjects_df['Subject']==mouse, 'Group'].values[0]
        dfibersniff_df = pd.read_csv(repo_path / f'{mouse}_{code}_fibersniff.csv',index_col=0)
        fibersniff_df = pd.read_csv(repo_path / f'{mouse}_{code}_fibersniffnotderived.csv',index_col=0)
        sr = pp.samplerate(dfibersniff_df)
        list_BOI = dfibersniff_df.columns[3:]
        
        diffmeanmaxdFF_list.append(sc.diffmeanmax_dFF(dfibersniff_df, list_BOI, mouse, group))
        mean_dFFs_list.append(sc.meandFF_sniffs(fibersniff_df, exp, session, mouse, group, joined=True))
        
#export to excel
diffdFFs_allmice = pd.concat(diffmeanmaxdFF_list)
diffdFFs_allmice.to_excel(groupanalysis_path / f'{exp}_{session}_diffmaxmeandFFs.xlsx')

#export to excel
meandFFs_allmice = pd.concat(mean_dFFs_list)
meandFFs_allmice.to_excel(groupanalysis_path / f'{exp}_{session}_globmeandFFs.xlsx')
    
#%% 3.5 - Plot PETH for all mice (each odor and all counts)

#------------------------------#
#for PETH, groups that will be plotted:
included_groups = ['CD','HFD']
event = 'onset' #or 'withdrawal'
timewindow = [6,15]
#------------------------------#
                
print('##########################################')
print(f'EXPERIMENT : {exp} - SESSION : {session}')
print('##########################################')
code = gp.session_code(session,exp)
repo_path = session_path /  f'length{EVENT_TIME_THRESHOLD}_interbout{THRESH_S}_o{ORDER}f{CUT_FREQ}'
groupanalysis_path = repo_path / 'Group analysis'
for type_event in ['Stim']:
    for odor in set(sniffs_df['Odor']):
        print(odor)
        group_list = []
        PETH_list = []
        for mouse in subjects_df['Subject']:
            if mouse in set(sniffs_df['Subject']):
                print("--------------")
                print(f'MOUSE : {mouse}')
                print("--------------")
                group = subjects_df.loc[subjects_df['Subject']==mouse, 'Group'].values[0]
                group_list.append(group)
                dfibersniff_df = pd.read_csv(repo_path / f'{mouse}_{code}_fibersniff.csv',index_col=0)
                sr = pp.samplerate(dfibersniff_df)
                #PETH data
                if type_event == 'Stim':
                    PETH_list.append(np.mean(plp.PETH_stim(dfibersniff_df, odor, 0, event, timewindow, mouse, sr, PRE_EVENT_TIME), axis=0))
                elif type_event == 'Sniff':
                    PETH_list.append(np.mean(plp.PETH_sniff(dfibersniff_df, odor, 0, event, timewindow, mouse, sr, PRE_EVENT_TIME), axis=0))

        
             
        #plot PETH
        PETH_array = np.stack(PETH_list)
        BOI = f'{type_event} {odor}'
        PETHarray_list=[]
        for group in included_groups:
            PETH_array_group = PETH_array
            list_todelete = []
            for (i,group_mouse) in enumerate(group_list): 
                if group not in group_mouse:
                    list_todelete.append(i)
            PETHarray_list.append(np.delete(PETH_array_group,(list_todelete),axis=0))
           # print(PETHarray_list,len(PETHarray_list))
        fig_PETHpooled = bp.plot_PETH_pooled(included_groups, PETHarray_list, BOI, event, timewindow, exp, session) 
        fig_PETHpooled.savefig(groupanalysis_path / f'{included_groups[0]}{included_groups[1]}{type_event}{odor}{event[0]}{timewindow[1]}_PETH.pdf')
        fig_PETHpooled.savefig(groupanalysis_path / f'{included_groups[0]}{included_groups[1]}{type_event}{odor}{event[0]}{timewindow[1]}_PETH.png')

                
#%% 3.6 - Plot PETH for all mice (each odor and each count)

#------------------------------#
#for PETH, groups that will be plotted:
included_groups = ['CD','HFD']
event = 'onset' #or 'withdrawal'
timewindow = [6,15]
#------------------------------#
                
print('##########################################')
print(f'EXPERIMENT : {exp} - SESSION : {session}')
print('##########################################')
code = gp.session_code(session,exp)
repo_path = session_path /  f'length{EVENT_TIME_THRESHOLD}_interbout{THRESH_S}_o{ORDER}f{CUT_FREQ}'
groupanalysis_path = repo_path / 'Group analysis'
for type_event in ['Stim']:
    for odor in set(sniffs_df['Odor']):
        for count in set(sniffs_df.loc[sniffs_df['Odor']==odor,'Count'].values):
            count+=(-1)
            group_list = []
            PETH_array = None
            for mouse in subjects_df['Subject']:
                if mouse in set(sniffs_df['Subject']):
                    print("--------------")
                    print(f'MOUSE : {mouse}')
                    print("--------------")
                    group = subjects_df.loc[subjects_df['Subject']==mouse, 'Group'].values[0]
                    group_list.append(group)
                    dfibersniff_df = pd.read_csv(repo_path / f'{mouse}_{code}_fibersniff.csv')
                    sr = pp.samplerate(dfibersniff_df)
                    #PETH data
                    if type_event == 'Stim':
                        PETH_list.append(np.mean(plp.PETH_stim(dfibersniff_df, odor, count, event, timewindow, mouse, sr, PRE_EVENT_TIME), axis=0))
                    elif type_event == 'Sniff':
                        PETH_list.append(np.mean(plp.PETH_sniff(dfibersniff_df, odor, count, event, timewindow, mouse, sr, PRE_EVENT_TIME), axis=0))
                        
            #plot PETH
            PETH_array = np.stack(PETH_list)
            BOI = f'{type_event} {odor} {count}'
            PETHarray_list=[]
            for group in included_groups:
                PETH_array_group = PETH_array
                list_todelete = []
                for (i,group_mouse) in enumerate(group_list): 
                    if group not in group_mouse:
                        list_todelete.append(i)
                PETHarray_list.append(np.delete(PETH_array_group,(list_todelete),axis=0))
            fig_PETHpooled = bp.plot_PETH_pooled(included_groups, PETHarray_list, BOI, event, timewindow, exp, session) 
            fig_PETHpooled.savefig(groupanalysis_path / f'{included_groups[0]}{included_groups[1]}{type_event}{odor}{count}{event[0]}{timewindow[1]}_PETH.pdf')
            fig_PETHpooled.savefig(groupanalysis_path / f'{included_groups[0]}{included_groups[1]}{type_event}{odor}{count}{event[0]}{timewindow[1]}_PETH.png')
