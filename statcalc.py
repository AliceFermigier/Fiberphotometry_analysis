# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 17:16:16 2023

@author: alice fermigier
"""

#%%
##########
#IMPORTED#
##########

import pandas as pd
import numpy as np
import math

#%%
###################
#DEFINED FUNCTIONS#
###################

def meandFF_sniffs(fiberbehav_df, exp, session, mouse, group, joined=True):
    """
    Calculates mean dFF during each behaviour
    Output : dataframe with 'Subject','Group','Mean dFF','Baseline', 'Post_baseline', and mean dFF for each behaviour
    """
    list_meandFF = []
    
    fiberbehavsnip_df = fiberbehav_df[fiberbehav_df['Time(s)'] > 15]
    
    if joined == True:
        #give same name to stims and sniffs with same odors
        for col in fiberbehavsnip_df.columns[3:]:
            fiberbehavsnip_df.rename(columns={col: col[:-2]}, inplace=True)
        fiberbehavsnip_df = fiberbehavsnip_df.groupby(level=0, axis=1, sort=False).sum()
        print(fiberbehavsnip_df.columns)
        #put mean dFF in list
        for behav in fiberbehavsnip_df.columns[3:]:
            meandFF_behav_df = fiberbehavsnip_df.groupby([behav], as_index=False).mean()
            print(meandFF_behav_df)
            list_meandFF.append(meandFF_behav_df.loc[meandFF_behav_df[behav]==1, 'Denoised dFF'].values[0])
                                
    if joined == False:
        #put mean dFF in list
        for behav in fiberbehavsnip_df.columns[3:]:
            meandFF_behav_df = fiberbehavsnip_df.groupby([behav], as_index=False).mean()
            print(meandFF_behav_df)
            list_meandFF.append(meandFF_behav_df.loc[meandFF_behav_df[behav]==1, 'Denoised dFF'].values[0])
    
    #create dataframe with values
    
    list_dFFs = [mouse, group]
    list_columns = ['Subject','Group']
    for (behav,meandFF) in zip(fiberbehavsnip_df.columns[3:],list_meandFF):
        list_dFFs.append(meandFF)
        list_columns.append(behav)
    meandFFs_df = pd.DataFrame(data=[list_dFFs], columns=list_columns)
    
    return(meandFFs_df)

def diffmeanmax_dFF_perbout(behavprocess_df, list_BOI, mouse, group):
    """
    Calculates dFF difference between beginning and end of bouts
    Also calculates mean and max dFF during each bout
    """
    list_behav_analyzed = []
    for behav in list_BOI:
        if behav not in ['Entry in arena','Gate opens','Shock']:
            list_behav_analyzed.append(behav)
    list_deltadFF = []
    list_meandFF = []
    list_maxdFF = []
    listind_behav = []
    listind_bouts = []
    
    for behav in list_behav_analyzed:
        list_starts = np.where(behavprocess_df[behav]==1)[0].tolist()
        list_stops = np.where(behavprocess_df[behav]==-1)[0].tolist()
        bout_n = 1
        for (start, stop) in zip(list_starts, list_stops):
            mean_start = behavprocess_df.loc[start-5:start+5, 'Denoised dFF'].mean()
            mean_stop = behavprocess_df.loc[stop-5:stop+5, 'Denoised dFF'].mean()
            delta = mean_stop-mean_start
            mean = behavprocess_df.loc[start:stop, 'Denoised dFF'].mean()
            maximum = behavprocess_df.loc[start:stop, 'Denoised dFF'].max()
            list_deltadFF.append(delta)
            list_meandFF.append(mean)
            list_maxdFF.append(maximum)
            listind_behav.append(behav)
            listind_bouts.append(bout_n)
            bout_n+=1
    
    diffdFF_df = pd.DataFrame(data = {'Subject':[mouse]*len(listind_behav), 'Group':[group]*len(listind_behav),
                                      'Behaviour':listind_behav, 'Bout':listind_bouts, 'Mean dFF':list_meandFF,
                                      'Max dFF' : list_maxdFF, 'Delta dFF':list_deltadFF})
    
    return(diffdFF_df)

def diffmeanmax_dFF(behavprocess_df, list_BOI, mouse, group):
    """
    Calculates dFF difference between beginning and end of bouts
    Also calculates mean dFF during each bout
    """
    list_behav_analyzed = []
    for behav in list_BOI:
        if behav not in ['Entry in arena','Gate opens','Shock']:
            list_behav_analyzed.append(behav)
    list_deltadFF = []
    list_meandFF = []
    list_maxdFF = []
    listind_behav = []
    
    for behav in list_behav_analyzed:
        list_starts = np.where(behavprocess_df[behav]==1)[0].tolist()
        list_stops = np.where(behavprocess_df[behav]==-1)[0].tolist()
        listind_behav.append(behav)
        list_deltadFF_behav = []
        list_meandFF_behav = []
        list_maxdFF_behav = []
        for (start, stop) in zip(list_starts, list_stops):
            mean_start = behavprocess_df.loc[start-5:start+5, 'Denoised dFF'].mean()
            mean_stop = behavprocess_df.loc[stop-5:stop+5, 'Denoised dFF'].mean()
            delta = mean_stop-mean_start
            mean = behavprocess_df.loc[start:stop, 'Denoised dFF'].mean()
            maximum = behavprocess_df.loc[start:stop, 'Denoised dFF'].max()
            list_deltadFF_behav.append(delta)
            list_meandFF_behav.append(mean)
            list_maxdFF_behav.append(maximum)
            
        list_deltadFF.append(np.mean(list_deltadFF_behav))
        list_meandFF.append(np.mean(list_meandFF_behav))
        list_maxdFF.append(np.mean(list_maxdFF_behav))
            
    
    diffdFF_df = pd.DataFrame(data = {'Subject':[mouse]*len(listind_behav), 'Group':[group]*len(listind_behav),
                                      'Behaviour':listind_behav, 'Mean dFF':list_meandFF,
                                      'Max dFF':list_maxdFF, 'Delta dFF':list_deltadFF})
    
    return(diffdFF_df)

def meanmax_dFF_stims(behavprocess_df, list_BOI, mouse, group, TIME_MEANMAX, sr):
    """
    Calculates dFF difference between beginning and end of bouts
    Also calculates mean dFF during each bout
    """
    list_behav_analyzed = []
    for behav in list_BOI:
        if behav not in ['Entry in arena','Gate opens','Shock']:
            list_behav_analyzed.append(behav)
    list_meandFF_before = []
    list_meandFF_after = []
    list_maxdFF_before = []
    list_maxdFF_after = []
    listind_behav = []
    
    for behav in list_behav_analyzed:
        start = np.where(behavprocess_df[behav]==1)[0][0]
        listind_behav.append(behav)
        
        mean_before = behavprocess_df.loc[start-TIME_MEANMAX*sr:start, 'Denoised dFF'].mean()
        maximum_before = behavprocess_df.loc[start-TIME_MEANMAX*sr:start, 'Denoised dFF'].max()
        mean_after = behavprocess_df.loc[start:start+TIME_MEANMAX*sr, 'Denoised dFF'].mean()
        maximum_after = behavprocess_df.loc[start:start+TIME_MEANMAX*sr, 'Denoised dFF'].max()
        list_meandFF_before.append(mean_before)
        list_maxdFF_before.append(maximum_before)
        list_meandFF_after.append(mean_after)
        list_maxdFF_after.append(maximum_after)
            
    
    meanmaxdFF_df = pd.DataFrame(data = {'Subject':[mouse]*len(listind_behav), 'Group':[group]*len(listind_behav),
                                      'Behaviour':listind_behav, 'Mean dFF before':list_meandFF_before,
                                      'Mean dFF during':list_meandFF_after, 
                                      'Max dFF before':list_maxdFF_before, 'Max dFF during':list_maxdFF_after})
    
    return(meanmaxdFF_df)

#%%
#####
#OLD#
#####

def meandFF_behav(list_BOI, fiberbehav_df, exp, session, mouse, group):
    """
    Calculates mean dFF during each behaviour
    Output : dataframe with 'Subject','Group','Mean dFF','Baseline', 'Post_baseline', and mean dFF for each behaviour
    """
    list_behav_analyzed = []
    list_meandFF = []
    
    fiberbehavsnip_df = fiberbehav_df[fiberbehav_df['Time(s)'] > 15]
    
    #get index of when the gate opens
    if exp == 'NewContext':
        ind_start_trial = fiberbehavsnip_df.index[fiberbehavsnip_df['New context'] == 1].tolist()[0]
    elif exp == 'Fear' and session == 'Test':
        ind_start_trial = fiberbehavsnip_df.index[fiberbehavsnip_df['Fear cage'] == 1].tolist()[0]
    elif exp == 'Fear' and session == 'Conditioning':
        ind_start_trial = fiberbehavsnip_df.index[fiberbehavsnip_df['Shock'] == 1].tolist()[0]
    elif exp == 'Consumption':
        ind_start_trial = 0
    else:
        ind_start_trial = fiberbehavsnip_df.index[fiberbehavsnip_df['Entry in arena'] == 1].tolist()[0]

        
    
    #create list of behaviours and list of corresponding mean dFFs
    for behav in list_BOI:
        if behav in fiberbehav_df.columns[2:].tolist():
            if fiberbehavsnip_df[behav].sum() > 2:
                list_behav_analyzed.append(behav)
                meandFF_behav_df = fiberbehavsnip_df.groupby([behav], as_index=False).mean()
                print(meandFF_behav_df)
                list_meandFF.append(meandFF_behav_df.loc[meandFF_behav_df[behav]==1, 'Denoised dFF'].values[0])
                                
    #calculate mean dFF during baseline
    
    meandFF_baseline = fiberbehavsnip_df.loc[:ind_start_trial, 'Denoised dFF'].mean()
    
    #calculate mean dFF when no exploration and after gate opens
    
    meandFF_postbaseline_df = fiberbehavsnip_df.loc[ind_start_trial:]
    for behav in list_behav_analyzed:
        meandFF_postbaseline_df = meandFF_postbaseline_df.loc[meandFF_postbaseline_df[behav]==0]
        
    meandFF_postbaseline = meandFF_postbaseline_df['Denoised dFF'].mean()
    
    #calculate mean dFF when no exploration on total trial
    meandFF_df = fiberbehavsnip_df
    for behav in list_behav_analyzed:
        meandFF_df = meandFF_df.loc[meandFF_df[behav]==0]
        
    meandFF = meandFF_df['Denoised dFF'].mean()
    
    #create dataframe with values
    
    list_dFFs = [mouse, group, meandFF,meandFF_baseline,meandFF_postbaseline]
    list_columns = ['Subject','Group','Mean dFF','Baseline', 'Post_baseline']
    for (behav,meandFF) in zip(list_behav_analyzed,list_meandFF):
        list_dFFs.append(meandFF)
        list_columns.append(behav)
    meandFFs_df = pd.DataFrame(data=[list_dFFs], columns=list_columns)
    
    return(meandFFs_df)