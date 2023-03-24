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
    else:
        ind_start_trial = fiberbehavsnip_df.index[fiberbehavsnip_df['Entry in arena'] == 1].tolist()[0]
    
    #create list of behaviours and list of corresponding mean dFFs
    for behav in list_BOI:
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