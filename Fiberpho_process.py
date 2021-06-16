# -*- coding: utf-8 -*-
"""
Created on Fri Mar 19 13:58:32 2021

Processing of fiberpho data
Full trace dFF according to Lerner et al. 2015
https://www.tdt.com/support/python-sdk/offline-analysis-examples/licking-bout-epoc-filtering/

@author: Alice
"""

#%%
##########
#IMPORTED#
##########

import pandas as pd
import numpy as np
import math
import os

#Custom
#put path to directory where python files are stored
sys.path.append('C:\\Users\\Alice\\Documents\\GitHub\\Fiberphotometry_analysis')

#%%
########
#LOADER#
########

from Fiberpho_loader import experiment_path, analysis_path, subjects_df, SAMPLERATE
from Fiberpho_loader import list_EVENT, list_TIMEWINDOW, PRE_EVENT_TIME, TIME_BEGIN

os.chdir(experiment_path)

#%%
###################
#DEFINED FUNCTIONS#
###################

def rem_artifacts(rawdata_df):
    """
    rawdata_df = pandas dataframe from csv file
    """
    
    
    
    return(data)

def rem_photobleach(data):
    """
    data = np array
    """
    
    return(filt_data)

def dFF(filt_data):
    """
    filt_data = np array
    """
    
    return(filt_dFF)


#%%
########
#SCRIPT#
########