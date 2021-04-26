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
import matplotlib.pyplot as plt
from pathlib import Path
import os

#%%
########
#LOADER#
########

# experiment Path mac : /Volumes/My Passport/Alice/Fiber/202103_CDHFDdvHPC/
# experiment Path windows : E:/Alice/Fiber/202103_CDHFDdvHPC/

#path to experiment and analysis folder
experiment_path = Path('E:/Alice/Fiber/202103_CDHFDdvHPC/')
analysis_path = experiment_path / '202103_Analysis/20210302_ORM/Test'
mouse_path = analysis_path / 'CD1/'

#import ID and groups of all mice
subjects_path = experiment_path / '202103_Analysis/subjects.xlsx'
subjects_df = pd.read_excel(subjects_path)

#all analysis files in Batch/Experiment/Analysis/Subject_Trial
#each trial of each subject has a specific directory called CD1_0, CD1_2, etc...
#for example, for ORM : CD1_0 = training of subject CD1 ; CD1_2 = test of subject CD1

rawdata_path = mouse_path / 'CD1_rawdata.csv'
rawdata_df = pd.read_csv(rawdata_path)

#to be created:
    
procdata_path = mouse_path / 'CD1_procdata.csv'

#%%
###################
#DEFINED FUNCTIONS#
###################

def dFF(rawdata_df):
    
    x = np.array(data.streams[ISOS].data)
    y = np.array(data.streams[GCAMP].data)
    bls = np.polyfit(x, y, 1)
    Y_fit_all = np.multiply(bls[0], x) + bls[1]
    Y_dF_all = y - Y_fit_all

    dFF = np.multiply(100, np.divide(Y_dF_all, Y_fit_all))
    std_dFF = np.std(dFF)
    
    return(dFF_df)

#%%
########
#SCRIPT#
########