#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  5 11:59:06 2024

@author: alice
"""

##########
#IMPORTED#
##########

import os as os

##########
#%LOADER#
##########

path = "/Users/alice/Desktop/FiberD1D2/Data"

#################################
#%DEFINED FUNCTIONS AND CLASSES#
#################################

##########
#%SCRIPT#
##########

for dirName, subdirList, fileList in os.walk(path):
    for f in fileList:
        for to_replace in ['_0','_2','_3']:
            if to_replace in f:
                old_path = f'{dirName}/{f}'
                new_name = f.replace(to_replace,'_1')
                new_path = f'{dirName}/{new_name}'
                os.rename(old_path, new_path)
                break