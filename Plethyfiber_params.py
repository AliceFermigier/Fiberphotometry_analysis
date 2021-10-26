# -*- coding: utf-8 -*-
"""
Created on Tue Oct 26 14:19:39 2021

Parameters for plethysmograph and fiberphotometry data

@author: afermigier
"""

############
#PARAMETERS#
############

#RecDuration = 80000;
Shift = 2000; #in ms
SamplingRate = 10000; #in Hz
LOW_CUTOFF_FREQ = 3;
HIGH_CUTOFF_FREQ = 30;
StimStart_list = [380, 440, 500, 560, 620]
StimDurTh = 30 #in s
PostStimDur = 10;
minpeakdist = 70; #quantity of datapoints in between 2 peaks
thresh = 0;
#minpeakheight = 0.005; % calculated below, based on 80th percentile of the signal
ArTfact = 0.2;
CROP = 500 #in ms

