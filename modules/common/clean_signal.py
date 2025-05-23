# -*- coding: utf-8 -*-
"""
Functions to clean fiberphotometry data
To be used after manual removal of big artifacts caused by patch cord movement (i.e. directly on dff signals)

@author: alice fermigier
"""

#%%
##########
#IMPORTED#
##########

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.signal import butter, filtfilt, detrend

import modules.common.preprocess as pp

#%%
###################
#DEFINED FUNCTIONS#
###################

def hampel_filter(data, window_size, n_sigmas=5):
    new_series = data.copy()
    k = 1.4826  # scaling factor for Gaussian distribution

    for i in range(window_size, len(data) - window_size):
        window = data.iloc[i - window_size:i + window_size + 1]
        median = np.nanmedian(window)
        mad = k * np.nanmedian(np.abs(window - median))
        
        if np.abs(data.iloc[i] - median) > n_sigmas * mad:
            new_series.iloc[i] = median

    return new_series

def clean_signal(rawdata_df, crop=[50,-10], detrending=False, apply_hampel=True, apply_filter=False):

    time = rawdata_df['Time(s)'][crop[0]:crop[1]]
    detrended_405 = rawdata_df['405 Deinterleaved'][crop[0]:crop[1]]
    detrended_465 = rawdata_df['465 Deinterleaved'][crop[0]:crop[1]]

    # --- Detrend ---
    if detrending:
        detrended_405 = detrend(detrended_405, type='linear')
        detrended_465 = detrend(detrended_465, type='linear')
        plt.plot(time, detrended_465, linewidth=1, color='deepskyblue', label='GCaMP')
        plt.plot(time, detrended_405, linewidth=1, color='blueviolet', label='ISOS')
        plt.legend()
        plt.title("Detrending")
        plt.show()

    # --- Hampel Filter ---
    if apply_hampel:
        detrended_hampel_405 = hampel_filter(detrended_405, window_size=5, n_sigmas=5)
        detrended_hampel_465 = hampel_filter(detrended_465, window_size=5, n_sigmas=5)
        plt.plot(time, detrended_hampel_465, linewidth=1, color='deepskyblue', label='GCaMP')
        plt.plot(time, detrended_hampel_405, linewidth=1, color='blueviolet', label='ISOS')
        plt.legend()
        plt.title("Hampel Filtering")
        plt.show()
        detrended_405 = detrended_hampel_405
        detrended_465 = detrended_hampel_465

    # --- High-pass Filter --- (not necessary if detrending is sufficient)
    if apply_filter:
        sr = pp.samplerate(rawdata_df)
        nyquist = sr / 2
        cutoff = 0.001 / nyquist
        b, a = butter(N=2, Wn=cutoff, btype='high')

        detrended_filtered_405 = filtfilt(b, a, detrended_405)
        detrended_filtered_465 = filtfilt(b, a, detrended_465)
        plt.plot(time, detrended_filtered_465, linewidth=1, color='deepskyblue', label='GCaMP')
        plt.plot(time, detrended_filtered_405, linewidth=1, color='blueviolet', label='ISOS')
        plt.legend()
        plt.title("High-pass Filtering (10Hz Cutoff)")
        plt.show()
        detrended_405 = detrended_filtered_405
        detrended_465 = detrended_filtered_465

    clean_deinterleaved_df = pd.DataFrame({
        'Time(s)': time,
        '405 Deinterleaved': detrended_405,
        '465 Deinterleaved': detrended_465
        })
    
    return clean_deinterleaved_df