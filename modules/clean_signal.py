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
from scipy.signal import butter, filtfilt, detrend

import nomenclature as nom

#%%
###################
#DEFINED FUNCTIONS#
###################

def hampel_filter(data, window_size=5, n_sigmas=5):
    """
    Hampel filter for outlier detection
    """
    n = len(data)
    new_series = data.copy()
    k = 1.4826  # scale factor for Gaussian distribution
    for i in range(window_size, n - window_size):
        window = data[i - window_size:i + window_size + 1]
        median = np.nanmedian(window)
        mad = k * np.nanmedian(np.abs(window - median))
        if np.abs(data[i] - median) > n_sigmas * mad:
            new_series[i] = median
    return new_series

def clean_signal(t, dff, raw_signals, apply_hampel=True, apply_filter=False):

    # --- Detrend ---
    dff = detrend(dff, type='linear')
    plt.plot(t, dff)
    plt.title("After Detrending")
    plt.show()

    # --- Hampel Filter ---
    if apply_hampel:
        dff1 = hampel_filter(dff, window_size=5, n_sigmas=5)
        plt.plot(t, dff, 'b', label='Original')
        plt.plot(t, dff1, 'r', label='Hampel')
        plt.legend()
        plt.title("Hampel Filtering")
        plt.show()
        dff = dff1

    # --- High-pass Filter --- (not necessary if detrend is sufficient)
    if apply_filter:
        b, a = butter(2, 10 / (raw_signals['sfreq'] / 2), btype='high')
        filt_dff = filtfilt(b, a, dff)
        plt.plot(t, dff, 'b', label='Before Filter')
        plt.plot(t, filt_dff, 'r', label='Filtered')
        plt.legend()
        plt.title("High-pass Filtering")
        plt.show()
        dff = filt_dff