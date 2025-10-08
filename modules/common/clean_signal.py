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
    k = 1.4826  # scaling factor for Gaussian distribution

    # Determine if input is a pandas Series or numpy array
    is_series = isinstance(data, pd.Series)
    original_data = data.values if is_series else data
    new_data = original_data.copy()

    for i in range(window_size, len(original_data) - window_size):
        window = original_data[i - window_size:i + window_size + 1]
        median = np.nanmedian(window)
        mad = k * np.nanmedian(np.abs(window - median))

        if np.abs(original_data[i] - median) > n_sigmas * mad:
            new_data[i] = median

    # Return result in the same format as input
    return pd.Series(new_data, index=data.index) if is_series else new_data

def highpass_filter(data_df, sr, cutoff=0.01, order=1):
    """
    High-pass filters the signal to remove slow trends.

    Parameters:
    - signal: 1D numpy array or list of your raw fluorescence values
    - cutoff: cutoff frequency in Hz (e.g., 0.01 Hz = 100 sec cycles)
    - fs: sampling rate in Hz (10 Hz in your case)
    - order: filter order (higher = sharper cutoff)

    Returns:
    - detrended signal as a NumPy array
    """
    nyq = 0.5 * sr
    norm_cutoff = cutoff / nyq
    b, a = butter(order, norm_cutoff, btype='high', analog=False)
    filtered_signal = filtfilt(b, a, data_df)
    return filtered_signal

def highpass_filter_with_padding(signal, sr, cutoff=0.01, order=3, pad_seconds=50):
    pad_len = int(sr * pad_seconds)
    pre_pad = signal[:pad_len][::-1] if pad_len < len(signal) else signal[::-1]
    post_pad = signal[-pad_len:][::-1] if pad_len < len(signal) else signal[::-1]
    padded = np.concatenate([pre_pad, signal, post_pad])

    nyq = 0.5 * sr
    norm_cutoff = cutoff / nyq
    b, a = butter(order, norm_cutoff, btype='high', analog=False)
    filtered = filtfilt(b, a, padded)

    return filtered[pad_len:-pad_len]

def clean_signal(rawdata_df, crop=[0,-10], detrending=False, apply_hampel=True):

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

    clean_deinterleaved_df = pd.DataFrame({
        'Time(s)': time,
        '405 Deinterleaved': detrended_405,
        '465 Deinterleaved': detrended_465
        })
    
    return clean_deinterleaved_df

def clean_signal_dualcolor(rawdata_df, crop=[10,-10], detrending=False, apply_hampel=True):
    time = rawdata_df['Time(s)'][crop[0]:crop[1]]
    detrended_405 = rawdata_df['405 Deinterleaved'][crop[0]:crop[1]]
    detrended_465 = rawdata_df['465 Deinterleaved'][crop[0]:crop[1]]
    detrended_560 = rawdata_df['560 Deinterleaved'][crop[0]:crop[1]]

    # --- Detrend ---
    if detrending:
        detrended_405 = detrend(detrended_405, type='linear')
        detrended_465 = detrend(detrended_465, type='linear')
        detrended_560 = detrend(detrended_560, type='linear')
        plt.plot(time, detrended_465, linewidth=1, color='deepskyblue', label='GCaMP')
        plt.plot(time, detrended_405, linewidth=1, color='blueviolet', label='ISOS')
        plt.plot(time, detrended_560, linewidth=1, color='orange', label='rGECO')
        plt.legend()
        plt.title("Detrending")
        plt.show()

    # --- Hampel Filter ---
    if apply_hampel:
        detrended_hampel_405 = hampel_filter(detrended_405, window_size=5, n_sigmas=5)
        detrended_hampel_465 = hampel_filter(detrended_465, window_size=5, n_sigmas=5)
        detrended_hampel_560 = hampel_filter(detrended_560, window_size=5, n_sigmas=5)
        plt.plot(time, detrended_hampel_465, linewidth=1, color='deepskyblue', label='GCaMP')
        plt.plot(time, detrended_hampel_405, linewidth=1, color='blueviolet', label='ISOS')
        plt.plot(time, detrended_hampel_560, linewidth=1, color='orange', label='rGECO')
        plt.legend()
        plt.title("Hampel Filtering")
        plt.show()
        detrended_405 = detrended_hampel_405
        detrended_465 = detrended_hampel_465
        detrended_560 = detrended_hampel_560

    clean_deinterleaved_df = pd.DataFrame({
        'Time(s)': time,
        '405 Deinterleaved': detrended_405,
        '465 Deinterleaved': detrended_465,
        '560 Deinterleaved' : detrended_560
        })
    
    return clean_deinterleaved_df

def highpass_filter_dff(dff, dualcolor = False):
    sr = pp.samplerate(dff)
    cutoff_freq = 0.01
    denoised_dff = dff['Denoised dFF']
    time = dff['Time(s)']

    if dualcolor == True:
        dff_560 = dff['Denoised 560 dFF']
        filtered_denoised_dff = highpass_filter_with_padding(
            denoised_dff, sr, cutoff=cutoff_freq, order=1, pad_seconds=50
        )
        filtered_560_denoised_dff = highpass_filter_with_padding(
            dff_560, sr, cutoff=cutoff_freq, order=1, pad_seconds=50
        )
        dff['Denoised 560 dFF'] = filtered_560_denoised_dff
        
    else:
        filtered_denoised_dff = highpass_filter_with_padding(
            denoised_dff, sr, cutoff=cutoff_freq, order=1, pad_seconds=50
        )

        # Plot settings
        fig, axs = plt.subplots(2, 1, figsize=(12, 6), sharex=True, gridspec_kw={'height_ratios': [1, 1]})
        
        # Unfiltered
        axs[0].plot(time, denoised_dff, color='black', linewidth=1)
        axs[0].set_title('Unfiltered dF/F')
        axs[0].set_ylabel('dF/F (%)')

        # Filtered
        axs[1].plot(time, filtered_denoised_dff, color='seagreen', linewidth=1)
        axs[1].set_title(f'Filtered dF/F (High-pass {cutoff_freq} Hz)')
        axs[1].set_xlabel('Time (s)')
        axs[1].set_ylabel('dF/F (%)')

        # Adjust layout
        plt.tight_layout()
        plt.show()

    dff['Denoised dFF'] = filtered_denoised_dff
    return dff