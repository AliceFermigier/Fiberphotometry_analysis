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
import plotly.express as px

import modules.common.preprocess as pp

#%%
###################
#DEFINED FUNCTIONS#
###################

def hampel_filter(data, window_size, n_sigmas=5):
    k = 1.4826  # scaling factor for Gaussian distribution
    
    original_data = data.values if isinstance(data, pd.Series) else data
    new_data = original_data.copy()

    artifact_idx = []

    for i in range(window_size, len(original_data) - window_size):
        window = original_data[i - window_size:i + window_size + 1]
        median = np.nanmedian(window)
        mad = k * np.nanmedian(np.abs(window - median))

        if np.abs(original_data[i] - median) > n_sigmas * mad:
            new_data[i] = median
            artifact_idx.append(i)

    filtered = pd.Series(new_data, index=data.index) if isinstance(data, pd.Series) else new_data

    return filtered, artifact_idx

def plot_hampel_results(time, raw_405, raw_465, filt_405, filt_465, art405, art465):

    df_plot = pd.DataFrame({
        "Time": time,
        "Raw 405": raw_405,
        "Filtered 405": filt_405,
        "Raw 465": raw_465,
        "Filtered 465": filt_465
    })

    df_long = df_plot.melt(id_vars="Time", var_name="Signal", value_name="Value")

    fig = px.line(
        df_long,
        x="Time",
        y="Value",
        color="Signal",
        title="Raw vs Hampel Filtered Signals",
        color_discrete_map={
        "405 Raw": "blueviolet",
        "465 Raw": "deepskyblue",
        "405 Hampel": "mediumpurple",
        "465 Hampel": "lightskyblue"}
    )

    # Artifact markers
    fig.add_scatter(
        x=time.iloc[art405],
        y=filt_405.iloc[art405],
        mode="markers",
        marker=dict(size=8,color="blueviolet",opacity=0.3),
        name="Artifacts 405",
        hovertemplate="Artifact<br>Time: %{x}<br>Value: %{y}"
    )

    fig.add_scatter(
        x=time.iloc[art465],
        y=filt_465.iloc[art465],
        mode="markers",
        marker=dict(size=8,color="deepskyblue",opacity=0.3),
        name="Artifacts 465",
        hovertemplate="Artifact<br>Time: %{x}<br>Value: %{y}"
    )

    fig.show()

def plot_hampel_results_dualcolor(time, raw_405, raw_465, raw_560,
                                   filt_405, filt_465, filt_560,
                                   art405, art465, art560):

    df_long = pd.DataFrame({
        "Time": time,
        "405 Raw":     raw_405.values if isinstance(raw_405, pd.Series) else raw_405,
        "465 Raw":     raw_465.values if isinstance(raw_465, pd.Series) else raw_465,
        "560 Raw":     raw_560.values if isinstance(raw_560, pd.Series) else raw_560,
        "405 Hampel":  filt_405.values if isinstance(filt_405, pd.Series) else filt_405,
        "465 Hampel":  filt_465.values if isinstance(filt_465, pd.Series) else filt_465,
        "560 Hampel":  filt_560.values if isinstance(filt_560, pd.Series) else filt_560,
    }).melt(id_vars="Time", var_name="Signal", value_name="Value")

    fig = px.line(
        df_long,
        x="Time",
        y="Value",
        color="Signal",
        title="Raw vs Hampel Filtered Signals",
        color_discrete_map={
            "405 Raw":    "blueviolet",
            "465 Raw":    "deepskyblue",
            "560 Raw":    "orangered",
            "405 Hampel": "mediumpurple",
            "465 Hampel": "lightskyblue",
            "560 Hampel": "lightsalmon",
        }
    )

    fig.update_layout(
        legend=dict(font=dict(color="black")),
        title_font=dict(color="black")
    )

    # --- Artifact markers ---
    fig.add_scatter(
        x=time.iloc[art405],
        y=filt_405.iloc[art405],
        mode="markers",
        marker=dict(size=8, color="blueviolet", opacity=0.5),
        name="Artifacts 405",
        hovertemplate="Artifact<br>Time: %{x}<br>Value: %{y}"
    )
    fig.add_scatter(
        x=time.iloc[art465],
        y=filt_465.iloc[art465],
        mode="markers",
        marker=dict(size=8, color="deepskyblue", opacity=0.5),
        name="Artifacts 465",
        hovertemplate="Artifact<br>Time: %{x}<br>Value: %{y}"
    )
    fig.add_scatter(
        x=time.iloc[art560],
        y=filt_560.iloc[art560],
        mode="markers",
        marker=dict(size=8, color="orangered", opacity=0.5),
        name="Artifacts 560",
        hovertemplate="Artifact<br>Time: %{x}<br>Value: %{y}"
    )

    fig.show()
    
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

def remove_high_artifacts(rawdata_df):

    time = rawdata_df['Time(s)']
    data_405 = rawdata_df['405 Deinterleaved']
    data_465 = rawdata_df['465 Deinterleaved']

    # --- Hampel Filter ---
    hampel_405, artifacts_405 = hampel_filter(data_405, window_size=5, n_sigmas=5)
    hampel_465, artifacts_465 = hampel_filter(data_465, window_size=5, n_sigmas=5)

    print(f"405 artifacts removed: {len(artifacts_405)}")
    print(f"465 artifacts removed: {len(artifacts_465)}")

    plot_hampel_results(
        time,
        data_405,
        data_465,
        hampel_405,
        hampel_465,
        artifacts_405,
        artifacts_465)

    clean_deinterleaved_df = pd.DataFrame({
        'Time(s)': time,
        '405 Deinterleaved': hampel_405,
        '465 Deinterleaved': hampel_465
        })
    
    return clean_deinterleaved_df

def remove_high_artifacts_dualcolor(rawdata_df):
    time = rawdata_df['Time(s)']
    data_405 = rawdata_df['405 Deinterleaved']
    data_465 = rawdata_df['465 Deinterleaved']
    data_560 = rawdata_df['560 Deinterleaved']

    # --- Hampel Filter ---
    hampel_405, artifacts_405 = hampel_filter(data_405, window_size=5, n_sigmas=5)
    hampel_465, artifacts_465 = hampel_filter(data_465, window_size=5, n_sigmas=5)
    hampel_560, artifacts_560 = hampel_filter(data_560, window_size=5, n_sigmas=5)

    print(f"405 artifacts removed: {len(artifacts_405)}")
    print(f"465 artifacts removed: {len(artifacts_465)}")
    print(f"560 artifacts removed: {len(artifacts_560)}")

    plot_hampel_results_dualcolor(
        time,
        data_405,
        data_465,
        data_560,
        hampel_405,
        hampel_465,
        hampel_560,
        artifacts_405,
        artifacts_465,
        artifacts_560)

    clean_deinterleaved_df = pd.DataFrame({
        'Time(s)': time,
        '405 Deinterleaved': hampel_405,
        '465 Deinterleaved': hampel_465,
        '560 Deinterleaved' : hampel_560
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