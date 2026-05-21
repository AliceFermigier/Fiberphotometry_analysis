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
from scipy.signal import butter, filtfilt
from scipy.optimize import curve_fit
import plotly.express as px

import modules.common.median_filtering as mf
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

def highpass_filter_with_padding(signal, sr, cutoff=0.01, order=3, pad_seconds=50):
    pad_n   = int(pad_seconds * sr)
    padded  = np.pad(signal, pad_n, mode='reflect')

    nyq = sr / 2
    norm_cutoff = cutoff / nyq
    b, a = butter(order, norm_cutoff, btype='high', analog=False)
    filtered = filtfilt(b, a, padded)

    return filtered[pad_n:-pad_n]

def lowpass_filter_with_padding(signal, sr, cutoff=0.01, order=1, pad_seconds=50):
    pad_n   = int(pad_seconds * sr)
    padded  = np.pad(signal, pad_n, mode='reflect')

    nyq     = sr / 2
    norm_cutoff = cutoff / nyq
    b, a = butter(order, norm_cutoff, btype='low', analog=False)
    filtered = filtfilt(b, a, padded)

    return filtered[pad_n:-pad_n]

def remove_high_artifacts(rawdata_df):

    time = rawdata_df['Time(s)']
    data_405 = rawdata_df['405 Deinterleaved'].copy()
    data_465 = rawdata_df['465 Deinterleaved'].copy()

    # --- Hampel Filter ---
    hampel_405, artifacts_405 = hampel_filter(data_405, window_size=10, n_sigmas=6)
    hampel_465, artifacts_465 = hampel_filter(data_465, window_size=10, n_sigmas=6)

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
    data_405 = rawdata_df['405 Deinterleaved'].copy()
    data_465 = rawdata_df['465 Deinterleaved'].copy()
    data_560 = rawdata_df['560 Deinterleaved'].copy()

    # --- Hampel Filter ---
    hampel_405, artifacts_405 = hampel_filter(data_405, window_size=10, n_sigmas=6)
    hampel_465, artifacts_465 = hampel_filter(data_465, window_size=10, n_sigmas=6)
    hampel_560, artifacts_560 = hampel_filter(data_560, window_size=10, n_sigmas=5)

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

def plot_highpass_filter_results(time, dff, filtered_dff, cutoff_freq, title = 'dF/F'):
    # Plot settings
    fig, axs = plt.subplots(2, 1, figsize=(12, 6), sharex=True, gridspec_kw={'height_ratios': [1, 1]})
    
    # Unfiltered
    axs[0].plot(time, dff, color='black', linewidth=1)
    axs[0].set_title(f'Unfiltered {title}')
    axs[0].set_ylabel('dF/F (%)')

    # Filtered
    axs[1].plot(time, filtered_dff, color='black', linewidth=1)
    axs[1].set_title(f'Filtered {title} (High-pass {cutoff_freq} Hz)')
    axs[1].set_xlabel('Time (s)')
    axs[1].set_ylabel('dF/F (%)')

    # Adjust layout
    plt.tight_layout()
    plt.show()

def plot_exponential_fit_results(time, dff, filtered_dff, title = 'dF/F'):
    # Plot settings
    fig, axs = plt.subplots(2, 1, figsize=(12, 6), sharex=True, gridspec_kw={'height_ratios': [1, 1]})

    # Unfiltered
    axs[0].plot(time, dff, color='black', linewidth=1)
    axs[0].set_title(f'Original {title}')
    axs[0].set_ylabel('dF/F (%)')

    # Filtered
    axs[1].plot(time, filtered_dff, color='black', linewidth=1)
    axs[1].set_title(f'Detrended {title} - Type : exponential')
    axs[1].set_xlabel('Time (s)')
    axs[1].set_ylabel('dF/F (%)')

    # Adjust layout
    plt.tight_layout()
    plt.show()

def highpass_filter_dff(dff_df, dualcolor=False, cutoff_freq=0.001,
                         use_smooth_baseline=True):
    sr   = pp.samplerate(dff_df)
    time = dff_df['Time(s)']

    def _transient_safe_highpass(df, column_name):
        raw_signal = df[column_name].values
        if use_smooth_baseline:
            # iterative_median_filter returns (result_df, best_win_s, fig)
            result_df, _, _ = mf.iterative_median_filter(df, column_name, step_size=1.0)
            smooth_signal   = result_df[column_name].values         
            baseline        = lowpass_filter_with_padding(
                smooth_signal, sr, cutoff=cutoff_freq, order=1, pad_seconds=50)
            return raw_signal - baseline
        else:
            return highpass_filter_with_padding(
                raw_signal, sr, cutoff=cutoff_freq, order=1, pad_seconds=50)

    # ── 465 dFF ───────────────────────────────────────────────────────────────
    raw_dff_values   = dff_df['dFF'].copy().values                   
    filtered_dff     = _transient_safe_highpass(dff_df.copy(), 'dFF')
    dff_df['dFF']    = filtered_dff
    plot_highpass_filter_results(time, raw_dff_values, filtered_dff, cutoff_freq)

    # ── 560 dFF (optional) ────────────────────────────────────────────────────
    if dualcolor:
        raw_560_values    = dff_df['560 dFF'].copy().values         
        filtered_560      = _transient_safe_highpass(dff_df.copy(), '560 dFF')
        dff_df['560 dFF'] = filtered_560
        plot_highpass_filter_results(time, raw_560_values, filtered_560, cutoff_freq, title='560 dF/F')

    return dff_df

def exp_func(t, A, tau, C):
    return A * np.exp(-t / tau) + C

def exponential_detrend(dff_df, dualcolor = False):

    dff_465 = dff_df['dFF'].copy()
    time = dff_df['Time(s)']

    # Initial parameter guesses
    p0 = [np.max(dff_465), np.mean(time), np.min(dff_465)]

    # Fit exponential decay
    params, _ = curve_fit(exp_func, time, dff_465, p0=p0)
    trend = exp_func(time, *params)
    detrended_dff_465 = dff_465 - trend

    plot_exponential_fit_results(time, dff_465, detrended_dff_465)

    dff_df['dFF'] = detrended_dff_465 

    if dualcolor == True:
        dff_560 = dff_df['560 dFF'].copy()
        # Initial parameter guesses
        p0 = [np.max(dff_560), np.mean(time), np.min(dff_560)]

        # Fit exponential decay
        params, _ = curve_fit(exp_func, time, dff_560, p0=p0)
        trend = exp_func(time, *params)
        detrended_dff_560 = dff_560 - trend

        plot_exponential_fit_results(time, dff_560, detrended_dff_560, title = '560 dF/F')

        dff_df['560 dFF'] = detrended_dff_560 

    return dff_df

def lowpass_dFF(dff, dual_color = False, order = 2, cut_freq = 10):

    sampling_rate = pp.samplerate(dff)
    raw_dff = dff['dFF']

    # Lowpass filter - zero phase filtering (with filtfilt) is used to avoid distorting the signal.
    b,a = butter(order, cut_freq, btype='low', fs=sampling_rate)
    dFF_lowpass = filtfilt(b,a, raw_dff)

    dff['dFF'] = dFF_lowpass

    if dual_color:
        raw_dff = dff['560 dFF']
        
        b,a = butter(order, cut_freq, btype='low', fs=sampling_rate)
        dFF_lowpass_560 = filtfilt(b,a, raw_dff)
        dff['560 dFF'] = dFF_lowpass_560

    return dff

def smoothing_moving_average(signal: np.ndarray, window_samples: int) -> np.ndarray:
    """
    Simple symmetric moving-average smoothing.
    Mirrors MATLAB's smooth(x, window) default behaviour.
    Edge values are computed with a progressively smaller window
    (same as np.convolve mode='same').
    """
    if window_samples < 1:
        return signal.copy()
    return np.convolve(signal, np.ones(window_samples) / window_samples, mode="same")