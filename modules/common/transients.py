# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 09:54:44 2024
@author: alice fermigier
"""
from scipy.signal import find_peaks
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize
from scipy.signal import butter, filtfilt 

import modules.common.preprocess as pp

def bandpass_filter(data, lowcut, highcut, order=3):
    sr = pp.samplerate(data)
    print('samplerate:', sr)
    print('lowcut:', lowcut)
    print('highcut:', highcut)
    signal = data['dFF'].values
    print('signal:', signal)
    nyquist = 0.5 * sr
    low = lowcut/nyquist
    high = highcut/nyquist
    b, a = butter(order, [low, high], btype='band')
    print("b:", b)
    print("a:", a)
    y = filtfilt(b, a, signal)
    return y
    
def plot_signal_and_spectrum(dfiber_df):
    plt.figure(figsize=(20, 5))
    time = dfiber_df['Time(s)']
    signal = dfiber_df['dFF']
    filtered_signal = dfiber_df['Filtered dFF']
    sr = pp.samplerate(dfiber_df)

    # Plot du signal original et filtré
    plt.subplot(1, 2, 1)
    plt.plot(time, signal, label='Original Signal', color='grey')
    plt.plot(time, filtered_signal, label='Filtered Signal', color='black')
    plt.xlabel('Time(s)')
    plt.ylabel('Amplitude')
    plt.title('Band-Pass Filtering')
    plt.legend()

    # Spectre de fréquence
    plt.subplot(1, 2, 2)
    freqs, psd = plt.psd(signal, Fs=sr, NFFT=1024, color='grey', label='Original Signal')
    plt.psd(filtered_signal, Fs=sr, NFFT=1024, color='black', label='Filtered Signal')
    plt.xlabel('Frequency(Hz)')
    plt.ylabel('Spectral Power')
    plt.legend()
    plt.tight_layout()
    plt.show()

def mad(data, axis=None):
    return np.median(np.abs(data - np.median(data, axis)), axis)

def transients(fiberpho_df, threshold='one_MAD'):
    # Détection des pics
    sr=pp.samplerate(fiberpho_df)
    fiberpeaks_df = fiberpho_df
    filtered_sig=fiberpeaks_df['Filtered dFF']
    if threshold == 'two_MAD':
        MAD_th = np.median(filtered_sig) + (2 * mad(filtered_sig))
    elif threshold == 'one_MAD':
        MAD_th = np.median(filtered_sig) + mad(filtered_sig)
    peaks, properties = find_peaks(filtered_sig, height=MAD_th, distance=sr*0.5)
    
    # Ajout des pics détectés au DataFrame
    fiberpeaks_df['Peaks'] = 0
    fiberpeaks_df.loc[peaks, 'Peaks'] = 1
    
    # Amplitudes des pics
    peak_amplitudes = properties['peak_heights']
    mean_peak_amplitudes = np.mean(peak_amplitudes)
    
    # Fréquence des pics (nombre de pics par unité de temps)
    peak_frequency = len(peaks) / (fiberpeaks_df['Time(s)'].iloc[-1] - fiberpeaks_df['Time(s)'].iloc[0])

    # Normalize amplitudes for colormap
    norm = Normalize(vmin=np.min(peak_amplitudes), vmax=np.max(peak_amplitudes))
    colors = cm.viridis(norm(peak_amplitudes))

    # Create figure and axis
    fig, ax = plt.subplots(figsize=(20, 5))

    # Plot main signal
    ax.plot(fiberpeaks_df['Time(s)'], filtered_sig, color='black', linewidth=.6, label='Filtered Signal')

    # Scatter peaks with color mapping
    sc = ax.scatter(
        fiberpeaks_df['Time(s)'].iloc[peaks],
        filtered_sig.iloc[peaks],
        c=colors,
        s=40,
        label='Peaks'
    )

    # Add colorbar and attach it to the axis
    cbar = fig.colorbar(cm.ScalarMappable(norm=norm, cmap='viridis'), ax=ax)
    cbar.set_label('Peak Amplitude')

    # Add labels and legend
    ax.set_xlabel('Time(s)')
    ax.set_ylabel('dFF')
    ax.set_title(f'Peaks in Photometry Signal - {threshold}')
    ax.legend()
    ax.margins(x=0.005)
    fig.subplots_adjust(right=0.3)
    plt.tight_layout()
    plt.show()
    
    return(fiberpeaks_df, peak_frequency, mean_peak_amplitudes, fig)