# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 09:54:44 2024

@author: alice
"""
from scipy.signal import find_peaks
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, filtfilt 

def bandpass_filter(data, lowcut, highcut, sr, order=3):
    nyquist = 0.5 * sr
    low = lowcut/nyquist
    high = highcut/nyquist
    b, a = butter(order, [low, high], btype='band')
    print("Coefficients du filtre (b):", b)
    print("Coefficients du filtre (a):", a)
    y = filtfilt(b, a, data)
    return y
    
def plot_signal_and_spectrum(time, signal, filtered_signal, fs):
    plt.figure(figsize=(20, 5))

    # Plot du signal original et filtré
    plt.subplot(1, 2, 1)
    plt.plot(time, signal, label='Signal Original', color='grey')
    plt.plot(time, filtered_signal, label='Signal Filtré', color='black')
    plt.xlabel('Temps (s)')
    plt.ylabel('Amplitude')
    plt.title('Filtrage Passe-bande du Signal de Photométrie')
    plt.legend()

    # Spectre de fréquence
    plt.subplot(1, 2, 2)
    freqs, psd = plt.psd(signal, Fs=fs, NFFT=1024, color='grey', label='Signal Original')
    plt.psd(filtered_signal, Fs=fs, NFFT=1024, color='black', label='Signal Filtré')
    plt.xlabel('Fréquence (Hz)')
    plt.ylabel('Puissance Spectrale')
    plt.legend()

    plt.tight_layout()
    plt.show()

def mad(data, axis=None):
    return np.median(np.abs(data - np.median(data, axis)), axis)

def transients(fiberpho_df, sr):
    # Détection des pics
    fiberpeaks_df = fiberpho_df
    filtered_sig=fiberpeaks_df['Filtered dFF']
    two_MAD_th = np.median(filtered_sig) + (2 * mad(filtered_sig))
    one_MAD_th = np.median(filtered_sig) + mad(filtered_sig)
    peaks, properties = find_peaks(filtered_sig, height=one_MAD_th, distance=sr*0.5)
    
    # Ajout des pics détectés au DataFrame
    fiberpeaks_df['Peaks'] = 0
    fiberpeaks_df.loc[peaks, 'Peaks'] = 1
    
    # Amplitudes des pics
    mean_peak_amplitudes = np.mean(properties['peak_heights'])
    
    # Fréquence des pics (nombre de pics par unité de temps)
    peak_frequency = len(peaks) / (fiberpeaks_df['Time(s)'].iloc[-1] - fiberpeaks_df['Time(s)'].iloc[0])
    
    #Plot transients
    plt.figure(figsize=(15, 5))
    plt.plot(fiberpeaks_df['Time(s)'], filtered_sig, label='Signal Filtré')
    plt.plot(fiberpeaks_df['Time(s)'].iloc[peaks],filtered_sig.iloc[peaks], 'ro', label='Pics')
    plt.xlabel('Temps (s)')
    plt.ylabel('Signal')
    plt.title('Détection des Pics dans le Signal de Photométrie')
    plt.legend()
    plt.show()
    
    return(fiberpeaks_df, peak_frequency, mean_peak_amplitudes)