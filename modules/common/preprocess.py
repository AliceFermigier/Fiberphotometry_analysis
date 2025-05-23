# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 15:51:56 2023

Functions for preprocessing fiberphotometry data

@author: alice fermigier
"""

#%%
##########
#IMPORTED#
##########

import pandas as pd
import numpy as np
from scipy import signal
from ast import literal_eval
import h5py

import modules.common.nomenclature as nom


#%%
###################
#DEFINED FUNCTIONS#
###################

def deinterleave(rawdata_df):
    """
    Deinterleaves photometry signals from raw data.
    
    --> Parameters:
        rawdata_df : pd.DataFrame
            A dataframe containing the following columns:
            - 'Time(s)' : Time of each sample
            - 'DI/O-1' : Digital input for 405nm signal
            - 'DI/O-2' : Digital input for 470nm signal
            - 'AIn-1' : Analog input (combined signal)
    
    --> Returns:
        deinterleaved_df : pd.DataFrame
            Dataframe with columns:
            - 'Time(s)' : Aligned time points for 405/470 signals
            - '405 Deinterleaved' : Deinterleaved signal for 405nm
            - '470 Deinterleaved' : Deinterleaved signal for 470nm
    """
    
    # Calculate the rising edges of DI/O-1 and DI/O-2
    derivative405 = rawdata_df['DI/O-1'].diff()
    derivative470 = rawdata_df['DI/O-2'].diff()
    
    # Find indices of rising edges (where derivative = 1)
    indices_405 = np.where(derivative405 == 1)[0]
    indices_470 = np.where(derivative470 == 1)[0]
    
    # Ensure indices do not exceed the length of the DataFrame when adding 250
    max_index = len(rawdata_df) - 1  # last valid index in the DataFrame
    indices_405 = indices_405[indices_405 + 250 <= max_index]
    indices_470 = indices_470[indices_470 + 250 <= max_index]
    
    # Extract the 405 and 470 signal values at (index + 250)
    list_405 = rawdata_df.loc[indices_405 + 250, 'AIn-1'].values
    list_470 = rawdata_df.loc[indices_470 + 250, 'AIn-1'].values

    # Make sure both lists have the same length
    min_length = min(len(list_405), len(list_470))
    list_405 = list_405[:min_length]
    list_470 = list_470[:min_length]

    # Create a time vector that corresponds to the deinterleaved signals
    timevector = np.linspace(0, rawdata_df['Time(s)'].max(), min_length)
    
    # Create a DataFrame to store deinterleaved signals
    deinterleaved_df = pd.DataFrame({
        'Time(s)': timevector,
        '405 Deinterleaved': list_405,
        '470 Deinterleaved': list_470
    })
    
    # Replace zeros with NaN (if necessary) to handle empty or missing signal values
    deinterleaved_df.replace(0, value=np.nan, inplace=True)
        
    return deinterleaved_df

def load_deinterleaved_doric(file_path):
    with h5py.File(file_path, 'r') as f:
        base = "DataAcquisition/FPConsole/Signals/Series0001/"
        
        t_405 = f[base + "AIN01xDIO01-Deinterleaved/Time"][:]
        sig_405 = f[base + "AIN01xDIO01-Deinterleaved/Values"][:]
        
        t_465 = f[base + "AIN01xDIO02-Deinterleaved/Time"][:]
        sig_465 = f[base + "AIN01xDIO02-Deinterleaved/Values"][:]

    # Truncate so that all signals have the same length
    min_len = min(len(sig_465), len(sig_405))
    sig_465 = sig_465[:min_len]
    sig_405 = sig_405[:min_len]
    t_465 = t_465[:min_len]
    t_405 = t_405[:min_len]

    deinterleaved_df = pd.DataFrame({
        'Time(s)': t_405,
        '405 Deinterleaved': sig_405,
        '465 Deinterleaved': sig_465
    })

    return deinterleaved_df

def samplerate(data_df):
    
    sr = len(data_df)/(data_df['Time(s)'].max()-data_df['Time(s)'].min())
    
    return sr

def update_artifacts_file(file_path, filecode, artifacts):
    """
    Update the artifacts file with a new entry or update an existing entry for the given filecode.
    
    Parameters:
    - file_path (Path): Path to the artifacts Excel file.
    - filecode (str): Unique identifier for the file (e.g., 'Experiment1_Session1_MouseA1').
    - artifacts (list of tuples): List of artifact intervals (e.g., [(start1, end1), (start2, end2)]).
    """
    # Load the existing file (or create a new one if it doesn't exist)
    df = nom.create_or_load_artifacts_file(file_path)
    
    # Convert the list of artifact intervals to a string for storage in Excel
    artifacts_str = str(artifacts)  # Store as a string representation of the list
    
    # Check if the filecode already exists in the DataFrame
    if filecode in df['Filecode'].values:
        print(f"Filecode '{filecode}' already exists. Updating artifact data.")
        df.loc[df['Filecode'] == filecode, 'Artifacts'] = artifacts_str
    else:
        print(f"Adding new entry for Filecode '{filecode}'.")
        new_entry = {'Filecode': filecode, 'Artifacts': artifacts_str}
        df = pd.concat([df, pd.DataFrame([new_entry])], ignore_index=True)
    
    # Save the updated DataFrame back to the Excel file
    df.to_excel(file_path, index=False)
    print(f"Updated Excel file at: {file_path}")

def linearfit_sklearn(sig_405, sig_465):
    # Fit 405 to 465 using linear regression
    from sklearn.linear_model import LinearRegression
    model = LinearRegression()
    
    if isinstance(sig_405, np.ndarray):
        sig_405 = sig_405.reshape(-1, 1)
    else:
        sig_405 = sig_405.to_numpy().reshape(-1, 1)
        sig_465 = sig_465.to_numpy()

    model.fit(sig_405, sig_465)
    fitted_405 = model.predict(sig_405)

    return fitted_405

def controlFit(control, signal): # Deprecated
    """
    Fits a linear model to control vs. signal and predicts signal using the fit.
    
    Parameters:
    - control (pd.Series or np.array): 405nm control signal
    - signal (pd.Series or np.array): 470nm signal to be fitted to the control
    
    Returns:
    - np.array: Fitted signal using a linear model
    """
    p = np.polyfit(control, signal, 1)  # Linear fit (y = p[0]*x + p[1])
    fitted_signal = (p[0] * control) + p[1]
    return fitted_signal

def remove_artifacts(data_df, artifact_intervals, col, method='mean'):
    """
    Helper function to remove artifacts from a specific column of the data.
    
    Parameters:
    - data_df (pd.DataFrame): Input data
    - artifact_intervals (list of tuples): List of artifact intervals as [(start, stop), ...]
    - col (str): Column to process ('405 Deinterleaved' or '470 Deinterleaved')
    - begin (int): Starting index for the segment
    - end (int): Ending index for the segment
    - sr (int): Sampling rate
    - method (str): 'mean' or 'fit' method to compute dFF
    
    Returns:
    - Tuple: Updated dFF segment, updated 'begin' index, and 'end' index
    """
    begin=0
    dFF_segment = np.full(len(data_df), np.nan)  # Create an array filled with NaNs of the same length as data_df
    artifact_intervals.append([len(data_df), 'End'])
    
    for x_start, x_stop in artifact_intervals:
        try:
            # Calculate 'end' as the first index where 'Time(s)' is greater than x_start, -1 to not overlap with artifact
            end = data_df.index[data_df['Time(s)'] < x_start][-1]
            
            # Extract the segment of data for processing
            segment = data_df.iloc[begin+1:end][col].values  # Use iloc for absolute indexing
            
            if method == 'mean':
                mean_fluorescence = np.nanmean(segment)
                dFF_values = ((segment - mean_fluorescence) / mean_fluorescence) * 100
                # Check for length match before assignment
                if len(dFF_values) == len(dFF_segment[begin+1:end]):
                    dFF_segment[begin+1:end] = dFF_values
                else:
                    print(f"Shape mismatch: dFF_values ({len(dFF_values)}) vs dFF_segment ({len(dFF_segment[begin+1:end])})")
            
            elif method == 'fit':
                dFF_values = linearfit_sklearn(data_df.iloc[begin+1:end]['405 Deinterleaved'].values, 
                                            data_df.iloc[begin+1:end]['465 Deinterleaved'].values)
                if len(dFF_values) == len(dFF_segment[begin+1:end]):
                    dFF_segment[begin+1:end] = dFF_values
                else:
                    print(f"Shape mismatch: dFF_values ({len(dFF_values)}) vs dFF_segment ({len(dFF_segment[begin+1:end])})")
            
            # Update 'begin' to the index just before the stop of the artifact
            if x_stop != 'End':
                begin = data_df.index[data_df['Time(s)'] > x_stop][0]
                
        except Exception as e:
            print(f"Error processing artifact interval ({x_start}, {x_stop}): {e}")
    
    return dFF_segment

def dFF(data_df, artifacts_df, filecode, method='mean'):
    """
    Calculates dFF (delta F over F) and removes artifacts from 405nm and 470nm photometry data.
    
    Parameters:
    - data_df (pd.DataFrame): Input photometry data containing 'Time(s)', '405 Deinterleaved', '470 Deinterleaved'
    - artifacts_df (pd.DataFrame): Dataframe containing artifact information
    - filecode (str): Unique identifier for the file being processed
    - sr (int): Sampling rate of the data
    - method (str): 'mean' or 'fit' method for calculating dFF
    
    Returns:
    - dFFdata_df (pd.DataFrame): DataFrame with 'Time(s)', '405 dFF', '470 dFF', and 'Denoised dFF'
    """
    dFFdata = np.full([3, len(data_df)], np.nan)  # [405 dFF, 470 dFF, Denoised dFF]

    if method == 'mean':
        for i, col in enumerate(['405 Deinterleaved', '465 Deinterleaved']):
            if filecode in artifacts_df['Filecode'].values:
                artifact_intervals = artifacts_df.loc[artifacts_df['Filecode'] == filecode, 'Artifacts'].values
                artifact_intervals = literal_eval(artifact_intervals[0])
                dFFdata[i] = remove_artifacts(data_df, artifact_intervals, col, method='mean')
            else:
                mean_fluorescence = np.nanmean(data_df[col])
                dFFdata[i] = ((data_df[col] - mean_fluorescence) / mean_fluorescence) * 100
    
    elif method == 'fit':
        if filecode in artifacts_df['Filecode'].values:
            artifact_intervals = artifacts_df.loc[artifacts_df['Filecode'] == filecode, 'Artifacts'].values
            artifact_intervals = literal_eval(artifact_intervals[0])
            dFFdata[0] = remove_artifacts(data_df, artifact_intervals, '405 Deinterleaved', method='fit')
            dFFdata[1] = data_df['465 Deinterleaved'].to_numpy()
        else:
            dFFdata[0] = linearfit_sklearn(data_df['405 Deinterleaved'], data_df['465 Deinterleaved'])
            dFFdata[1] = data_df['465 Deinterleaved'].to_numpy()

        # Calculate Denoised dFF
        dFFdata[2] = ((dFFdata[1] - dFFdata[0]) / dFFdata[0]) * 100

    dFFdata_df = pd.DataFrame({
        'Time(s)': data_df['Time(s)'],
        '405 dFF': dFFdata[0],
        '465 dFF': dFFdata[1],
        'Denoised dFF': dFFdata[2]
    })

    return dFFdata_df

def interpolate_dFFdata(data_df, method='linear'):
    """
    Interpolates or fills NaN values in dFF data.
    
    Parameters:
    -----------
    data_df : pd.DataFrame
        DataFrame containing columns for '405 dFF', '470 dFF', 'Denoised dFF' with NaN values if artifacts were removed
        
    method : str, optional
        The method to fill NaN values. Options:
        - 'linear' : Linearly interpolates NaN values between valid data points.
        - 'mean pad' : Replaces all NaN values with the mean of the respective column.
        
    Returns:
    --------
    data_df : pd.DataFrame
        The same DataFrame, but with NaN values filled.
    """
    
    if method not in ['linear', 'mean pad']:
        raise ValueError(f"Unsupported method: '{method}'. Choose 'linear' or 'mean pad'.")
    
    # Process only the dFF columns (ignoring 'Time(s)')
    dff_columns = data_df.columns[1:]  # Exclude 'Time(s)' column
    
    if method == 'linear':
        # Use interpolate with 'both' to ensure it fills NaNs at the beginning and end
        data_df[dff_columns] = data_df[dff_columns].interpolate(method='linear', limit_direction='both')
    
    elif method == 'mean pad':
        # Compute the mean of each column and fill NaNs with this mean
        col_means = data_df[dff_columns].mean(skipna=True)  # Mean of each dFF column, ignoring NaNs
        data_df[dff_columns] = data_df[dff_columns].fillna(col_means)
        
    return data_df

def butterfilt(data_df, order, cut_freq):
    """
    Applies a Butterworth low-pass filter to all columns (except 'Time(s)') in the input DataFrame.
    
    Parameters:
    -----------
    data_df : pd.DataFrame
        DataFrame containing 'Time(s)' and one or more data columns to be filtered.
        
    order : int
        The order of the Butterworth filter (e.g., 3, 4, etc.).
        
    cut_freq : float
        The cutoff frequency (in Hz) for the low-pass Butterworth filter.
    
    Returns:
    --------
    filtered_df : pd.DataFrame
        DataFrame containing the filtered data for each column (excluding 'Time(s)').
    """
    # Calculate sampling rate using the first two time points to avoid index issues
    sampling_rate = 1000 / (data_df['Time(s)'].iloc[1000] - data_df['Time(s)'].iloc[0])
    
    # Create the Butterworth filter coefficients
    sos = signal.butter(order, cut_freq, btype='low', output='sos', fs=sampling_rate)
    
    # Create a copy of the DataFrame to avoid modifying the original data
    filtered_df = data_df.copy()
    
    # Apply filter to all columns except 'Time(s)'
    data_columns = data_df.columns[1:]  # Exclude 'Time(s)'
    filtered_df[data_columns] = filtered_df[data_columns].apply(lambda col: signal.sosfilt(sos, col), axis=0)
    
    return filtered_df

def smoothing_SMA(data_df,win_size):
    """
    Simple moving average (SMA) smoothing
    win_size : size of the moving window
    """
    for col in data_df.columns[1:]:
        data_df[col] = data_df[col].rolling(win_size, min_periods=1).mean()
        
    return data_df