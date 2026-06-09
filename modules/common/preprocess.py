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
import warnings
from ast import literal_eval
import h5py
from sklearn.linear_model import LinearRegression, HuberRegressor
from statsmodels.nonparametric.smoothers_lowess import lowess

import modules.common.nomenclature as nom
import modules.common.median_filtering as mf


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

        try:
            # Data structure with laptop (Julien's stup)
            t_405 = f[base + "DeinterleavedDIO01/Time"][:]
            sig_405 = f[base + "DeinterleavedDIO01/AIN01"][:]

            t_465 = f[base + "DeinterleavedDIO02/Time"][:]
            sig_465 = f[base + "DeinterleavedDIO02/AIN01"][:]
        except KeyError:
            # Data structure in ephys room
            t_405 = f[base + "AIN01xDIO01-Deinterleaved/Time"][:]
            sig_405 = f[base + "AIN01xDIO01-Deinterleaved/Values"][:]

            t_465 = f[base + "AIN01xDIO02-Deinterleaved/Time"][:]
            sig_465 = f[base + "AIN01xDIO02-Deinterleaved/Values"][:]

    # Ensure all signals have same length
    min_len = min(len(sig_465), len(sig_405), len(t_465), len(t_405))
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

def load_lockin_dualcolor_doric(file_path):
    with h5py.File(file_path, 'r') as f:
        base = "DataAcquisition/FPConsole/Signals/Series0001/"

        try:
            # Data structure with laptop (Julien's setup)
            t_405 = f[base + "LockInAOUT01/Time"][:]
            sig_405 = f[base + "LockInAOUT01/AIN01"][:]

            t_465 = f[base + "LockInAOUT02/Time"][:]
            sig_465 = f[base + "LockInAOUT02/AIN01"][:]

            t_560 = f[base + "LockInAOUT03/Time"][:]
            sig_560 = f[base + "LockInAOUT03/AIN02"][:]

        except KeyError as e:
            raise RuntimeError(f"Could not load .doric file due to missing dataset(s): {e}")

    # Ensure all signals have same length
    min_len = min(len(sig_465), len(sig_405), len(sig_560), len(t_465), len(t_405), len(t_560))
    sig_465 = sig_465[:min_len]
    sig_405 = sig_405[:min_len]
    sig_560 = sig_560[:min_len]
    t_465 = t_465[:min_len]
    t_405 = t_405[:min_len]
    t_560 = t_560[:min_len]

    lockin_df = pd.DataFrame({
        'Time(s)': t_405,
        '405 Deinterleaved': sig_405,
        '465 Deinterleaved': sig_465,
        '560 Deinterleaved': sig_560
    })

    return lockin_df

def samplerate(data_df):
    
    time_diffs = np.diff(data_df['Time(s)'].dropna().values)
    time_diffs = time_diffs[time_diffs > 0]  # guard against duplicates
    sr = 1.0 / np.median(time_diffs)  # median is robust to outliers
    
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

def linearfit_sklearn(sig_405, sig_465, filt_405, filt_465, trim=[10, -10], filtered_405=False, model_name = 'linear'):
    '''
    Fit the 405nm isosbestic signal to a target channel (465nm or 560nm) using
    linear regression, and return the fitted baseline for dFF computation.

    The model is always trained on trimmed, median-filtered signals to avoid
    contamination from edge artifacts and transient events. Prediction can be
    made on either the raw or the filtered 405nm signal (see `filtered_405`).

    Parameters
    ----------
    sig_405 : array-like
        Raw 405nm isosbestic signal.
    sig_465 : array-like
        Raw target signal (465nm or 560nm).
    filt_405 : array-like
        Median-filtered 405nm signal, used for model training.
    filt_465 : array-like
        Median-filtered target signal, used for model training.
    trim : list of int, optional
        [start, stop] indices used to trim both filtered signals before fitting,
        removing edge artifacts. Default is [10, -10].
    filtered_405 : bool, optional
        If True, predict on filt_405 (smooth fitted baseline — recommended for
        weak signals such as 560nm, where 405nm noise would otherwise dominate).
        If False, predict on sig_405 (fitted baseline retains raw 405nm
        noise — acceptable for 465nm where signal amplitude is large).
        Default is False.
    smooth_405_sigma_s : smoothing option for 405nm signal
    model_name : str, optional
        Regression model to use:
        - 'linear' : ordinary least-squares (sklearn LinearRegression).
        - 'huber'  : robust Huber regression (sklearn HuberRegressor,
                     epsilon=1.35), which down-weights large residuals and is
                     less sensitive to remaining transients in the fitting window.
        Default is 'linear'.

    Returns
    -------
    fitted_405 : np.ndarray
        Fitted baseline at the scale of the target signal, same length as
        sig_405. Subtract from the raw target signal to isolate dFF:
        dFF = (target - fitted_405) / fitted_405 * 100.
    '''

    if model_name == 'linear':
        model = LinearRegression()
    elif model_name == 'huber':
        model = HuberRegressor(epsilon=1.35, max_iter=300)

    # Convert all inputs to numpy
    if not isinstance(sig_405,  np.ndarray): sig_405  = sig_405.to_numpy()
    if not isinstance(sig_465,  np.ndarray): sig_465  = sig_465.to_numpy()
    if not isinstance(filt_405, np.ndarray): filt_405 = filt_405.to_numpy()
    if not isinstance(filt_465, np.ndarray): filt_465 = filt_465.to_numpy()

    # Fit on trimmed filtered signals, predict on full raw signal
    model.fit(filt_405[trim[0]:trim[1]].reshape(-1, 1),
              filt_465[trim[0]:trim[1]])
    
    if filtered_405:
        fitted_405 = model.predict(filt_405.reshape(-1, 1))
    else:
        fitted_405 = model.predict(sig_405.reshape(-1, 1))

    return fitted_405

def remove_artifacts(data_df, filtered_data_df, artifact_intervals, col,
                     method='fit', filtered_405=False,
                     stitch_window=15, lowess_frac=0.1):
    """
    Remove disconnection artifacts from a signal column, segment by segment.

    Parameters
    ----------
    data_df : pd.DataFrame
    filtered_data_df : pd.DataFrame
        Median-filtered version of data_df (used only for method='fit').
    artifact_intervals : list of [start, stop]
    col : str
        Column to process.
    method : str
        'fit'    – regress against 405nm isosbestic (existing behaviour).
        'mean'   – normalise each segment by its own mean, then stitch.
        'lowess' – fit a LOWESS baseline per segment, then stitch.
    filtered_405 : bool
        Passed through to linearfit_sklearn (method='fit' only).
    stitch_window : int
        Samples at each segment edge used for boundary-median stitching
        (methods 'mean' and 'lowess').
    lowess_frac : float
        Fraction of data used for each local LOWESS fit (method='lowess').

    Returns
    -------
    np.ndarray  dFF signal with NaN over artifact gaps.
    """
    begin           = 0
    dFF_out     = np.full(len(data_df), np.nan)
    lowess_out  = np.full(len(data_df), np.nan)
    prev_end_median = None                            
    intervals       = artifact_intervals + [[len(data_df), 'End']]

    for x_start, x_stop in intervals:
        try:
            end     = data_df.index[data_df['Time(s)'] < x_start][-1]
            segment = data_df.iloc[begin+1:end][col].values

            if len(segment) == 0:
                if x_stop != 'End':
                    begin = data_df.index[data_df['Time(s)'] > x_stop][0]
                continue

            # ── Trim slice for baseline estimation ───────────────────────
            if begin == 0:
                trim = slice(10, None)
            elif x_stop == 'End':
                trim = slice(None, -10)
            else:
                trim = slice(None, None)

            F0 = np.full(len(segment), np.nan)
            # ── Compute dFF values for this segment ───────────────────────
            if method == 'mean':
                mean_f = np.nanmean(segment[trim])
                if mean_f == 0:
                    if x_stop != 'End':
                        begin = data_df.index[data_df['Time(s)'] > x_stop][0]
                    continue
                dFF_values = (segment - mean_f) / mean_f * 100

            elif method == 'lowess':
                time_seg = data_df.iloc[begin+1:end]['Time(s)'].values
                good     = ~np.isnan(segment)
                if good.sum() < 10:
                    if x_stop != 'End':
                        begin = data_df.index[data_df['Time(s)'] > x_stop][0]
                    continue
                # Fit LOWESS on trimmed data to avoid edge contamination,
                # then interpolate back to the full segment time axis
                t_trim   = time_seg[trim][~np.isnan(segment[trim])]
                f_trim   = segment[trim][~np.isnan(segment[trim])]
                fitted   = lowess(f_trim, t_trim,
                                  frac=lowess_frac, return_sorted=False)
                F0       = np.interp(time_seg, t_trim, fitted)
                with np.errstate(invalid='ignore', divide='ignore'):
                    dFF_values = np.where(F0 != 0,
                                          (segment - F0) / F0 * 100,
                                          np.nan)

            elif method == 'fit':
                # linearfit_sklearn expects a list trim, not a slice
                if begin == 0:
                    trim_list = [10, -1]
                elif x_stop == 'End':
                    trim_list = [0, -10]
                else:
                    trim_list = [0, -1]
                dFF_values = linearfit_sklearn(
                    data_df.iloc[begin+1:end]['405 Deinterleaved'].values,
                    data_df.iloc[begin+1:end][col].values,
                    filtered_data_df.iloc[begin+1:end]['405 Deinterleaved'].values,
                    filtered_data_df.iloc[begin+1:end][col].values,
                    trim=trim_list,
                    filtered_405=filtered_405)

            # ── Boundary stitching (mean and lowess only) ─────────────────
            if method in ('mean', 'lowess'):
                k = max(1, min(stitch_window, len(dFF_values) // 6))
                if prev_end_median is not None:
                    offset     = prev_end_median - np.nanmedian(dFF_values[:k])
                    dFF_values = dFF_values + offset
                prev_end_median = np.nanmedian(dFF_values[-k:])

            # ── Write to output ───────────────────────────────────────────
            target = dFF_out[begin+1:end]
            if len(dFF_values) == len(target):
                dFF_out[begin+1:end] = dFF_values
                lowess_out[begin+1:end] = F0
            else:
                print(f"Shape mismatch at ({x_start}, {x_stop}): "
                      f"{len(dFF_values)} vs {len(target)}")

            if x_stop != 'End':
                begin = data_df.index[data_df['Time(s)'] > x_stop][0]

        except Exception as e:
            print(f"Error processing artifact interval ({x_start}, {x_stop}): {e}")

    return dFF_out, lowess_out

def dFF(data_df, artifacts_df, filecode, method='fit', apply_median_filter = True):
    """
    Calculates dFF (delta F over F) and removes artifacts from 405nm and 465nm photometry data.
    
    Parameters:
    - data_df (pd.DataFrame): Input photometry data containing 'Time(s)', '405 Deinterleaved', '465 Deinterleaved'
    - artifacts_df (pd.DataFrame): Dataframe containing artifact information
    - filecode (str): Unique identifier for the file being processed
    - sr (int): Sampling rate of the data
    - method (str): 'mean' or 'fit' method for calculating dFF
    
    Returns:
    - dFFdata_df (pd.DataFrame): DataFrame with 'Time(s)', '405 Fitted', '465 Fitted', and 'dFF'
    """
    dFFdata = np.full([3, len(data_df)], np.nan)

    if method == 'mean':
        for i, col in enumerate(['405 Deinterleaved', '465 Deinterleaved']):
            if filecode in artifacts_df['Filecode'].values:
                artifact_intervals = artifacts_df.loc[artifacts_df['Filecode'] == filecode, 'Artifacts'].values
                artifact_intervals = literal_eval(artifact_intervals[0])
                dFFdata[i],_ = remove_artifacts(data_df, data_df, artifact_intervals, col, method='mean')
            else:
                mean_fluorescence = np.nanmean(data_df[col])
                dFFdata[i] = ((data_df[col] - mean_fluorescence) / mean_fluorescence) * 100
    
    elif method == 'fit':
        if apply_median_filter == True:
            # find best window from 465 nm and filter
            result_df, best_win_s, _ = mf.iterative_median_filter(data_df, '465 Deinterleaved',verbose=True)
            filtered_465 = result_df['465 Deinterleaved']
            # filter 405 nm with the same window
            filtered_405 = mf.median_filter_dff(data_df, '405 Deinterleaved', best_win_s)['405 Deinterleaved']
            filtered_data_df = pd.DataFrame({'Time(s)'           : data_df['Time(s)'].values,
                                                '465 Deinterleaved' : filtered_465.values,
                                                '405 Deinterleaved' : filtered_405.values})
        else:
            filtered_data_df = data_df.copy()

        if filecode in artifacts_df['Filecode'].values:
            artifact_intervals = artifacts_df.loc[artifacts_df['Filecode'] == filecode, 'Artifacts'].values
            artifact_intervals = literal_eval(artifact_intervals[0])
            dFFdata[0],_ = remove_artifacts(data_df, filtered_data_df, artifact_intervals, '465 Deinterleaved', method='fit')
            dFFdata[1] = data_df['465 Deinterleaved'].to_numpy()
        else:
            dFFdata[0] = linearfit_sklearn(data_df['405 Deinterleaved'], data_df['465 Deinterleaved'],
                                            filtered_data_df['405 Deinterleaved'], filtered_data_df['465 Deinterleaved'])
            dFFdata[1] = data_df['465 Deinterleaved'].to_numpy()

        # Calculate Denoised dFF
        dFFdata[2] = ((dFFdata[1] - dFFdata[0]) / dFFdata[0]) * 100

        # Replace first and last 10 frames with 1st quartile 
        q1 = np.nanpercentile(dFFdata[2], 25)
        dFFdata[2][:10]  = q1
        dFFdata[2][-10:] = q1

    dFFdata_df = pd.DataFrame({
        'Time(s)': data_df['Time(s)'],
        '405 Fitted': dFFdata[0],
        '465 Fitted': dFFdata[1],
        'dFF': dFFdata[2]
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

def dff_lowess_560(data_df, col, lowess_frac=0.1):
    """
    Compute dF/F for 560nm using a global LOWESS baseline with
    artifact masking.

    Strategy
    --------
    1. Mask artifact intervals as NaN in the raw signal.
    2. Fit LOWESS on the non-NaN points only → slow-varying F0.
    3. Interpolate F0 linearly across the masked gaps.
    4. dF/F = (F - F0) / F0 * 100 everywhere outside artifacts.

    Because F0 is a single continuous curve fitted over the whole
    recording, segment boundaries are naturally stitched — no DC
    jump, no per-segment centering.

    Parameters
    ----------
    lowess_frac : float
        Fraction of data used for each local fit. Larger = smoother.
        Start around 0.05–0.15 and tune visually.
    stitch : bool
        If True, also apply boundary-median offset correction as a
        safety net (useful when frac is small and F0 wiggles near edges).
    """
    
    time  = data_df['Time(s)'].values
    raw   = data_df[col].values.copy().astype(float)

    # ── 1. Mask artifacts (skipped entirely if list is empty) ────────
    mask = np.ones(len(raw), dtype=bool)

    # ── 2. LOWESS on non-NaN, non-artifact points ────────────────────
    good    = mask & ~np.isnan(raw)
    fitted  = lowess(raw[good], time[good], frac=lowess_frac, return_sorted=False)
    F0_full = np.interp(time, time[good], fitted)

    # ── 3. dF/F ──────────────────────────────────────────────────────
    with np.errstate(invalid='ignore', divide='ignore'):
        dFF = np.where(mask, (raw - F0_full) / F0_full * 100, np.nan)

    return dFF, F0_full

def dFF_dualcolor(data_df, artifacts_df, filecode, method_560='lowess', apply_median_filter=True):
    """
    Calculates dFF for dual-color fiber photometry (465nm + 560nm).

    Parameters
    ----------
    data_df : pd.DataFrame
        Must contain 'Time(s)', '405 Deinterleaved', '465 Deinterleaved', '560 Deinterleaved'.
    artifacts_df : pd.DataFrame
        Artifact information.
    filecode : str
        Unique file identifier.
    method_560 : str

        'lowess' : LOcally WEighted Scatterplot Smoothing
    apply_median_filter : bool
        If True, apply iterative hybrid median filter before fitting to suppress
        transient contamination of the isosbestic fit.
    """

    # ── Optional median filtering ─────────────────────────────────────────────
    if apply_median_filter:
        # Find best window from 465nm
        result_df, best_win_s, _ = mf.iterative_median_filter(
            data_df, '465 Deinterleaved', verbose=True)
        filtered_465 = result_df['465 Deinterleaved']

        # Filter 405nm with the same window
        filtered_405 = mf.median_filter_dff(
            data_df, '405 Deinterleaved', best_win_s)['405 Deinterleaved']

        filtered_data_df = pd.DataFrame({
            'Time(s)'           : data_df['Time(s)'].values,
            '405 Deinterleaved' : filtered_405.values,
            '465 Deinterleaved' : filtered_465.values,
        })

        if method_560 == 'fit':
            # Filter 560nm with the same window
            filtered_560 = mf.median_filter_dff(
                data_df, '560 Deinterleaved', best_win_s)['560 Deinterleaved']
            filtered_data_df['560 Deinterleaved'] = filtered_560.values
    else:
        filtered_data_df = data_df.copy()

    # ── fitted560 branch ──────────────────────────────────────────────────────
    if method_560 == 'fit':
        dFFdata = np.full([6, len(data_df)], np.nan)

        if filecode in artifacts_df['Filecode'].values:
            artifact_intervals = artifacts_df.loc[
                artifacts_df['Filecode'] == filecode, 'Artifacts'].values
            artifact_intervals = literal_eval(artifact_intervals[0])
            dFFdata[0],_ = remove_artifacts(data_df, filtered_data_df, artifact_intervals,
                                          '465 Deinterleaved', method='fit')
            dFFdata[1] = data_df['465 Deinterleaved'].to_numpy()
            dFFdata[2],_ = remove_artifacts(data_df, filtered_data_df, artifact_intervals,
                                          '560 Deinterleaved', method='fit', filtered_405=True)
            dFFdata[3] = data_df['560 Deinterleaved'].to_numpy()
        else:
            dFFdata[0] = linearfit_sklearn(
                data_df['405 Deinterleaved'],      data_df['465 Deinterleaved'],
                filtered_data_df['405 Deinterleaved'], filtered_data_df['465 Deinterleaved'])
            dFFdata[1] = data_df['465 Deinterleaved'].to_numpy()
            dFFdata[2] = linearfit_sklearn(
                data_df['405 Deinterleaved'],      data_df['560 Deinterleaved'],
                filtered_data_df['405 Deinterleaved'], filtered_data_df['560 Deinterleaved'],
                filtered_405=True)
            dFFdata[3] = data_df['560 Deinterleaved'].to_numpy()

        dFFdata[4] = ((dFFdata[1] - dFFdata[0]) / dFFdata[0]) * 100
        dFFdata[5] = ((dFFdata[3] - dFFdata[2]) / dFFdata[2]) * 100

        for row, label in [(4, '465'), (5, '560')]:
            q1 = np.nanpercentile(dFFdata[row], 25)
            dFFdata[row][:10]  = q1
            dFFdata[row][-10:] = q1

        dFFdata_df = pd.DataFrame({
            'Time(s)'      : data_df['Time(s)'],
            '405 Fitted'   : dFFdata[0],
            '465 Fitted'   : dFFdata[1],
            'dFF'          : dFFdata[4],
            '405 Fitted 560': dFFdata[2],
            '560 Fitted'   : dFFdata[3],
            '560 dFF'      : dFFdata[5],
        })

    # ── unfitted560 branch (mean) ────────────────────────────────────────────────────
    elif method_560 == 'mean':
        dFFdata = np.full([4, len(data_df)], np.nan)

        if filecode in artifacts_df['Filecode'].values:
            artifact_intervals = artifacts_df.loc[
                artifacts_df['Filecode'] == filecode, 'Artifacts'].values
            artifact_intervals = literal_eval(artifact_intervals[0])
            dFFdata[0],_ = remove_artifacts(data_df, filtered_data_df, artifact_intervals,
                                          '465 Deinterleaved', method='fit')
            dFFdata[1] = data_df['465 Deinterleaved'].to_numpy()
            dFFdata[2],_ = remove_artifacts(data_df, filtered_data_df, artifact_intervals,
                                          '560 Deinterleaved', method='mean')
        else:
            # 465: fit isosbestic (median-filtered if requested)
            dFFdata[0] = linearfit_sklearn(
                data_df['405 Deinterleaved'],          data_df['465 Deinterleaved'],
                filtered_data_df['405 Deinterleaved'], filtered_data_df['465 Deinterleaved'])
            dFFdata[1] = data_df['465 Deinterleaved'].to_numpy()
            # 560: normalize on its own mean (no fitting, median filter not applied)
            mean_fluorescence_560 = np.nanmean(data_df['560 Deinterleaved'][10:-10])
            dFFdata[2] = ((data_df['560 Deinterleaved'] - mean_fluorescence_560)
                          / mean_fluorescence_560) * 100

        dFFdata[3] = ((dFFdata[1] - dFFdata[0]) / dFFdata[0]) * 100

        for row in [3, 2]:
            q1 = np.nanpercentile(dFFdata[row], 25)
            dFFdata[row][:10]  = q1
            dFFdata[row][-10:] = q1

        dFFdata_df = pd.DataFrame({
            'Time(s)'    : data_df['Time(s)'],
            '405 Fitted' : dFFdata[0],
            '465 Fitted' : dFFdata[1],
            '560 dFF'    : dFFdata[2],
            'dFF'        : dFFdata[3],
        })

    elif method_560 == 'lowess':
        dFFdata = np.full([6, len(data_df)], np.nan)

        if filecode in artifacts_df['Filecode'].values:
            artifact_intervals = artifacts_df.loc[
                artifacts_df['Filecode'] == filecode, 'Artifacts'].values
            artifact_intervals = literal_eval(artifact_intervals[0])
            dFFdata[0],_ = remove_artifacts(data_df, filtered_data_df, artifact_intervals,
                                          '465 Deinterleaved', method='fit')
            dFFdata[1] = data_df['465 Deinterleaved'].to_numpy()
            dFFdata[2], dFFdata[3] = remove_artifacts(data_df, filtered_data_df, artifact_intervals, 
                                                      '560 Deinterleaved', method='lowess', 
                                                      stitch_window=15, lowess_frac=0.05)
        else:
            # 465: fit isosbestic (median-filtered if requested)
            dFFdata[0] = linearfit_sklearn(
                data_df['405 Deinterleaved'],          data_df['465 Deinterleaved'],
                filtered_data_df['405 Deinterleaved'], filtered_data_df['465 Deinterleaved'])
            dFFdata[1] = data_df['465 Deinterleaved'].to_numpy()
            # 560: lowess normalization
            dFFdata[2], dFFdata[3] = dff_lowess_560(data_df, '560 Deinterleaved', lowess_frac=0.05)

        dFFdata[4] = ((dFFdata[1] - dFFdata[0]) / dFFdata[0]) * 100

        for row in [2, 4]:
            q1 = np.nanpercentile(dFFdata[row], 25)
            dFFdata[row][:10]  = q1
            dFFdata[row][-10:] = q1

        dFFdata_df = pd.DataFrame({
            'Time(s)'      : data_df['Time(s)'],
            '405 Fitted'   : dFFdata[0],
            '465 Fitted'   : dFFdata[1],
            'dFF'          : dFFdata[4],
            '405 Fitted 560': dFFdata[3],
            '560 Fitted'   : data_df['560 Deinterleaved'].to_numpy(),
            '560 dFF'      : dFFdata[2],
        })

    return dFFdata_df

def downsample(rawdata_df, target_frequency=20):
    """
    Downsample a fiberphotometry DataFrame using an anti-aliased decimation filter.

    Parameters
    ----------
    rawdata_df : pd.DataFrame
        Must contain a 'Time(s)' column plus one or more signal columns.
    target_frequency : float
        Desired output sampling frequency in Hz.

    Returns
    -------
    pd.DataFrame
        Downsampled DataFrame with the same column layout.
    """
    # 1. Calculate current frequency from the 'Time' column
    # We use median diff to be robust against occasional dropped frames
    dt = np.median(np.diff(rawdata_df['Time(s)'].values))
    current_fs = 1.0 / dt
    print('Downsampling data')
    print(f'Current frequency : {current_fs}Hz. Target frequency : {target_frequency}Hz')

    # 2. Calculate integer downsampling factor (q)
    # decimate requires an integer; we round to the nearest whole number
    # if the rounded downsampling factor is too far from the exact factor, raises warning
    exact_factor = current_fs / target_frequency
    downsampling_factor = int(round(exact_factor))

    if abs(exact_factor - downsampling_factor) > 0.05:
        warnings.warn(
            f"Downsampling factor {exact_factor:.2f} rounded to {downsampling_factor}. "
            f"Effective output frequency: {current_fs / downsampling_factor:.2f} Hz "
            f"(target was {target_frequency} Hz)."
        )

    if downsampling_factor <= 1:
        print("Target frequency is higher than or equal to current frequency. Returning original.")
        return rawdata_df
    
    downsampled_data = {}
    # 3. Apply decimate to all columns except 'Time'
    # We iterate through all columns and skip 'Time' specifically
    signal_cols = [c for c in rawdata_df.columns if c != 'Time(s)']
    for col in signal_cols:
        # decimate applies an anti-aliasing low-pass filter before downsampling
        # we use zero-phase to avoid phase distortion
        downsampled_data[col] = signal.decimate(rawdata_df[col].values, 
                                               downsampling_factor, 
                                               ftype='fir',
                                               zero_phase=True)

    # 4. Reconstruct DataFrame
    downsampled_df = pd.DataFrame(downsampled_data)

    # 5. Re-create the Time column
    downsampled_df['Time(s)'] = rawdata_df['Time(s)'].values[::downsampling_factor][:len(downsampled_df)]

    # Reorder columns to put 'Time' first (standard for Doric/Fiber data)
    cols = ['Time(s)'] + [c for c in downsampled_df.columns if c != 'Time(s)']
    return downsampled_df[cols]

def apply_excluded_regions(df,signal_col,exclusion_df,filecode):
    """
    Replace manually excluded regions by NaNs.

    Parameters
    ----------
    df : pd.DataFrame

    signal_col : str
        Signal column to modify.

    exclusion_df : pd.DataFrame
        Must contain:
            - Filecode
            - Excluded

    filecode : str

    Returns
    -------
    df : pd.DataFrame
        Modified dataframe with:
            - signal replaced by NaNs
            - exclusion mask column
    """

    signal = df[signal_col].to_numpy().copy()

    exclusion_mask = np.zeros(len(signal), dtype=bool)

    # NO EXCLUSIONS
    if filecode not in exclusion_df['Filecode'].values:

        df[f'{signal_col} ExclusionMask'] = exclusion_mask.astype(int)

        return df

    # LOAD INTERVALS
    exclusion_intervals = exclusion_df.loc[
        exclusion_df['Filecode'] == filecode,
        'Artifacts'
    ].values[0]

    exclusion_intervals = literal_eval(exclusion_intervals)

    # BUILD MASK
    for start, stop in exclusion_intervals:
        idx = (
            (df['Time(s)'] >= start) &
            (df['Time(s)'] <= stop)
        )
        exclusion_mask[idx] = True
    # REPLACE WITH NaNs
    signal[exclusion_mask] = np.nan

    df[signal_col] = signal
    df[f'{signal_col} ExclusionMask'] = (
        exclusion_mask.astype(int)
    )
    return df