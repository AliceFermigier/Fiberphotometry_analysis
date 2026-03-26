##########
#IMPORTED#
##########

import pandas as pd
import numpy as np

import modules.common.transients as tr
import modules.common.preprocess as pp

def extract_dff_summary(fiberbehav_df, mouse, batch, group,
                         zone_cols, behav_cols,
                         merge_into=None,
                         dff_col='dFF',
                         fps=20,
                         use_zscore=False):
    """
    Extract mean dFF and AUC dFF during each zone/behaviour for one animal.

    Parameters
    ----------
    fiberbehav_df : pd.DataFrame
        Combined fiber + behavior dataframe (fiberbehavnotderived.csv).
    mouse, batch, group : str
        Animal identifiers.
    zone_cols : list of str
        Binary zone columns. Time spent in each → mean dFF + AUC dFF.
        E.g. ['Closed arm', 'Open arm', 'Center'].
    behav_cols : list of str
        Binary behaviour columns. Same metrics computed independently.
        E.g. ['Head dipping'].
    merge_into : dict, optional
        {behaviour: zone} pairs where frames of a behaviour should ALSO be
        included in the zone's signal. 
        E.g. {'Head dipping': 'Open arm'} ensures head dipping frames
        are counted in Open arm dFF in addition to their own column.
    dff_col : str
        Column name for raw dFF signal.
    fps : float
        Frame rate, used to convert frame count to seconds for AUC 
        (AUC = sum(dFF) / fps, i.e. dFF integrated over time in seconds).
    use_zscore : bool
        If True, z-score the full dFF trace before extracting metrics.

    Returns
    -------
    record : dict
        Flat dictionary with mean dFF and AUC dFF per zone/behaviour.
    """
    merge_into = merge_into or {}

    # ── Signal preparation ────────────────────────────────────────────────────
    if dff_col not in fiberbehav_df.columns:
        raise ValueError(f"Column '{dff_col}' not found in dataframe.")

    signal = fiberbehav_df[dff_col].values.copy()

    if use_zscore:
        signal = (signal - np.mean(signal)) / np.std(signal, ddof=1)

    record = {
        'Mouse' : mouse,
        'Batch' : batch,
        'Group' : group,
    }

    # ── Helper: compute metrics for a boolean mask ────────────────────────────
    def _metrics(mask, label):
        if mask.sum() == 0:
            print(f"  Warning: no frames found for '{label}', filling with NaN.")
            record[f'{label} mean dFF']  = float('nan')
            record[f'{label} AUC dFF']   = float('nan')
            return
        sig_masked             = signal[mask]
        record[f'{label} mean dFF'] = round(float(np.mean(sig_masked)), 6)
        # AUC: integral over time = sum(dFF * dt) where dt = 1/fps
        record[f'{label} AUC dFF']  = round(float(np.sum(sig_masked) / fps), 6)

    # ── Zone metrics ──────────────────────────────────────────────────────────
    # Build zone masks, expanding with merged behaviours where specified
    for col in zone_cols:
        if col not in fiberbehav_df.columns:
            print(f"  Warning: zone column '{col}' not found, filling with NaN.")
            record[f'{col} mean dFF'] = float('nan')
            record[f'{col} AUC dFF']  = float('nan')
            continue

        mask = fiberbehav_df[col].values.astype(bool)

        # Merge any behaviours that should count toward this zone
        # e.g. head dipping frames → also included in Open arm
        for behav, target_zone in merge_into.items():
            if target_zone == col and behav in fiberbehav_df.columns:
                behav_mask = fiberbehav_df[behav].values.astype(bool)
                mask = mask | behav_mask  # union

        _metrics(mask, col)

    # ── Behaviour metrics (independent, no merging) ───────────────────────────
    for col in behav_cols:
        if col not in fiberbehav_df.columns:
            print(f"  Warning: behaviour column '{col}' not found, filling with NaN.")
            record[f'{col} mean dFF'] = float('nan')
            record[f'{col} AUC dFF']  = float('nan')
            continue

        mask = fiberbehav_df[col].values.astype(bool)
        _metrics(mask, col)

    return record

def variance_transients(fiberbehav_df, list_BOI, mouse, group, exp, batch, threshold):
    """
    Calculates variance, transient frequency, and amplitude during whole trace and during behaviours
    
    Output:
    - A dataframe with variance, transient frequency and amplitude.
    """
    
    # Calculate variance during whole trace, baseline and post-baseline periods
    variance = np.var(fiberbehav_df['dFF'])
    
    # Calculate transients for whole trace, baseline and post-baseline periods
    peaks_df, peak_frequency, peak_amplitude, transients_fig = tr.transients(fiberbehav_df, threshold)

    # Store the results in a dataframe
    results_df = pd.DataFrame({
        'Batch': batch,
        'Subject': mouse,
        'Group': group,
        'Variance': variance,
        'Transients Frequency': peak_frequency,
        'Transients Amplitude': peak_amplitude,
    }, index=[0])

    if list_BOI != []:
        for behavior in list_BOI:
            if behavior not in fiberbehav_df.columns:
                print(f"Warning: Behavior '{behavior}' not found in DataFrame.")
                continue

            # Find rows where behavior is active
            behavior_mask = fiberbehav_df[behavior] == 1

            # Peaks during that behavior
            behavior_peaks = peaks_df[(peaks_df['Peaks'] == 1) & behavior_mask]

            # Duration of behavior period in seconds
            behavior_time = fiberbehav_df.loc[behavior_mask, 'Time(s)']
            if behavior_time.empty:
                duration = np.nan
            else:
                duration = behavior_time.iloc[-1] - behavior_time.iloc[0]

            # Calculate frequency and mean amplitude during behavior
            if duration and len(behavior_peaks) > 0:
                freq = len(behavior_peaks) / duration
                amp = behavior_peaks['Filtered dFF'].mean()
            else:
                freq = np.nan
                amp = np.nan

            # Add to results
            results_df[f'{behavior} Transients Frequency'] = freq
            results_df[f'{behavior} Transients Amplitude'] = amp
    
    return results_df, transients_fig