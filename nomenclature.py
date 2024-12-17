# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 00:48:26 2024

Functions for nomenclature

@author: alice
"""

import os

def create_directory(path):
    """Create a directory if it doesn't already exist."""
    os.makedirs(path, exist_ok=True)


def get_experiment_sessions(proto_df, exp):
    """
    Get a list of session names for a given experiment from the protocol DataFrame.
    
    Parameters:
    - proto_df (pd.DataFrame): DataFrame containing protocol information.
    - exp (str): Name of the experiment.
    
    Returns:
    - list: List of session names for the experiment.
    """
    try:
        sessions = proto_df.loc[proto_df['Task'] == exp, 'Sessions'].values[0]
        return sessions.split()
    except (KeyError, IndexError):
        raise ValueError(f"Sessions not found for experiment '{exp}' in the protocol DataFrame.")


def get_experiment_data_path(proto_df, data_path, exp):
    """
    Get the raw data path for the experiment from the protocol DataFrame.
    
    Parameters:
    - proto_df (pd.DataFrame): DataFrame containing protocol information.
    - data_path (Path): Base path to data storage.
    - exp (str): Name of the experiment.
    
    Returns:
    - Path: Full path to the experiment's data directory.
    """
    try:
        relative_path = proto_df.loc[proto_df['Task'] == exp, 'Data_path'].values[0]
        return data_path / relative_path
    except (KeyError, IndexError):
        raise ValueError(f"Data path not found for experiment '{exp}' in the protocol DataFrame.")


def setup_experiment_directory(analysis_path, exp, session_names):
    """
    Set up the main experiment directory and create subdirectories for each session.
    
    Parameters:
    - analysis_path (Path): Base path for analysis storage.
    - exp (str): Name of the experiment.
    - session_names (list): List of session names for the experiment.
    
    Returns:
    - Path: Full path to the experiment directory.
    """
    exp_path = analysis_path / exp
    create_directory(exp_path)
    
    for session_name in session_names:
        session_path = exp_path / session_name
        create_directory(session_path)
    
    return exp_path


def setup_preprocessing_directory(data_path_exp):
    """
    Create a preprocessing directory inside the raw data path for the experiment.
    
    Parameters:
    - data_path_exp (Path): Path to the experiment's raw data directory.
    
    Returns:
    - Path: Full path to the preprocessing directory.
    """
    pp_path = data_path_exp / 'Preprocessing'
    create_directory(pp_path)
    return pp_path
