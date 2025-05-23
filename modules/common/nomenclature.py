# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 00:48:26 2024

Functions for nomenclature

@author: alice
"""

import os
import pandas as pd
import shutil

def create_directory(path):
    """Create a directory if it doesn't already exist."""
    os.makedirs(path, exist_ok=True)

def get_experiment_sessions(batches, proto_df, exp):
    """
    Get a list of session names for a given experiment from the protocol DataFrame.
    
    Parameters:
    - batches (list): List of included batches of mice.
    - proto_df (pd.DataFrame): DataFrame containing protocol information.
    - exp (str): Name of the experiment.
    
    Returns:
    - dict: Dictionary mapping batch names to lists of session names.
    """
    sessions_list = []
    
    for B in batches:
        try:
            filtered_df = proto_df[(proto_df['Task'] == exp) & (proto_df['Batch'] == B)]
            
            if filtered_df.empty:
                raise ValueError(f"No sessions found for experiment '{exp}' in batch '{B}'.")
            
            sessions_list.append(filtered_df['Sessions'].values[0])  # Extract session value
        
        except KeyError as e:
            raise ValueError(f"Missing column in DataFrame: {e}")
        except IndexError:
            raise ValueError(f"Sessions not found for experiment '{exp}' in batch '{B}'.")
    
    return sessions_list

def get_experiment_data_path(batches, proto_df, data_path, exp):
    """
    Get the raw data path for the experiment from the protocol DataFrame.
    
    Parameters:
    - batches : included batches of mice
    - proto_df (pd.DataFrame): DataFrame containing protocol information.
    - data_path (Path): Base path to data storage.
    - exp (str): Name of the experiment.
    
    Returns:
    - dictionary of paths to the experiment's data directories
    """
    datapathexp_dict = {B: None for B in batches}
    
    for B in batches:
        try:
            filtered_df = proto_df[(proto_df['Task'] == exp) & (proto_df['Batch'] == B)]
            
            if filtered_df.empty:
                raise ValueError(f"No sessions found for experiment '{exp}' in batch '{B}'.")
            
            relative_path = filtered_df['Data_path'].values[0]
            
            datapathexp_dict[B] = data_path / relative_path

        except (KeyError, IndexError):
            raise ValueError(f"Data path not found for experiment '{exp}' in the protocol DataFrame.")
            
    return datapathexp_dict


def setup_experiment_directory(analysis_path, exp):
    """
    Set up the main experiment directory.
    
    Parameters:
    - analysis_path (Path): Base path for analysis storage.
    - exp (str): Name of the experiment.
    
    Returns:
    - Path: Full path to the experiment directory.
    """
    exp_path = analysis_path / exp
    create_directory(exp_path)
    
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

def create_or_load_artifacts_file(file_path, option='create_or_load'):
    """
    Create an Excel file with 'Filecode' and 'Artifacts' columns if it doesn't exist.
    If it exists, load the file as a DataFrame.
    
    Parameters:
    - file_path (Path): Path to the artifacts Excel file.
    
    Returns:
    - pd.DataFrame: DataFrame containing the existing or newly created file.
    """
    if not file_path.exists():
        # Create a new DataFrame with required columns
        df = pd.DataFrame(columns=['Filecode', 'Artifacts'])
        df.to_excel(file_path, index=False)
    elif option=='create_or_load':
        df = pd.read_excel(file_path)
    else:
        df = None
    return df

def create_or_load_sniffs_file(file_path, option='create_or_load'):
    """
    Create an Excel file with ['Subject', 'Group', 'Odor', 'Count', 'Stim', 'Sniffs'] columns if it doesn't exist.
    If it exists, load the file as a DataFrame.
    
    Parameters:
    - file_path (Path): Path to the artifacts Excel file.
    
    Returns:
    - pd.DataFrame: DataFrame containing the existing or newly created file.
    """
    if not file_path.exists():
        # Create a new DataFrame with required columns
        df = pd.DataFrame(columns=['Subject', 'Group', 'Odor', 'Count', 'Stim', 'Sniffs'])
        df.to_excel(file_path, index=False)
    elif option=='create_or_load':
        df = pd.read_excel(file_path)
    else:
        df = None
    return df
            
def move_behav_files(data_path, behav_path):
    # Créer le dossier cible s'il n'existe pas
    os.makedirs(behav_path, exist_ok=True)

    for filename in os.listdir(data_path):
        file_path = os.path.join(data_path, filename)

        # Vérifier si le fichier correspond au modèle "behav_*"
        if os.path.isfile(file_path) and filename.startswith("behav_"):
            # Déplacer le fichier vers le dossier cible
            shutil.move(file_path, os.path.join(behav_path, filename))            

def correction_behav_files(directory):
    for filename in os.listdir(directory):
        file_path = os.path.join(directory, filename)

        # Ignorer les dossiers et les fichiers déjà nommés correctement
        if os.path.isfile(file_path) and not filename.startswith("behav_"):
            # Extraire la partie centrale du nom du fichier
            filename_base = os.path.splitext(filename)[0]
            parts = filename_base.split('_')
            
            # Trouver le coden ou générer le code si oubli
            if parts[-3] not in ['0','1','2','3']:
                code = '0'
            else:
                code = parts[-3]

            # Vérifier si le nom de la souris est répété
            if parts[-1] != parts[-2]:
                # Supprimer le fichier non correspondant
                os.remove(file_path)   
            elif parts[-1] == parts[-2]:
                # Changer le nom du fichier
                new_name = f"behav_{code}_{parts[-1]}.csv"
                new_path = os.path.join(directory, new_name)
                os.rename(file_path, new_path)
                
