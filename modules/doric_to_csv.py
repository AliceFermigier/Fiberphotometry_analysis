import h5py
import pandas as pd
import numpy as np

def doric_to_csv(doric_path, csv_path):

    # Open .doric (HDF5-based format)
    with h5py.File(doric_path, 'r') as f:
        # Navigate into the main group (usually 'Data')
        group = f['Data'] if 'Data' in f else f[list(f.keys())[0]]
        signals = {}
        time = None

        # Loop through available signals
        for key in group:
            dataset = group[key]
            if 'Values' in dataset:
                values = dataset['Values'][()]
                label = dataset.attrs.get('Label', key).decode() if isinstance(dataset.attrs.get('Label', key), bytes) else dataset.attrs.get('Label', key)
                signals[label] = values

            if time is None and 'Time' in dataset:
                time = dataset['Time'][()]

        if time is None:
            raise ValueError("Time vector not found in .doric file.")

        # Build DataFrame
        df = pd.DataFrame({'Time': time})
        for label, values in signals.items():
            df[label] = values

    # Save to CSV
    df.to_csv(csv_path, index=False)
    print(f"Saved: {csv_path}")
    return