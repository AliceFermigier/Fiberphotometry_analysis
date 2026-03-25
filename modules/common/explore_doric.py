"""
explore_doric.py
Prints the full nomenclature (groups, datasets, attributes) of a .doric file.
Usage: python explore_doric.py path/to/your_file.doric
"""

import sys
import h5py
import numpy as np

def print_structure(name, obj):
    """Visitor callback: prints every group and dataset."""
    indent = "  " * name.count("/")
    if isinstance(obj, h5py.Group):
        print(f"{indent}📁 {name}/")
    elif isinstance(obj, h5py.Dataset):
        shape = obj.shape
        dtype = obj.dtype
        print(f"{indent}📄 {name}  |  shape={shape}  dtype={dtype}")

    # Print attributes if any
    if obj.attrs:
        for key, val in obj.attrs.items():
            # Truncate long array attributes for readability
            if isinstance(val, np.ndarray) and val.size > 6:
                display = f"{val[:3]} ... {val[-3:]} (len={val.size})"
            else:
                display = val
            print(f"{indent}    🏷  {key}: {display}")


def explore(filepath):
    print(f"\n{'='*60}")
    print(f"  File: {filepath}")
    print(f"{'='*60}\n")

    with h5py.File(filepath, "r") as f:
        # Root-level attributes
        if f.attrs:
            print("📌 Root attributes:")
            for key, val in f.attrs.items():
                print(f"    🏷  {key}: {val}")
            print()

        # Walk the full tree
        f.visititems(print_structure)

    print(f"\n{'='*60}")
    print("  Done.")
    print(f"{'='*60}\n")