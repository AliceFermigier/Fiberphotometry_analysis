import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns

def compute_behavior_metrics(df, bin_size):
    """
    df: fiberbehav aligned dataframe
    """

    # Create time bins
    df["time_bin"] = (df["Time(s)"] // bin_size).astype(int)

    # Aggregate behavior metrics inside bins
    grouped = df.groupby("time_bin").agg(
        airpuffs_count = ("Airpuffs", "sum"),
        licking_time = ("Licks_filtered", lambda x: x.sum() * (df["Time(s)"].diff().median())),  # convert frames to ms
        nose_in_airport_time = ("Nose_in_any_airport", lambda x: x.sum() * df["Time(s)"].diff().median())
    ).reset_index()

    return grouped

def plot_behavior_metrics(metrics_df, mouse, bin_size, save_dir=None):
    """
    Plot behavioral metrics and save PNG/PDF if save_dir is provided.
    """
    # Create save directory
    if save_dir is not None:
        save_dir.mkdir(parents=True, exist_ok=True)

    # ---- Plot ----
    plt.figure(figsize=(12,4))
    plt.plot(metrics_df["time_bin"] * bin_size, metrics_df["airpuffs_count"], label="Airpuffs")
    plt.plot(metrics_df["time_bin"] * bin_size, metrics_df["licking_time"], label="Licking (s)")
    plt.plot(metrics_df["time_bin"] * bin_size, metrics_df["nose_in_airport_time"], label="Nose in Airport (s)")
    
    plt.title(f"Behavior Metrics — Mouse {mouse}")
    plt.xlabel("Time (s)")
    plt.legend()
    plt.tight_layout()

    # ---- Save ----
    if save_dir is not None:
        png_path = save_dir / f"{mouse}_behavior_metrics.png"
        pdf_path = save_dir / f"{mouse}_behavior_metrics.pdf"
        plt.savefig(png_path, dpi=300)
        plt.savefig(pdf_path)
        print(f"Saved behavioral metrics plot for {mouse} to:\n  {png_path}\n  {pdf_path}")

    plt.show()

def compute_and_plot_heatmap(df, mouse, bins=(50,50), save_dir=None):
    """
    Compute + plot occupancy heatmap, save PNG/PDF if save_dir given.
    """
    if save_dir is not None:
        save_dir.mkdir(parents=True, exist_ok=True)

    x = df["center_x"]
    y = df["center_y"]

    heatmap, xedges, yedges = np.histogram2d(x, y, bins=bins)

    plt.figure(figsize=(6,5))
    sns.heatmap(heatmap.T, cmap="inferno", cbar=True)
    plt.title(f"Occupancy Heatmap — Mouse {mouse}")
    plt.xlabel("X bin")
    plt.ylabel("Y bin")
    plt.tight_layout()

    # ---- Save ----
    if save_dir is not None:
        png_path = save_dir / f"{mouse}_heatmap.png"
        pdf_path = save_dir / f"{mouse}_heatmap.pdf"
        plt.savefig(png_path, dpi=300)
        plt.savefig(pdf_path)
        print(f"Saved occupancy heatmap for {mouse} to:\n  {png_path}\n  {pdf_path}")

    plt.show()
