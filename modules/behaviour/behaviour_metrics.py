import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from scripts.loader import project_root

def ensure_behavior_columns(df, required_behaviors):
    """
    Ensures all required behavior columns exist in df.
    If missing, they are added and filled with zeros.
    """
    for col in required_behaviors:
        if col not in df.columns:
            print(f"[INFO] Column '{col}' missing — creating empty column.")
            df[col] = 0
    return df

def compute_behavior_metrics(df, bin_size):
    # Make sure all expected behaviors exist
    required_cols = ["Airpuffs", "Licks_filtered", "Nose_in_any_airport"]
    df = ensure_behavior_columns(df, required_cols)

    df["time_bin"] = (df["Time(s)"] // bin_size).astype(int)

    metrics = df.groupby("time_bin").agg(
        airpuffs_count=("Airpuffs", "sum"),
        licking_time=("Licks_filtered",
                      lambda x: x.sum() * df["Time(s)"].diff().median()),
        nose_in_airport_time=("Nose_in_any_airport",
                              lambda x: x.sum() * df["Time(s)"].diff().median())
    ).reset_index()

    return metrics

def plot_behavior_raster(df, mouse, batch, behaviors, save_dir=None):
    """
    Create a raster plot where each behavior (binary column) is a row.
    A tick is drawn at the time of each bout (value==1).
    """
    import json
    behavior_colors_path = Path(project_root) / "modules/behaviour/behaviour_colors.json"

    # Load colors from JSON
    with open(behavior_colors_path, "r") as f:
        behavior_colors = json.load(f)

    # Create save directory
    if save_dir is not None:
        save_dir.mkdir(parents=True, exist_ok=True)

    # Prepare figure
    plt.figure(figsize=(14, 0.7 * len(behaviors) + 2))
    ax = plt.gca()

    time = df["Time(s)"].values

    # Raster plot: each behavior draws ticks at y = behavior index
    for i, behav in enumerate(behaviors):
        if behav not in df.columns:
            print(f"[!] Behavior '{behav}' not found in dataframe — skipping.")
            continue

        # Extract the 1's (bouts)
        bout_idx = np.where(df[behav].values == 1)[0]
        bout_times = time[bout_idx]

        # Color from json, fallback to grey
        color, alpha = behavior_colors.get(behav, ("grey", 0.8))

        # Plot ticks (vertical lines)
        ax.vlines(
            bout_times,
            i + 0.1,      # small height for visual clarity
            i + 0.9,
            color=color,
            alpha=alpha,
            linewidth=2
        )

    # Format axes
    ax.set_yticks(np.arange(len(behaviors)) + 0.5)
    ax.set_yticklabels(behaviors)

    ax.set_xlabel("Time (s)")
    ax.set_title(f"Behavior Raster Plot — Mouse {mouse} - Batch {batch}")

    plt.tight_layout()

    # Save if required
    if save_dir is not None:
        png_path = save_dir / f"{batch}_{mouse}_behavior_raster.png"
        pdf_path = save_dir / f"{batch}_{mouse}_behavior_raster.pdf"
        plt.savefig(png_path, dpi=300)
        plt.savefig(pdf_path)
        print(f"Saved raster plot for {mouse} to:\n  {png_path}\n  {pdf_path}")

    plt.show()

def plot_behavior_metrics(metrics_df, mouse, batch, bin_size, save_dir=None):
    """
    Plot behavioral metrics and save PNG/PDF if save_dir is provided.
    """
    # Ensure expected columns
    expected_cols = ["airpuffs_count", "licking_time", "nose_in_airport_time"]
    for col in expected_cols:
        if col not in metrics_df.columns:
            print(f"[INFO] Metrics column '{col}' missing — filling with zeros.")
            metrics_df[col] = 0

    # Create save directory
    if save_dir is not None:
        save_dir.mkdir(parents=True, exist_ok=True)

    # ---- Plot ----
    plt.figure(figsize=(12,4))
    plt.plot(metrics_df["time_bin"] * bin_size, metrics_df["airpuffs_count"], label="Airpuffs")
    plt.plot(metrics_df["time_bin"] * bin_size, metrics_df["licking_time"], label="Licking (s)")
    plt.plot(metrics_df["time_bin"] * bin_size, metrics_df["nose_in_airport_time"], label="Nose in Airport (s)")
    
    plt.title(f"Behavior Metrics — Mouse {mouse} - Batch {batch}")
    plt.xlabel("Time (s)")
    plt.legend()
    plt.tight_layout()

    # ---- Save ----
    if save_dir is not None:
        png_path = save_dir / f"{batch}_{mouse}_behavior_metrics.png"
        pdf_path = save_dir / f"{batch}_{mouse}_behavior_metrics.pdf"
        plt.savefig(png_path, dpi=300)
        plt.savefig(pdf_path)
        print(f"Saved behavioral metrics plot for {mouse} to:\n  {png_path}\n  {pdf_path}")

    plt.show()

def compute_and_plot_heatmap(df, mouse, batch, bins=(50,50), n_bins=1, save_dir=None):
    """
    Compute + plot occupancy heatmap, save PNG/PDF if save_dir given.
    """
    if save_dir is not None:
        save_dir.mkdir(parents=True, exist_ok=True)

    # -------------------------
    # Split session into equal bins
    # -------------------------
    total_len = len(df)
    bin_len = total_len // n_bins
    dfs = [df.iloc[i*bin_len : (i+1)*bin_len] for i in range(n_bins)]
    dfs[-1] = df.iloc[(n_bins-1)*bin_len :]  # include leftovers

    # -------------------------
    # Figure layout: 2 × n_bins
    # -------------------------
    fig, axes = plt.subplots(
        2, n_bins,
        figsize=(4*n_bins, 8),
        gridspec_kw={'height_ratios': [1, 4]}
    )

    for i, subdf in enumerate(dfs):
        x = subdf["center_x"].values
        y = subdf["center_y"].values

        # -------------------------
        # 1) Trajectory
        # -------------------------
        ax_traj = axes[0, i]
        ax_traj.plot(x, y, color="black", linewidth=1)
        ax_traj.set_title(f"Bin {i+1}", fontsize=10)

        # Remove ticks & labels
        ax_traj.set_xticks([])
        ax_traj.set_yticks([])
        ax_traj.set_xlabel("")
        ax_traj.set_ylabel("")

        # Maintain exact aspect ratio
        ax_traj.set_aspect('equal', adjustable='box')

        # -------------------------
        # 2) Heatmap
        # -------------------------
        ax_hm = axes[1, i]
        heatmap, _, _ = np.histogram2d(x, y, bins=bins)

        sns.heatmap(
            heatmap.T,
            cmap="inferno",
            ax=ax_hm,
            cbar=True,
            square=False,
            vmin=0,
            vmax=100,
            cbar_kws={"shrink": 0.35, "pad": 0.02}
        )

        # Remove ticks & labels
        ax_hm.set_xticks([])
        ax_hm.set_yticks([])
        ax_hm.set_xlabel("")
        ax_hm.set_ylabel("")

        # Maintain equal aspect ratio for heatmap
        ax_hm.set_aspect('equal', adjustable='box')

    plt.suptitle(f"Occupancy Heatmaps — Mouse {mouse}", fontsize=16)
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    # ---- Save ----
    if save_dir is not None:
        png_path = save_dir / f"{batch}_{mouse}_heatmap_multibin.png"
        pdf_path = save_dir / f"{batch}_{mouse}_heatmap_multibin.pdf"
        fig.savefig(png_path, dpi=300)
        fig.savefig(pdf_path)
        print(f"[✔] Multi-bin heatmap saved for {mouse}:\n  {png_path}\n  {pdf_path}")

    plt.show()
