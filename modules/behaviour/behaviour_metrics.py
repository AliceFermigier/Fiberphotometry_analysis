import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import json
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
        png_path = save_dir / f"{batch}_{mouse}_behavior_metrics_binsize{bin_size}s.png"
        pdf_path = save_dir / f"{batch}_{mouse}_behavior_metrics_binsize{bin_size}s.pdf"
        plt.savefig(png_path, dpi=300)
        plt.savefig(pdf_path)
        print(f"Saved behavioral metrics plot for {mouse} to:\n  {png_path}\n  {pdf_path}")

    plt.show()

def compute_and_plot_heatmap(df, mouse, batch, ports_json, arena_json,
                             n_bins=1, bins=(50, 50), save_dir=None):

    PAD = 100  # px padding around arena

    # ----------------------------------------------
    # Load ports & arena rectangle
    # ----------------------------------------------
    with open(ports_json, "r") as f:
        ports = json.load(f)

    with open(arena_json, "r") as f:
        arena = json.load(f)

    rect = arena["Arena_rectangle_px"]
    x1 = rect["x1"] - PAD
    y1 = rect["y1"] - PAD
    x2 = rect["x2"] + PAD
    y2 = rect["y2"] + PAD

    # Rectangle outline
    arena_x = [x1, x2, x2, x1, x1]
    arena_y = [y1, y1, y2, y2, y1]

    # Limits for both plots
    x_min, x_max = x1, x2
    y_min, y_max = y1, y2

    # ----------------------------------------------
    # Split session into bins
    # ----------------------------------------------
    total_len = len(df)
    bin_len = total_len // n_bins
    dfs = [df.iloc[i*bin_len : (i+1)*bin_len] for i in range(n_bins)]
    dfs[-1] = df.iloc[(n_bins-1)*bin_len:]   # last bin takes remainder

    # ----------------------------------------------
    # Create figure
    # ----------------------------------------------
    fig, axes = plt.subplots(
        2, n_bins,
        figsize=(4*n_bins, 8),
        gridspec_kw={"height_ratios": [1, 4]}
    )

    for i, subdf in enumerate(dfs):

        x = subdf["center_x"].values
        y = subdf["center_y"].values

        # ------------------------------------------------------
        # Trajectory panel
        # ------------------------------------------------------
        ax_t = axes[0, i]
        ax_t.plot(x, y, color="black", linewidth=1)

        # Arena outline
        ax_t.plot(arena_x, arena_y, color="whitesmoke", linewidth=1)

        # Ports
        ax_t.scatter(ports["lick_port"]["x"], ports["lick_port"]["y"], 
                     c="lime", s=30)
        ax_t.scatter(ports["airpuff_left"]["x"], ports["airpuff_left"]["y"],
                     c="red", s=30)
        ax_t.scatter(ports["airpuff_right"]["x"], ports["airpuff_right"]["y"],
                     c="red", s=30)

        ax_t.set_xlim(x_min, x_max)
        ax_t.set_ylim(y_min, y_max)
        ax_t.set_aspect("equal")
        ax_t.set_xticks([])
        ax_t.set_yticks([])
        ax_t.set_title(f"Bin {i+1}")

        # ------------------------------------------------------
        # Heatmap panel
        # ------------------------------------------------------

        # 🔥 CRUCIAL FIX: histogram must use arena rectangle range
        heatmap, xedges, yedges = np.histogram2d(
            x, y,
            bins=bins,
            range=[[x_min, x_max], [y_min, y_max]]
        )

        ax_h = axes[1, i]

        # Construct extent manually (sns does not handle this)
        extent = [x_min, x_max, y_min, y_max]

        ax_h.imshow(
            heatmap.T,
            origin="lower",
            cmap="jet",
            vmin=0,
            vmax=60,
            extent=extent,
            interpolation="nearest",
            aspect="equal"
        )

        # Arena outline + ports
        ax_h.plot(arena_x, arena_y, color="black", linewidth=1)
        ax_h.scatter(ports["lick_port"]["x"], ports["lick_port"]["y"], 
                     c="lime", s=30)
        ax_h.scatter(ports["airpuff_left"]["x"], ports["airpuff_left"]["y"],
                     c="red", s=30)
        ax_h.scatter(ports["airpuff_right"]["x"], ports["airpuff_right"]["y"],
                     c="red", s=30)

        ax_h.set_xlim(x_min, x_max)
        ax_h.set_ylim(y_min, y_max)
        ax_h.set_xticks([])
        ax_h.set_yticks([])

    plt.suptitle(f"Occupancy Heatmaps — Mouse {mouse}", fontsize=16)
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    # ----------------------------------------------
    # Saving
    # ----------------------------------------------
    if save_dir is not None:
        save_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_dir / f"{batch}_{mouse}_heatmaps_{n_bins}bins.png", dpi=300)
        fig.savefig(save_dir / f"{batch}_{mouse}_heatmaps_{n_bins}bins.pdf")

    plt.show()

def extract_behav_summary(behav_df, mouse, batch, group,
                           fps,
                           behav_cols=[],
                           speed_col='Speed',
                           immobility_threshold=0.1):
    """
    Extract a summary of behavioral metrics from a behav_df into a flat dict,
    suitable for building a cross-animal Excel summary.

    The function is experiment-agnostic: any binary columns can be passed via 
    behav_cols, that get time + percentage columns; Locomotion metrics are always computed.

    Parameters
    ----------
    behav_df : pd.DataFrame
        Output of analyze_mouse_position(), one row per frame.
    mouse : str
        Mouse identifier.
    batch : str
        Batch identifier.
    group : str
        Experimental group.
    fps : float
        Video frame rate (frames per second). Used to convert frames → seconds.
    scale : float
        Scale in cm per pixel. Used for distance computation.
    behav_cols : list of str, optional
        Binary columns — time + percentage columns are
        generated. E.g. ['Closed arm', 'Open arm', 'Center'].
    speed_col : str
        Column name for instantaneous speed (cm/s).
    immobility_threshold : float
        Speed (cm/s) below which the animal is considered immobile.

    Returns
    -------
    record : dict
        Flat dictionary with all summary metrics for one animal.
    """
    record = {
        'Mouse' : mouse,
        'Batch' : batch,
        'Group' : group,
    }

    # ── Session duration ──────────────────────────────────────────────────────
    n_frames      = len(behav_df)
    total_time_s  = n_frames / fps
    record['Total time (s)'] = round(total_time_s, 2)

    # ── Zone times + percentages ──────────────────────────────────────────────
    for col in behav_cols:
        if col not in behav_df.columns:
            print(f"  Warning: zone column '{col}' not found, filling with NaN.")
            record[f'{col} time (s)']  = float('nan')
            record[f'{col} (%)']       = float('nan')
            continue
        frames_in_behav         = behav_df[col].sum()
        time_in_behav            = frames_in_behav  / fps
        pct                     = (frames_in_behav  / n_frames) * 100
        record[f'{col} time (s)'] = round(time_in_behav, 2)
        record[f'{col} (%)']      = round(pct, 2)

    # ── Locomotion metrics ────────────────────────────────────────────────────
    if speed_col in behav_df.columns:
        speed = behav_df[speed_col].values

        # Immobility: frames below threshold → seconds
        immobile_frames          = (speed < immobility_threshold).sum()
        record['Immobility (s)'] = round(immobile_frames / fps, 2)

        # Total distance: sum of per-frame displacement (speed / fps * scale)
        # Speed is already in cm/s, so distance per frame = speed / fps
        total_distance           = (speed / fps).sum()
        record['Total distance (cm)'] = round(total_distance, 2)

        record['Mean speed (cm/s)'] = round(float(speed.mean()), 3)
    else:
        print(f"  Warning: speed column '{speed_col}' not found.")
        record['Immobility (s)']      = float('nan')
        record['Total distance (cm)'] = float('nan')
        record['Mean speed (cm/s)']   = float('nan')

    return record