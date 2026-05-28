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

    metrics["airpuffs_count_cumulative"]      = metrics["airpuffs_count"].cumsum()
    metrics["licking_time_cumulative"]        = metrics["licking_time"].cumsum()
    metrics["nose_in_airport_time_cumulative"]= metrics["nose_in_airport_time"].cumsum()

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

def _arena_bounds(arena_json):
    with open(arena_json) as f:
        r = json.load(f)["Arena_rectangle_px"]
    return r["x1"], r["y1"], r["x2"], r["y2"]

def load_ports_px(ports_json):
    """Load port positions as {name: np.array([x, y])}."""
    with open(ports_json) as f:
        ports = json.load(f)
    return {name: np.array([pt["x"], pt["y"]]) for name, pt in ports.items()}

def compute_reference_ports(all_ports_dict):
    """
    Compute reference port positions as the mean across all mice.
    No single mouse is privileged; total distortion is minimised.

    Parameters
    ----------
    all_ports_dict : dict {mouse: {port_name: np.array([x, y])}}
    """
    port_names = list(next(iter(all_ports_dict.values())).keys())
    return {
        name: np.mean(
            [ports[name] for ports in all_ports_dict.values() if name in ports],
            axis=0
        )
        for name in port_names
    }

def estimate_port_transform(src_ports, ref_ports):
    """
    Estimate the affine transform that maps src port positions to ref port positions.
    With 3 ports the system is exactly determined; with more it is solved via least squares.

    Parameters
    ----------
    src_ports, ref_ports : dict {name: np.array([x, y])}

    Returns
    -------
    M : np.ndarray, shape (2, 3)
        Affine matrix such that [x', y'] = M @ [x, y, 1].
    """
    common = [n for n in src_ports if n in ref_ports]
    if len(common) < 2:
        raise ValueError(f"Need ≥ 2 common ports, found {len(common)}: {common}")

    src_pts = np.array([src_ports[n] for n in common])
    ref_pts = np.array([ref_ports[n] for n in common])

    N = len(src_pts)
    A = np.zeros((2 * N, 6))
    b = np.zeros(2 * N)
    for i, ((x, y), (xr, yr)) in enumerate(zip(src_pts, ref_pts)):
        A[2*i]   = [x, y, 1, 0, 0, 0]
        A[2*i+1] = [0, 0, 0, x, y, 1]
        b[2*i],  b[2*i+1] = xr, yr

    params = np.linalg.lstsq(A, b, rcond=None)[0]
    return params.reshape(2, 3)


def apply_transform(x, y, M):
    """Apply a 2×3 affine matrix to coordinate arrays x, y."""
    pts = np.stack([np.asarray(x, float),
                    np.asarray(y, float),
                    np.ones(len(x))])
    result = M @ pts
    return result[0], result[1]


def get_aligned_arena_bounds(arena_json, M, pad=100):
    """
    Apply affine transform to the arena rectangle corners and return
    (x_min, y_min, x_max, y_max) of the transformed arena with padding.
    """
    x1, y1, x2, y2 = _arena_bounds(arena_json)
    cx = np.array([x1 - pad, x2 + pad, x2 + pad, x1 - pad])
    cy = np.array([y1 - pad, y1 - pad, y2 + pad, y2 + pad])
    tx, ty = apply_transform(cx, cy, M)
    return tx.min(), ty.min(), tx.max(), ty.max()

def plot_group_heatmap(aligned_positions_list, mouse_list,
                       ref_ports=None, arena_bounds=None,
                       bins=(50, 50), n_bins=1,
                       cmap="jet", vmax=None, label='Group', save_dir=None):
    n_mice = len(aligned_positions_list)

    if arena_bounds is not None:
        x_min, y_min, x_max, y_max = arena_bounds
    else:
        all_x = np.concatenate([x for x, y in aligned_positions_list])
        all_y = np.concatenate([y for x, y in aligned_positions_list])
        x_min, x_max = all_x.min() - 50, all_x.max() + 50
        y_min, y_max = all_y.min() - 50, all_y.max() + 50

    hist_range = [[x_min, x_max], [y_min, y_max]]
    extent     = [x_min, x_max, y_min, y_max]
    outline_x  = [x_min, x_max, x_max, x_min, x_min]
    outline_y  = [y_min, y_min, y_max, y_max, y_min]

    fig, axes = plt.subplots(2, n_bins, figsize=(4 * n_bins, 8),
                              gridspec_kw={"height_ratios": [1, 4]})
    if n_bins == 1:
        axes = axes[:, np.newaxis]

    def _draw_ports(ax):
        if ref_ports is None:
            return
        for pname, pos in ref_ports.items():
            c = "lime" if "lick" in pname.lower() else "red"
            ax.scatter(pos[0], pos[1], c=c, s=40, zorder=5)

    for bin_idx in range(n_bins):
        heatmaps, traj_xs, traj_ys = [], [], []

        for x, y in aligned_positions_list:
            n = len(x)
            start = int(bin_idx * n / n_bins)
            stop  = n if bin_idx == n_bins - 1 else int((bin_idx + 1) * n / n_bins)
            xb, yb = x[start:stop], y[start:stop]
            traj_xs.append(xb); traj_ys.append(yb)
            hm, _, _ = np.histogram2d(xb, yb, bins=bins, range=hist_range)
            hm = hm / hm.sum() if hm.sum() > 0 else hm
            heatmaps.append(hm)

        group_map = np.mean(heatmaps, axis=0)

        ax_t = axes[0, bin_idx]
        for xb, yb in zip(traj_xs, traj_ys):
            ax_t.plot(xb, yb, color="black", linewidth=0.4, alpha=0.3)
        ax_t.plot(outline_x, outline_y, color="grey", linewidth=1) 
        _draw_ports(ax_t)
        ax_t.set(xlim=(x_min, x_max), ylim=(y_min, y_max),
                 xticks=[], yticks=[], aspect="equal", title=f"Bin {bin_idx + 1}")

        ax_h = axes[1, bin_idx]
        nonzero = group_map[group_map > 0]
        _vmax = vmax or (float(np.nanpercentile(nonzero, 99)) if nonzero.size else 1)
        ax_h.imshow(group_map.T, origin="lower", cmap=cmap,
                    vmin=0, vmax=_vmax, extent=extent,
                    interpolation="bilinear", aspect="equal")
        ax_h.plot(outline_x, outline_y, color="black", linewidth=1)
        _draw_ports(ax_h)
        ax_h.set(xlim=(x_min, x_max), ylim=(y_min, y_max), xticks=[], yticks=[])

    plt.suptitle(f"Group Occupancy Heatmap — {label} (n={n_mice})", fontsize=14)
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    if save_dir is not None:
        save_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_dir / f"group_heatmap_{n_bins}bins.png", dpi=300)
        fig.savefig(save_dir / f"group_heatmap_{n_bins}bins.pdf")
    plt.show()
    return fig

def _get_bouts(binary_col: pd.Series, time_col: pd.Series) -> list[dict]:
    """
    Detect contiguous runs of 1 in a binary column.
 
    Returns
    -------
    list of dict with keys:
        start_idx, end_idx   – positional indices in the original DataFrame
        start_time, end_time – values from time_col
        duration_s           – end_time - start_time
    """
    arr = binary_col.to_numpy()
    t   = time_col.to_numpy()
 
    bouts = []
    in_bout = False
    start = 0
 
    for i, v in enumerate(arr):
        if v == 1 and not in_bout:
            in_bout = True
            start   = i
        elif v != 1 and in_bout:
            in_bout = False
            bouts.append(dict(start_idx=start, end_idx=i - 1,
                               start_time=t[start], end_time=t[i - 1],
                               duration_s=t[i - 1] - t[start]))
    if in_bout:                           # bout that reaches the last sample
        bouts.append(dict(start_idx=start, end_idx=len(arr) - 1,
                           start_time=t[start], end_time=t[-1],
                           duration_s=t[-1] - t[start]))
    return bouts
 
 
def _freezing_stats(freezing_col: pd.Series,
                    time_col:     pd.Series,
                    start_idx:    int,
                    end_idx:      int) -> dict:
    """
    Compute freezing time (s) and percentage for a window [start_idx, end_idx].
    Uses the actual per-sample time-step so uneven sampling is handled correctly.
    """
    seg_freeze = freezing_col.iloc[start_idx : end_idx + 1].to_numpy()
    seg_time   = time_col.iloc[start_idx : end_idx + 1].to_numpy()
 
    if len(seg_time) < 2:
        return dict(freezing_time_s=0.0, freezing_pct=0.0)
 
    # dt for each sample = half the gap to the previous + half the gap to the next
    dt = np.diff(seg_time, prepend=seg_time[0], append=seg_time[-1])
    dt = (dt[:-1] + dt[1:]) / 2          # central differences → same length as seg
 
    freezing_time = float(np.sum(seg_freeze * dt))
    total_time    = float(np.sum(dt))
    freezing_pct  = 100.0 * freezing_time / total_time if total_time > 0 else 0.0
 
    return dict(freezing_time_s=round(freezing_time, 4),
                freezing_pct=round(freezing_pct, 2))
 
 
def _pre_cs_window(start_idx: int, n_samples: int, time_col: pd.Series) -> tuple[int, int]:
    """
    Return (pre_start_idx, pre_end_idx) for a baseline window of *n_samples*
    immediately before *start_idx*, clipped to the beginning of the recording.
    """
    pre_end   = start_idx - 1
    pre_start = max(0, start_idx - n_samples)
    return pre_start, pre_end
 
 
# ─────────────────────────────────────────────────────────────────────────────
# Phase detection
# ─────────────────────────────────────────────────────────────────────────────
 
# Expected CS counts per phase  (cs_plus_count, cs_minus_count)
_PHASE_PROFILES = {
    "Habituation":   (4,  4),
    "Conditioning":  (6,  6),
    "Extinction":    (12, 4),
}
 
def _detect_phase(n_cs_plus: int, n_cs_minus: int) -> str:
    """
    Heuristic phase label based on the number of CS+ and CS- presentations.
    Returns 'Unknown' if no profile matches.
    """
    for phase, (np_, nm_) in _PHASE_PROFILES.items():
        if n_cs_plus == np_ and n_cs_minus == nm_:
            return phase
    return f"Unknown (CS+={n_cs_plus}, CS-={n_cs_minus})"
 
 
# ─────────────────────────────────────────────────────────────────────────────
# Main public function
# ─────────────────────────────────────────────────────────────────────────────
 
def compute_behavior_metrics_FC(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute per-CS-presentation behavioural metrics for a Fear Conditioning session.
 
    Parameters
    ----------
    df : pd.DataFrame
        fiberbehav_notderived_df with at minimum:
        'Time(s)', 'CS+', 'CS-', 'Freezing'
        Optional but used if present: 'Shock', 'Speed'
 
    Returns
    -------
    pd.DataFrame
        One row per CS presentation, sorted chronologically.
        Columns
        -------
        Phase               – Habituation / Conditioning / Extinction / Unknown
        CS_type             – 'CS+' or 'CS-'
        CS_number           – 1-based index within CS type
        Start_time_s        – onset of the CS
        End_time_s          – offset of the CS
        CS_duration_s       – CS duration
        Freezing_time_s     – seconds spent freezing during CS
        Freezing_pct        – % of CS spent freezing
        PreCS_freezing_time_s  – freezing in the pre-CS baseline window
        PreCS_freezing_pct     – % of baseline window spent freezing
        PreCS_duration_s    – duration of the pre-CS baseline window actually used
        n_shocks            – number of shock bouts during this CS  (0 if no Shock col)
        Mean_speed          – mean speed during CS  (NaN if no Speed col)
    """
 
    required = {'Time(s)', 'CS+', 'CS-', 'Freezing'}
    missing  = required - set(df.columns)
    if missing:
        raise ValueError(f"DataFrame is missing required columns: {missing}")
 
    time     = df['Time(s)'].reset_index(drop=True)
    freezing = df['Freezing'].reset_index(drop=True)
    cs_plus  = df['CS+'].reset_index(drop=True)
    cs_minus = df['CS-'].reset_index(drop=True)
 
    has_shock = 'Shock' in df.columns
    has_speed = 'Speed' in df.columns
    if has_shock:
        shock = df['Shock'].reset_index(drop=True)
    if has_speed:
        speed = df['Speed'].reset_index(drop=True)
 
    # ── Detect CS bouts ──────────────────────────────────────────────────────
    cs_plus_bouts  = _get_bouts(cs_plus,  time)
    cs_minus_bouts = _get_bouts(cs_minus, time)
 
    phase = _detect_phase(len(cs_plus_bouts), len(cs_minus_bouts))
 
    rows = []
 
    # ── Helper to build one row ───────────────────────────────────────────────
    def _make_row(cs_type, cs_num, bout):
        si, ei = bout['start_idx'], bout['end_idx']
        n_cs_samples = ei - si + 1
 
        # CS-period freezing
        freeze_stats = _freezing_stats(freezing, time, si, ei)
 
        # Pre-CS baseline (same number of samples as the CS)
        pre_si, pre_ei = _pre_cs_window(si, n_cs_samples, time)
        if pre_ei >= pre_si:
            pre_stats = _freezing_stats(freezing, time, pre_si, pre_ei)
            pre_dur   = float(time.iloc[pre_ei] - time.iloc[pre_si])
        else:
            pre_stats = dict(freezing_time_s=np.nan, freezing_pct=np.nan)
            pre_dur   = 0.0
 
        # Shocks
        n_shocks = 0
        if has_shock:
            shock_bouts = _get_bouts(shock.iloc[si:ei + 1].reset_index(drop=True),
                                     time.iloc[si:ei + 1].reset_index(drop=True))
            n_shocks = len(shock_bouts)
 
        # Speed
        mean_speed = float(speed.iloc[si:ei + 1].mean()) if has_speed else np.nan
 
        return {
            'Phase':                   phase,
            'CS_type':                 cs_type,
            'CS_number':               cs_num,
            'Start_time_s':            round(bout['start_time'], 4),
            'End_time_s':              round(bout['end_time'],   4),
            'CS_duration_s':           round(bout['duration_s'], 4),
            'Freezing_time_s':         freeze_stats['freezing_time_s'],
            'Freezing_pct':            freeze_stats['freezing_pct'],
            'PreCS_freezing_time_s':   pre_stats['freezing_time_s'],
            'PreCS_freezing_pct':      pre_stats['freezing_pct'],
            'PreCS_duration_s':        round(pre_dur, 4),
            'n_shocks':                n_shocks,
            'Mean_speed':              round(mean_speed, 4) if not np.isnan(mean_speed) else np.nan,
        }
 
    for i, bout in enumerate(cs_plus_bouts):
        rows.append(_make_row('CS+', i + 1, bout))
 
    for i, bout in enumerate(cs_minus_bouts):
        rows.append(_make_row('CS-', i + 1, bout))
 
    metrics_df = (pd.DataFrame(rows)
                    .sort_values('Start_time_s')
                    .reset_index(drop=True))
 
    return metrics_df
