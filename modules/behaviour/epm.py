import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import LinearSegmentedColormap
import pathlib

import modules.behaviour.mouse_position as mp
import modules.behaviour.behaviour_metrics as bm

def mask_outside_arena(coords, epm_coordinates, dist_scale, buffer_cm=2):
    """
    Set all body part coordinates to NaN for frames where the head is detected
    outside the arena (open arm + buffer, closed arm, center), then interpolate.

    Parameters:
    -----------
    coords : pd.DataFrame
        DeepLabCut coordinates with columns like 'head_x', 'head_y', etc.
    epm_coordinates : dict
        EPM zone bounding boxes.
    dist_scale : float
        Scale in cm/px.
    buffer_cm : float
        Buffer around open arm in cm (default 6).

    Returns:
    --------
    coords_clean : pd.DataFrame
        Coordinates with outside-arena frames interpolated out.
    """
    buffer_px = buffer_cm / dist_scale

    open_box_buffered = {
        'xL':   epm_coordinates['open_xL']   - buffer_px,
        'xR':   epm_coordinates['open_xR']   + buffer_px,
        'yBot': epm_coordinates['open_yBot'] - buffer_px,
        'yTop': epm_coordinates['open_yTop'] + buffer_px,
    }
    closed_box = {
        'xL': epm_coordinates['closed_xL'], 'xR': epm_coordinates['closed_xR'],
        'yBot': epm_coordinates['closed_yBot'], 'yTop': epm_coordinates['closed_yTop']
    }
    center_box = {
        'xL': epm_coordinates['center_xL'], 'xR': epm_coordinates['center_xR'],
        'yBot': epm_coordinates['center_yBot'], 'yTop': epm_coordinates['center_yTop']
    }

    center_x = coords['center_x']
    center_y = coords['center_y']

    in_arena = np.array([
        is_in_zone(x, y, open_box_buffered)
        or is_in_zone(x, y, closed_box)
        or is_in_zone(x, y, center_box)
        for x, y in zip(center_x, center_y)
    ])

    coords_clean = coords.copy()
    coords_clean[~in_arena] = np.nan
    coords_clean = coords_clean.interpolate(method='linear', limit_direction='both')

    return coords_clean, ~in_arena
    
def analyze_mouse_position(coords, epm_coordinates, arena_scale, bodypart='center'):
    """
    Analyze mouse position in elevated plus maze using absolute coordinates.
    
    Parameters:
    -----------
    coords : pd.DataFrame
        DeepLabCut coordinates with columns like 'head_x', 'head_y', etc.
    epm_coordinates : dict
        Dictionary with keys: 'open_xL', 'open_xR', 'open_yBot', 'open_yTop',
                              'closed_xL', 'closed_xR', 'closed_yBot', 'closed_yTop',
                              'center_xL', 'center_xR', 'center_yBot', 'center_yTop'
    bodypart : str
        Body part to analyze ('head', 'center', 'tailbase', etc.)
    
    Returns:
    --------
    behav_df : pd.DataFrame
        DataFrame with original coordinates, zone classification, and speed
    """
    dist_scale = arena_scale['Scale_cm_per_px']
    video_fps  = arena_scale['Video_fps']

    # Step 1: clean coordinates before anything else
    coords_clean, outside_mask = mask_outside_arena(coords, epm_coordinates, dist_scale)
    n_excluded = outside_mask.sum()
    if n_excluded > 0:
        print(f"  {n_excluded} frames ({n_excluded / len(coords) * 100:.1f}%) removed and interpolated (outside arena).")

    # Step 2: classify on cleaned coordinates
    coords_x = coords_clean[f'{bodypart}_x']
    coords_y = coords_clean[f'{bodypart}_y']
    coords_byzone_df = classify_position(coords_x, coords_y, epm_coordinates)
    coords_byzone_df = split_center_by_origin(coords_byzone_df)

    # Step 3: speed on cleaned coordinates
    speed_df = mp.compute_speed(coords_clean, video_fps, dist_scale)

    behav_df = pd.concat([coords_clean, coords_byzone_df, speed_df], axis=1)
    return behav_df

def is_in_zone(x, y, zone_box):
    """
    Check if coordinates (x, y) fall within a zone's bounding box.
    
    Parameters:
    -----------
    x, y : float
        Coordinates to check
    zone_box : dict
        Dictionary with keys 'xL', 'xR', 'yBot', 'yTop'
    
    Returns:
    --------
    bool
        True if coordinates are within the zone
    """
    x_min, x_max = min(zone_box['xL'], zone_box['xR']), max(zone_box['xL'], zone_box['xR'])
    y_min, y_max = min(zone_box['yBot'], zone_box['yTop']), max(zone_box['yBot'], zone_box['yTop'])
    
    return (x_min <= x <= x_max) and (y_min <= y <= y_max)

def split_center_by_origin(coords_byzone_df):
    """
    Split 'Center' frames into 'Closed arm -> Center' and 'Open arm -> Center'
    based on the last non-center zone the mouse was in before entering center.

    Frames where the mouse enters center with no prior zone history default
    to 'Closed arm -> Center'.

    Parameters
    ----------
    coords_byzone_df : pd.DataFrame
        Output of classify_position(), must contain 'Center', 'Open arm', 'Closed arm'.

    Returns
    -------
    pd.DataFrame with 'Center' replaced by two columns:
        'Closed arm -> Center' and 'Open arm -> Center'
    """
    n = len(coords_byzone_df)
    closed_to_center = np.zeros(n, dtype=int)
    open_to_center   = np.zeros(n, dtype=int)

    # Track the last zone seen before center (default: closed arm)
    last_zone = 'Closed arm'

    for i in range(n):
        if coords_byzone_df.at[i, 'Open arm'] == 1:
            last_zone = 'Open arm'
        elif coords_byzone_df.at[i, 'Closed arm'] == 1:
            last_zone = 'Closed arm'
        elif coords_byzone_df.at[i, 'Center'] == 1:
            if last_zone == 'Open arm':
                open_to_center[i] = 1
            else:
                closed_to_center[i] = 1
        # Head dipping and undefined frames don't update last_zone

    result = coords_byzone_df.drop(columns=['Center']).copy()
    result['Closed arm to Center'] = closed_to_center
    result['Open arm to Center']   = open_to_center

    return result

def classify_position(coords_x, coords_y, epm_coordinates):
    """
    Classify mouse position into zones using absolute coordinates.
    
    Parameters:
    -----------
    coords_x : array-like
        X coordinates of mouse position
    coords_y : array-like
        Y coordinates of mouse position
    epm_coordinates : dict
        Dictionary containing bounding box coordinates for each zone
    
    Returns:
    --------
    coords_byzone_df : pd.DataFrame
        DataFrame with binary columns for each zone (Closed arm, Open arm, Center)
    """
    n_frames = len(coords_x)
    closed_arm = np.zeros(n_frames, dtype=int)
    open_arm = np.zeros(n_frames, dtype=int)
    center = np.zeros(n_frames, dtype=int)
    head_dips = np.zeros(n_frames, dtype=int)
    
    # Extract coordinates for each zone
    open_box = {
        'xL': epm_coordinates['open_xL'],
        'xR': epm_coordinates['open_xR'],
        'yBot': epm_coordinates['open_yBot'],
        'yTop': epm_coordinates['open_yTop']
    }
    
    closed_box = {
        'xL': epm_coordinates['closed_xL'],
        'xR': epm_coordinates['closed_xR'],
        'yBot': epm_coordinates['closed_yBot'],
        'yTop': epm_coordinates['closed_yTop']
    }
    
    center_box = {
        'xL': epm_coordinates['center_xL'],
        'xR': epm_coordinates['center_xR'],
        'yBot': epm_coordinates['center_yBot'],
        'yTop': epm_coordinates['center_yTop']
    }
    
    # Classify each frame
    for i, (x, y) in enumerate(zip(coords_x, coords_y)):
        if is_in_zone(x, y, center_box):
            center[i] = 1
        elif is_in_zone(x, y, closed_box):
            closed_arm[i] = 1
        elif is_in_zone(x, y, open_box):
            open_arm[i] = 1
        else:
            head_dips[i] = 1

        # If not in any zone, all remain 0 (edge/undefined area)
    
    coords_byzone_df = pd.DataFrame({
        'Closed arm': closed_arm,
        'Open arm': open_arm,
        'Center': center,
        'Head dipping' : head_dips
    })
    
    return coords_byzone_df

def plot_epm_behavior(behav_df, epm_coordinates, mouse, batch, bodypart='nose', 
                      n_bins=1, bins=(50, 50), save_dir=None, figsize=None):
    """
    Create behavioral plots for EPM test: pie chart and heatmap(s).
    
    Parameters:
    -----------
    behav_df : pd.DataFrame
        Behavioral dataframe from analyze_mouse_position() with columns:
        'head_x', 'head_y', 'Closed arm', 'Open arm', 'Center', 'Speed (cm/s)'
    epm_coordinates : dict
        Dictionary with zone coordinates (open_xL, open_xR, etc.)
    bodypart : str
        Body part used for analysis ('head', 'center', 'tailbase')
    n_bins : int
        Number of time bins to split the session (default: 1)
    bins : tuple
        Number of bins for heatmap (nrows, ncols)
    save_dir : str or Path
        Directory to save figures (optional)
    figsize : tuple
        Figure size for the heatmap plot
    
    Returns:
    --------
    fig_pie : matplotlib figure
        Pie chart figure
    fig_heatmap : matplotlib figure
        Heatmap figure(s)
    """
    
    # =====================================================
    # Plot 1: Pie chart of time spent in each zone
    # =====================================================
    fig_pie, ax_pie = plt.subplots(figsize=(8, 6))
    
    closed_time = behav_df['Closed arm'].sum()
    open_time = behav_df['Open arm'].sum() + behav_df['Head dipping'].sum()
    center_time = behav_df['Center'].sum()
    
    # Handle cases where mouse doesn't visit all zones
    total_time = closed_time + open_time + center_time
    if total_time == 0:
        print("Warning: No position data classified into zones")
        return fig_pie, None
    
    sizes = [closed_time, open_time, center_time]
    labels = [
        f'Closed Arm\n{closed_time/total_time*100:.1f}%',
        f'Open Arm\n{open_time/total_time*100:.1f}%',
        f'Center\n{center_time/total_time*100:.1f}%'
    ]
    colors = ['#FF6B6B', '#4ECDC4', '#FFE66D']
    explode = (0.05, 0.05, 0.05)
    
    ax_pie.pie(sizes, labels=labels, colors=colors, explode=explode,
               autopct='', startangle=90, textprops={'fontsize': 12})
    ax_pie.set_title('Time Distribution in EPM Zones', fontsize=14, fontweight='bold')
    
    if save_dir is not None:
        save_dir = pathlib.Path(save_dir)
        save_dir.mkdir(parents=True, exist_ok=True)
        fig_pie.savefig(save_dir / f'{batch}_{mouse}_epm_pie_chart.png', dpi=300, bbox_inches='tight')
        fig_pie.savefig(save_dir / f'{batch}_{mouse}_epm_pie_chart.pdf', bbox_inches='tight')
    
    # =====================================================
    # Plot 2: Heatmap(s) of trajectory and occupancy
    # =====================================================
    fig_heatmap = plot_epm_heatmap(
        behav_df, epm_coordinates, mouse, batch, bodypart=bodypart,
        n_bins=n_bins, bins=bins, save_dir=save_dir, figsize=figsize
    )
    
    return fig_pie, fig_heatmap

def plot_epm_heatmap(behav_df, epm_coordinates, mouse, batch, bodypart='nose',
                     n_bins=1, bins=(50, 50), save_dir=None, figsize=None):
    """
    Plot mouse trajectory and occupancy heatmap with EPM zone outlines.
    
    Parameters:
    -----------
    behav_df : pd.DataFrame
        Behavioral dataframe with position data
    epm_coordinates : dict
        EPM zone coordinates
    bodypart : str
        Body part to plot ('head', 'center', 'tailbase')
    n_bins : int
        Number of time bins
    bins : tuple
        Histogram bin resolution
    save_dir : str or Path
        Save directory (optional)
    figsize : tuple
        Figure size
    
    Returns:
    --------
    fig : matplotlib figure
        Heatmap figure
    """
    
    if figsize is None:
        figsize = (5*n_bins, 10)
    
    # Define arena limits (with padding)
    x_min = min(epm_coordinates['open_xL'], epm_coordinates['open_xR'],
                epm_coordinates['closed_xL'], epm_coordinates['closed_xR'],
                epm_coordinates['center_xL'], epm_coordinates['center_xR']) - 50
    x_max = max(epm_coordinates['open_xL'], epm_coordinates['open_xR'],
                epm_coordinates['closed_xL'], epm_coordinates['closed_xR'],
                epm_coordinates['center_xL'], epm_coordinates['center_xR']) + 50
    
    y_min = min(epm_coordinates['open_yBot'], epm_coordinates['open_yTop'],
                epm_coordinates['closed_yBot'], epm_coordinates['closed_yTop'],
                epm_coordinates['center_yBot'], epm_coordinates['center_yTop']) - 50
    y_max = max(epm_coordinates['open_yBot'], epm_coordinates['open_yTop'],
                epm_coordinates['closed_yBot'], epm_coordinates['closed_yTop'],
                epm_coordinates['center_yBot'], epm_coordinates['center_yTop']) + 50
    
    # Split into time bins
    total_len = len(behav_df)
    bin_len = total_len // n_bins
    dfs = [behav_df.iloc[i*bin_len : (i+1)*bin_len] for i in range(n_bins)]
    dfs[-1] = behav_df.iloc[(n_bins-1)*bin_len:]  # last bin takes remainder
    
    # Create figure
    fig, axes = plt.subplots(
        2, n_bins,
        figsize=figsize,
        gridspec_kw={"height_ratios": [1, 1.2]}
    )
    
    # Ensure axes is 2D even if n_bins=1
    if n_bins == 1:
        axes = axes.reshape(2, 1)
    
    for bin_idx, subdf in enumerate(dfs):
        x = subdf[f'{bodypart}_x'].values
        y = subdf[f'{bodypart}_y'].values
        
        # ===== Trajectory Panel =====
        ax_t = axes[0, bin_idx]
        ax_t.plot(x, y, color='black', linewidth=0.5, alpha=0.7)
        
        # Draw EPM zone outlines
        _draw_epm_zones(ax_t, epm_coordinates)
        
        ax_t.set_xlim(x_min, x_max)
        ax_t.set_ylim(y_min, y_max)
        ax_t.set_aspect('equal')
        ax_t.set_xticks([])
        ax_t.set_yticks([])
        ax_t.set_title(f'Trajectory - Bin {bin_idx+1}', fontsize=11, fontweight='bold')
        
        # ===== Heatmap Panel =====
        heatmap, xedges, yedges = np.histogram2d(
            x, y,
            bins=bins,
            range=[[x_min, x_max], [y_min, y_max]]
        )
        
        ax_h = axes[1, bin_idx]
        extent = [x_min, x_max, y_min, y_max]
        
        im = ax_h.imshow(
            heatmap.T,
            origin='lower',
            cmap='hot',
            extent=extent,
            interpolation='bilinear',
            aspect='equal'
        )
        
        # Draw EPM zone outlines
        _draw_epm_zones(ax_h, epm_coordinates, linewidth=2)
        
        ax_h.set_xlim(x_min, x_max)
        ax_h.set_ylim(y_min, y_max)
        ax_h.set_xticks([])
        ax_h.set_yticks([])
        ax_h.set_title(f'Occupancy Heatmap - Bin {bin_idx+1}', fontsize=11, fontweight='bold')
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax_h, fraction=0.046, pad=0.04)
        cbar.set_label('Frames', fontsize=9)
    
    plt.suptitle('Mouse Trajectory and Occupancy in EPM', fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    
    if save_dir is not None:
        import pathlib
        save_dir = pathlib.Path(save_dir)
        save_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_dir / f'{batch}_{mouse}_epm_heatmap_{n_bins}bins.png', dpi=300, bbox_inches='tight')
        fig.savefig(save_dir / f'{batch}_{mouse}_epm_heatmap_{n_bins}bins.pdf', bbox_inches='tight')
    
    return fig
 
def _draw_epm_zones(ax, epm_coordinates, linewidth=1.5, alpha=0.8):
    """
    Draw EPM zone rectangles on an axis.
    
    Parameters:
    -----------
    ax : matplotlib axis
        Axis to draw on
    epm_coordinates : dict
        EPM zone coordinates
    linewidth : float
        Line width for rectangles
    alpha : float
        Line alpha (transparency)
    """
    
    # Open arm
    open_x_min = min(epm_coordinates['open_xL'], epm_coordinates['open_xR'])
    open_x_max = max(epm_coordinates['open_xL'], epm_coordinates['open_xR'])
    open_y_min = min(epm_coordinates['open_yBot'], epm_coordinates['open_yTop'])
    open_y_max = max(epm_coordinates['open_yBot'], epm_coordinates['open_yTop'])
    
    open_rect = patches.Rectangle(
        (open_x_min, open_y_min),
        open_x_max - open_x_min, open_y_max - open_y_min,
        linewidth=linewidth, edgecolor='#4ECDC4', facecolor='none',
        alpha=alpha, label='Open Arm'
    )
    ax.add_patch(open_rect)
    
    # Closed arm
    closed_x_min = min(epm_coordinates['closed_xL'], epm_coordinates['closed_xR'])
    closed_x_max = max(epm_coordinates['closed_xL'], epm_coordinates['closed_xR'])
    closed_y_min = min(epm_coordinates['closed_yBot'], epm_coordinates['closed_yTop'])
    closed_y_max = max(epm_coordinates['closed_yBot'], epm_coordinates['closed_yTop'])
    
    closed_rect = patches.Rectangle(
        (closed_x_min, closed_y_min),
        closed_x_max - closed_x_min, closed_y_max - closed_y_min,
        linewidth=linewidth, edgecolor='#FF6B6B', facecolor='none',
        alpha=alpha, label='Closed Arm'
    )
    ax.add_patch(closed_rect)
    
    # Center
    center_x_min = min(epm_coordinates['center_xL'], epm_coordinates['center_xR'])
    center_x_max = max(epm_coordinates['center_xL'], epm_coordinates['center_xR'])
    center_y_min = min(epm_coordinates['center_yBot'], epm_coordinates['center_yTop'])
    center_y_max = max(epm_coordinates['center_yBot'], epm_coordinates['center_yTop'])
    
    center_rect = patches.Rectangle(
        (center_x_min, center_y_min),
        center_x_max - center_x_min, center_y_max - center_y_min,
        linewidth=linewidth, edgecolor='#FFE66D', facecolor='none',
        alpha=alpha, label='Center'
    )
    ax.add_patch(center_rect)
 
def plot_epm_dff_heatmap_grouped(x_list, y_list, dFF_list, subject_list, closed_arm_list,
                                  epm_coordinates, group, bodypart='nose',
                                  bins=(50, 50), show_individual=False,
                                  cmap='RdBu_r', vmin=None, vmax=None,
                                  use_zscore=False,
                                  signal_name='465nm',     
                                  save_dir=None, figsize=None):
    """
    Plot a group-averaged dFF spatial heatmap over EPM zones, with optional
    per-subject panels alongside the group average.

    Strategy
    --------
    For each subject, a per-subject mean-dFF map is computed (sum_dFF / count).
    The group average is then the mean of those per-subject maps (NaN-safe),
    so every mouse contributes equally regardless of session length.

    Parameters
    ----------
    x_list, y_list, dFF_list, closed_arm_list : list of np.ndarray
        Position, dFF and behaviour (closed arm) arrays, one entry per subject.
        closed_arm_list entries are 1d binary arrays the same length as the
        corresponding dFF array, where 1 = animal in the closed arm. Used as a
        mask to compute the z-score baseline (mean/std) when use_zscore=True.
    subject_list : list of str
        Subject labels matching the above lists.
    epm_coordinates : dict
        EPM zone coordinates.
    group : str
        Group label (used in title and filename).
    bodypart : str
        Body part label (used in title and filename).
    bins : tuple
        (nx, ny) spatial resolution.
    show_individual : bool
        If True, show one panel per subject + one group-average panel.
        If False, show only the group-average panel.
    cmap : str
        Colormap (diverging recommended, e.g. 'RdBu_r').
    vmin, vmax : float, optional
        Shared color scale. If None, symmetric limits are derived from the data.
    use_zscore : bool, optional
        If True, z-score each subject's dFF using the mean/std computed only
        from closed-arm timepoints (per closed_arm_list), then apply that
        baseline to the entire trace.
    signal_name : str, optional
        Wavelength of led. Will appear in the plot title  
    save_dir : str or Path, optional
        Save directory for PNG / PDF output.
    figsize : tuple, optional
        Figure size override.

    Returns
    -------
    fig : matplotlib.figure.Figure
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import pathlib

    n_subjects = len(subject_list)

    # Z-score each subject's dFF using only closed-arm samples as the baseline
    # (mean/std), then apply that baseline to the full trace
    if use_zscore:
        zscored_dFF_list = []
        for dff, closed_mask in zip(dFF_list, closed_arm_list):
            closed_mask = np.asarray(closed_mask).astype(bool)
            closed_dff = dff[closed_mask]
            zscored_dFF_list.append((dff - np.median(closed_dff)) / np.std(closed_dff, ddof=1))
        dFF_list = zscored_dFF_list

    signal_label = f'Mean {"z-scored (closed-arm baseline) " if use_zscore else ""}dFF ({signal_name})'

    # ── Arena limits ──────────────────────────────────────────────────────────
    c = epm_coordinates
    x_min = min(c['open_xL'], c['open_xR'], c['closed_xL'], c['closed_xR'],
                c['center_xL'], c['center_xR']) - 50
    x_max = max(c['open_xL'], c['open_xR'], c['closed_xL'], c['closed_xR'],
                c['center_xL'], c['center_xR']) + 50
    y_min = min(c['open_yBot'], c['open_yTop'], c['closed_yBot'], c['closed_yTop'],
                c['center_yBot'], c['center_yTop']) - 50
    y_max = max(c['open_yBot'], c['open_yTop'], c['closed_yBot'], c['closed_yTop'],
                c['center_yBot'], c['center_yTop']) + 50
    extent = [x_min, x_max, y_min, y_max]
    hist_range = [[x_min, x_max], [y_min, y_max]]

    # ── Per-subject mean-dFF maps ─────────────────────────────────────────────
    per_subject_maps = []

    for x, y, dff in zip(x_list, y_list, dFF_list):
        dff_sum, _, _ = np.histogram2d(x, y, bins=bins, range=hist_range, weights=dff)
        counts, _, _  = np.histogram2d(x, y, bins=bins, range=hist_range)
        with np.errstate(invalid='ignore'):
            subj_map = np.where(counts > 0, dff_sum / counts, np.nan)
        per_subject_maps.append(subj_map)

    # ── Group average (mean across subjects, NaN-safe) ────────────────────────
    stack     = np.stack(per_subject_maps, axis=0)          # (n_subjects, nx, ny)
    group_map = np.nanmean(stack, axis=0)                   # (nx, ny)
    group_sem = np.nanstd(stack, axis=0) / np.sqrt(         # optional SEM map
                    np.sum(~np.isnan(stack), axis=0).clip(min=1))

    # ── Shared colour scale ───────────────────────────────────────────────────
    if vmin is None or vmax is None:
        abs_max    = np.nanmax(np.abs(group_map))
        _vmin, _vmax = -abs_max, abs_max
    else:
        _vmin, _vmax = vmin, vmax

    # ── Layout ───────────────────────────────────────────────────────────────
    n_panels = (n_subjects + 1) if show_individual else 1
    if figsize is None:
        figsize = (4 * n_panels, 5)

    fig, axes = plt.subplots(1, n_panels, figsize=figsize)
    if n_panels == 1:
        axes = [axes]

    def _plot_panel(ax, dff_map, title, add_cbar=False):
        im = ax.imshow(
            dff_map.T,
            origin='lower',
            cmap=cmap,
            extent=extent,
            interpolation='bilinear',
            aspect='equal',
            vmin=_vmin,
            vmax=_vmax,
        )
        _draw_epm_zones(ax, epm_coordinates, linewidth=0.8, alpha=0.4)
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(title, fontsize=10, fontweight='bold')
        if add_cbar:
            cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
            cbar.set_label(signal_label, fontsize=9)
        return im

    # Individual panels
    if show_individual:
        for idx, (subj_map, subj_label) in enumerate(zip(per_subject_maps, subject_list)):
            _plot_panel(axes[idx], subj_map, title=subj_label, add_cbar=False)

    # Group average panel (always last, always has colorbar)
    avg_title = f'Group Average\n(n={n_subjects})'
    im = _plot_panel(axes[-1], group_map, title=avg_title, add_cbar=True)

    # Shared colorbar for individual panels
    if show_individual and n_subjects > 0:
        fig.colorbar(im, ax=axes[:-1], fraction=0.02, pad=0.02, label=signal_label)

    plt.suptitle(
        f'dFF Heatmap — {group} | {bodypart} | {signal_name}',
        fontsize=13, fontweight='bold'
    )
    plt.tight_layout(rect=[0, 0, 1, 0.93])

    # ── Save ─────────────────────────────────────────────────────────────────
    if save_dir is not None:
        save_dir = pathlib.Path(save_dir)
        save_dir.mkdir(parents=True, exist_ok=True)
        tag = group.replace(' ', '_')
        if use_zscore:
            tag = f'zscored_{tag}'
        if show_individual:
            tag = f'indiv_{tag}'
        stem = f'{tag}_{bodypart}_{signal_name}_dFF_heatmap' 
        fig.savefig(save_dir / f'{stem}.png', dpi=300, bbox_inches='tight')
        fig.savefig(save_dir / f'{stem}.pdf', bbox_inches='tight')

    return fig

def plot_epm_dff_heatmap_grouped(x_list, y_list, dFF_list, subject_list, closed_arm_list,
                                  epm_coordinates, group, bodypart='nose',
                                  bins=(50, 50), show_individual=False,
                                  cmap='RdBu_r', vmin=None, vmax=None,
                                  use_zscore=False,
                                  signal_name='465nm',     
                                  save_dir=None, figsize=None):
    """
    Plot a group-averaged dFF spatial heatmap over EPM zones, with optional
    per-subject panels alongside the group average.

    Strategy
    --------
    For each subject, a per-subject mean-dFF map is computed (sum_dFF / count).
    The group average is then the mean of those per-subject maps (NaN-safe),
    so every mouse contributes equally regardless of session length.

    Parameters
    ----------
    x_list, y_list, dFF_list, closed_arm_list : list of np.ndarray
        Position, dFF and behaviour (closed arm) arrays, one entry per subject.
    subject_list : list of str
        Subject labels matching the above lists.
    epm_coordinates : dict
        EPM zone coordinates.
    group : str
        Group label (used in title and filename).
    bodypart : str
        Body part label (used in title and filename).
    bins : tuple
        (nx, ny) spatial resolution.
    show_individual : bool
        If True, show one panel per subject + one group-average panel.
        If False, show only the group-average panel.
    cmap : str
        Colormap (diverging recommended, e.g. 'RdBu_r').
    vmin, vmax : float, optional
        Shared color scale. If None, symmetric limits are derived from the data.
    use_zscore : bool, optional
        Says if data is zcored or not before plotting
    signal_name : str, optional
        Wavelength of led. Will appear in the plot title  
    save_dir : str or Path, optional
        Save directory for PNG / PDF output.
    figsize : tuple, optional
        Figure size override.

    Returns
    -------
    fig : matplotlib.figure.Figure
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import pathlib

    n_subjects = len(subject_list)

    # Z-score each subject's dFF independently before mapping is zscoring specified
    if use_zscore:
        dFF_list = [
            (dff - np.mean(dff)) / np.std(dff, ddof=1)
            for dff in dFF_list
        ]

    signal_label = f'Mean {"z-scored " if use_zscore else ""}dFF ({signal_name})'

    # ── Arena limits ──────────────────────────────────────────────────────────
    c = epm_coordinates
    x_min = min(c['open_xL'], c['open_xR'], c['closed_xL'], c['closed_xR'],
                c['center_xL'], c['center_xR']) - 50
    x_max = max(c['open_xL'], c['open_xR'], c['closed_xL'], c['closed_xR'],
                c['center_xL'], c['center_xR']) + 50
    y_min = min(c['open_yBot'], c['open_yTop'], c['closed_yBot'], c['closed_yTop'],
                c['center_yBot'], c['center_yTop']) - 50
    y_max = max(c['open_yBot'], c['open_yTop'], c['closed_yBot'], c['closed_yTop'],
                c['center_yBot'], c['center_yTop']) + 50
    extent = [x_min, x_max, y_min, y_max]
    hist_range = [[x_min, x_max], [y_min, y_max]]

    # ── Per-subject mean-dFF maps ─────────────────────────────────────────────
    per_subject_maps = []

    for x, y, dff in zip(x_list, y_list, dFF_list):
        dff_sum, _, _ = np.histogram2d(x, y, bins=bins, range=hist_range, weights=dff)
        counts, _, _  = np.histogram2d(x, y, bins=bins, range=hist_range)
        with np.errstate(invalid='ignore'):
            subj_map = np.where(counts > 0, dff_sum / counts, np.nan)
        per_subject_maps.append(subj_map)

    # ── Group average (mean across subjects, NaN-safe) ────────────────────────
    stack     = np.stack(per_subject_maps, axis=0)          # (n_subjects, nx, ny)
    group_map = np.nanmean(stack, axis=0)                   # (nx, ny)
    group_sem = np.nanstd(stack, axis=0) / np.sqrt(         # optional SEM map
                    np.sum(~np.isnan(stack), axis=0).clip(min=1))

    # ── Shared colour scale ───────────────────────────────────────────────────
    if vmin is None or vmax is None:
        abs_max    = np.nanmax(np.abs(group_map))
        _vmin, _vmax = -abs_max, abs_max
    else:
        _vmin, _vmax = vmin, vmax

    # ── Layout ───────────────────────────────────────────────────────────────
    n_panels = (n_subjects + 1) if show_individual else 1
    if figsize is None:
        figsize = (4 * n_panels, 5)

    fig, axes = plt.subplots(1, n_panels, figsize=figsize)
    if n_panels == 1:
        axes = [axes]

    def _plot_panel(ax, dff_map, title, add_cbar=False):
        im = ax.imshow(
            dff_map.T,
            origin='lower',
            cmap=cmap,
            extent=extent,
            interpolation='bilinear',
            aspect='equal',
            vmin=_vmin,
            vmax=_vmax,
        )
        _draw_epm_zones(ax, epm_coordinates, linewidth=0.8, alpha=0.4)
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(title, fontsize=10, fontweight='bold')
        if add_cbar:
            cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
            cbar.set_label(signal_label, fontsize=9)
        return im

    # Individual panels
    if show_individual:
        for idx, (subj_map, subj_label) in enumerate(zip(per_subject_maps, subject_list)):
            _plot_panel(axes[idx], subj_map, title=subj_label, add_cbar=False)

    # Group average panel (always last, always has colorbar)
    avg_title = f'Group Average\n(n={n_subjects})'
    im = _plot_panel(axes[-1], group_map, title=avg_title, add_cbar=True)

    # Shared colorbar for individual panels
    if show_individual and n_subjects > 0:
        fig.colorbar(im, ax=axes[:-1], fraction=0.02, pad=0.02, label=signal_label)

    plt.suptitle(
        f'dFF Heatmap — {group} | {bodypart} | {signal_name}',
        fontsize=13, fontweight='bold'
    )
    plt.tight_layout(rect=[0, 0, 1, 0.93])

    # ── Save ─────────────────────────────────────────────────────────────────
    if save_dir is not None:
        save_dir = pathlib.Path(save_dir)
        save_dir.mkdir(parents=True, exist_ok=True)
        tag = group.replace(' ', '_')
        if use_zscore:
            tag = f'zscored_{tag}'
        if show_individual:
            tag = f'indiv_{tag}'
        stem = f'{tag}_{bodypart}_{signal_name}_dFF_heatmap' 
        fig.savefig(save_dir / f'{stem}.png', dpi=300, bbox_inches='tight')
        fig.savefig(save_dir / f'{stem}.pdf', bbox_inches='tight')

    return fig

# ─────────────────────────────────────────────────────────────────────────────
# EPM GROUP ALIGNMENT & PLOTTING
# ─────────────────────────────────────────────────────────────────────────────

def load_epm_landmarks(epm_coordinates):
    """
    Extract the 4 corners of the open arm and 4 corners of the closed arm
    as named landmark points for affine alignment.

    Parameters
    ----------
    epm_coordinates : dict

    Returns
    -------
    dict {landmark_name: np.array([x, y])}
    """
    return {
        'open_TL':   np.array([epm_coordinates['open_xL'],  epm_coordinates['open_yTop']]),
        'open_TR':   np.array([epm_coordinates['open_xR'],  epm_coordinates['open_yTop']]),
        'open_BL':   np.array([epm_coordinates['open_xL'],  epm_coordinates['open_yBot']]),
        'open_BR':   np.array([epm_coordinates['open_xR'],  epm_coordinates['open_yBot']]),
        'closed_TL': np.array([epm_coordinates['closed_xL'], epm_coordinates['closed_yTop']]),
        'closed_TR': np.array([epm_coordinates['closed_xR'], epm_coordinates['closed_yTop']]),
        'closed_BL': np.array([epm_coordinates['closed_xL'], epm_coordinates['closed_yBot']]),
        'closed_BR': np.array([epm_coordinates['closed_xR'], epm_coordinates['closed_yBot']]),
    }

def compute_reference_epm(all_landmarks_dict):
    """
    Compute mean landmark positions across all mice.
    Mirrors compute_reference_ports().

    Parameters
    ----------
    all_landmarks_dict : dict {mouse: {landmark_name: np.array([x, y])}}

    Returns
    -------
    dict {landmark_name: np.array([x, y])}
    """
    landmark_names = list(next(iter(all_landmarks_dict.values())).keys())
    return {
        name: np.mean(
            [lm[name] for lm in all_landmarks_dict.values() if name in lm],
            axis=0
        )
        for name in landmark_names
    }

def get_aligned_epm_coordinates(epm_coordinates, M):
    """
    Apply affine transform M to all epm_coordinates corner points,
    returning a transformed epm_coordinates dict in the reference space.

    Parameters
    ----------
    epm_coordinates : dict
    M : np.ndarray, shape (2, 3)

    Returns
    -------
    dict  (same keys as epm_coordinates)
    """
    # Collect all named points
    keys_x = ['open_xL', 'open_xR', 'closed_xL', 'closed_xR', 'center_xL', 'center_xR']
    keys_y = ['open_yBot', 'open_yTop', 'closed_yBot', 'closed_yTop', 'center_yBot', 'center_yTop']

    # Transform the 6 x-coords paired with a neutral y (0) and vice-versa is wrong —
    # we need to transform actual corner pairs, then unpack.
    corners = {
        'open_xL':    (epm_coordinates['open_xL'],    epm_coordinates['open_yBot']),
        'open_xR':    (epm_coordinates['open_xR'],    epm_coordinates['open_yTop']),
        'closed_xL':  (epm_coordinates['closed_xL'],  epm_coordinates['closed_yBot']),
        'closed_xR':  (epm_coordinates['closed_xR'],  epm_coordinates['closed_yTop']),
        'center_xL':  (epm_coordinates['center_xL'],  epm_coordinates['center_yBot']),
        'center_xR':  (epm_coordinates['center_xR'],  epm_coordinates['center_yTop']),
    }

    xs = np.array([v[0] for v in corners.values()])
    ys = np.array([v[1] for v in corners.values()])
    tx, ty = bm.apply_transform(xs, ys, M)
    t = dict(zip(corners.keys(), zip(tx, ty)))

    return {
        'open_xL':    t['open_xL'][0],    'open_xR':    t['open_xR'][0],
        'open_yBot':  t['open_xL'][1],    'open_yTop':  t['open_xR'][1],
        'closed_xL':  t['closed_xL'][0],  'closed_xR':  t['closed_xR'][0],
        'closed_yBot':t['closed_xL'][1],  'closed_yTop':t['closed_xR'][1],
        'center_xL':  t['center_xL'][0],  'center_xR':  t['center_xR'][0],
        'center_yBot':t['center_xL'][1],  'center_yTop':t['center_xR'][1],
    }

def plot_group_epm_heatmap(aligned_positions_list, ref_epm_coordinates,
                            bins=(50, 50), n_bins=1,
                            cmap='jet', vmax=None, label='Group', save_dir=None):
    """
    Plot group-level occupancy heatmap for EPM, aligned across mice.

    Parameters
    ----------
    aligned_positions_list : list of (x_array, y_array)
    ref_epm_coordinates : dict
        EPM coordinates in the reference (aligned) space, used for zone outlines.
    bins : tuple
    n_bins : int
        Number of time bins.
    cmap, vmax, label, save_dir : as in plot_group_heatmap
    """
    n_mice = len(aligned_positions_list)

    all_x = np.concatenate([x for x, y in aligned_positions_list])
    all_y = np.concatenate([y for x, y in aligned_positions_list])
    x_min, x_max = all_x.min() - 50, all_x.max() + 50
    y_min, y_max = all_y.min() - 50, all_y.max() + 50

    hist_range = [[x_min, x_max], [y_min, y_max]]
    extent     = [x_min, x_max, y_min, y_max]

    fig, axes = plt.subplots(2, n_bins, figsize=(4 * n_bins, 8),
                              gridspec_kw={"height_ratios": [1, 4]})
    if n_bins == 1:
        axes = axes[:, np.newaxis]

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

        # ── Trajectory panel ──────────────────────────────────────────────────
        ax_t = axes[0, bin_idx]
        for xb, yb in zip(traj_xs, traj_ys):
            ax_t.plot(xb, yb, color='black', linewidth=0.4, alpha=0.3)
        _draw_epm_zones(ax_t, ref_epm_coordinates, linewidth=1.5)
        ax_t.set(xlim=(x_min, x_max), ylim=(y_min, y_max),
                 xticks=[], yticks=[], aspect='equal',
                 title=f'Bin {bin_idx + 1}')

        # ── Heatmap panel ─────────────────────────────────────────────────────
        ax_h = axes[1, bin_idx]
        nonzero = group_map[group_map > 0]
        _vmax = vmax or (float(np.nanpercentile(nonzero, 99)) if nonzero.size else 1)
        ax_h.imshow(group_map.T, origin='lower', cmap=cmap,
                    vmin=0, vmax=_vmax, extent=extent,
                    interpolation='bilinear', aspect='equal')
        _draw_epm_zones(ax_h, ref_epm_coordinates, linewidth=2)
        ax_h.set(xlim=(x_min, x_max), ylim=(y_min, y_max), xticks=[], yticks=[])

    plt.suptitle(f'Group Occupancy Heatmap — {label} (n={n_mice})', fontsize=14)
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    if save_dir is not None:
        save_dir = pathlib.Path(save_dir)
        save_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_dir / f'group_epm_heatmap_{n_bins}bins.png', dpi=300)
        fig.savefig(save_dir / f'group_epm_heatmap_{n_bins}bins.pdf')
    plt.show()
    return fig

def plot_group_epm_pie(all_behav_dfs, all_groups, label='Group', save_dir=None):
    """
    Plot mean zone occupancy as a pie chart for each group, with per-mouse
    values shown as scatter dots on a companion bar chart.

    Parameters
    ----------
    all_behav_dfs : list of pd.DataFrame
        One behav_df per mouse (output of analyze_mouse_position).
    all_groups : list of str/int
        Group label for each mouse (same order as all_behav_dfs).
    label : str
    save_dir : Path or None
    """
    groups = sorted(set(all_groups))
    colors = ['#FF6B6B', '#4ECDC4', '#FFE66D']
    zone_labels = ['Closed Arm', 'Open Arm', 'Center']

    fig, axes = plt.subplots(1, len(groups), figsize=(5 * len(groups), 5))
    if len(groups) == 1:
        axes = [axes]

    for ax, group in zip(axes, groups):
        indices = [i for i, g in enumerate(all_groups) if g == group]
        fractions = []  # shape: (n_mice, 3)

        for i in indices:
            df = all_behav_dfs[i]
            closed = df['Closed arm'].sum()
            open_  = df['Open arm'].sum() + df['Head dipping'].sum()
            center = df['Center'].sum()
            total  = closed + open_ + center
            if total > 0:
                fractions.append([closed / total, open_ / total, center / total])

        fractions = np.array(fractions)  # (n_mice, 3)
        mean_fracs = fractions.mean(axis=0)

        wedges, _ = ax.pie(
            mean_fracs,
            colors=colors,
            startangle=90,
            wedgeprops=dict(width=0.6),   # donut style — cleaner for group plots
        )
        # Annotate wedges with mean %
        for wedge, frac, zlabel in zip(wedges, mean_fracs, zone_labels):
            angle = (wedge.theta1 + wedge.theta2) / 2
            x = 0.75 * np.cos(np.radians(angle))
            y = 0.75 * np.sin(np.radians(angle))
            ax.text(x, y, f'{frac*100:.1f}%', ha='center', va='center', fontsize=10)

        ax.legend(wedges, zone_labels, loc='lower center',
                  bbox_to_anchor=(0.5, -0.15), fontsize=9, frameon=False)
        n = len(fractions)
        ax.set_title(f'Group {group}\n(n={n})', fontsize=12, fontweight='bold')

    plt.suptitle(f'Zone Occupancy — {label}', fontsize=14, fontweight='bold')
    plt.tight_layout()

    if save_dir is not None:
        save_dir = pathlib.Path(save_dir)
        save_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_dir / 'group_epm_pie.png', dpi=300, bbox_inches='tight')
        fig.savefig(save_dir / 'group_epm_pie.pdf', bbox_inches='tight')
    plt.show()
    return fig