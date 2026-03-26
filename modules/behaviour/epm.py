import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import LinearSegmentedColormap
import pathlib

import modules.behaviour.mouse_position as mp
    
def analyze_mouse_position(coords, epm_coordinates, arena_scale, bodypart='head'):
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
    coords_x = coords[f'{bodypart}_x']
    coords_y = coords[f'{bodypart}_y']

    dist_scale = arena_scale['Scale_cm_per_px']
    video_fps = arena_scale['Video_fps']
    
    coords_byzone_df = classify_position(coords_x, coords_y, epm_coordinates)
    speed_df = mp.compute_speed(coords, video_fps, dist_scale)
    
    behav_df = pd.concat([coords, coords_byzone_df, speed_df], axis=1)
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
 
def plot_epm_dff_heatmap(fiberbehav_df, epm_coordinates, mouse, batch, bodypart='nose',
                          n_bins=1, bins=(50, 50), save_dir=None, figsize=None,
                          cmap='RdBu_r', vmin=None, vmax=None):
    """
    Plot fiber photometry dFF signal as a spatial heatmap over EPM zones.
    Each spatial bin displays the mean dFF value of all frames where the
    specified bodypart was located in that bin.

    Parameters:
    -----------
    fiberbehav_df : pd.DataFrame
        Combined fiber + behavior dataframe containing dFF and bodypart coordinates.
    epm_coordinates : dict
        EPM zone coordinates (same format as used in plot_epm_heatmap).
    mouse : str
        Mouse identifier (used for saving).
    batch : str
        Batch identifier (used for saving).
    bodypart : str
        Body part column prefix to use for x/y position (e.g. 'nose', 'center').
    n_bins : int
        Number of time bins to split the session into.
    bins : tuple
        (nx, ny) resolution of the spatial histogram.
    save_dir : str or Path, optional
        Directory to save PNG and PDF outputs.
    figsize : tuple, optional
        Figure size. Defaults to (5*n_bins, 5).
    cmap : str
        Colormap. 'RdBu_r' works well for signed dFF (negative=blue, positive=red).
    vmin, vmax : float, optional
        Color scale limits. If None, uses symmetric limits based on the data max.

    Returns:
    --------
    fig : matplotlib.figure.Figure
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    import pathlib

    if figsize is None:
        figsize = (5 * n_bins, 5)

    # --- Arena limits (with padding) ---
    coord_vals = epm_coordinates
    x_min = min(coord_vals['open_xL'], coord_vals['open_xR'],
                coord_vals['closed_xL'], coord_vals['closed_xR'],
                coord_vals['center_xL'], coord_vals['center_xR']) - 50
    x_max = max(coord_vals['open_xL'], coord_vals['open_xR'],
                coord_vals['closed_xL'], coord_vals['closed_xR'],
                coord_vals['center_xL'], coord_vals['center_xR']) + 50

    y_min = min(coord_vals['open_yBot'], coord_vals['open_yTop'],
                coord_vals['closed_yBot'], coord_vals['closed_yTop'],
                coord_vals['center_yBot'], coord_vals['center_yTop']) - 50
    y_max = max(coord_vals['open_yBot'], coord_vals['open_yTop'],
                coord_vals['closed_yBot'], coord_vals['closed_yTop'],
                coord_vals['center_yBot'], coord_vals['center_yTop']) + 50

    # --- Time bins ---
    total_len = len(fiberbehav_df)
    bin_len = total_len // n_bins
    dfs = [fiberbehav_df.iloc[i * bin_len:(i + 1) * bin_len] for i in range(n_bins)]
    dfs[-1] = fiberbehav_df.iloc[(n_bins - 1) * bin_len:]  # last bin keeps remainder

    # --- Figure ---
    fig, axes = plt.subplots(1, n_bins, figsize=figsize)
    if n_bins == 1:
        axes = [axes]

    for bin_idx, subdf in enumerate(dfs):
        x = subdf[f'{bodypart}_x'].values
        y = subdf[f'{bodypart}_y'].values
        dff = subdf['dFF'].values

        # Build mean-dFF map:
        # 1) Sum of dFF per bin
        dff_sum, xedges, yedges = np.histogram2d(
            x, y,
            bins=bins,
            range=[[x_min, x_max], [y_min, y_max]],
            weights=dff
        )
        # 2) Frame count per bin (for averaging)
        counts, _, _ = np.histogram2d(
            x, y,
            bins=bins,
            range=[[x_min, x_max], [y_min, y_max]]
        )

        # Mean dFF (NaN where no data)
        with np.errstate(invalid='ignore'):
            dff_map = np.where(counts > 0, dff_sum / counts, np.nan)

        # Symmetric color scale around 0 if not specified
        if vmin is None or vmax is None:
            abs_max = np.nanmax(np.abs(dff_map))
            _vmin, _vmax = -abs_max, abs_max
        else:
            _vmin, _vmax = vmin, vmax

        ax = axes[bin_idx]
        extent = [x_min, x_max, y_min, y_max]

        im = ax.imshow(
            dff_map.T,
            origin='lower',
            cmap=cmap,
            extent=extent,
            interpolation='bilinear',
            aspect='equal',
            vmin=_vmin,
            vmax=_vmax
        )

        _draw_epm_zones(ax, epm_coordinates, linewidth=2)

        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        ax.set_xticks([])
        ax.set_yticks([])

        title = f'dFF Heatmap - Bin {bin_idx + 1}' if n_bins > 1 else 'dFF Heatmap'
        ax.set_title(title, fontsize=11, fontweight='bold')

        cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label('Mean dFF', fontsize=9)

    plt.suptitle(
        f'Fiber Photometry dFF — {bodypart} position\n{batch} | {mouse}',
        fontsize=13, fontweight='bold'
    )
    plt.tight_layout(rect=[0, 0, 1, 0.93])

    if save_dir is not None:
        save_dir = pathlib.Path(save_dir)
        save_dir.mkdir(parents=True, exist_ok=True)
        stem = f'{batch}_{mouse}_dFF_heatmap_{bodypart}_{n_bins}bins'
        fig.savefig(save_dir / f'{stem}.png', dpi=300, bbox_inches='tight')
        fig.savefig(save_dir / f'{stem}.pdf', bbox_inches='tight')

    return fig

def plot_epm_dff_heatmap_grouped(x_list, y_list, dFF_list, subject_list,
                                  epm_coordinates, group, bodypart='nose',
                                  bins=(50, 50), show_individual=False,
                                  cmap='RdBu_r', vmin=None, vmax=None,
                                  use_zscore=False,  
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
    x_list, y_list, dFF_list : list of np.ndarray
        Position and dFF arrays, one entry per subject.
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

    signal_label = 'Mean z-scored dFF' if use_zscore else 'Mean dFF'

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
        f'dFF Heatmap — {group} | {bodypart}',
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
        stem = f'{tag}_{bodypart}_dFF_heatmap'
        fig.savefig(save_dir / f'{stem}.png', dpi=300, bbox_inches='tight')
        fig.savefig(save_dir / f'{stem}.pdf', bbox_inches='tight')

    return fig