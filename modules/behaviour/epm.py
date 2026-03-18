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
    
    coords_byzone_df = classify_position(coords_x, coords_y, epm_coordinates)
    speed_df = mp.compute_speed(coords, dist_scale, bodypart)
    
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
        if is_in_zone(x, y, closed_box):
            closed_arm[i] = 1
        elif is_in_zone(x, y, open_box):
            open_arm[i] = 1
        elif is_in_zone(x, y, center_box):
            center[i] = 1
        # If not in any zone, all remain 0 (edge/undefined area)
    
    coords_byzone_df = pd.DataFrame({
        'Closed arm': closed_arm,
        'Open arm': open_arm,
        'Center': center
    })
    
    return coords_byzone_df

def plot_epm_behavior(behav_df, epm_coordinates, bodypart='head', 
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
    open_time = behav_df['Open arm'].sum()
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
        fig_pie.savefig(save_dir / 'epm_pie_chart.png', dpi=300, bbox_inches='tight')
        fig_pie.savefig(save_dir / 'epm_pie_chart.pdf', bbox_inches='tight')
    
    plt.show()
    
    # =====================================================
    # Plot 2: Heatmap(s) of trajectory and occupancy
    # =====================================================
    fig_heatmap = plot_epm_heatmap(
        behav_df, epm_coordinates, bodypart=bodypart,
        n_bins=n_bins, bins=bins, save_dir=save_dir, figsize=figsize
    )
    
    return fig_pie, fig_heatmap

def plot_epm_heatmap(behav_df, epm_coordinates, bodypart='head',
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
        fig.savefig(save_dir / f'epm_heatmap_{n_bins}bins.png', dpi=300, bbox_inches='tight')
        fig.savefig(save_dir / f'epm_heatmap_{n_bins}bins.pdf', bbox_inches='tight')
    
    plt.show()
    
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
 