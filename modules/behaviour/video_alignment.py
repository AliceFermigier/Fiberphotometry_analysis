#%%
import cv2
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from pathlib import Path
import importlib
import json
import math
from tqdm import tqdm

import modules.common.preprocess as pp
importlib.reload(pp)
import modules.behaviour.camera_processing as cp

from scripts.loader import project_root

#%%

def create_overlay_frame(index, fiberbehav_df, behavior_cols, window, figsize=(10, 4)):
    """
    Generates a matplotlib plot as image for a specific time window.
    """
    sr = len(fiberbehav_df) / (fiberbehav_df['Time(s)'].max() - fiberbehav_df['Time(s)'].min())
    idx_center = index
    half_window = int(window * sr / 2)
    idx_start = max(0, idx_center - half_window)
    idx_end = min(len(fiberbehav_df), idx_center + half_window)
    window_df = fiberbehav_df.iloc[idx_start:idx_end]
    
    # Setup figure
    n_behavior = len(behavior_cols)
    has_speed = 'Speed' in fiberbehav_df.columns
    height_ratios = [5] + [1]*n_behavior
    if has_speed:
        height_ratios += [2]
    
    fig, axs = plt.subplots(
        len(height_ratios),
        1,
        figsize=figsize,
        sharex=True,
        gridspec_kw={'height_ratios': height_ratios}
    )
    fig.subplots_adjust(hspace=0.2)
    if len(height_ratios) == 1:
        axs = [axs]  # ensure list

    t = window_df['Time(s)']

    # Plot fiber signal
    axs[0].plot(t, window_df['Denoised dFF'], color='black')
    axs[0].set_ylabel('Denoised dFF')
    axs[0].axvline(window_df['Time(s)'].iloc[idx_center - idx_start], color='red', linestyle='--')
    axs[0].set_ylim(fiberbehav_df['Denoised dFF'].min(),window_df['Denoised dFF'].max())

    # Plot behaviors
    behavior_colors_path = Path(project_root) / "modules/behaviour/behaviour_colors.json"
    with open(behavior_colors_path, "r") as f:
        behaviors_to_plot = json.load(f)
    
    for i, behavior in enumerate(behavior_cols):
        color, alpha = behaviors_to_plot.get(behavior, ('grey',0.05))
        axs[i + 1].fill_between(
            t,
            0,
            window_df[behavior],
            step='post',
            color=color,
            alpha=alpha
        )
        axs[i + 1].set_yticks([])
        axs[i + 1].set_ylabel(behavior, rotation=0, labelpad=35, va='center')

    # Optional: plot speed
    if has_speed:
        axs[-1].plot(t, window_df['Speed'], color='black')
        axs[-1].set_ylabel('Speed')
        axs[-1].set_ylim(fiberbehav_df['Speed'].min(),window_df['Speed'].max())

    # Hide x-axis labels and bottom spines for all but the last axis
    for ax in axs:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(True)

    for ax in axs[:-1]:
        ax.tick_params(labelbottom=False)
        ax.spines['bottom'].set_visible(False)
    
    axs[-1].set_xlabel('Time (s)')

    # Convert to image
    canvas = FigureCanvas(fig)
    canvas.draw()
    img = np.frombuffer(canvas.buffer_rgba(), dtype='uint8')  
    img = img.reshape(canvas.get_width_height()[::-1] + (4,))  
    img = cv2.cvtColor(img, cv2.COLOR_RGBA2RGB)  
    plt.close(fig)
    return img

def get_video_time(video_path, file_path):
    '''
    Get timestamps of video frames, in seconds
    '''
    camera_df = cp.get_camera_flashes(file_path)
    camera_times = camera_df['Time(s)'].values
    cap = cv2.VideoCapture(str(video_path))
    n_frames = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
    cap.release()

    # camera_times is expected to be a list or array of camera flash timestamps (length == n_frames)
    # If not the same length, interpolate linearly
    if len(camera_times) != n_frames:
        video_time = np.linspace(camera_times[0], camera_times[-1], n_frames)
    else:
        video_time = np.array(camera_times)

    return video_time

def align_fiber_to_video(fiber_df, video_time):
    """
    For each video frame timestamp, find the closest fiberphotometry time index.
    
    --> Parameters:
        fiber_df : pd.DataFrame
            Must contain a 'Time(s)' column.
        video_time : np.ndarray
            Array of timestamps for each video frame.
    
    --> Returns:
        fiber_indices : np.ndarray
            Array of indices into fiber_df for each video frame.
    """
    fiber_time = fiber_df['Time(s)'].values

    # Interpolate: video_time → index of closest fiber_time
    fiber_indices = np.interp(video_time, fiber_time, np.arange(len(fiber_time)))

    # Round and clip to valid indices
    fiber_indices = np.clip(np.round(fiber_indices).astype(int), 0, len(fiber_time) - 1)

    return fiber_indices

def make_combined_video(video_path, fiberbehav_df, output_path, fiber_indices=None, window=10, verbose=True, test=True):
    """
    Function to align video with fiberphotometry signal and behavior

    Example output :

    +------------------------+          <- Behavior Video
    |                        |
    |      mouse.avi         |
    |                        |
    +------------------------+
    |   ΔF/F trace plot (real-time line w/ history)   |
    +------------------------+
    |  Gantt chart showing behaviors as horizontal bars  |
    +------------------------+

    Rolling Gantt logic : For each video frame, Use a fixed window (e.g. last 10 seconds)
    For each behavior, draw a colored bar (1 = filled, 0 = empty) like this:
    Walk:   ███     ████
    Groom:      ████    
    Rear:   █     █     

    --> Parameters:
        fiberbehav_df : Fiberphotometry signal dataframe
        video_path : Path to mouse video (avi)
        output_path : Path to save the combined video (mp4)
        window : Time window in seconds for trace and gantt display (default: 10 seconds)
    """
    # Check if the video file exists
    if not os.path.exists(video_path):
        raise FileNotFoundError(f"Video file not found: {video_path}")
    
    # Check if the DataFrame is empty or None
    if fiberbehav_df is None or fiberbehav_df.empty:
        raise ValueError("Provided fiberbehav_df is empty or None.")
    
    # Check output file extension
    valid_extensions = ['.mp4', '.avi', '.mov', '.mkv']
    _, ext = os.path.splitext(output_path)
    if ext.lower() not in valid_extensions:
        raise ValueError(f"Unsupported video extension '{ext}'. Must be one of: {valid_extensions}")
    
    if verbose:
        print("All input validations passed. Proceeding with video generation...")

    # Extract list of behaviors
    behavior_cols = [
        col for col in fiberbehav_df.columns
        if set(fiberbehav_df[col].dropna().unique()).issubset({0, 1})
    ]
    if verbose:
        print(f"Behaviors : {behavior_cols}")

    cap = cv2.VideoCapture(str(video_path))
    fps = cap.get(cv2.CAP_PROP_FPS)
    if test:
        n_frames = int(30 * fps)
    else:
        n_frames = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
    width = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH))
    height = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))

    fourcc = cv2.VideoWriter_fourcc(*'mp4v')
    out = cv2.VideoWriter(str(output_path), fourcc, fps, (width, height + 300))
    
    if verbose:
        print(f"[INFO] Processing video: {video_path}")
        print(f"[INFO] Total frames: {n_frames}, FPS: {fps}")
        print(f"[INFO] Output video will be saved to: {output_path}\n")

    for i in tqdm(range(n_frames), desc="Rendering frames", unit="frame", leave=False):
        ret, frame = cap.read()
        if not ret:
            print(f"[WARNING] Failed to read frame {i}. Skipping.")
            continue
        
        # Get corresponding fiber_df row
        if fiber_indices is not None and i < len(fiber_indices):
            data_idx = fiber_indices[i]
        else:
            data_idx = min(i, len(fiberbehav_df) - 1)

        time_val = fiberbehav_df['Time(s)'].iloc[data_idx]
        overlay_img = create_overlay_frame(data_idx, fiberbehav_df, behavior_cols, window=window)

        # Resize overlay to match width
        overlay_img = cv2.resize(overlay_img, (width, 300))
        combined_frame = np.vstack((frame, overlay_img))

        out.write(combined_frame)

        if verbose and i % int(fps) == 0:  # Every ~1 second
            print(f"  ⏱ Frame {i}/{n_frames} → Time {time_val:.2f}s → Data idx {data_idx}")

    cap.release()
    out.release()
    print(f"\n✅ Combined video saved to {output_path}")

#%%
exp_path = Path(r'E:\FiberPhotometry\202504_OptoFluidACh\Data\20250512_EPM')
analysis_path = Path(r'E:\FiberPhotometry\202504_OptoFluidACh\Analysis\EPM_2\length0_interbout0_o4f1')
video_path = exp_path / '768.avi'
raw_file_path = exp_path / '768.doric'
fiberbehav_df = pd.read_csv(analysis_path / '1_768_fiberbehavnotderived.csv')
output_path = exp_path / '768_combined.mp4'

#%%
video_time = get_video_time(video_path, raw_file_path)
print(video_time)
#%%
fiber_indices = align_fiber_to_video(fiberbehav_df, video_time)
print(fiber_indices)
#%%
make_combined_video(video_path, fiberbehav_df, output_path, fiber_indices=None, window=10, verbose=True)
# %%
