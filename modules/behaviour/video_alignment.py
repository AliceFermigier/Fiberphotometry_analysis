
import cv2
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from pathlib import Path
from PIL import Image
import importlib
import json

import modules.common.preprocess as pp
importlib.reload(pp)
import modules.behaviour.camera_processing as cp

from scripts.loader import project_root

def create_overlay_frame(timepoint, fiberbehav_df, window, figsize=(10, 4)):
    """
    Generates a matplotlib plot as image for a specific time window.
    """
    sr = len(fiberbehav_df) / (fiberbehav_df['Time(s)'].max() - fiberbehav_df['Time(s)'].min())
    idx_center = np.argmin(np.abs(fiberbehav_df['Time(s)'] - timepoint))
    half_window = int(window * sr / 2)
    idx_start = max(0, idx_center - half_window)
    idx_end = min(len(fiberbehav_df), idx_center + half_window)
    window_df = fiberbehav_df.iloc[idx_start:idx_end]
    
    # Setup figure
    has_speed = 'Speed' in window_df.columns
    behavior_cols = [col for col in window_df.columns if set(window_df[col].dropna().unique()).issubset({0,1})]
    fig_height = 1.5 + len(behavior_cols) + (1 if has_speed else 0)
    
    fig, axs = plt.subplots(fig_height, 1, figsize=figsize, sharex=True)
    fig.subplots_adjust(hspace=0.2)
    if fig_height == 1:
        axs = [axs]  # ensure list

    t = window_df['Time(s)']

    # Plot fiber signal
    axs[0].plot(t, window_df['Denoised dFF'], color='black')
    axs[0].set_ylabel('Denoised dFF')

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
    
    axs[-1].set_xlabel('Time (s)')
    axs[0].set_title(f'Time: {timepoint:.2f}s')

    # Convert to image
    canvas = FigureCanvas(fig)
    canvas.draw()
    img = np.frombuffer(canvas.tostring_rgb(), dtype='uint8')
    img = img.reshape(canvas.get_width_height()[::-1] + (3,))
    plt.close(fig)
    return img

def get_video_start(video_path, file_path):
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

def align_video_to_fiber(fiber_df, video_time):
    fiber_time = fiber_df['Time(s)'].values

    # Get frame numbers corresponding to fiber time points
    frame_indices = np.interp(fiber_time, video_time, np.arange(len(video_time)))

    # Round to nearest frame
    frame_indices = np.clip(np.round(frame_indices).astype(int), 0, len(video_time) - 1)

    return frame_indices

def make_combined_video(video_path, fiberbehav_df, output_path, frame_indices=None, window=10):
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
    cap = cv2.VideoCapture(str(video_path))
    fps = cap.get(cv2.CAP_PROP_FPS)
    n_frames = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
    width = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH))
    height = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))

    fourcc = cv2.VideoWriter_fourcc(*'mp4v')
    out = cv2.VideoWriter(str(output_path), fourcc, fps, (width, height + 300))  # Add space for overlay

    for i in range(n_frames):
        ret, frame = cap.read()
        if not ret:
            break
        
        # Get corresponding fiber_df row
        if frame_indices is not None and i < len(frame_indices):
            data_idx = frame_indices[i]
        else:
            data_idx = min(i, len(fiber_df) - 1)

        time_sec = i / fps
        overlay_img = create_overlay_frame(data_idx, fiberbehav_df, window=window)

        # Resize overlay to match width
        overlay_img = cv2.resize(overlay_img, (width, 300))
        combined_frame = np.vstack((frame, overlay_img))

        out.write(combined_frame)

    cap.release()
    out.release()
    print(f"Video saved to {output_path}")


