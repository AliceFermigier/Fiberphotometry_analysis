#%%
import cv2
import sys
import os
import shutil
import numpy as np
import pandas as pd
import gc
import subprocess
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['text.antialiased'] = True
plt.rcParams['lines.antialiased'] = True
plt.rcParams.update({
    'font.size': 35,   
    'axes.titlesize': 35,
    'axes.labelsize': 35,
    'xtick.labelsize': 32,
    'ytick.labelsize': 32
})
plt.rcParams['font.family'] = 'Arial'
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from pathlib import Path
import importlib
import json
from tqdm import tqdm

project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

import modules.common.preprocess as pp
importlib.reload(pp)
import modules.behaviour.camera_processing as cp
importlib.reload(cp)

ffmpeg_path = shutil.which("ffmpeg")
if not ffmpeg_path:
    ffmpeg_path = r"C:\Users\afermigier\Programmes\ffmpeg-7.1.1-essentials_build\bin\ffmpeg.exe"

if not os.path.isfile(ffmpeg_path):
    raise FileNotFoundError(f"ffmpeg not found at: {ffmpeg_path}")

#%%

def make_even(val):
    return val if val % 2 == 0 else val - 1

def create_overlay_frame(index, fiberbehav_df, behavior_cols, window):
    """
    Generates a matplotlib plot as image for a specific time window.
    """
    sr = len(fiberbehav_df) / (fiberbehav_df['Time(s)'].max() - fiberbehav_df['Time(s)'].min())
    quarter_window = int(window * sr / 4)
    idx_center = index
    idx_start = max(0, idx_center - (3*quarter_window))
    idx_end = min(len(fiberbehav_df), idx_center + quarter_window)
    window_df = fiberbehav_df.iloc[idx_start:idx_end]
    
    n_behavior = len(behavior_cols)
    has_speed = 'Speed' in fiberbehav_df.columns
    has_560 = '560 dFF' in fiberbehav_df.columns
    height_ratios = [8]
    if has_560:
        height_ratios += [8]
    height_ratios += [1]*n_behavior
    if has_speed:
        height_ratios += [3]
    
    center_time = window_df['Time(s)'].iloc[idx_center - idx_start]
    end_time = window_df['Time(s)'].iloc[-1]

    height = sum(height_ratios)/2
    size = (15 , height)
    
    fig, axs = plt.subplots(
        len(height_ratios),
        1,
        dpi=72,
        figsize=size,
        sharex=True,
        gridspec_kw={'height_ratios': height_ratios}
    )

    if len(height_ratios) == 1:
        axs = [axs]

    t = window_df['Time(s)']

    axs[0].plot(t, window_df['dFF'], color='darkgreen')
    axs[0].set_ylabel('465 dFF')
    axs[0].set_ylim(max(-1.5, fiberbehav_df['dFF'].min()), min(13, fiberbehav_df['dFF'].max()))
    axs[0].axvspan(center_time, end_time, color='white', zorder=10)

    if has_560:
        axs[1].plot(t, window_df['560 dFF'], color='darkorange')
        axs[1].set_ylabel('560 dFF')
        axs[1].set_ylim(max(-0.3, fiberbehav_df['560 dFF'].min()), min(0.8, fiberbehav_df['560 dFF'].max()))
        axs[1].axvspan(center_time, end_time, color='white', zorder=10)
    
    behavior_colors_path = Path(project_root) / "modules/behaviour/behaviour_colors.json"
    with open(behavior_colors_path, "r") as f:
        behaviors_to_plot = json.load(f)

    pad = 1
    if has_560:
        pad = 2
    
    for i, behavior in enumerate(behavior_cols):
        color, alpha = behaviors_to_plot.get(behavior, ('grey', 0.05))
        behavior_data = window_df[behavior].values
        if np.any(behavior_data > 0):
            axs[i + pad].fill_between(
                t, 0, behavior_data,
                step='post', color=color, alpha=alpha
            )
        axs[i + pad].set_yticks([])
        axs[i + pad].set_ylabel(behavior, rotation=0, labelpad=10, va='center', ha='right')
        axs[i + pad].spines['left'].set_visible(False)
        axs[i + pad].set_ylim(0.1, 1.1)
        axs[i + pad].axvspan(center_time, end_time, color='white', linewidth=3, zorder=10)

    if has_speed:
        global_min = fiberbehav_df['Speed'].dropna().min()
        global_max = fiberbehav_df['Speed'].dropna().max()
        axs[-1].set_ylim(global_min, global_max)
        axs[-1].plot(t, window_df['Speed'], color='black')
        axs[-1].set_ylabel('Speed')
        axs[-1].axvspan(center_time, end_time, color='white', zorder=10)
        
    for ax in axs:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.set_autoscale_on(False)

    for ax in axs[:-1]:
        ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
        
    axs[-1].set_xlabel('Time (s)')

    fig.subplots_adjust(left=0.22, right=0.97, bottom=0.15, top=0.99)

    canvas = FigureCanvas(fig)
    canvas.draw()
    img = np.frombuffer(canvas.buffer_rgba(), dtype='uint8')  
    img = img.reshape(canvas.get_width_height()[::-1] + (4,))  
    img = cv2.cvtColor(img, cv2.COLOR_RGBA2BGR)  
    plt.close(fig)
    del fig, axs, canvas
    gc.collect()
    return img, size

def get_video_time(video_path, file_path, csv_path=None, automated_alignment=False, 
                   bonsai_setup=True, slope=1.0, intercept=0.0):
    print("🔍 Checking video path:", video_path)
    if not os.path.exists(video_path):
        print(f"❌ Video path does not exist: {video_path}")
        return None

    cap = cv2.VideoCapture(str(video_path))
    if not cap.isOpened():
        print(f"❌ Failed to open video: {video_path}")
        return None

    n_frames = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
    print(f'{n_frames} frames to process')
    cap.release()

    if bonsai_setup:
        camera_df = cp.get_timestamps_from_bonsai_csv(csv_path)
        camera_times = camera_df['Time(s)'].values
    elif automated_alignment:
        camera_df = cp.get_camera_flashes_from_csv(csv_path)
        camera_times = camera_df['Time(s)'].values
    else:
        camera_df = cp.get_camera_flashes(file_path)
        camera_times = camera_df['Time(s)'].values

    camera_times_doric = (camera_times - intercept) / slope

    if len(camera_times_doric) != n_frames:
        video_time = np.linspace(camera_times_doric[0], camera_times_doric[-1], n_frames)
    else:
        video_time = camera_times_doric

    return video_time

def align_fiber_to_video(fiber_df, video_time):
    fiber_time = fiber_df['Time(s)'].values
    fiber_indices = np.interp(video_time, fiber_time, np.arange(len(fiber_time)))
    fiber_indices = np.clip(np.round(fiber_indices).astype(int), 0, len(fiber_time) - 1)
    return fiber_indices

def make_combined_video(video_path, 
                        fiberbehav_df, 
                        output_path, 
                        fiber_indices=None, 
                        window=10, 
                        verbose=True, 
                        test=False,
                        start_frame=0,
                        end_frame=None):
    """
    Combines video with fiberphotometry signal and behavior plots, side by side:

    +----------------+-------------------------------+
    |                |   ΔF/F trace                  |
    |   mouse.avi    |   Behavior bars               |
    |                |   Speed (optional)            |
    +----------------+-------------------------------+

    Rolling Gantt logic : For each video frame, Use a fixed window (e.g. last 10 seconds)
    For each behavior, draw a colored bar (1 = filled, 0 = empty) like this:
    Walk:   ███     ████
    Groom:      ████    
    Rear:   █     █     
    """
    if not os.path.exists(video_path):
        raise FileNotFoundError(f"Video file not found: {video_path}")
    
    if fiberbehav_df is None or fiberbehav_df.empty:
        raise ValueError("Provided fiberbehav_df is empty or None.")
    
    if verbose:
        print("All input validations passed. Proceeding with video generation...")

    behavior_cols = [
        col for col in fiberbehav_df.columns
        if set(fiberbehav_df[col].dropna().unique()).issubset({0, 1})
    ]
    if verbose:
        print(f"Behaviors : {behavior_cols}")

    if fiber_indices is None or len(fiber_indices) < 1:
        raise ValueError("fiber_indices must be provided and non-empty.")

    cap = cv2.VideoCapture(str(video_path))
    fps = cap.get(cv2.CAP_PROP_FPS)
    cap.set(cv2.CAP_PROP_POS_FRAMES, start_frame)

    total_frames = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
    n_frames = 200 if test else total_frames
    end_frame = end_frame if end_frame is not None else total_frames

    # Probe overlay size using a dummy frame index
    _, size = create_overlay_frame(10, fiberbehav_df, behavior_cols, window=window)
    width_img, height_img = size  # figsize in inches (width, height)

    vid_width  = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH))
    vid_height = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))

    # Scale overlay so its height matches the video height
    coeff = vid_height / height_img
    overlay_w = make_even(int(width_img * coeff))
    overlay_h = make_even(vid_height)

    combined_w = make_even(vid_width + overlay_w)
    combined_h = make_even(vid_height)

    output_path = str(output_path)
    if not output_path.endswith('.mp4'):
        output_path += '.mp4'

    out = cv2.VideoWriter(
        output_path,
        cv2.VideoWriter_fourcc(*'mp4v'),
        fps,
        (combined_w, combined_h)
    )

    if verbose:
        print(f"🎞️ Video initialized. Size: {combined_w}x{combined_h}, FPS: {fps}, Frames: {n_frames}")

    for i in tqdm(range(start_frame, end_frame), desc="Rendering frames", unit="frame", leave=False):
        ret, frame = cap.read()
        if not ret:
            print(f"[WARNING] Skipping unreadable frame {i}")
            continue

        if i >= len(fiber_indices):
            print(f"[WARNING] No fiber index for frame {i}")
            continue

        data_idx = fiber_indices[i]
        overlay_img, _ = create_overlay_frame(data_idx, fiberbehav_df, behavior_cols, window=window)

        # Resize overlay to (overlay_w, overlay_h) — height matches video
        overlay_img = cv2.resize(overlay_img, (overlay_w, overlay_h))

        combined = np.hstack((frame, overlay_img))
        out.write(combined)

    cap.release()
    out.release()
    print(f"✅ Combined video saved to: {output_path}")

def export_behavior_videos(
    video_path,
    fiberbehav_df,
    fiber_indices,
    behavior_col,
    output_dir,
    window=10,
    pre_time=5,
    post_time=5,
    verbose=True
):
    """
    Export short behavioral clips centered on each behavioral event.

    Parameters
    ----------
    video_path : str
        Path to the raw mouse video (avi/mp4).

    fiberbehav_df : pd.DataFrame
        Fiber + behavior dataframe with 'Time(s)'.

    fiber_indices : ndarray
        Output of align_fiber_to_video().

    behavior_col : str
        Column name of the behavior (binary 0/1).

    output_dir : Path or str
        Directory where clips will be saved.

    window : int
        Window size (passed to create_overlay_frame).

    pre_time : float
        Time before behavior onset (seconds).

    post_time : float
        Time after behavior onset (seconds).
    """
    os.makedirs(output_dir, exist_ok=True)

    # Sampling rate
    sr = pp.samplerate(fiberbehav_df)

    # Detect event onsets
    behavior = fiberbehav_df[behavior_col].values
    onsets = np.where((behavior[:-1] == 0) & (behavior[1:] == 1))[0] + 1

    if verbose:
        print(f"Found {len(onsets)} {behavior_col} onsets.")

    fiber_time = fiberbehav_df["Time(s)"].values
    video_fps = cv2.VideoCapture(video_path).get(cv2.CAP_PROP_FPS)

    for i, onset_idx in enumerate(onsets):
        onset_time = fiber_time[onset_idx]

        # Define clip time window
        t_start = onset_time - pre_time
        t_end   = onset_time + post_time

        # Convert to video frames
        frame_start = np.searchsorted(fiber_time, t_start)
        frame_end   = np.searchsorted(fiber_time, t_end)

        # Map to video frames via fiber_indices
        video_start = np.searchsorted(fiber_indices, frame_start)
        video_end   = np.searchsorted(fiber_indices, frame_end)

        # Make sure valid
        if video_end <= video_start:
            print(f"⚠️ Skipping bout {i}: invalid frame range")
            continue

        # Build output filename
        clip_path = Path(output_dir) / f"{behavior_col}_bout_{i+1}.mp4"

        if verbose:
            print(f"➡️ Exporting bout {i+1} → {clip_path.name} "
                  f"({video_end - video_start} frames)")

        # Call your existing video generator
        make_combined_video(
            video_path=video_path,
            fiberbehav_df=fiberbehav_df,
            output_path=str(clip_path),
            fiber_indices=fiber_indices,
            window=window,
            start_frame=int(video_start),
            end_frame=int(video_end),
            verbose=False
        )

    print("Done exporting all behavioral bout videos.")

def concatenate_videos(video_parts_dir: Path, base_name: str, output_path: Path, delete_temp=False):
    """
    Concatenate multiple video chunks into a single final video using ffmpeg.

    Parameters:
    - video_parts_dir: Path where chunked videos are stored.
    - base_name: Common prefix of chunked video files
    - output_path: Full path to the final output video
    - delete_temp: If True, deletes intermediate part files after concatenation.
    """
    # Find matching files
    video_parts = sorted(video_parts_dir.glob(f"{base_name}*.mp4"))
    if not video_parts:
        raise FileNotFoundError(f"No video chunks found matching: {base_name}*.mp4")

    print(f"Found {len(video_parts)} chunks to concatenate.")

    # Create a temporary text file for ffmpeg
    list_file = video_parts_dir / "concat_list.txt"
    with open(list_file, "w", encoding='utf-8') as f:
        for part in video_parts:
            f.write(f"file '{part.as_posix()}'\n")

    print(f"File list written to: {list_file}")

    # Run ffmpeg to concatenate
    cmd = [
        ffmpeg_path,
        "-f", "concat",
        "-safe", "0",
        "-i", str(list_file),
        "-c", "copy",
        str(output_path)
    ]

    print(f"Running ffmpeg to generate: {output_path}")
    subprocess.run(cmd, check=True)
    print(f"✅ Final video saved to: {output_path}")

    # Clean up
    list_file.unlink()
    if delete_temp:
        print("Deleting chunk files...")
        for part in video_parts:
            part.unlink()

#%%

if __name__ == "__main__":
        ##### TO BE CHANGED #####
    batch = 1
    for mouse in ['991','1014']:    
        exp='RewardAirpuffs'
        behavior = "Licks_filtered"
        data_path = Path(r'F:\202603_FiberGRABACh-FlexRGECO\Data\20260409_RewardAirpuff')
        analysis_path = Path(r'F:\202603_FiberGRABACh-FlexRGECO\Analysis\RewardAirpuffs\length0_interbout3_o4fNone')
        #########################

        video_name = f'{mouse}.avi'
        pp_path = data_path / 'Preprocessing'
        video_path = data_path / 'Behaviour' / f'{video_name}'
        raw_doric_path = data_path / f'{mouse}_0000.doric'
        fiberbehav_df = pd.read_csv(analysis_path / f'{batch}_{mouse}_fiberbehavnotderived.csv')
        output_path = analysis_path / f'Videos_{behavior}/{batch}_{mouse}' / f'{video_name[:-4]}_combined'
        camera_csv_path = data_path / f'camera_flashes_{mouse}.csv'
        led_flashes_path = data_path / f'miniscope_sync_{mouse}.csv'


        print(f"{mouse}")
        # Extract sync channel from Doric raw data and Bonsai corresponding sync data
        ttl_sync_df = cp.extract_sync_channel(raw_doric_path, sync_channel = "DIO04")
        led_df = cp.get_timestamps_from_bonsai_csv(led_flashes_path) 
        # Compute linear regression to correct for differences between clocks
        print(f"Syncing Doric and Bonsai clocks")
        slope, intercept = cp.time_mapping(ttl_sync_df, led_df)

        video_time = get_video_time(video_path,
                                    raw_doric_path,
                                    csv_path=camera_csv_path,
                                    automated_alignment=False, 
                                    bonsai_setup=True,
                                    slope=slope,
                                    intercept=intercept)
        
        # Drop frames with no corresponding fiber signal
        fiber_start_time = fiberbehav_df['Time(s)'].iloc[0]
        valid_frame_indices = np.where(video_time >= fiber_start_time)[0]
        video_time_trimmed = video_time[valid_frame_indices]

        # Align fiber data to trimmed video timestamps
        fiber_indices = align_fiber_to_video(fiberbehav_df, video_time_trimmed)
        
        export_behavior_videos(
            video_path,
            fiberbehav_df,
            fiber_indices,
            behavior_col = behavior,
            output_dir = analysis_path / f'Videos_{behavior}/{batch}_{mouse}',
            window = 10,
            pre_time = 5,
            post_time = 10
        )
        '''

        chunk_size = int(15 * 20)  # 30 seconds * 30 FPS
        total_frames = len(video_time)
        test = False

        if test :
            start = 7500
            end = 8500
            print(f"Processing frames {start} to {end}...")
            
            chunk_output = output_path.parent / f"{output_path.stem}_test"
            
            make_combined_video(
                video_path=video_path,
                fiberbehav_df=fiberbehav_df,
                output_path=chunk_output,
                fiber_indices=fiber_indices,
                window=10,
                verbose=True,
                test=test,
                start_frame=start,
                end_frame=end
            )

        else:
            for start in range(0, total_frames, chunk_size):
                end = min(start + chunk_size, total_frames)
                print(f"Processing frames {start} to {end}...")
                
                chunk_output = output_path.parent / f"{output_path.stem}_part{start//chunk_size}"
                
                make_combined_video(
                    video_path=video_path,
                    fiberbehav_df=fiberbehav_df,
                    output_path=chunk_output,
                    fiber_indices=fiber_indices,
                    window=10,
                    verbose=False,
                    test=test,
                    start_frame=start,
                    end_frame=end
                )
    
            # Concatenate videos 
            video_parts_dir = exp_path / 'Videos'
            base_name = f"{mouse}_combined_part"
            output_path = video_parts_dir / f"{mouse}_combined_full.mp4"

            concatenate_videos(video_parts_dir, base_name, output_path, delete_temp=True)
    '''
    # %%
