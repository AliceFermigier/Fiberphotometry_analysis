import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cv2
from scipy.ndimage import gaussian_filter
from scipy.signal import savgol_filter
import importlib

import modules.common.clean_signal as cs
importlib.reload(cs)

def get_dlc_data(data_path, threshold=0.99, interpolate=True):
    '''
    Gets DLC data and filters it based on likelihood (default = 0.99)
    Creates new filtered csv with 
    '''
    dlc_df = pd.read_csv(data_path, header=[1])
    coords_array=[]
    for bodypart in dlc_df.columns.tolist()[1::3]:
        x_filtered = []
        y_filtered = []
        for x,y,likelihood in zip(dlc_df[f'{bodypart}'][1:], dlc_df[f'{bodypart}.1'][1:], dlc_df[f'{bodypart}.2'][1:]):
            if float(likelihood)>=threshold:
                x_filtered.append(x)
                y_filtered.append(y)
            else:
                x_filtered.append(np.nan)
                y_filtered.append(np.nan)

        bodypart_df = pd.DataFrame({
            f'{bodypart}_x': pd.to_numeric(x_filtered, errors='coerce'),
            f'{bodypart}_y': pd.to_numeric(y_filtered, errors='coerce')
        })

        if interpolate:
            bodypart_df = bodypart_df.interpolate(method='linear', limit_direction='both')
        
        coords_array.append(bodypart_df)

    coordinates_df = pd.concat(coords_array, axis=1)
        
    return coordinates_df

def load_video_frame(video_path):
    cap = cv2.VideoCapture(video_path)
    ret, frame = cap.read()
    cap.release()
    if ret:
        return cv2.cvtColor(frame, cv2.COLOR_BGR2RGB)
    else:
        raise FileNotFoundError("Video frame could not be read.")   
    
def compute_speed(coordinates_df, dist_scale=0.1322, frame_rate=19, bodypart='center'):
    '''
    dist_scale in cm/px
    frame_rate in fps
    '''
    dx = np.diff(coordinates_df[f'{bodypart}_x'])
    dy = np.diff(coordinates_df[f'{bodypart}_y'])
    distance = dist_scale * np.sqrt(dx**2 + dy**2)
    speed = distance * frame_rate
    speed, _ = cs.hampel_filter(speed, window_size=15) # remove big artifacts
    speed = savgol_filter(speed, 5, 2)  # smoothing

    speed = np.concatenate([[0], speed])  # pad first frame with 0
    speed_df = pd.DataFrame({'Speed': speed}, index=coordinates_df.index)

    return speed_df
 
def create_heatmap(coords_x, coords_y, dff, minx, maxx, miny, maxy):
    heat_map = np.zeros((101, 101))
    counts = np.zeros((101, 101))

    for i in range(len(coords_x)):
        col = int((coords_x[i] - minx) * 100 / (maxx - minx))
        row = int((coords_y[i] - miny) * 100 / (maxy - miny))
        if 0 <= row < 101 and 0 <= col < 101:
            heat_map[row, col] += dff[i]
            counts[row, col] += 1

    with np.errstate(invalid='ignore'):
        heat_map = np.divide(heat_map, counts, where=counts != 0)
    heat_map[heat_map == 0] = np.nan
    filtered_map = gaussian_filter(heat_map, sigma=3)
    return filtered_map