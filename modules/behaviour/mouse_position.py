import numpy as np
import matplotlib.pyplot as plt
import cv2
from scipy.ndimage import gaussian_filter
from scipy.signal import savgol_filter

def load_video_frame(video_path):
    cap = cv2.VideoCapture(video_path)
    ret, frame = cap.read()
    cap.release()
    if ret:
        return cv2.cvtColor(frame, cv2.COLOR_BGR2RGB)
    else:
        raise FileNotFoundError("Video frame could not be read.")
    
def compute_speed(x, y, dist_scale=30, frame_rate=20):
    dx = np.diff(x)
    dy = np.diff(y)
    distance = dist_scale * np.sqrt(dx**2 + dy**2)
    speed = distance / frame_rate
    speed = savgol_filter(speed, 5, 2)  # smoothing
    speed[speed > 1] = 1
    return speed

def create_heatmap(x, y, dff, minx, maxx, miny, maxy):
    heat_map = np.zeros((101, 101))
    counts = np.zeros((101, 101))

    for i in range(len(x)):
        col = int((x[i] - minx) * 100 / (maxx - minx))
        row = int((y[i] - miny) * 100 / (maxy - miny))
        if 0 <= row < 101 and 0 <= col < 101:
            heat_map[row, col] += dff[i]
            counts[row, col] += 1

    with np.errstate(invalid='ignore'):
        heat_map = np.divide(heat_map, counts, where=counts != 0)
    heat_map[heat_map == 0] = np.nan
    filtered_map = gaussian_filter(heat_map, sigma=3)
    return filtered_map