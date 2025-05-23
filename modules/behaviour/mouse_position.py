import numpy as np
import matplotlib.pyplot as plt
import cv2
import os
import pickle
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

def classify_position(x, y, x1, x2, y1, y2):
    closed_arm = ((y <= y1) | (y >= y2)) & (x >= x1) & (x <= x2)
    open_arm = (~closed_arm) & ((x <= x1) | (x >= x2))
    center = (~closed_arm) & (~open_arm)
    return closed_arm.astype(int), open_arm.astype(int), center.astype(int)

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

def plot_results(t, dff, closed_arm, open_arm, center, x, y, x1, x2, y1, y2, heatmap):
    fig, axs = plt.subplots(2, 4, figsize=(20, 8))

    # ΔF/F trace
    axs[0, 0].plot(t, dff, 'k')
    axs[0, 0].fill_between(t, -100, 100, where=closed_arm==1, color='gray', alpha=0.5)
    axs[0, 0].fill_between(t, -100, 100, where=open_arm==1, color='goldenrod', alpha=0.5)
    axs[0, 0].set_title('ΔF/F Trace')

    # XY path
    axs[1, 0].plot(x, y, 'k')
    axs[1, 0].axvline(x1); axs[1, 0].axvline(x2)
    axs[1, 0].axhline(y1); axs[1, 0].axhline(y2)
    axs[1, 0].set_title('XY Path')

    # Time spent pie
    sizes = [np.mean(open_arm) * 100, np.mean(closed_arm) * 100, np.mean(center) * 100]
    axs[1, 1].pie(sizes, labels=['Open Arm', 'Closed Arm', 'Center'], colors=['goldenrod', 'gray', 'white'])
    axs[1, 1].set_title('Time Spent')

    # Heatmap
    im = axs[1, 3].imshow(heatmap, cmap='viridis', interpolation='none')
    plt.colorbar(im, ax=axs[1, 3])
    axs[1, 3].set_title('Heatmap')
    
    plt.tight_layout()
    plt.show()

    def analyze_mouse_position(dff, t, coords, video_path):
    x = coords[:, 1]
    y = coords[:, 2]

    frame = load_video_frame(video_path)

    # You would replace this with real UI input in practice
    x1, x2 = 100, 300
    y1, y2 = 150, 350
    minx, maxx, miny, maxy = 0, 400, 0, 400

    closed_arm, open_arm, center = classify_position(x, y, x1, x2, y1, y2)
    speed = compute_speed(x, y)
    heatmap = create_heatmap(x, y, dff, minx, maxx, miny, maxy)
    plot_results(t, dff, closed_arm, open_arm, center, x, y, x1, x2, y1, y2, heatmap)