import cv2
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import json
import os
from pathlib import Path

def get_click_coordinates(image, n_points=2, title='Click to select points'):
    plt.imshow(cv2.cvtColor(image, cv2.COLOR_BGR2RGB))
    plt.title(title)
    plt.axis('on')
    coords = plt.ginput(n_points, timeout=-1)
    plt.close()
    return coords

def define_epm_boundaries(video_path):
    # Read 10th frame
    cap = cv2.VideoCapture(video_path)
    cap.set(cv2.CAP_PROP_POS_FRAMES, 9)  
    ret, frame = cap.read()
    cap.release()
    
    if not ret:
        raise ValueError("Could not read frame from video.")

    # Ask if orientation is correct
    plt.imshow(cv2.cvtColor(frame, cv2.COLOR_BGR2RGB))
    plt.title("Is the EPM orientation correct? (Open up/down, closed right/left) (Y=Yes, N=No)")
    plt.axis("off")
    plt.show(block=True)
    plt.pause(0.5)
    answer = input("Is the EPM in the correct orientation? (Open up/down, closed right/left) (Y/N): ").strip().upper()
    plt.close('all')

    # Select open_arms (bottom-left and top-right corners)
    open_arms = get_click_coordinates(frame, 2, "Click open arms bottom-left, then top-right")
    open_xL, open_yBot = open_arms[0]
    open_xR, open_yTop = open_arms[1]
    open_xL, open_xR = sorted([open_xL, open_xR])
    open_yBot, open_yTop = sorted([open_yBot, open_yTop])

    # Select closed_arms (bottom-left and top-right corners)
    closed_arms = get_click_coordinates(frame, 2, "Click closed arms bottom-left, then top-right")
    closed_xL, closed_yBot = closed_arms[0]
    closed_xR, closed_yTop = closed_arms[1]
    closed_xL, closed_xR = sorted([closed_xL, closed_xR])
    closed_yBot, closed_yTop = sorted([closed_yBot, closed_yTop])

    # Center coordinates
    if answer == 'N': #closed arm top/down
        center_xL = closed_xL
        center_xR = closed_xR
        center_yBot = open_yBot
        center_yTop = open_yTop
    
    else:
        center_xL = open_xL
        center_xR = open_xR
        center_yBot = closed_yBot
        center_yTop = closed_yTop

    # Show confirmation
    fig, ax = plt.subplots()
    ax.imshow(cv2.cvtColor(frame, cv2.COLOR_BGR2RGB))

    # Open arms rectangle (green)
    open_rect = patches.Rectangle(
        (open_xL, open_yBot),
        open_xR - open_xL,
        open_yTop - open_yBot,
        linewidth=2,
        edgecolor='yellow',
        facecolor='yellow',
        alpha=0.3,
        label='Open Arms'
    )
    ax.add_patch(open_rect)

    # Closed arms rectangle (red)
    closed_rect = patches.Rectangle(
        (closed_xL, closed_yBot),
        closed_xR - closed_xL,
        closed_yTop - closed_yBot,
        linewidth=2,
        edgecolor='red',
        facecolor='red',
        alpha=0.3,
        label='Closed Arms'
    )
    ax.add_patch(closed_rect)

    # Center rectangle (blue)
    center_rect = patches.Rectangle(
        (center_xL, center_yBot),
        center_xR - center_xL,
        center_yTop - center_yBot,
        linewidth=2,
        edgecolor='green',
        facecolor='green',
        alpha=0.3,
        label='Center'
    )
    ax.add_patch(center_rect)

    ax.legend()
    plt.title("EPM Boundary Confirmation")
    plt.axis("off")
    plt.show(block=True)

    epm_coordinates = {
        'open_xL': open_xL, 'open_xR': open_xR,
        'open_yBot': open_yBot, 'open_yTop': open_yTop,
        'closed_xL': closed_xL, 'closed_xR': closed_xR,
        'closed_yBot': closed_yBot, 'closed_yTop': closed_yTop,
        'center_xL': center_xL, 'center_xR': center_xR,
        'center_yBot': center_yBot, 'center_yTop': center_yTop
    }

    return epm_coordinates

def save_boundaries_to_json(boundaries, output_path):
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(boundaries, f, indent=4)
    print(f"Boundaries saved to {output_path}")

"""
if __name__ == "__main__":
    video_folder = r'E:\FiberPhotometry\202404_DualColourGRABAChxFlexGECO\Data\20240410_EPM\Videos'
    video_name = '466'
    video_path = f'{video_folder}\{video_name}.avi'
    output_json = f'{video_folder}\{video_name}_epm_boundaries.json'

    boundaries = define_epm_boundaries(video_path)
    print("\nEPM Boundaries:")
    for key, value in boundaries.items():
        print(f"{key}: {value:.2f}")
    
    save_boundaries_to_json(boundaries, output_json)
"""
