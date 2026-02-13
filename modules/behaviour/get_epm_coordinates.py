import cv2
import numpy as np
import matplotlib.pyplot as plt
import json
import os

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
    plt.show()
    answer = input("Is the EPM in the correct orientation? (Open up/down, closed right/left) (Y/N): ").strip().upper()

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
    plt.imshow(cv2.cvtColor(frame, cv2.COLOR_BGR2RGB))
    plt.axvline(x1, color='g', label='x1/x2 (Center)')
    plt.axvline(x2, color='g')
    plt.axhline(y1, color='g')
    plt.axhline(y2, color='g')
    plt.axvline(minx, color='r', linestyle='--', label='min/max x (Maze)')
    plt.axvline(maxx, color='r', linestyle='--')
    plt.axhline(miny, color='r', linestyle='--', label='min/max y (Maze)')
    plt.axhline(maxy, color='r', linestyle='--')
    plt.title("EPM Boundaries")
    plt.legend()
    plt.axis('on')
    plt.show()

    epm_coordinates = {
        'x1': x1, 'x2': x2,
        'y1': y1, 'y2': y2,
        'minx': minx, 'maxx': maxx,
        'miny': miny, 'maxy': maxy,
        'rotation angle' : rotation_angle
    }

    return epm_coordinates

def save_boundaries_to_json(boundaries, output_path):
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(boundaries, f, indent=4)
    print(f"Boundaries saved to {output_path}")

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
