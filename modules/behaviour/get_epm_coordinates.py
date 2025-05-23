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
    plt.title("Is the EPM orientation correct? (Y=Yes, N=No)")
    plt.axis("off")
    plt.show()
    answer = input("Is the EPM in the correct orientation? (Y/N): ").strip().upper()

    if answer == 'N':
        frame = cv2.rotate(frame, cv2.ROTATE_90_CLOCKWISE)  # Adjust rotation if needed

    # Select center zone (bottom-left and top-right corners)
    center_points = get_click_coordinates(frame, 2, "Click center bottom-left, then top-right")
    x1, y2 = center_points[0]
    x2, y1 = center_points[1]
    x1, x2 = sorted([x1, x2])
    y1, y2 = sorted([y1, y2])

    # Select entire EPM bounding box
    area_points = get_click_coordinates(frame, 2, "Click maze bottom-left, then top-right")
    minx, miny = area_points[0]
    maxx, maxy = area_points[1]
    minx, maxx = sorted([minx, maxx])
    miny, maxy = sorted([miny, maxy])

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
        'miny': miny, 'maxy': maxy
    }

    return epm_coordinates

def save_boundaries_to_json(boundaries, output_path):
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(boundaries, f, indent=4)
    print(f"Boundaries saved to {output_path}")

if __name__ == "__main__":
    video_path = "your_video_file.avi"  # <-- Replace this with your video filename
    output_json = "epm_boundaries.json" # <-- Replace or update path as needed

    boundaries = define_epm_boundaries(video_path)
    print("\nEPM Boundaries:")
    for key, value in boundaries.items():
        print(f"{key}: {value:.2f}")
    
    save_boundaries_to_json(boundaries, output_json)
