import cv2
import numpy as np
import matplotlib.pyplot as plt
import json
import os

def get_click_coordinates(image, n_points=1, title='Click to select point'):
    plt.imshow(cv2.cvtColor(image, cv2.COLOR_BGR2RGB))
    plt.title(title)
    plt.axis('on')
    coords = plt.ginput(n_points, timeout=-1)
    plt.close()
    return coords

def define_ports(video_path):
    """
    Click:
    1) licking port
    2) airpuff port left
    3) airpuff port right
    """
    cap = cv2.VideoCapture(video_path)
    cap.set(cv2.CAP_PROP_POS_FRAMES, 0)
    ret, frame = cap.read()
    cap.release()

    if not ret:
        raise ValueError("Could not read frame from video.")

    # ---- LICK PORT ----
    lick_port = get_click_coordinates(frame, 1,"Click the LICKING PORT")[0]

    # ---- AIRPUFF PORT LEFT ----
    air_left = get_click_coordinates(frame, 1,"Click the LEFT AIRPUFF PORT")[0]

    # ---- AIRPUFF PORT RIGHT ----
    air_right = get_click_coordinates(frame, 1,"Click the RIGHT AIRPUFF PORT")[0]

    # Show confirmation
    plt.imshow(cv2.cvtColor(frame, cv2.COLOR_BGR2RGB))
    plt.scatter([lick_port[0]], [lick_port[1]], c='g', label='Lick Port')
    plt.scatter([air_left[0]], [air_left[1]], c='r', label='Air Left')
    plt.scatter([air_right[0]], [air_right[1]], c='b', label='Air Right')
    plt.legend()
    plt.title("Port locations")
    plt.axis('on')
    plt.show(block=True)

    ports = {
        "lick_port": {"x": lick_port[0], "y": lick_port[1]},
        "airpuff_left": {"x": air_left[0], "y": air_left[1]},
        "airpuff_right": {"x": air_right[0], "y": air_right[1]},
    }

    return ports

def save_ports_to_json(ports, output_path):
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(ports, f, indent=4)
    print(f"Port coordinates saved to {output_path}")
