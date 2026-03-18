import cv2
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import json
import os

class RectangleSelector:
    """Handles click-drag-release rectangle selection on a Matplotlib figure."""
    def __init__(self, ax):
        self.ax = ax
        self.rect = None
        self.x0 = None
        self.y0 = None
        self.x1 = None
        self.y1 = None
        self.press_cid = ax.figure.canvas.mpl_connect("button_press_event", self.on_press)
        self.motion_cid = ax.figure.canvas.mpl_connect("motion_notify_event", self.on_motion)
        self.release_cid = ax.figure.canvas.mpl_connect("button_release_event", self.on_release)
        self.finished = False

    def on_press(self, event):
        if event.inaxes != self.ax:
            return
        self.x0, self.y0 = event.xdata, event.ydata

        if self.rect is None:
            self.rect = Rectangle((self.x0, self.y0), 0, 0,
                                  fill=False, edgecolor="red", linewidth=2)
            self.ax.add_patch(self.rect)

    def on_motion(self, event):
        if self.rect is None or event.inaxes != self.ax or self.x0 is None:
            return

        x1, y1 = event.xdata, event.ydata
        self.rect.set_width(x1 - self.x0)
        self.rect.set_height(y1 - self.y0)
        self.rect.set_xy((self.x0, self.y0))
        self.ax.figure.canvas.draw_idle()

    def on_release(self, event):
        if event.inaxes != self.ax:
            return
        self.x1, self.y1 = event.xdata, event.ydata
        self.finished = True

        plt.close()  # close figure when done

    def get_rectangle(self):
        if not self.finished:
            return None
        return {
            "x1": float(self.x0),
            "y1": float(self.y0),
            "x2": float(self.x1),
            "y2": float(self.y1)
        }

def get_scale_and_arena_rect(video_path, real_world_distance_cm, real_world_distance_name="arena long side", frame_number=1000):
    """
    First selects two points for scale, then lets the user draw a rectangle
    by click-drag-release for arena boundaries.
    Real word distance corresponds to a known distance
    """

    # ---------------------------
    # Load frame
    # ---------------------------
    cap = cv2.VideoCapture(video_path)
    cap.set(cv2.CAP_PROP_POS_FRAMES, frame_number)
    success, frame = cap.read()
    cap.release()
    if not success:
        raise ValueError(f"Cannot read frame {frame_number} from {video_path}")

    frame_rgb = cv2.cvtColor(frame, cv2.COLOR_BGR2RGB)

    # ---------------------------
    # 1️⃣ SCALE (2 clicks)
    # ---------------------------
    plt.imshow(frame_rgb)
    plt.title(f"Click two extremities of {real_world_distance_name} with known distance ({real_world_distance_cm} cm)")
    scale_points = plt.ginput(2, timeout=0)
    plt.close()

    if len(scale_points) != 2:
        raise ValueError("You must click exactly 2 points for scale.")

    (x1, y1), (x2, y2) = scale_points
    pixel_dist = np.hypot(x2 - x1, y2 - y1)
    scale = real_world_distance_cm / pixel_dist

    # ---------------------------
    # 2️⃣ ARENA BOUNDARY (rectangle drag)
    # ---------------------------
    fig, ax = plt.subplots()
    ax.imshow(frame_rgb)
    ax.set_title("Drag to draw arena rectangle")

    selector = RectangleSelector(ax)
    plt.show(block=True)  # waits until user draws

    rect = selector.get_rectangle()
    if rect is None:
        raise RuntimeError("Rectangle was not drawn.")

    # ---------------------------
    # Prepare output
    # ---------------------------
    output = {
        "Scale_cm_per_px": float(scale),
        "Scale_points_px": {
            "P1": {"x": float(x1), "y": float(y1)},
            "P2": {"x": float(x2), "y": float(y2)}
        },
        "Arena_rectangle_px": rect
    }

    return output

def save_to_json(input_dict, output_path):
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(input_dict, f, indent=4)
    print(f"Video scale and arena coordinates saved to {output_path}")