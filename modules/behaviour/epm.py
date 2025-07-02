import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import modules.behaviour.mouse_position as mp

def plot_pie(fiberbehav_df):

    open_arm = fiberbehav_df['Open arm'].sum()
    closed_arm = fiberbehav_df['Closed arm'].sum()
    center = fiberbehav_df['Center'].sum()
    all = sum(open_arm, closed_arm, center)

    fig = plt.plots(figsize=(8, 8))

    # Time spent pie
    sizes = [open_arm/all * 100, closed_arm/all * 100, center/all * 100]
    fig = plt.pie(sizes, labels=['Open Arm', 'Closed Arm', 'Center'], colors=['goldenrod', 'gray', 'white'])
    fig.set_title('Time Spent')

    return fig

def plot_heatmap(t, dff, closed_arm, open_arm, center, x, y, x1, x2, y1, y2, heatmap):
     
    return
    
def analyze_mouse_position(coords, epm_coordinates, bodypart='head'):
    coords_x = coords[f'{bodypart}_x']
    coords_y = coords[f'{bodypart}_y']

    x1, x2 = epm_coordinates['x1'], epm_coordinates['x2']
    y1, y2 = epm_coordinates['y1'], epm_coordinates['y2']
    rotation_angle = epm_coordinates['rotation angle']

    coords_byzone_df = classify_position(coords_x, coords_y, x1, x2, y1, y2, rotation_angle)
    speed_df = mp.compute_speed(coords)

    behav_df = pd.concat([coords, coords_byzone_df, speed_df], axis=1)

    return behav_df

def classify_position(coords_x, coords_y, x1, x2, y1, y2, rotation_angle):
    closed_arm = np.zeros(len(coords_x))
    open_arm = np.zeros(len(coords_x))
    center = np.zeros(len(coords_x))
    if rotation_angle == 0:
        for i, (x, y) in enumerate(zip(coords_x, coords_y)): 
            if ((x <= x1) or (x >= x2)) and (y2 >= y >= y1):
                closed_arm[i]=1
            elif (x1 < x < x2) and (y1 <= y <= y2):
                center[i]=1
            else:
                open_arm[i]=1
    elif rotation_angle == 90:
            if ((y <= y1) or (y >= y2)) and (x2 >= x >= x1):
                closed_arm[i] = 1
            elif (y1 < y < y2) and (x1 <= x <= x2):
                center[i] = 1
            else:
                open_arm[i] = 1
    else:
            raise ValueError("Unsupported rotation angle. Use 0 or 90.")

    coords_byzone_df = pd.DataFrame(
         {'Closed arm':closed_arm,
          'Open arm':open_arm,
          'Center':center}
    )
    return coords_byzone_df
