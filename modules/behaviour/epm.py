import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import modules.behaviour.mouse_position as mp

def plot_pie(coords_byzone_df):
    fig, axs = plt.subplots(2, 4, figsize=(20, 8))

    # Time spent pie
    sizes = [np.mean(open_arm) * 100, np.mean(closed_arm) * 100, np.mean(center) * 100]
    axs[1, 1].pie(sizes, labels=['Open Arm', 'Closed Arm', 'Center'], colors=['goldenrod', 'gray', 'white'])
    axs[1, 1].set_title('Time Spent')

def plot_heatmap(t, dff, closed_arm, open_arm, center, x, y, x1, x2, y1, y2, heatmap):

    # Heatmap
    im = axs[1, 3].imshow(heatmap, cmap='viridis', interpolation='none')
    plt.colorbar(im, ax=axs[1, 3])
    axs[1, 3].set_title('Heatmap')
    
def analyze_mouse_position(coords, epm_coordinates, bodypart='head'):
    coords_x = coords[f'{bodypart}_x']
    coords_y = coords[f'{bodypart}_y']

    x1, x2 = epm_coordinates['x1'], epm_coordinates['x2']
    y1, y2 = epm_coordinates['y1'], epm_coordinates['y2']

    coords_byzone_df = classify_position(coords_x, coords_y, x1, x2, y1, y2)
    speed_df = mp.compute_speed(coords_x, coords_y)

    behav_df = pd.concat([coords, coords_byzone_df, speed_df], axis=1)

    return behav_df

def classify_position(coords_x, coords_y, x1, x2, y1, y2):
    closed_arm = np.zeros(len(coords_x))
    open_arm = np.zeros(len(coords_x))
    center = np.zeros(len(coords_x))
    for i,x,y in enumerate(zip(coords_x, coords_y)): 
            if ((x<=x1) or (x>=x2)) and (y2>=y>=y1):
                closed_arm[i]=1
            elif ((y2>=y>=y1) and (y2>=y>=y1)):
                center[i]=1
            else:
                open_arm[i]=1
    coords_byzone_df = pd.DataFrame(
         {'Closed arm':closed_arm,
          'Open arm':open_arm,
          'Center':center}
    )
    return coords_byzone_df
