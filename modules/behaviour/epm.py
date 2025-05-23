import numpy as np
import matplotlib.pyplot as plt

import modules.behaviour.mouse_position as mp

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

def analyze_mouse_position(dff, t, coords, epm_coordinates):
    x = coords[:, 1]
    y = coords[:, 2]

    x1, x2 = epm_coordinates['x1'], epm_coordinates['x2']
    y1, y2 = epm_coordinates['y1'], epm_coordinates['y2']
    minx, maxx = epm_coordinates['minx'], epm_coordinates['maxx']
    miny, maxy = epm_coordinates['miny'], epm_coordinates['maxy']

    closed_arm, open_arm, center = classify_position(x, y, x1, x2, y1, y2)
    speed = mp.compute_speed(x, y)
    heatmap = mp.create_heatmap(x, y, dff, minx, maxx, miny, maxy)
    plot_results(t, dff, closed_arm, open_arm, center, x, y, x1, x2, y1, y2, heatmap)

    return fibertracking_df

def classify_position(x, y, x1, x2, y1, y2):
    closed_arm = ((y <= y1) | (y >= y2)) & (x >= x1) & (x <= x2)
    open_arm = (~closed_arm) & ((x <= x1) | (x >= x2))
    center = (~closed_arm) & (~open_arm)
    return closed_arm.astype(int), open_arm.astype(int), center.astype(int)
