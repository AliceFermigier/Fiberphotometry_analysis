import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import plotly.express as px
import modules.common.preprocess as pp

def txt_to_df(capacitance_txt_path):
    # Load all columns
    df = pd.read_csv(
        capacitance_txt_path,
        header=None,
        names=["time(ms)", "Capacitance", "recording"],
        low_memory=False
    )

    # Find indices for recording start and stop
    try:
        start_idx = df.index[df["recording"] == "RECORDING_START"][0]
        stop_idx = df.index[df["recording"] == "RECORDING_STOP"][0]
    except IndexError:
        raise ValueError("RECORDING_START or RECORDING_STOP not found in file")

    # Slice the dataframe to include only the recording period
    df = df.loc[start_idx:stop_idx]

    # Reset time to start at 0 and convert to seconds
    df["Time(s)"] = (df["time(ms)"] - df["time(ms)"].iloc[0]) / 1000

    # Drop the original millisecond column
    df = df[["Time(s)", "Capacitance"]]

    # Convert capacitance to integer
    df["Capacitance"] = df["Capacitance"].astype(int)

    return df

def extract_lick_bouts(capacitance_df, threshold):
    licks_df = (
        capacitance_df[capacitance_df["Capacitance"] > threshold][["Time(s)"]]
        .reset_index(drop=True)
    )
    return licks_df

def plot_licks_and_threshold(capacitance_df, threshold=500, downsample=100):

    capacitance_df_downsampled = pp.downsample(capacitance_df, target_frequency=downsample)

    fig = px.line(
        capacitance_df_downsampled,
        x="Time(s)",
        y="Capacitance",
        title="Capacitance over time",
        labels={"Time(s)": "Time (s)", "Capacitance": "Capacitance"},
    )

    # Threshold line
    fig.add_hline(
        y=threshold,
        line_dash="dash",
        line_color="red",
        annotation_text=f"Threshold ({threshold})",
        annotation_position="top right",
    )

    fig.update_layout(
        xaxis_title="Time (s)",
        yaxis_title="Capacitance",
        hovermode="x unified",   # shows all values at a given x on hover
    )

    fig.show()

def compute_distance(df, port, scale, nose_x='nose_x', nose_y='nose_y'):
    """
    Compute distance from nose to port in centimeters.
    
    Parameters:
    -----------
    df : DataFrame
        Contains nose position columns (in pixels)
    port : dict
        Port coordinates with 'x' and 'y' keys (in pixels)
    scale : float
        Scale factor from scale_and_coords["Scale_cm_per_px"]
    nose_x, nose_y : str
        Column names for nose coordinates (in pixels)
    
    Returns:
    --------
    Distance in centimeters
    """
    px = port["x"]
    py = port["y"]
    pixel_dist = np.sqrt((df[nose_x] - px)**2 + (df[nose_y] - py)**2)
    return pixel_dist * scale  # Convert pixels to cm

def filter_licking(df, ports, scale_and_coords, lick_col='Licks', lick_radius_cm=0.5):
    """
    Filter licks based on distance from lick port.
    
    Parameters:
    -----------
    lick_radius_cm : float
        Radius in centimeters within which a lick is valid 
    """
    video_scale = scale_and_coords["Scale_cm_per_px"]
    dist = compute_distance(df, ports["lick_port"], video_scale)
    true_lick = (df[lick_col] == 1) & (dist < lick_radius_cm)
    df["Licks_filtered"] = true_lick.astype(int)
    return df

def detect_airpuff_entry(df, ports, scale_and_coords, radius_cm=3.0):
    """
    Detect when nose enters airpuff ports.
    
    Parameters:
    -----------
    radius_cm : float
        Radius in centimeters defining port entry zone
    """
    video_scale = scale_and_coords["Scale_cm_per_px"]
    distL = compute_distance(df, ports["airpuff_left"], video_scale)
    distR = compute_distance(df, ports["airpuff_right"], video_scale)
    df["Nose_in_any_airport"] = ((distL < radius_cm) | (distR < radius_cm)).astype(int)
    return df