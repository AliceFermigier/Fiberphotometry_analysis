import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def txt_to_df(capacitance_txt_path):
    # Load all columns
    df = pd.read_csv(
        capacitance_txt_path,
        header=None,
        names=["time(ms)", "capacitance", "recording"],
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
    df["time(s)"] = (df["time(ms)"] - df["time(ms)"].iloc[0]) / 1000

    # Drop the original millisecond column
    df = df[["time(s)", "capacitance"]]

    # Convert capacitance to integer
    df["capacitance"] = df["capacitance"].astype(int)

    return df

def load_mouse_data(mouse, batch, datapath_exp_dict):
    data_path_exp = datapath_exp_dict[batch]
    capacitance_txt_path = data_path_exp / f"{mouse}.txt"

    df = txt_to_df(capacitance_txt_path)
    return df

def extract_lick_bouts(licks_df, threshold):
    lick_bouts_df = pd.Dataframe(data = {"time(ms)":licks_df["time(ms)"],"licks":np.zeros(len(licks_df["time(ms)"]))})

    return lick_bouts_df

def plot_licks_and_threshold(licks_df, threshold):

    plt.plot(
        licks_df["time(ms)"],
        licks_df["capacitance"],
        linestyle="-",
        label="Capacitance"
    )

    # Add horizontal threshold line
    plt.axhline(
        y=threshold,
        color="red",
        linestyle="--",
        linewidth=1.5,
        label=f"Threshold = {threshold}"
    )

    plt.xlabel("Time (ms)")
    plt.ylabel("Capacitance")
    plt.grid(True)

    # Clean legend (avoid duplicates)
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys())

    plt.show()

def compute_distance(df, port, nose_x='nose_x', nose_y='nose_y'):
    px = port["x"]
    py = port["y"]
    return np.sqrt((df[nose_x] - px)**2 + (df[nose_y] - py)**2)

def filter_licking(df, ports, lick_col='Licks', lick_radius=20):
    dist = compute_distance(df, ports["lick_port"])
    true_lick = (df[lick_col] == 1) & (dist < lick_radius)
    df["Licks_filtered"] = true_lick.astype(int)
    return df

def detect_airpuff_entry(df, ports, radius=30):
    distL = compute_distance(df, ports["airpuff_left"])
    distR = compute_distance(df, ports["airpuff_right"])

    df["Nose_in_any_airport"]   = ((distL < radius) | (distR < radius)).astype(int)

    return df