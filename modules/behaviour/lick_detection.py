import pandas as pd
import matplotlib.pyplot as plt
import os
import plotly.graph_objects as go
from ipywidgets import interact, FloatSlider

def txt_to_df(capacitance_txt_path):

    licks_df = pd.read_csv(capacitance_txt_path, header=None, names=["time(ms)", "capacitance", "recording"])

    return licks_df

def plot_licks_and_define_threshold(licks_df, threshold=None):
    """
    Plots capacitance vs time with Plotly.
    Allows defining a threshold for lick detection.
    Exports binary series (lick=1, no lick=0) as CSV.
    """

    plt.plot(licks_df["time(ms)"], licks_df["capacitance"], linestyle="-")
    plt.xlabel("Time (ms)")
    plt.ylabel("Capacitance")
    plt.grid(True)

    # Add vertical lines where timestamp is not empty
    for _, row in licks_df.dropna().iterrows():
        if row["recording"] == "RECORDING_START":
            plt.axvline(x=row["time(ms)"], color="green", linestyle="--", label="Start")
        elif row["recording"] == "RECORDING_STOP":
            plt.axvline(x=row["time(ms)"], color="red", linestyle="--", label="Stop")

    # Avoid duplicate labels in legend
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys())

    plt.show()

    return