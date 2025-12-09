#!/usr/bin/env python3

import argparse
import glob
import logging
import os
import sys
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.signal import savgol_filter

logging.getLogger('matplotlib.font_manager').setLevel(logging.WARNING)

def load_file(filepath):
    """Loads a CSV file and ensures numeric conversion for relevant columns."""
    df = pd.read_csv(filepath, skiprows=1)  # Skip header row
    df.drop(0, axis=0, inplace=True)        # Drop redundant index row
    df = df.apply(pd.to_numeric, errors='coerce')  # Convert all data to numeric
    df.dropna(subset=["Time", "Stroke", "Force"], inplace=True)
    return df

def smooth_data(series, window=5, poly=2):
    """Applies Savitzky-Golay smoothing to a given series."""
    return savgol_filter(series, window_length=window, polyorder=poly, mode="nearest")

def truncate_after_force_zero(data):
    """Truncates data after force returns to zero following a rise."""
    force_values = data["Smoothed Force"].values
    for i in range(1, len(force_values)):
        if force_values[i] <= 0 and force_values[i-1] > 0:
            return data.iloc[:i]
    return data

def process_single_file(filepath):
    """Processes a single CSV file and returns a smoothed, truncated DataFrame."""
    df = load_file(filepath)

    df["Smoothed Distance"] = smooth_data(df["Stroke"])
    df["Smoothed Force"] = smooth_data(df["Force"])

    df = df[["Smoothed Distance", "Smoothed Force"]]
    df = truncate_after_force_zero(df)
    
    return df

def plot_all_trials(data_list, labels):
    """Plots each trial's force vs distance curve with unique colors."""
    plt.figure(figsize=(16, 9))
    plt.subplots_adjust(right=0.65)
    plt.rcParams.update({'font.size': 25})

    # Extended color palette
    colors = [
        'b', 'g', 'r', 'c', 'm', 'y', 'orange', 'purple',
        'brown', 'pink', 'olive', 'teal', 'navy', 'maroon'
    ]

    max_x = 0  # track maximum x across all data

    for i, data in enumerate(data_list):
        color = colors[i % len(colors)]

        plt.plot(data["Smoothed Distance"], data["Smoothed Force"],
                 label=f'{labels[i]}',
                 linewidth=3, color=color, zorder=2)

        max_idx = data["Smoothed Force"].idxmax()
        max_force = data.loc[max_idx, "Smoothed Force"]
        max_dist = data.loc[max_idx, "Smoothed Distance"]

        plt.scatter(max_dist, max_force, color=color, marker="o", s=200,
                    label=f'{labels[i]} Max {max_force:.2f} N', zorder=3)

        # update global max x
        if data["Smoothed Distance"].max() > max_x:
            max_x = data["Smoothed Distance"].max()

    plt.xlabel("Stroke (mm)", fontsize=35)
    plt.ylabel("Force (N)", fontsize=35)

    # ✅ X axis: 0 → max, ticks every 5
    plt.xlim(0, max_x)
    plt.xticks(np.arange(0, max_x + 1, 5))

    # ✅ Y axis autoscale (no buffer)
    plt.autoscale(enable=True, axis="y", tight=True)

    plt.suptitle("Retention Force vs Stroke for Each Trial", fontsize=35, fontweight='bold')
    plt.grid()

    legend = plt.legend(loc='upper left', bbox_to_anchor=(1.05, 1), fontsize=15)
    legend.get_frame().set_facecolor('lightgrey')
    legend.get_frame().set_alpha(0.8)

    output_folder = Path("out")
    output_folder.mkdir(exist_ok=True)
    output_path = output_folder / "trials_force_vs_distance.png"
    plt.savefig(output_path, dpi=600)
    print(f"Plot saved to {output_path}")
    plt.show()

def main():
    parser = argparse.ArgumentParser(description="Plot force vs distance for each individual trial in given folders.")
    parser.add_argument("-d", "--directories", required=True, nargs="+",
                        help="Folders containing CSV files.")
    args = parser.parse_args()

    all_data = []
    all_labels = []

    for directory in args.directories:
        csv_files = glob.glob(os.path.join(directory, "*.csv"))
        if not csv_files:
            print(f"No CSV files found in {directory}")
            continue

        for file_path in sorted(csv_files):
            trial_data = process_single_file(file_path)
            all_data.append(trial_data)
            all_labels.append(os.path.splitext(os.path.basename(file_path))[0])  # file name without extension

    if all_data:
        plot_all_trials(all_data, all_labels)
    else:
        print("No valid data found to plot.")

if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    sys.exit(main())

