#!/usr/bin/env python3
#to run go into downloads, bridge, src
#python all.py -d /Users/kirstenpeterson/Downloads/BRIDGE_Connection/data/V1 /Users/kirstenpeterson/Downloads/BRIDGE_Connection/data/V2 /Users/kirstenpeterson/Downloads/BRIDGE_Connection/data/V3 /Users/kirstenpeterson/Downloads/BRIDGE_Connection/data/V4 /Users/kirstenpeterson/Downloads/BRIDGE_Connection/data/Friction_Full

import argparse
import glob
import logging
import os
import sys
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
from scipy.signal import savgol_filter

logging.getLogger('matplotlib.font_manager').setLevel(logging.WARNING)

def load_file(filepath):
    """Loads a CSV file and ensures numeric conversion for relevant columns."""
    df = pd.read_csv(filepath, skiprows=1)  # Skip header row
    df.drop(0, axis=0, inplace=True)  # Drop redundant index if necessary
    df = df.apply(pd.to_numeric, errors='coerce')  # Convert all data to numeric
    df.dropna(subset=["Time", "Stroke", "Force"], inplace=True)  # Remove invalid rows
    return df

def smooth_data(series, window=5, poly=2):
    """Applies Savitzky-Golay smoothing to a given series."""
    return savgol_filter(series, window_length=window, polyorder=poly, mode="nearest")

def truncate_after_force_zero(avg_data):
    """Truncates data after force returns to zero following an increase."""
    force_values = avg_data["Average Force"].values
    for i in range(1, len(force_values)):
        if force_values[i] <= 0 and force_values[i-1] > 0:
            return avg_data.iloc[:i]  # Truncate data up to this point
    return avg_data  # Return full data if force never returns to zero

def compute_average_force_and_stroke(data_paths):
    """Computes the average force and stroke across multiple CSV files."""
    all_data = []

    for filepath in data_paths:
        df = load_file(filepath)
        all_data.append(df[["Time", "Stroke", "Force"]])  # Ensure correct columns

    # Merge all dataframes on Time, keeping all time positions
    combined_df = pd.concat(all_data)

    # Compute mean per row
    avg_data = combined_df.groupby("Time").agg({
        "Stroke": "mean",
        "Force": "mean"
    }).reset_index()
    
    avg_data.rename(columns={"Stroke": "Average Distance", "Force": "Average Force"}, inplace=True)
    
    # Apply smoothing
    avg_data["Smoothed Distance"] = smooth_data(avg_data["Average Distance"])
    avg_data["Smoothed Force"] = smooth_data(avg_data["Average Force"])
    
    # Truncate data after force returns to zero
    avg_data = truncate_after_force_zero(avg_data)

    return avg_data

def plot_combined_force_vs_distance(avg_data_list, labels):
    """Plots smoothed average force against smoothed average distance for multiple datasets, keeping max force markers."""
    plt.figure(figsize=(16, 9))
    plt.subplots_adjust(right=0.65)  # Adjust subplot to make space for the legend
    plt.rcParams.update({'font.size': 25})  # Set all text to font size 25

    colors = ['b', 'g', 'r', 'c', 'm', 'y']
    marks = ["*", "d", "o", "X", "^", "s"]
    line_styles = ['-', '--', '-.', ':', (0, (3, 1, 1, 1)), (0, (5, 2, 2, 2))]

    for i, avg_data in enumerate(avg_data_list):
        color = colors[i % len(colors)]  # Cycle through colors if more than available
        mark = marks[i % len(marks)]  # Cycle through markers
        line_style = line_styles[i % len(line_styles)]  # Cycle through line styles
        plt.plot(avg_data["Smoothed Distance"], avg_data["Smoothed Force"],
             label=f'{labels[i]} Retention Force (N)',
             linewidth=4, color=color, linestyle=line_style, zorder=2)

        # Compute average max force and its corresponding distance
        max_force_idx = avg_data["Smoothed Force"].idxmax()
        max_force_value = avg_data.loc[max_force_idx, "Smoothed Force"]
        max_force_distance = avg_data.loc[max_force_idx, "Smoothed Distance"]

        closest_idx = (avg_data["Smoothed Distance"] - 30).abs().idxmin()
        force_at_closest = avg_data.loc[closest_idx, "Smoothed Force"]

        print(f"Force at closest distance ({avg_data.loc[closest_idx, 'Smoothed Distance']}): {force_at_closest}")

        # Add star marker at the max force location
        plt.scatter(max_force_distance, max_force_value,
                    color=color, marker=mark, s=350,
                    label=f"{labels[i]} Max Force {max_force_value:.2f} N", linewidth=1, zorder=3)

    # Create legend with opaque grey background
    #legend = plt.legend(loc='lower right', frameon=True, fontsize=20,ncol=2)  # Increase legend font size
    legend = plt.legend(loc='upper left', bbox_to_anchor=(1.05, 1), frameon=True, fontsize=15, ncol=1)
    frame = legend.get_frame()
    frame.set_facecolor('lightgrey')
    frame.set_alpha(0.8)

    plt.xlabel("Stroke (mm)", fontsize=35)   # Increase x-axis label font size
    plt.ylabel("Force (N)", fontsize=35)     # Increase y-axis label font size
    #plt.xlim(0,20)
    #plt.yscale('log')
    plt.xlim(10**-1.5, 20)
    plt.xscale('log')
    plt.xticks([0.1, 0.3, 1, 3, 10, 20], ['0.1', '0.3', '1', '3', '10', '20'])
    #plt.title("Average Retention Force vs Stroke Results", fontsize=35, fontweight='bold')
    #plt.title("Average Retention Force vs Stroke Results", fontsize=35, fontweight='bold', loc='center')
    plt.suptitle("Average Retention Force vs Stroke Results", fontsize=35, fontweight='bold')

    plt.grid()

    output_folder = Path("out")
    output_folder.mkdir(exist_ok=True)
    output_path = output_folder / "combined_force_vs_distance.png"
    plt.savefig(output_path, dpi=600)
    print(f"Plot saved to {output_path}")
    plt.show()

def main():
    parser = argparse.ArgumentParser(description="Compute and plot smoothed average force vs average distance for multiple datasets.")
    parser.add_argument("-d", "--directories", required=True, nargs="+",
                        help="Directories containing CSV files.")
    args = parser.parse_args()

    avg_data_list = []
    labels = []

    for directory in args.directories:
        input_paths = glob.glob(os.path.join(directory, "*.csv"))
        if input_paths:
            avg_data = compute_average_force_and_stroke(input_paths)
            avg_data_list.append(avg_data)
            labels.append(os.path.basename(directory))
        else:
            print(f"No CSV files found in directory {directory}")

    if avg_data_list:
        plot_combined_force_vs_distance(avg_data_list, labels)
    else:
        print("No data to plot.")

if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    sys.exit(main())
