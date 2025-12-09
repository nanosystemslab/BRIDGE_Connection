#!/usr/bin/env python3

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
    df.dropna(subset=["Time", "Stroke"], inplace=True)  # Remove invalid rows
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

    # Compute mean and standard deviation per row
    avg_data = combined_df.groupby("Time").agg({
        "Stroke": "mean",
        "Force": "mean"
    }).reset_index()
    std_data = combined_df.groupby("Time").agg({
        "Stroke": "std",
        "Force": "std"
    }).reset_index()
    
    avg_data["Stroke Std"] = std_data["Stroke"]
    avg_data["Force Std"] = std_data["Force"]
    avg_data.rename(columns={"Stroke": "Average Distance", "Force": "Average Force"}, inplace=True)
    
    # Apply smoothing
    avg_data["Smoothed Distance"] = smooth_data(avg_data["Average Distance"])
    avg_data["Smoothed Force"] = smooth_data(avg_data["Average Force"])
    
    # Truncate data after force returns to zero
    avg_data = truncate_after_force_zero(avg_data)

    print("Computed Standard Deviations:")
    print(avg_data[["Time", "Force Std", "Stroke Std"]])

    return avg_data
def plot_average_force_vs_distance(avg_data):
    """Plots the smoothed average force against the smoothed average distance with std deviation."""
    plt.figure(figsize=(16, 9))
    plt.rcParams.update({'font.size': 25})  # Set all text to font size 25
    
    # Compute average max force and its corresponding distance
    max_force_idx = avg_data["Smoothed Force"].idxmax()
    max_force_value = avg_data.loc[max_force_idx, "Smoothed Force"]
    max_force_distance = avg_data.loc[max_force_idx, "Smoothed Distance"]

    # Plot smoothed force vs smoothed distance with increased line thickness
    plt.plot(avg_data["Smoothed Distance"], avg_data["Smoothed Force"], 
             color='b', label='Retention Force (N)', linewidth=4, zorder=2)  # Increased linewidth to 4

    # Apply std deviation shading (unsmoothed)
    plt.fill_between(avg_data["Smoothed Distance"], 
                     avg_data["Smoothed Force"] - avg_data["Force Std"], 
                     avg_data["Smoothed Force"] + avg_data["Force Std"], 
                     color='b', alpha=0.2, label='Standard Deviation', zorder=1)

    # Add orange star marker at the max force location with high zorder
    plt.scatter(max_force_distance, max_force_value, 
                color='orange', marker='*', s=350, 
                label=f"Max Force {max_force_value:.2f} N", linewidth=1, zorder=3)

    # Create legend with opaque grey background
    legend = plt.legend(loc='upper right', frameon=True, fontsize=25)  # Increase legend font size
    frame = legend.get_frame()
    frame.set_facecolor('lightgrey')
    frame.set_alpha(0.8)
    
    plt.xlabel("Stroke (mm)", fontsize=35)   # Increase x-axis label font size
    plt.ylabel("Force (N)", fontsize=35)     # Increase y-axis label font size
    plt.title("Average Retention Force versus Stroke Results", fontsize=35, fontweight='bold')
    plt.grid()

    output_folder = Path("out")
    output_folder.mkdir(exist_ok=True)
    output_path = output_folder / "smoothed_force_vs_distance.png"
    plt.savefig(output_path, dpi=600)
    print(f"Plot saved to {output_path}")
    plt.show()






def main():
    parser = argparse.ArgumentParser(description="Compute and plot smoothed average force vs average distance.")
    parser.add_argument("-i", "--input", required=True, nargs="+",
                        help="Input CSV file paths or directory containing CSV files.")
    args = parser.parse_args()
    
    input_paths = []
    for path in args.input:
        if os.path.isdir(path):
            input_paths.extend(glob.glob(os.path.join(path, "*.csv")))
        else:
            input_paths.append(path)
    
    avg_data = compute_average_force_and_stroke(input_paths)
    plot_average_force_vs_distance(avg_data)

if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    sys.exit(main())
