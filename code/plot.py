import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse


def plotting(chr_file, locations, homozygosity, output):

    # Read data from CSV files
    chromosomes = pd.read_csv(chr_file, sep=";")
    data = pd.read_csv(locations, sep=";")
    homozygosity = pd.read_csv(homozygosity, sep=";")

    scale = 1000000

    pdf_width = 11.69
    pdf_height = 8.27
    plt.figure(figsize=(pdf_width, pdf_height))

    line_width = 10
    fig, ax = plt.subplots()
    ax.bar(chromosomes['Chr'], chromosomes['length'] / scale, width=0.4, color='grey', label='Chromosomes')

    # Inverted y-axis
    ax.invert_yaxis()
    ax.xaxis.tick_top()

    # Print blue regions from the test locations
    for run in range(1, 30):
        subset = data[data['CHR'] == run]
        if len(subset) > 0:
            for _, row in subset.iterrows():
                start_position = row['BP1'] / scale
                end_position = row['BP2'] / scale
                ax.plot([run, run], [start_position, end_position], lw=line_width, color='blue')
        else:
            print(f"No variant on chromosome {run}")

    # Print homozygosity regions in red
    for run in range(1, 30):
        subset_h = homozygosity[homozygosity['CHR'] == run]
        if len(subset_h) > 0:
            for _, row in subset_h.iterrows():
                start_position_h = row['BP1'] / scale
                end_position_h = row['BP2'] / scale
                ax.plot([run, run], [start_position_h, end_position_h], lw=line_width, color='red', alpha=0.5)
        else:
            print(f"No homozygosity on chromosome {run}")

    ax.set_ylabel("Length in Mb")
    ax.set_xlabel("Chromosomes")
    ax.set_xticks(chromosomes['Chr'])
    ax.set_xticklabels(chromosomes['Chr'])
    ax.grid(False)

    plt.savefig(output, format="pdf")
    plt.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Plot output with chromosomes",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--chr", required=True, help=".csv file")
    parser.add_argument("--locations", required=True, help=".csv file")
    parser.add_argument("--homozygosity", required=True, help=".csv file")
    parser.add_argument("--plot", required=True, help=".csv file")
    args = parser.parse_args()

    plotting(args.chr, args.locations, args.homozygosity, args.plot)