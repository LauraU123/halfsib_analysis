import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
from matplotlib.ticker import FuncFormatter


def plotting(chr_file, chr_nr, locations, homozyg, output):

    chromosomes = pd.read_csv(chr_file, sep=";")
    data = pd.read_csv(locations, sep=";")
    homozygosity = pd.read_csv(homozyg, sep=";")
    scale = 1000000
    pdf_width = 17.69
    pdf_height = 8.27
    line_width = 5
    plt.figure(figsize=(pdf_width, pdf_height))

    fig, ax = plt.subplots()
    ax.bar(chromosomes['Chr'], chromosomes['length'] / scale, width=0.4, color='grey')

    ax.invert_yaxis()
    ax.xaxis.tick_top()
    # add common haplotype locations in blue
    for run in range(1, int(chr_nr)+1):
        subset = data[data['CHR'] == run]
        if len(subset) > 0:
            for _, row in subset.iterrows():
                start_position = row['BP1'] / scale
                end_position = row['BP2'] / scale
                ax.plot([run, run], [start_position, end_position], lw=line_width, color='blue')

        else:
            print(f"No variant on chromosome {run}")
    # add homozygosity in paternal haplotypes in red
    for run in range(1, 30):
        subset = homozygosity[homozygosity['CHR'] == run]
        if len(subset) > 0:
            for _, row in subset.iterrows():
                start_position = row['BP1'] / scale
                end_position = row['BP2'] / scale
                ax.plot([run, run], [start_position, end_position], lw=line_width, color='red', alpha=0.5)

        else:
            print(f"No homozygosity on chromosome {run}")

    def scale_format(x, pos):
        return f"{int(x*scale):,}"
    ax.plot([], [], lw=line_width, color='blue', label='Linked Haplotypes')
    ax.plot([], [], lw=line_width, color='red', label='Homozygosity')
    ax.plot([], [], lw=line_width, color='purple', label='Overlap')
    # y axis ticks and scaling 
    y_ticks = ax.get_yticks()
    extra_lines = []
    for i in range(len(y_ticks)-1):
        extra_lines.append(y_ticks[i] + y_ticks[i+1]/2)
        extra_lines.append(y_ticks[i] - y_ticks[i+1]/2)

    all_lines = np.sort(np.concatenate((y_ticks, extra_lines)))
    for y in all_lines:
            ax.axhline(y=y, color='black', linestyle='--', linewidth=0.5)
    ax.yaxis.set_major_formatter(FuncFormatter(scale_format))

    ax.legend(loc='lower right')
    ax.set_ylabel("Length in Mb")
    ax.set_xlabel("Chromosomes")
    ax.set_xticks(chromosomes['Chr'])
    ax.set_xticklabels(chromosomes['Chr'])
    ax.tick_params(axis='x', labelsize=6)
    max_length_mb = chromosomes['length'].max()/scale
    ax.set_ylim(0,max_length_mb)
    ax.invert_yaxis()
    ax.grid(False)

    plt.savefig(output, format="pdf")
    plt.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Plot output with chromosomes",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--chr_file", required=True, help="number of chromosomes")
    parser.add_argument("--linked", required=True, help=".csv file")
    parser.add_argument("--chr_nr", required=True, help=".csv file")
    parser.add_argument("--homozyg", required=True, help=".csv file with ")
    parser.add_argument("--plot", required=True, help="output plot in pdf format")
    args = parser.parse_args()

    plotting(args.chr_file, args.chr_nr, args.linked, args.homozyg, args.plot)