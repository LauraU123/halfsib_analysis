import pandas as pd
import argparse
import numpy as np


"""def variant_only(var_file, homozygosity_file, output_file):
    This script returns regions which are variants but not homozygous. 
    variants = pd.read_csv(var_file, sep=";")
    homozygosity = pd.read_csv(homozygosity_file, sep=";")
    non_overlapping = []

    for _, row1 in variants.iterrows():
        chr1, start1, end1 = row1['CHR'], row1['BP1'], row1['BP2']

        # Find any overlapping rows in df2
        overlaps = homozygosity[(homozygosity['CHR'] == chr1) & (homozygosity['BP1'] < end1) & (homozygosity['BP2'] > start1)]

        if overlaps.empty:
            non_overlapping.append((chr1, start1, end1))
        else:
            for _, row2 in overlaps.iterrows():
                start2, end2 = row2['BP1'], row2['BP2']
                if start1 < start2:
                    non_overlapping.append((chr1, start1, start2))  # Before overlap
                if end1 > end2:
                    non_overlapping.append((chr1, end2, end1))  # After overlap
        for i in non_overlapping:
            print(i[0], i[1], i[2])

    non_overlapping = pd.DataFrame(non_overlapping, columns=['CHR', 'BP1', 'BP2'])
    non_overlapping.to_csv(output_file, sep=';', index=False)"""


def variant_only(var_file, homozygosity_file, output_file):
    """This script returns regions which are variants but not homozygous. """
    variants = pd.read_csv(var_file, sep=";")
    homozygosity = pd.read_csv(homozygosity_file, sep=";")
    non_overlapping = []

    for _, row1 in variants.iterrows():
        chr1, start1, end1 = row1['CHR'], row1['BP1'], row1['BP2']
        out = []
        # Find any overlapping rows in df2
        overlaps = homozygosity[(homozygosity['CHR'] == chr1) & (homozygosity['BP1'] < end1) & (homozygosity['BP2'] > start1)]

        if overlaps.empty:
            non_overlapping.append((chr1, start1, end1))
        else:
            for _, row2 in overlaps.iterrows():
                start2, end2 = row2['BP1'], row2['BP2']
                if start1 < start2:
                    non_overlapping.append((chr1, start1, start2))  # Before overlap
                if end1 > end2:
                    non_overlapping.append((chr1, end2, end1))  # After overlap
        with open(output_file, "w") as f:
            for i in non_overlapping:
                f.write(f"chromosome : {int(i[0])} location : {i[1]} - {i[2]} length : {i[2]-i[1]}\n")


    #non_overlapping = pd.DataFrame(non_overlapping, columns=['CHR', 'BP1', 'BP2'])
    #non_overlapping.to_csv(output_file, sep=';', index=False)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Find variants which are not homozygous in the father",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--variants", required=True, help=".csv file")
    parser.add_argument("--homozygosity", required=True, help=".csv file")
    parser.add_argument("--output_csv", required=True, help=".csv file")
    args = parser.parse_args()

    variant_only(args.variants, args.homozygosity, args.output_csv)
