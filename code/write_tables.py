import pandas as pd
import argparse
import numpy as np

def writing_to_table(var_file, homozygosity_file, output_file_common, output_file_all):
    """locations which are common but not homozygous in the father"""
    lst = []
    lst_all = []
    common = pd.read_csv(var_file, sep=";")
    homozygosity = pd.read_csv(homozygosity_file, sep=";")
    non_overlapping = []
    all_ = []
    
    for _, row1 in common.iterrows():
        chr1, start1, end1 = row1['CHR'], row1['BP1'], row1['BP2']
        out = []
        overlaps = homozygosity[(homozygosity['CHR'] == chr1) & (homozygosity['BP1'] < end1) & (homozygosity['BP2'] > start1)]
        all_.append((chr1, start1, end1))

        if overlaps.empty:
            non_overlapping.append((chr1, start1, end1))
        else:
            for _, row2 in overlaps.iterrows():
                start2, end2 = row2['BP1'], row2['BP2']
                if start1 < start2:
                    non_overlapping.append((chr1, start1, start2)) 
                if end1 > end2:
                    non_overlapping.append((chr1, end2, end1))
    for i in non_overlapping:
        lst.append([int(i[0]),  i[1], i[2],  i[2]-i[1]])
    for i in all_:
        lst_all.append([int(i[0]),  i[1],  i[2],  i[2]-i[1]])
    with open(output_file_common, "w") as f:
        df = pd.DataFrame(lst, columns=["Chromosome", "Start", "End", "Length"])
        df.to_csv(output_file_common, index=False, sep=" ")
    
    with open(output_file_all, "w") as file_:
        df_1 = pd.DataFrame(lst_all, columns=["Chromosome", "Start", "End", "Length"])
        df_1.to_csv(output_file_all, index=False, sep=" ")

    #non_overlapping = pd.DataFrame(non_overlapping, columns=['CHR', 'BP1', 'BP2'])
    #non_overlapping.to_csv(output_file, sep=';', index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Find common regions which are not homozygous in the father. Output two tables - one with all common regions, one with only non-homozygous.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--linked", required=True, help=".csv file")
    parser.add_argument("--homozygosity", required=True, help=".csv file")
    parser.add_argument("--common", required=True, help=".csv file")
    parser.add_argument("--with_homozyg", required=True, help=".csv file")
    args = parser.parse_args()

    writing_to_table(args.linked, args.homozygosity, args.common, args.with_homozyg)
