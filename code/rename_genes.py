import pandas as pd
import numpy as np
import argparse

"""
This script renames each gene to a standard name in the format chr_position. 
"""


def rename_genes(input_, output):
    print(input_)
    df = pd.read_csv(input_, sep="\t", header=None)
    print(df)
    df.columns = ["chr", "name", "0", "pos", "col1", "col2"]
    df["new_column"] = np.nan
    list_ = []

    for chromosome, position in zip(df["chr"], df["pos"]):
        """adding the chromosome number"""
        chromosome = str(chromosome)
        position = str(position)
        if len(str(chromosome)) == 1: name =  f"0{chromosome}_"
        elif len(str(chromosome)) == 2: name = f"{chromosome}_"
        
        """adding the position"""
        if len(str(position)) == 6: name = name + f"000{position}"
        elif len(str(position)) == 5: name = name + f"0000{position}"
        elif len(str(position)) == 7: name = name + f"00{position}"
        elif len(str(position)) == 8: name = name + f"0{position}"
        elif len(str(position)) == 9: name = name + f"{position}"
        list_.append(name)
    df["new_column"] = list_
    print("done!")
    df.to_csv(output, header=None)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Standardise gene nomenclature from .map file to .csv file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--input_", required=True, help=".map file")
    parser.add_argument("--output", required=True, help=".csv file")
    args = parser.parse_args()
    rename_genes(args.input_, args.output)