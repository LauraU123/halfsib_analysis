# write homozygosity locations
import pandas as pd
import argparse


def find_homozygosity(inputfile, outputfile):
    df = pd.read_csv(inputfile, delim_whitespace=True)
    with open(outputfile, "w") as f:
        f.write("CHR;BP1;BP2\n")
        for chr, start, end in zip(df["CHR"], df["POS1"], df["POS2"]):
            f.write(f"{chr};{start/1000000};{end/1000000}\n") 

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Write homozygosity file in correct format",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--input_", required=True, help=".hom file")
    parser.add_argument("--output", required=True, help=".csv file")
    args = parser.parse_args()

    find_homozygosity(args.input_, args.output)