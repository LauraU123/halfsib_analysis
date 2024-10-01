import pandas as pd
import numpy as np
import argparse


def map_file(filepath):
    """Outputs list of marker names from tab-delim. file."""
    with open(filepath) as f:
        return [line.split(",")[-1].strip() for line in f]

def markers_and_trios(filepath):
    recode = {
        "2 2": 2, "1 1": 0, "1 2": 1, "2 1": 1, "0 0": 3,
        "A A": 0, "B B": 2, "A B": 1, "B A": 1
    }
    dictionary = {}
    trios = []
    with open(filepath) as f:
        for line in f:
            parts = line.split("\t")
            trio = parts[1:4]
            trios.append(trio)
            dictionary[parts[1]] = [recode.get(marker.strip(), None) for marker in parts[6:]] 
    return dictionary, trios

def process_haplotypes(chromosome, dictionary, trios, loc_list, outputfile):
    """processing haplotype trios"""
    haplo = {}
    print(f"Processing haplotype for chromosome {chromosome}...")

    for trio in trios:
        if set(trio).issubset(dictionary):
            trio_data = np.array([dictionary[t] for t in trio])
            df = pd.DataFrame(data=trio_data, index=["offspring", "father", "mother"], columns=loc_list).T
            sequence = []

            for loc, (offspring, father, mother) in df[df.index.str.startswith(chromosome)].iterrows():
                sequence.append(determine_haplotype(offspring, father, mother))
           
            haplo[f"{trio[0]}_1"] = ''.join(sequence)
   
    save_to_csv(haplo, outputfile)

def determine_haplotype(offspring, father, mother):
    """determining haplotype based on offspring, father, mother genotypes"""
    if offspring == 2:
        return 'B'
    elif offspring == 0:
        return 'A'
    elif offspring == 1:
        if father == 2 or (father == 1 and mother == 0):
            return 'B'
        elif father == 0 or (father == 1 and mother == 2):
            return 'A'
        else:
            return 'N'
    else:
        return 'B' if father == 2 else 'A' if father == 0 else 'N'

def save_to_csv(data, outputfile):
    """Saving dictionary to csv"""
    df_output = pd.DataFrame.from_dict(data, orient='index')
    print(df_output)
    df_output.to_csv(outputfile, header=False)
    print(f"Output written to {outputfile}.")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Find paternal haplotypes from input data",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--ped", required=True, help="input .ped file")
    parser.add_argument("--markers", required=True, help=".csv file with gene map")
    parser.add_argument("--output", required=False, help=".csv file")
    parser.add_argument("--chr", required=True, help="number of chromosomes")
    args = parser.parse_args()
    
    print(f"Processing chromosome {args.chr}...")
    dictionary, trios = markers_and_trios(args.ped)
    loc_list = map_file(args.markers)
    process_haplotypes(args.chr, dictionary, trios, loc_list, args.output)
