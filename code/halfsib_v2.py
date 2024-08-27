import pandas as pd
import numpy as np

chromosomes = ["01","02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29" ]

def map_file(filepath):
    """Output:list of  marker names from tab-delimited file."""
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


def process_haplotypes(chromosome, dictionary, trios, loc_list, outputfile, mode="haplotype"):
    """Processes either haplotypes or identical loci based on the mode."""
    haplo = {}
    print(f"Processing {mode} for chromosome {chromosome}...")

    for trio in trios:
        if set(trio).issubset(dictionary):
            trio_data = np.array([dictionary[t] for t in trio])
            df = pd.DataFrame(data=trio_data, index=["offspring", "father", "mother"], columns=loc_list).T
            sequence = []

            for loc, (offspring, father, mother) in df[df.index.str.startswith(chromosome)].iterrows():
                if mode == "haplotype":
                    sequence.append(determine_haplotype(offspring, father, mother))
                elif mode == "identical":
                    sequence.append(determine_identical(offspring, father, mother))
           
            haplo[f"{trio[0]}_1"] = ''.join(sequence)
   
    save_to_csv(haplo, outputfile)

def determine_haplotype(offspring, father, mother):
    """Determines the haplotype based on offspring, father, and mother genotypes."""
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

def determine_identical(offspring, father, mother):
    """Determines if the loci are identical in mother and father."""
    if offspring == 1:
        return 'U'
    elif offspring == 3:
        return 'I' if father == mother else 'N'
    else:
        return 'I'

def save_to_csv(data, outputfile):
    """Saves the data dictionary to a CSV file."""
    df_output = pd.DataFrame.from_dict(data, orient='index')
    df_output.to_csv(outputfile, header=False)
    print(f"Output written to {outputfile}.")

for chr in chromosomes:
    print(f"Processing chromosome {chr}...")
    dictionary, trios = markers_and_trios("example1/data/new_filtered.ped")
    loc_list = map_file("example1/data/new_filtered.csv")
    process_haplotypes(chr, dictionary, trios, loc_list, f"example1/output/{chr}_output.csv", mode="haplotype")
    process_haplotypes(chr, dictionary, trios, loc_list, f"example1/output/{chr}_comparison.csv", mode="identical")

#for chr in chromosomes:
#    print(chr)
#    dataframe_chr = write_csv(chr, "example1/new_3.ped", "example1/new_4.csv", f"example1/{chr}_output.csv")
#    dataframe_identical = write_identical(chr, "example1/new_3.ped", "example1/new_4.csv", f"example1/{chr}_comparison.csv")

# for example 3
"""
for chr in chromosomes:
    print(f"Processing chromosome {chr}...")
    dictionary, trios = markers_and_trios("example3/new.ped")
    loc_list = map_file("example3/new.csv")
    process_haplotypes(chr, dictionary, trios, loc_list, f"example3/{chr}_output.csv", mode="haplotype")
    process_haplotypes(chr, dictionary, trios, loc_list, f"example3/{chr}_comparison.csv", mode="identical")"""