
import pandas as pd
import argparse

chromosomes = ["01","02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29"]

def map_file(filepath):
    """
    Outputs a dictionary mapping marker names to their respective positions.
    """
    loc_dict = {}
    with open(filepath) as location:
        for loc in location:
            split_loc = loc.split(",")
            key = split_loc[-1].strip()
            entry = split_loc[4]
            loc_dict[key] = entry
    return loc_dict


def specific_chr_map(dict_, chr):
    """Return map of specific chromosome"""
    output_dict = dict()
    n = 0
    for marker_name, position in dict_.items():
        if int(marker_name[:2]) == chr:
            output_dict[n] = position
            n+=1
    return output_dict


def fewest_or_most_of_char(lst, min_or_max, char):
    """
    Finds the index of the list entry with the fewest or most occurrences of a specific character.
    """
    counts = {i: entry.count(char) for i, entry in enumerate(lst)}
    return min(counts, key=counts.get) if min_or_max == "min" else max(counts, key=counts.get)

def character_match(char1, char2):
    """
    Checks if two characters match, considering 'N' as a wildcard.
    """
    ignored_mismatches = {"N": {"A", "B"}, "A": {"N"}, "B": {"N"}}
    return char1 == char2 or (char1 in ignored_mismatches and char2 in ignored_mismatches[char1])

def find_common_subsequences(input_file, min_common_length=5, fuse_adjacent=True):
    """
    Identifies common subsequences across all haplotypes in the input file.
    Optionally fuses adjacent subsequences separated by one or two positions.
    """
    haplotypes = pd.read_csv(input_file, header=None)[1].tolist()
    chr_length = len(haplotypes[0])
    haplotype_least_N = fewest_or_most_of_char(haplotypes, "min", "N")
    common_subsequences = {}
    start = 0
    prev_start, prev_end = None, None

    while start < chr_length:
        max_subseq = ""
        end = start + min_common_length

        while end <= chr_length:
            subsequence = haplotypes[haplotype_least_N][start:end]
            if all(character_match(subsequence[i], haplo[start + i]) for haplo in haplotypes for i in range(len(subsequence))):
                if len(subsequence) > len(max_subseq):
                    max_subseq = subsequence
                end += 1
            else:
                break

        if max_subseq:
            if fuse_adjacent and prev_start is not None and (start - prev_end) in [1]:
                if prev_start in common_subsequences:
                    prev_seq = common_subsequences.pop(prev_start)
                    fused_seq = prev_seq + 'N' * (start - prev_end) + max_subseq
                    common_subsequences[prev_start] = fused_seq
            else:
                common_subsequences[start] = max_subseq
            prev_start, prev_end = start, start + len(max_subseq)
            start += len(max_subseq)
        else:
            start += 1
    return common_subsequences

def find_where_different(comparison_file, common_sequences, max_identical_ratio=0.97):
    """
    Filters common sequences to exclude those that are too identical in the comparison file.
    """
    haplotypes = pd.read_csv(comparison_file, header=None)[1].tolist()
    filtered_sequences = {}

    for loc, sequence in common_sequences.items():
        short_haplotypes = [haplo[loc:loc + len(sequence)] for haplo in haplotypes]
        most_I = fewest_or_most_of_char(short_haplotypes, "max", "I")
        entry_similarity = short_haplotypes[most_I]
        nr_of_unique = entry_similarity.count("U")

        if len(entry_similarity) > 1 and nr_of_unique>2:
            identical_count = entry_similarity.count("I")
            if identical_count != 0:
                identical_ratio = identical_count/ len(entry_similarity)
            elif identical_count == 0:
                identical_ratio = 0
            if identical_ratio <= max_identical_ratio:
                filtered_sequences[loc] = sequence
    return filtered_sequences

### Output Function

def write_output_top_seq(filtered_sequences, output_filename):
    """
    Writes the filtered sequences to an output file, sorted by sequence length.
    """
    with open(output_filename, 'w') as f:
        f.write("Longest common paternal sequences in the half-sibs\n")
        for loc, seq in sorted(filtered_sequences.items(), key=lambda x: len(x[1]), reverse=True):
            f.write(f"{seq} at position {loc}\n")


def write_common_locations(filtered_sequences, output_filename, locations, chromosome, min_len = 1277230):
    
    """Writes the filtered sequences and their genomic locations to an output file."""
    
    with open(output_filename, "a") as location_file:
        #output_file.write("Longest common paternal sequences in the half-sibs\n")
        for loc, seq in sorted(filtered_sequences.items(), key=lambda x: len(x[1]), reverse=True):
            position_in_chr = int(locations.get(loc, -1))
            end_pos = loc + len(seq) - 1
            if end_pos in locations:
                length = int(locations[end_pos]) - position_in_chr
                loc_len = int(locations[end_pos]) - int(locations[loc])
                n_fraction = seq.count("N") / len(seq)
                if n_fraction < 0.3 and loc_len>= int(min_len):
                    
                    chr_label = chromosome.lstrip('0')
                    location_file.write(f"{chr_label};{position_in_chr/1000000};{int(locations[end_pos])/1000000}\n")
                    #output_file.write(f"{length}:{int(locations[end_pos]) - position_in_chr}: {seq} at position {position_in_chr} \n")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Find paternal haplotypes from input data",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--map", required=True, help="input .ped file")
    parser.add_argument("--locations", required=True, help="input .ped file")
    parser.add_argument("--top", required=True, help=".csv file")
    parser.add_argument("--comparison", required=False, help=".csv file")
    parser.add_argument("--length", required=True, help="minimum length of the common subsequence in bp")
    args = parser.parse_args()

    input_ = map_file(args.map)
    with open(args.locations, "a") as f:
        f.write("CHR;BP1;BP2\n")

    for chr in chromosomes:
        print(f"Processing chromosome {chr}...")
        all_common_subsequences = find_common_subsequences(f"results/example1/{chr}_output.csv", fuse_adjacent=True)
        #filtered_common = find_where_different(args.comparison, all_common_subsequences, max_identical_ratio=0.7)
        map_ = specific_chr_map(input_, int(chr))
        write_common_locations(all_common_subsequences, args.top, map_, chr, args.length)

# example 3
"""
input_ = map_file("example3/new.csv")
with open("example3/locations.csv", "a") as f:
    f.write("CHR;BP1;BP2\n")
for chr in chromosomes:
    print(f"Processing chromosome {chr}...")
    all_common_subsequences = find_common_subsequences(f"example3/{chr}_output.csv", fuse_adjacent=True)
    filtered_common = find_where_different(f"example3/{chr}_comparison.csv", all_common_subsequences, max_identical_ratio=0.7)
    map_ = specific_chr_map(input_, int(chr))
    write_common_locations(filtered_common, f"example3/output/top_sequences_{chr}.txt", map_, chr)

"""
