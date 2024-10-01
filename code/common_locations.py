
import pandas as pd
import argparse

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
    """ find  index of list with fewest or most occurrences of a specific character."""
    counts = {i: entry.count(char) for i, entry in enumerate(lst)}
    return min(counts, key=counts.get) if min_or_max == "min" else max(counts, key=counts.get)

def character_match(char1, char2):
    """
    Check if characters match. A or B can only match themselves or N. N can match both.
    """
    ignored_mismatches = {"N": {"A", "B"}, "A": {"N"}, "B": {"N"}}
    return char1 == char2 or (char1 in ignored_mismatches and char2 in ignored_mismatches[char1])

def find_common_subsequences(input_file, min_common_length=5, fuse_nr = 1, fuse_adjacent=True):
    """Identify common subsequences across all haplotypes in the input - optionally fuses adjacent subsequences separated by one or two positions.
    """
    haplotypes = pd.read_csv(input_file, header=None)[1].tolist()
    chr_length = len(haplotypes[0])
    haplotype_least_N = fewest_or_most_of_char(haplotypes, "min", "N")
    common_subsequences = {}
    start = 0
    prev_start, prev_end = None, None

    while start < chr_length:
        max_subseq = ""
        end = start + int(min_common_length)

        while end <= chr_length:
            subsequence = haplotypes[haplotype_least_N][start:end]
            if all(character_match(subsequence[i], haplo[start + i]) for haplo in haplotypes for i in range(len(subsequence))):
                if len(subsequence) > len(max_subseq):
                    max_subseq = subsequence
                end += 1
            else:
                break

        if max_subseq:
            if fuse_adjacent and prev_start is not None and (start - prev_end) in [fuse_nr]:
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


### Output

def write_common_locations(filtered_sequences, output_filename, locations, chromosome, n_fraction_max= 0.3, min_len = 1000000):
    
    """Writing filtered sequences and genomic locations to outputfile"""

    with open(output_filename, "a") as location_file:
        for loc, seq in sorted(filtered_sequences.items(), key=lambda x: len(x[1]), reverse=True):
            position_in_chr = int(locations.get(loc, -1))
            end_pos = loc + len(seq) - 1
            if end_pos in locations:
                length = int(locations[end_pos]) - position_in_chr
                loc_len = int(locations[end_pos]) - int(locations[loc])
                n_fraction = seq.count("N") / len(seq)
                if n_fraction < float(n_fraction_max) and loc_len>= int(min_len):
                    chr_label = chromosome.lstrip('0')
                    location_file.write(f"{chr_label};{position_in_chr/1000000};{int(locations[end_pos])/1000000}\n")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Find common locations from haplotypes paternal haplotypes",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--map", required=True, help="input map file")
    parser.add_argument("--output", required=True, help="output file")
    parser.add_argument("--hapl", required=True, help="input haplotype file")
    parser.add_argument("--min_markers", required=True, help="minimum length of markers to be included in the analysis")
    parser.add_argument("--length", required=True, help="minimum length of the common subsequence in bp")
    parser.add_argument("--fuse_adjacent", required=True, help="Should neighbouring variants with marker missing in between be fused")
    parser.add_argument("--fuse_adjacent_nr", required=False, help="if neighbouring variants with marker missing are accepted, how many gaps are accepted. Default 1")
    parser.add_argument("--n_fraction_max", required=True, help="maximum undetermined char fraction in the variants. e.g in NNB it is 2/3, so it is not included.")
    parser.add_argument("--chr", required=True, help="number of chromosomes")
    args = parser.parse_args()

    input_ = map_file(args.map)
    print(f"Processing chromosome {args.chr}...")
    all_common_subsequences = find_common_subsequences(args.hapl, args.min_markers, args.fuse_adjacent_nr, args.fuse_adjacent)
    map_ = specific_chr_map(input_, int(args.chr))
    write_common_locations(all_common_subsequences, args.output, map_, args.chr, args.n_fraction_max, args.length)
    print(f"Chromosome {args.chr} written to {args.output}")