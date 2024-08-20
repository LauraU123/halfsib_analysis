
import pandas as pd

chromosomes = ["01","02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29"]

#note to self

"""
1.  add ratio options of longest I to length of string.  
"""
def map_file(filepath):
    """
    Outputs a dictionary mapping marker names to their respective positions.
    """
    loc_dict = {}
    with open(filepath) as location:
        for loc in location:
            split_loc = loc.split(",")
            key = split_loc[5].strip()
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

def find_common_subsequences(input_file, min_common_length=10, fuse_adjacent=True):
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
        #print(loc-1, len(sequence))
        short_haplotypes = [haplo[loc:loc + len(sequence)] for haplo in haplotypes]
        most_I = fewest_or_most_of_char(short_haplotypes, "max", "I")
        entry_similarity = short_haplotypes[most_I]
        nr_of_unique = entry_similarity.count("U")
        
        #if entry_similarity.startswith("N"):
        #    entry_similarity.replace("N", "I", 1)
       
        #if entry_similarity.startswith("I"):
        #    entry_similarity = entry_similarity.lstrip("I")
        #    sequence = sequence[len(sequence) - len(entry_similarity):]
        #    loc += len(sequence) - len(entry_similarity)
       
        #if entry_similarity.endswith("I"):
        #    entry_similarity = entry_similarity.rstrip("I")
        #    sequence = sequence[:len(entry_similarity)]
        if len(entry_similarity) > 1 and nr_of_unique>2:
            identical_count = entry_similarity.count("I")
            if identical_count != 0:
                identical_ratio = identical_count/ len(entry_similarity)
            elif identical_count == 0:
                identical_ratio = 0
            #print(loc, sequence, entry_similarity, identical_ratio)
            if identical_ratio <= 0.96:
                print(loc, sequence, entry_similarity, identical_ratio)
                filtered_sequences[loc] = sequence
                #print(loc, sequence, entry_similarity, identical_ratio, uniq_ratio)

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

"""
def write_common_locations(dictionary_subsequences, output_filename, locations, chr):
    with open("example3/locations.csv", "a") as f:
        with open(output_filename, "w") as outputfile:
            outputfile.write("Longest common paternal sequences in the half-sibs\n")
            for location, common_string in sorted(dictionary_subsequences.items(), key = lambda x: len(x[1]), reverse=True):
                positioninchromosome = int(locations[location])
                end_pos = location+len(common_string)-1
                if end_pos in locations:
                    length = int(locations[end_pos]) - positioninchromosome
                    number_of_N = common_string.count("N")
                    n_fraction = number_of_N/len(common_string)
                    if n_fraction < 0.3: 
                        if int(chr[0]) == 0: 
                            f.write(f"{chr[1]};{positioninchromosome/1000000};{int(locations.get(end_pos))/1000000}\n")
                        else: 

                            f.write(f"{chr};{positioninchromosome/1000000};{int(locations.get(end_pos))/1000000}\n")
                        outputfile.write(f"{length}:{int(end_pos)-int(locations[location])}: {common_string} at position {positioninchromosome} \n")"""

def write_common_locations(filtered_sequences, output_filename, locations, chromosome, min_len = 1277230):
    
    """Writes the filtered sequences and their genomic locations to an output file."""
    
    with open("example3/locations.csv", "a") as location_file, open(output_filename, "w") as output_file:
        output_file.write("Longest common paternal sequences in the half-sibs\n")
        for loc, seq in sorted(filtered_sequences.items(), key=lambda x: len(x[1]), reverse=True):
            
            position_in_chr = int(locations.get(loc, -1))
            end_pos = loc + len(seq) - 1
            #print(loc, end_pos, locations[loc], locations[end_pos])
            if end_pos in locations:
                length = int(locations[end_pos]) - position_in_chr
                #print(loc, end_pos, locations[loc], locations[end_pos])
                loc_len = int(locations[end_pos]) - int(locations[loc])
                print(loc_len)
                n_fraction = seq.count("N") / len(seq)
                #print(n_fraction)
                if n_fraction < 0.3 and loc_len>= min_len:
                    
                    chr_label = chromosome.lstrip('0')
                    location_file.write(f"{chr_label};{position_in_chr/1000000};{int(locations[end_pos])/1000000}\n")
                    output_file.write(f"{length}:{int(locations[end_pos]) - position_in_chr}: {seq} at position {position_in_chr} \n")

### Running the Script 

# Example 1


input_ = map_file("example1/new_4.csv")
with open("example1/locations.csv", "a") as f:
    f.write("CHR;BP1;BP2\n")


for chr in chromosomes:
    print(f"Processing chromosome {chr}...")
    all_common_subsequences = find_common_subsequences(f"example1/{chr}_output.csv", fuse_adjacent=True)
    filtered_common = find_where_different(f"example1/{chr}_comparison.csv", all_common_subsequences, max_identical_ratio=0.7)
    map_ = specific_chr_map(input_, int(chr))
    write_common_locations(filtered_common, f"example1/output/top_sequences_{chr}.txt", map_, chr)

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
