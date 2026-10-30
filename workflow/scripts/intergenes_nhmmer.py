from Bio import SeqIO
import os
import argparse
import pandas as pd
import logging
from naming import warn


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Get the intergenic sequences of an annotated assembly."
    )
    parser.add_argument(
        "--mitos_dir", help="MITOS2 output directory of the assembly", nargs="?"
    )
    parser.add_argument(
        "--genbank", help="GenBank file of the assembly (CPGAVAS2 or MitoHiFi)", nargs="?"
    )
    parser.add_argument(
        "--output_dir", help="Directory for intergenes_filter.fasta", required=True
    )

    args = parser.parse_args()

    return args


# Function to parse Fasta file
def parse_fasta(fasta_file):
    sequences = {}
    current_seq = None

    with open(fasta_file, "r") as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                current_seq = line[1:]
                sequences[current_seq] = ""
            else:
                sequences[current_seq] += line

    return sequences


def extract_intergenes_sequences(
    genome_file, position_list, output_file, positions_dict
):
    logging.info(f"Extracting intergene regions")
    ranges = []
    positions = set()
    sequences = {}

    # Convert position list to numeric ranges
    for entry in position_list:
        start, end = map(int, entry.split("-"))
        ranges.append((start, end))

    # Parse the genome fasta file
    genome = parse_fasta(genome_file)

    for seq_id, sequence in genome.items():
        seq_length = len(sequence)

        # Mark all positions in ranges as occupied
        for start, end in ranges:
            positions.update(range(start + 10, end - 10))

        # Extract non-overlapping sequences
        current_seq = ""
        range_index = 1
        current_range_start = None  # Initialize the variable
        for pos in range(seq_length):
            if pos not in positions:
                current_seq += sequence[pos]
            else:
                if current_seq:
                    if seq_id not in sequences:
                        sequences[seq_id] = []
                    range_header = f"{current_range_start}-{pos}"
                    sequences[seq_id].append((range_header, current_seq))
                    current_seq = ""
                    range_index += 1
                current_range_start = pos + 1

        # Add the last sequence if it extends to the end of the sequence
        if current_seq:
            if seq_id not in sequences:
                sequences[seq_id] = []
            range_header = f"{current_range_start}-{seq_length}"
            sequences[seq_id].append((range_header, current_seq))

    for entry in positions_dict:
        # Split the position string into two parts
        start, end = map(int, entry["position"].split("-"))
        # Update the positions
        start += 10
        end -= 10
        if end <= 0:
            end += 20

        # Update the position in the dictionary
        entry["position"] = f"{start}-{end}"
    print("positions_dict_update", positions_dict)

    # Write the non-overlapping sequences to a new fasta file
    with open(output_file, "w") as file:
        for seq_id, seq_list in sequences.items():
            for range_header, seq in seq_list:
                if range_header.split("-")[0] == "None":
                    start_range = 0
                else:
                    start_range = int(range_header.split("-")[0])
                end_range = int(range_header.split("-")[1])
                names_in_range = []
                for item in positions_dict:
                    # Split the position string into start and end
                    start_pos, end_pos = map(int, item["position"].split("-"))

                    # Check if the range overlaps with the current position
                    if start_pos <= end_range and end_pos >= start_range:
                        names_in_range.append(item["name"])
                if len(names_in_range) > 1:
                    file.write(f">{names_in_range[0]}-{names_in_range[1]}\n")
                    file.write(f"{seq}\n")
                else:
                    file.write(f">before_after-{names_in_range[0]}\n")
                    file.write(f"{seq}\n")


# Function to extract the gene positions in the genome
def get_positions(input_file):
    df_bed = pd.read_csv(
        input_file,
        header=None,
        delimiter="\t",
    )
    positions = []
    positions_list = []  # Initialize as a list
    for idx, row in df_bed.iterrows():
        position = f"{row[1]}-{row[2]}"
        positions.append(position)

        position_dict = {}
        position_dict["position"] = position
        position_dict["name"] = row[3]
        positions_list.append(position_dict)

    return positions, positions_list


def filter_fasta_by_size(input_file, output_file, min_size, max_size=10000):
    with open(input_file, "r") as input_handle, open(output_file, "w") as output_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            sequence_size = len(record.seq)
            if sequence_size >= min_size and sequence_size <= max_size:
                SeqIO.write(record, output_handle, "fasta")


def extract_unannotated_regions_gb(genbank_file, padding=10):
    unannotated_regions = []

    # Parse the GenBank file
    for record in SeqIO.parse(genbank_file, "genbank"):
        genome_sequence = record.seq  # Full sequence of the genome
        annotated_features = []

        # Extract all annotated feature positions with their names
        for feature in record.features:
            try:
                if "location" in dir(feature):
                    start = int(feature.location.start)
                    end = int(feature.location.end)
                    name = None
                    if "gene" in feature.qualifiers:
                        name = feature.qualifiers["gene"][0]
                    elif "product" in feature.qualifiers:
                        name = feature.qualifiers["product"][0]
                    else:
                        name = "unknown"
                    annotated_features.append((start, end, name))
            except:
                pass

        # Sort the annotated features
        annotated_features = sorted(annotated_features)

        # Find gaps (unannotated regions) and neighboring features
        previous_end = 0
        previous_name = "start_of_sequence"
        for start, end, name in annotated_features:
            if name != "unknown":
                if previous_end < start:  # There's a gap
                    # Include padding and ensure indices stay within bounds
                    left = max(0, previous_end - padding)
                    right = min(len(genome_sequence), start + padding)
                    unannotated_seq = genome_sequence[left:right]
                    unannotated_regions.append(
                        (previous_name, name, str(unannotated_seq))
                    )
                previous_end = max(previous_end, end)
                previous_name = name

        # Handle sequence end, if there's unannotated sequence after the last feature
        if previous_end < len(genome_sequence):
            left = max(0, previous_end - padding)
            right = len(genome_sequence)
            unannotated_seq = genome_sequence[left:right]
            unannotated_regions.append(
                (previous_name, "end_of_sequence", str(unannotated_seq))
            )

    return unannotated_regions


if __name__ == "__main__":
    args = parse_arguments()

    intergenes_file = f"{args.output_dir}/intergenes.fasta"
    filtered_file = f"{args.output_dir}/intergenes_filter.fasta"
    os.makedirs(args.output_dir, exist_ok=True)

    if args.mitos_dir and os.path.exists(f"{args.mitos_dir}/result.bed"):
        positions, positions_list = get_positions(f"{args.mitos_dir}/result.bed")
        extract_intergenes_sequences(
            f"{args.mitos_dir}/mitos2.fasta",
            positions,
            intergenes_file,
            positions_list,
        )
        filter_fasta_by_size(intergenes_file, filtered_file, 30, 10000)

    elif args.genbank and os.path.exists(args.genbank):
        unannotated_regions = extract_unannotated_regions_gb(args.genbank)

        with open(intergenes_file, "w") as output:
            for left_gene, right_gene, seq in unannotated_regions:
                output.write(f">{left_gene}-{right_gene}\n{seq}\n")

        filter_fasta_by_size(intergenes_file, filtered_file, 30)

    else:
        warn(f"No annotation found for {args.output_dir}, so there are no intergenic sequences to search.")
        open(filtered_file, "w").close()
