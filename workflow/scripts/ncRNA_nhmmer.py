from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
import argparse
import os
import logging
import sys
from naming import warn

FORMAT = "%(asctime)s [%(levelname)s] %(message)s"

logging.basicConfig(
    level=logging.INFO,
    stream=sys.stdout,
    format=FORMAT,
    datefmt="%Y-%m-%d %H:%M:%S",
)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Get the rRNA and tRNA sequences of an annotated assembly."
    )
    parser.add_argument(
        "--mitos_dir", help="MITOS2 output directory of the assembly", nargs="?"
    )
    parser.add_argument(
        "--genbank", help="GenBank file of the assembly (CPGAVAS2 or MitoHiFi)", nargs="?"
    )
    parser.add_argument(
        "--output_dir", help="Directory for rRNA-tRNA.fasta", required=True
    )

    args = parser.parse_args()

    return args


# Edit header for easy parse
def edit_header(input, output):

    logging.info(f"Editing result.fas Headers")
    with open(input, "r") as input:
        output = open(output, "w+")
        for line in input:
            if line.startswith(">"):
                header = line.split("; ")[-1]
                output.write(f">{header}")
            else:
                output.write(line)
    output.close()


# Function to get the rRNA and tRNA fasta sequences to run it with nhmmer
def select_sequences(input_file, output_file, search_list):
    selected_sequences = []
    with open(output_file, "w") as outfile:
        for record in SeqIO.parse(input_file, "fasta"):
            header = str(record.description)
            if any(search_string in header for search_string in search_list):
                SeqIO.write(record, outfile, "fasta")
                selected_sequences.append(record)

    return selected_sequences


def extract_ncRNA_to_fasta(genbank_file, output_file):
    sequences = []
    seen_names = {}

    # Parse the GenBank file
    with open(genbank_file, "r") as file:
        for record in SeqIO.parse(file, "genbank"):
            for feature in record.features:
                if feature.type in ["rRNA", "tRNA"]:  # Check for rRNA and tRNA
                    try:
                        if feature.location is None or feature.location.start < 0:
                            print(f"Fixing feature with invalid location: {feature}")
                            start = (
                                max(0, feature.location.start)
                                if feature.location
                                else 0
                            )
                            end = (
                                max(0, feature.location.end) if feature.location else 0
                            )
                            feature.location = FeatureLocation(start, end)

                        # Extract the name or identifier if available
                        name = feature.qualifiers.get("product", ["Unknown"])[0]

                        # Extract the sequence
                        sequence = feature.extract(record.seq)

                        # Genes with the same product (e.g. tRNA-Leu) keep their own entries
                        seen_names[name] = seen_names.get(name, 0) + 1
                        if seen_names[name] > 1:
                            name = f"{name}_{seen_names[name]}"

                        sequences.append((name, str(sequence)))
                    except KeyError:
                        pass  # Skip features with missing data

    # Write sequences to a single FASTA file
    with open(output_file, "w") as output:
        for name, seq in sequences:
            output.write(f">{name}\n{seq}\n")


if __name__ == "__main__":
    args = parse_arguments()

    output_file = f"{args.output_dir}/rRNA-tRNA.fasta"
    os.makedirs(args.output_dir, exist_ok=True)

    if args.mitos_dir and os.path.exists(f"{args.mitos_dir}/result.fas"):
        edit_header(
            f"{args.mitos_dir}/result.fas",
            f"{args.output_dir}/mitos2_annotation_edited.fas",
        )

        search_list = ["trn", "rrn"]
        select_sequences(
            f"{args.output_dir}/mitos2_annotation_edited.fas",
            output_file,
            search_list,
        )
    elif args.genbank and os.path.exists(args.genbank):
        extract_ncRNA_to_fasta(args.genbank, output_file)
    else:
        warn(f"No annotation found for {args.output_dir}, so there are no ncRNA sequences to search.")
        open(output_file, "w").close()
