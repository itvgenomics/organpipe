from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation
import warnings
import os
import argparse
from timeit import default_timer as timer
import pandas as pd
import sys
import subprocess
import logging
import glob


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Uses NHMMER to recheck the MITOS2 annotation"
    )
    parser.add_argument(
        "--organelle", help="Organelle type: mito or chloro", required=True
    )
    parser.add_argument(
        "--sequencing_type",
        help="Equipment used in the sequencing (e.g. illumina or pacbio)",
        required=True,
    )
    parser.add_argument(
        "--sample",
        help="Sample dir",
        required=False,
        nargs="?",
    )
    parser.add_argument(
        "--kmer",
        help="kmer used",
        nargs="?",
    )
    parser.add_argument(
        "--seed",
        help="Seed used",
        nargs="?",
    )
    parser.add_argument(
        "--software",
        help="Software used (e.g cpgavas2 or chloe)",
        nargs="?",
    )

    args = parser.parse_args()

    return args


def list_directories(path):
    # List to hold directories
    directories = []

    # Iterate through the items in the specified path
    for item in os.listdir(path):
        # Create the full path
        full_path = os.path.join(path, item)
        # Check if it is a directory
        if os.path.isdir(full_path):
            directories.append(item)

    return directories


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


# Function to create the HMMER output directory
def create_dirs(outpath):
    if not os.path.exists(f"{outpath}/nhmmer"):
        logging.info(f"Creating {outpath}/nhmmer")
        os.makedirs(f"{outpath}/nhmmer")
    else:
        logging.info(f"{outpath}/nhmmer exists")


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
    sequences = {}

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

                        # Extract the feature type (rRNA or tRNA)
                        feature_type = feature.type

                        # Extract the name or identifier if available
                        name = feature.qualifiers.get("product", ["Unknown"])[0]

                        # Extract the sequence
                        sequence = feature.extract(record.seq)

                        # Save to the dictionary with a combined identifier
                        sequences[name] = str(sequence)
                    except KeyError:
                        pass  # Skip features with missing data

    # Write sequences to a single FASTA file
    with open(output_file, "w") as output:
        for name, seq in sequences.items():
            output.write(f">{name}\n{seq}\n")


if __name__ == "__main__":
    args = parse_arguments()

    if args.organelle == "mito":
        create_dirs(f"results/{args.sample}")

        if args.sequencing_type == "short":
            annotations_path = (
                f"results/{args.sample}/mitos2/{args.seed}_kmer{args.kmer}"
            )
            annotations_dir = list_directories(annotations_path)

            for annotation in annotations_dir:
                output_path = f"results/{args.sample}/nhmmer/{args.seed}_kmer{args.kmer}/{annotation}"

                if not os.path.exists(output_path):
                    os.makedirs(output_path)

                if os.path.exists(f"{annotations_path}/{annotation}/result.fas"):
                    edit_header(
                        f"{annotations_path}/{annotation}/result.fas",
                        f"{output_path}/mitos2_annotation_edited.fas",
                    )

                    search_list = ["trn", "rrn"]
                    select_sequences(
                        f"{output_path}/mitos2_annotation_edited.fas",
                        f"{output_path}/rRNA-tRNA.fasta",
                        search_list,
                    )

        else:
            annotations_path = f"results/{args.sample}/mitohifi/"
            annotations_dir = list_directories(annotations_path)

            for annotation in annotations_dir:
                output_path = f"results/{args.sample}/nhmmer/{args.seed}"

                if not os.path.exists(output_path):
                    os.makedirs(output_path)

                if os.path.exists(
                    f"{annotations_path}/{annotation}/final_mitogenome.gb"
                ):
                    extract_ncRNA_to_fasta(
                        f"{annotations_path}/{annotation}/final_mitogenome.gb",
                        f"{output_path}/rRNA-tRNA.fasta",
                    )

    if args.organelle == "chloro":
        create_dirs(f"results/{args.sample}")

        if args.sequencing_type == "short":
            annotations_path = f"results/{args.sample}/genbanks"

            for annotation in os.listdir(annotations_path):
                output_path = f"results/{args.sample}/nhmmer/{args.seed}_kmer{args.kmer}/{annotation.replace(".cpgavas2.gb", "")}"

                if annotation.endswith(".cpgavas2.gb") and (
                    "Circularized_assembly_1" in annotation or "Option_1" in annotation
                ) and args.kmer in annotation and args.seed in annotation:
                    
                    if not os.path.exists(output_path):
                        os.makedirs(output_path)
                    
                    extract_ncRNA_to_fasta(
                        os.path.join(annotations_path, annotation),
                        f"{output_path}/rRNA-tRNA.fasta",
                    )