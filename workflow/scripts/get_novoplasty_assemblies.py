import argparse
import os
import re
from Bio import SeqIO
from Bio.Seq import Seq
from naming import warn


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Run and parse mitochondrial genome assembly data"
    )
    parser.add_argument(
        "--organelle",
        help="Organelle type: mito or chloro",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--genome_range",
        help="Genome Range",
        default="12000-22000",
        nargs="?",
    )
    parser.add_argument(
        "--sample",
        help="Sample dir",
        required=True,
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
    args = parser.parse_args()

    return args


def clean_sequence(sequence, name):
    """Replace every base other than ACGT with N, so that the coordinates do not change."""
    sequence = str(sequence).upper()
    cleaned, replaced = re.subn(r"[^ACGT]", "N", sequence)
    if replaced:
        warn(f"{name}: {replaced} base(s) other than A, C, G or T were replaced with N.")
    return Seq(cleaned)


if __name__ == "__main__":
    args = parse_arguments()

    sample = args.sample
    seed = args.seed
    kmer = args.kmer
    genome_range_min = str(args.genome_range).split("-")[0]
    genome_range_max = str(args.genome_range).split("-")[1]

    asm_dir = f"results/{sample}/novoplasty/{seed}/kmer{kmer}/"
    fasta_file = f"results/{sample}/assemblies/novoplasty/{seed}_kmer{kmer}.fasta"

    os.makedirs(os.path.dirname(fasta_file), exist_ok=True)

    all_sequences = []
    files = os.listdir(asm_dir)

    # Check the files to get the right ones
    circularized = any(file.startswith("Circularized_") for file in files)
    options = any(file.startswith("Option_") for file in files)

    for file in files:
        full_file_path = os.path.join(asm_dir, file)

        if file.startswith("Circularized_"):
            with open(full_file_path, "r") as input_handle:
                for record in SeqIO.parse(input_handle, "fasta"):
                    header = str(f"{file.split('.fasta')[0]}_{seed}_{kmer}").replace(
                        "Circularized_assembly", "CA"
                    )
                    # Clean the sequence
                    record.seq = clean_sequence(record.seq, header)
                    record.id = header
                    record.description = header
                    all_sequences.append(record)  # Append the record to the list

        elif file.startswith("Option_") and not circularized:
            with open(full_file_path, "r") as input_handle:
                for record in SeqIO.parse(input_handle, "fasta"):
                    if (
                        int(genome_range_min)
                        <= len(record.seq)
                        <= int(genome_range_max)
                    ):
                        header = str(
                            f"{file.split('.fasta')[0]}_{seed}_{kmer}"
                        ).replace("Option_", "OPT_")
                        record.seq = clean_sequence(record.seq, header)
                        record.id = header
                        record.description = header
                        all_sequences.append(record)  # Append the record to the list

        elif file.startswith("Contigs_") and not circularized and not options:
            with open(full_file_path, "r") as input_handle:
                for record in SeqIO.parse(input_handle, "fasta"):
                    if len(record.seq) >= int(genome_range_min) and len(
                        record.seq
                    ) <= int(genome_range_max):
                        header = str(
                            f"{file.split('.fasta')[0]}_{seed}_{kmer}"
                        ).replace("Contigs_", "CTG_")
                        record.seq = clean_sequence(record.seq, header)
                        record.id = header
                        record.description = header
                        all_sequences.append(record)  # Append the record to the list

    if not all_sequences:
        warn(f"Seed {seed} with k-mer {kmer} produced no assembly within {args.genome_range} bp.")

    # Now write all sequences to the output file in one go
    with open(fasta_file, "w") as output_handle:
        SeqIO.write(all_sequences, output_handle, "fasta")
