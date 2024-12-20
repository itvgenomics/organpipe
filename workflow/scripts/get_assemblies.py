import argparse
import os
import subprocess
from Bio import SeqIO


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


def clean_sequence(sequence):
    """Remove non-ATCG characters from the sequence."""
    return ''.join(filter(lambda x: x in 'ATCG', str(sequence)))


if __name__ == "__main__":
    args = parse_arguments()

    sample = args.sample
    seed = args.seed
    kmer = args.kmer
    genome_range_min = str(args.genome_range).split("-")[0]
    genome_range_max = str(args.genome_range).split("-")[1]

    os.makedirs(f"results/{sample}/assemblies/", exist_ok=True)

    asm_dir = f"results/{sample}/novoplasty/{seed}/kmer{kmer}/"
    fasta_file = f"results/{sample}/assemblies/{seed}_kmer{kmer}.fasta"

    all_sequences = []
    circularized = False
    options = False

    # Check the files to get the right ones
    for file in os.listdir(asm_dir):
        if file.startswith("Circularized_"):
            circularized = True
            break
        elif file.startswith("Option_"):
            options = True
            break  


    for file in os.listdir(asm_dir):
        full_file_path = os.path.join(asm_dir, file)

        if file.startswith("Circularized_"):
            with open(full_file_path, "r") as input_handle:
                for record in SeqIO.parse(input_handle, "fasta"):
                    # Clean the sequence
                    cleaned_seq = clean_sequence(record.seq)
                    header = str(f"{file.split('.fasta')[0]}_{seed}_{kmer}")
                    record.id = header
                    record.description = header
                    record.seq = cleaned_seq  # Update the record with the cleaned sequence
                    all_sequences.append(record)  # Append the record to the list

        elif file.startswith("Option_") and not circularized:
            with open(full_file_path, "r") as input_handle:
                for record in SeqIO.parse(input_handle, "fasta"):
                    cleaned_seq = clean_sequence(record.seq)
                    if int(genome_range_min) <= len(cleaned_seq) <= int(genome_range_max):
                        header = str(f"{file.split('.fasta')[0]}_{seed}_{kmer}")
                        record.id = header
                        record.description = header
                        record.seq = cleaned_seq
                        all_sequences.append(record)  # Append the record to the list

        elif file.startswith("Contigs_") and not circularized and not options:
            with open(full_file_path, "r") as input_handle:
                for record in SeqIO.parse(input_handle, "fasta"):
                    cleaned_seq = clean_sequence(record.seq)
                    if len(cleaned_seq) >= 1000 and len(cleaned_seq) <= int(genome_range_max):
                        header = str(f"{file.split('.fasta')[0]}_{seed}_{kmer}")
                        record.id = header
                        record.description = header
                        record.seq = cleaned_seq
                        all_sequences.append(record)  # Append the record to the list

    # Now write all sequences to the output file in one go
    if len(all_sequences) > 0:
        with open(fasta_file, "w") as output_handle:
            SeqIO.write(all_sequences, output_handle, "fasta")
    else:
        with open(fasta_file, "w") as output_handle:
            output_handle.write(">INVALIDSEED")
        
