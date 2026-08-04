import argparse
import os
from pathlib import Path


def split_novoplasty_fasta(input_file, output_dir):
    # Ensure the output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    with open(input_file, "r") as fasta_file:
        header = ""
        sequence = ""
        index = 1  # Initialize an index counter

        for line in fasta_file:
            line = line.strip()
            if line.startswith(">"):  # This line is a header
                if header:  # If we have a previous header, write it to a file
                    write_sequence_to_file(header, sequence, output_dir, index)
                    index += 1  # Increment the index for the next sequence
                header = line[1:]  # Remove '>' from the header
                sequence = ""  # Reset sequence for the new header
            else:
                sequence += line  # Append the sequence line

        # Write the last sequence to a file
        if header:
            write_sequence_to_file(header, sequence, output_dir, index)
            index += 1  # Increment the index for the next sequence

def write_sequence_to_file(header, sequence, output_dir, index):
    # Create a safe filename by replacing invalid characters
    safe_header = header.replace("/", "_").replace("\\", "_").replace(" ", "_").replace("+", "_")
    filename = os.path.join(output_dir, "{}_{}.fasta".format(safe_header, index))

    with open(filename, "w") as output_file:
        output_file.write(">{}_{}\n".format(header, index))  # Write the header
        output_file.write("{}\n".format(sequence))  # Write the sequence

def split_getorganelle_fasta(input_fasta, output_dir):

    input_fasta = Path(input_fasta)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    out = None

    with open(input_fasta) as fin:
        for line in fin:
            if line.startswith(">"):
                if out:
                    out.close()

                header = line[1:].strip().split()[0]
                filename = output_dir / f"{header}.fasta"
                out = open(filename, "w")
                out.write(line)
            else:
                out.write(line)

    if out:
        out.close()


def main():
    parser = argparse.ArgumentParser(
        description="Split a FASTA file into multiple files, one for each sequence."
    )
    parser.add_argument(
        "--fasta_file", required=True, help="Path to the input FASTA file."
    )
    parser.add_argument(
        "--output_dir",
        default=".",
        help="Directory to save the output FASTA files. Default is the current directory.",
    )
    parser.add_argument(
        "--assembler", required=True, help="Assembler used for the FASTA file."
    )

    args = parser.parse_args()

    if args.assembler.lower() == "novoplasty":
        split_novoplasty_fasta(args.fasta_file, args.output_dir)
    elif args.assembler.lower() == "getorganelle":
        split_getorganelle_fasta(args.fasta_file, args.output_dir)

if __name__ == "__main__":
    main()
