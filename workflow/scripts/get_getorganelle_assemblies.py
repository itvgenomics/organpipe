import argparse
from pathlib import Path
from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser(
        description="Merge GetOrganelle FASTA files, adding database and status to sequence headers."
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input directory containing results/{sample}/getorganelle/"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output merged FASTA file"
    )
    args = parser.parse_args()

    input_dir = Path(args.input)

    fasta_files = sorted(input_dir.glob("*.fasta"))

    counter = 1
    records = []

    for fasta in fasta_files:
        parts = fasta.stem.split(".")

        if len(parts) < 3:
            print(f"Skipping unexpected filename: {fasta.name}")
            continue

        database = parts[0]
        status = parts[2]

        for record in SeqIO.parse(fasta, "fasta"):
            record.id = f"{database}_{status}_{counter}"
            record.name = record.id
            record.description = ""
            records.append(record)
            counter += 1

    output_path = Path(args.output)

    # Create output directory if necessary
    output_path.parent.mkdir(parents=True, exist_ok=True)

    SeqIO.write(records, output_path, "fasta")

    print(f"Wrote {len(records)} sequences to {output_path}")


if __name__ == "__main__":
    main()
