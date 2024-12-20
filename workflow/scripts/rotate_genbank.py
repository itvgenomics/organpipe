import argparse
import os
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Rotate a GenBank file to start at a specified gene."
    )
    parser.add_argument(
        "--start_gene", type=str, help="Name of the gene to set as the starting point"
    )
    parser.add_argument(
        "--organelle", type=str, help="Type of organelle (e.g., 'mito', 'chloro')"
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

    return parser.parse_args()


def rotate_genbank(input_file, output_file, start_gene, organelle):
    # Parse the GenBank file
    record = SeqIO.read(input_file, "genbank")

    # Locate the start position using the start_gene
    start_pos = None

    if organelle == "mito":
        for feature in record.features:
            if feature.type == "rRNA" and "product" in feature.qualifiers:
                if start_gene in feature.qualifiers["product"]:
                    start_pos = feature.location.start
                    break

        if start_pos is None:
            raise ValueError(f"Product '{start_gene}' not found in the GenBank file.")

    elif organelle == "chloro":
        for feature in record.features:
            if feature.type == "gene" and "gene" in feature.qualifiers:
                if start_gene in feature.qualifiers["gene"]:
                    start_pos = feature.location.start
                    break

        if start_pos is None:
            raise ValueError(f"Gene '{start_gene}' not found in the GenBank file.")

    # Rotate the sequence
    sequence = record.seq
    rotated_sequence = sequence[start_pos:] + sequence[:start_pos]
    record.seq = rotated_sequence

    # Update feature locations
    new_features = []
    seq_length = len(sequence)

    def adjust_position(pos):
        """Adjust a single position based on the new sequence start."""
        if start_pos <= pos:
            return (pos - start_pos) % seq_length
        else:
            return (pos - start_pos + seq_length) % seq_length

    def adjust_location(location):
        """Adjust a FeatureLocation or CompoundLocation."""
        if location.parts:
            # It's a CompoundLocation (e.g., join())
            new_parts = []
            for part in location.parts:
                new_start = adjust_position(part.start)
                new_end = adjust_position(part.end)
                if new_start <= new_end:
                    new_parts.append(
                        FeatureLocation(new_start, new_end, strand=part.strand)
                    )
                else:  # Handle crossing sequence boundaries
                    new_parts.append(
                        FeatureLocation(
                            new_start, new_end + seq_length, strand=part.strand
                        )
                    )
            return sum(
                new_parts[1:], new_parts[0]
            )  # Combine parts into a new CompoundLocation
        else:
            # It's a single FeatureLocation
            new_start = adjust_position(location.start)
            new_end = adjust_position(location.end)
            if new_start <= new_end:
                return FeatureLocation(new_start, new_end, strand=location.strand)
            else:  # Handle crossing sequence boundaries
                return FeatureLocation(
                    new_start, new_end + seq_length, strand=location.strand
                )

    for feature in record.features:
        new_location = adjust_location(feature.location)
        new_features.append(
            SeqFeature(
                location=new_location, type=feature.type, qualifiers=feature.qualifiers
            )
        )

    # Replace the features with updated ones
    record.features = new_features

    # Write the rotated sequence to a new GenBank file
    with open(output_file, "w") as out_handle:
        SeqIO.write(record, out_handle, "genbank")


if __name__ == "__main__":
    args = parse_arguments()
    start_gene = args.start_gene
    organelle = args.organelle
    sample = args.sample
    seed = args.seed
    kmer = args.kmer
    sofware = args.software

    gb_dir = f"results/{sample}/genbanks"

    if organelle == "mito":
        for file in os.listdir(gb_dir):
            if seed in str(file) and kmer in str(file) and ".rotated." not in str(file):
                output_file = str(file).replace(".gb", ".rotated.gb")
                try:
                    rotate_genbank(
                        input_file=os.path.join(gb_dir, file),
                        output_file=os.path.join(gb_dir, output_file),
                        start_gene=start_gene,
                        organelle=organelle,
                    )
                except:
                    pass

    elif organelle == "chloro":
        for file in os.listdir(gb_dir):
            if (
                seed in str(file)
                and kmer in str(file)
                and sofware in str(file)
                and ".rotated." not in str(file)
            ):
                output_file = str(file).replace(".gb", ".rotated.gb")
                try:
                    rotate_genbank(
                        input_file=os.path.join(gb_dir, file),
                        output_file=os.path.join(gb_dir, output_file),
                        start_gene=start_gene,
                        organelle=organelle,
                    )
                except:
                    pass
