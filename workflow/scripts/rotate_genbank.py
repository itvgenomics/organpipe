import argparse
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from naming import warn

# Strand of the start gene after rotation, so all assemblies of a sample have the same
# orientation. Plastomes usually have psbA on the reverse strand, next to the start.
START_GENE_STRAND = {"mito": 1, "chloro": -1}


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
    parser.add_argument("--input", help="GenBank file to rotate", required=True)
    parser.add_argument("--output", help="Rotated GenBank file", required=True)

    return parser.parse_args()


def find_start_feature(record, start_gene, organelle):
    if organelle == "mito":
        for feature in record.features:
            if feature.type == "rRNA" and "product" in feature.qualifiers:
                if start_gene in feature.qualifiers["product"]:
                    return feature

        raise ValueError(f"Product '{start_gene}' not found in the GenBank file.")

    for feature in record.features:
        if feature.type == "gene" and "gene" in feature.qualifiers:
            if start_gene in feature.qualifiers["gene"]:
                return feature

    raise ValueError(f"Gene '{start_gene}' not found in the GenBank file.")


def orient_record(record, start_gene, organelle):
    """Reverse complement the record if the start gene is on the other strand."""
    strand = find_start_feature(record, start_gene, organelle).location.strand
    if strand is None or strand == START_GENE_STRAND[organelle]:
        return record

    return record.reverse_complement(
        id=True, name=True, description=True, features=True, annotations=True, dbxrefs=True
    )


def rotate_genbank(input_file, output_file, start_gene, organelle):
    # Parse the GenBank file
    record = orient_record(SeqIO.read(input_file, "genbank"), start_gene, organelle)

    # Locate the start position using the start_gene
    start_pos = find_start_feature(record, start_gene, organelle).location.start

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

    try:
        rotate_genbank(
            input_file=args.input,
            output_file=args.output,
            start_gene=args.start_gene,
            organelle=args.organelle,
        )
    except Exception as e:
        warn(f"Could not rotate {args.input}: {e}")
