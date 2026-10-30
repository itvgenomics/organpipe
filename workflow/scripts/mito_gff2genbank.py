from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.SeqRecord import SeqRecord
import argparse
from Bio.Seq import Seq
import os
import logging
import sys
from naming import warn


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Create GenBank file from the circularized mtGenome"
    )
    parser.add_argument("--code", help="Genetic code", required=True)
    parser.add_argument(
        "--annotation_dir", help="MITOS2 output directory of the assembly", required=True
    )
    parser.add_argument("--output", help="GenBank file to write", required=True)

    args = parser.parse_args()

    return args


def create_genbank(fasta_file, gff_file, output_file, code):
    logging.info(f"Writing GenBank file")

    # Read the FASTA sequence
    fasta_seq = SeqIO.read(fasta_file, "fasta")

    # sequence = Seq(str(fasta_seq.seq), generic_dna)
    sequence = Seq(str(fasta_seq.seq))

    # Create a SeqRecord object with the FASTA sequence
    seq_record = SeqRecord(
        sequence,
        id="Circularized_Contig",
        description="",
        annotations={"molecule_type": "DNA", "topology": "circular"},
    )

    # Read the GFF file
    with open(gff_file) as gff_handle:
        # Parse the GFF file and add features to the SeqRecord
        for line in gff_handle:
            if line.startswith("#"):
                continue

            parts = line.strip().split("\t")
            seq_id, source, feature_type, start, end, _, strand, _, attributes = parts

            strand_dict = {"-": -1, "+": +1}

            # Create a SeqFeature object for each feature in the GFF
            if feature_type in ["rRNA", "tRNA"]:
                if int(start) > int(end):
                    f1 = FeatureLocation(int(start) - 1, len(str(fasta_seq.seq)))
                    f2 = FeatureLocation(0, int(end))
                    seq_feature = SeqFeature(
                        CompoundLocation([f1, f2]), type=feature_type
                    )
                else:
                    seq_feature = SeqFeature(
                        FeatureLocation(
                            int(start) - 1, int(end), strand=strand_dict[strand]
                        ),
                        type=feature_type,
                    )
                # Parse the attributes and add relevant qualifiers to the SeqFeature
                attr_parts = attributes.split(";")
                for attr in attr_parts:
                    attr = attr.strip()
                    if "gene_id" in attr:
                        key, value = attr.split("=")
                        seq_feature.qualifiers["product"] = value

                # Add the SeqFeature to the SeqRecord
                seq_record.features.append(seq_feature)

            elif feature_type == "gene" and source == "mitos":
                if int(start) > int(end):
                    f1 = FeatureLocation(int(start) - 1, len(str(fasta_seq.seq)))
                    f2 = FeatureLocation(0, int(end))
                    seq_feature = SeqFeature(CompoundLocation([f1, f2]), type="gene")
                else:
                    seq_feature = SeqFeature(
                        FeatureLocation(
                            int(start) - 1, int(end), strand=strand_dict[strand]
                        ),
                        type="gene",
                    )
                # Parse the attributes and add relevant qualifiers to the SeqFeature
                attr_parts = attributes.split(";")
                for attr in attr_parts:
                    attr = attr.strip()
                    if "gene_id" in attr and ("trn" not in attr or "rrn" not in attr):
                        key, value = attr.split("=")
                        seq_feature.qualifiers["gene"] = value

                # Add the SeqFeature to the SeqRecord
                seq_record.features.append(seq_feature)

            elif feature_type == "exon":
                attr_parts = attributes.split(";")
                gene = attr_parts[1].replace("Name=", "")
                # tRNA and rRNA features are written above; every other exon is a coding gene
                if not gene.startswith(("trn", "rrn")):
                    if int(start) > int(end):
                        f1 = FeatureLocation(int(start) - 1, len(str(fasta_seq.seq)))
                        f2 = FeatureLocation(0, int(end))
                        seq_feature = SeqFeature(CompoundLocation([f1, f2]), type="CDS")
                    else:
                        seq_feature = SeqFeature(
                            FeatureLocation(
                                int(start) - 1, int(end), strand=strand_dict[strand]
                            ),
                            type="CDS",
                        )

                    # Parse the attributes and add relevant qualifiers to the SeqFeature
                    attr_parts = attributes.split(";")
                    for attr in attr_parts:
                        attr = attr.strip()
                        if "Name=" in attr:
                            key, value = attr.split("=")
                            if strand == "+":
                                seq_feature.qualifiers["translation"] = sequence[
                                    int(start) - 1 : int(end)
                                ].translate(table=code)
                                seq_feature.qualifiers["gene"] = value

                            elif strand == "-":
                                seq_feature.qualifiers["translation"] = (
                                    sequence[int(start) - 1 : int(end)]
                                    .reverse_complement()
                                    .translate(table=code)
                                )
                                seq_feature.qualifiers["gene"] = value
                else:
                    continue

                # Add the SeqFeature to the SeqRecord
                seq_record.features.append(seq_feature)

    # Write the SeqRecord to a GenBank file
    SeqIO.write(seq_record, output_file, "genbank")


if __name__ == "__main__":
    args = parse_arguments()

    FORMAT = "%(asctime)s [%(levelname)s] %(message)s"

    logging.basicConfig(
        level=logging.INFO,
        stream=sys.stdout,
        format=FORMAT,
        datefmt="%d-%m-%Y %H:%M:%S",
    )

    logging.info("Running genbank.py")

    logging.info(
        f"""
-   GENETIC CODE: {args.code}
-   ANNOTATION: {args.annotation_dir}
"""
    )

    fasta_file = f"{args.annotation_dir}/mitos2.fasta"
    gff_file = f"{args.annotation_dir}/result.gff"

    if os.path.exists(fasta_file) and os.path.exists(gff_file):
        create_genbank(
            fasta_file=fasta_file,
            gff_file=gff_file,
            output_file=args.output,
            code=args.code,
        )
    else:
        warn(f"MITOS2 produced no annotation in {args.annotation_dir}, so there is no GenBank file.")
