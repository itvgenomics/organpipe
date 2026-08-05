from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.SeqRecord import SeqRecord
import argparse
from Bio.Seq import Seq
import os
import re
import logging
import sys
import subprocess


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Create GenBank file from the circularized mtGenome"
    )
    parser.add_argument("--code", help="Genetic code", required=True)
    parser.add_argument("--sample", help="Sample code", required=True)
    parser.add_argument(
        "--seed", help="Seed used to assemble the mitogenome", required=False, nargs="?"
    )
    parser.add_argument(
        "--kmer", help="kmer used to assemble the mitogenome", required=False, nargs="?"
    )
    parser.add_argument(
        "--assembler", help="Assembler used to assemble the mitogenome", required=True
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
                if gene in ["nad1", "nad2","cox1", "cox2", "atp8", "atp6", "cox3", "nad3", "nad4l", "nad4", "nad5", "nad6", "cob"]:
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
-   SAMPLE: {args.sample}
-   SEED: {args.seed}
-   KMER: {args.kmer}
"""
    )

    sample = args.sample
    seed = args.seed
    kmer = args.kmer
    code = args.code

    if args.assembler == "novoplasty":
        annotations_path = f"results/{sample}/mitos2/novoplasty/{seed}_kmer{kmer}"

        if not os.path.exists(f"results/{sample}/genbanks/novoplasty"):
            os.makedirs(f"results/{sample}/genbanks/novoplasty")

        annotations_dir = list_directories(annotations_path)
        for annotation in annotations_dir:
            gb_output_file = f"results/{sample}/genbanks/novoplasty/{annotation}.gb"
            fasta_output_file = f"results/{sample}/genbanks/novoplasty/{annotation}.fasta"

            if os.path.exists(
                f"{annotations_path}/{annotation}/mitos2.fasta"
            ) and os.path.exists(f"{annotations_path}/{annotation}/result.gff"):

                create_genbank(
                    fasta_file=f"{annotations_path}/{annotation}/mitos2.fasta",
                    gff_file=f"{annotations_path}/{annotation}/result.gff",
                    output_file=gb_output_file,
                    code=code,
                )

    elif args.assembler == "getorganelle":
        annotations_path = f"results/{sample}/mitos2/getorganelle"

        if not os.path.exists(f"results/{sample}/genbanks/getorganelle"):
            os.makedirs(f"results/{sample}/genbanks/getorganelle")

        annotations_dir = list_directories(annotations_path)
        for annotation in annotations_dir:
            gb_output_file = f"results/{sample}/genbanks/getorganelle/{annotation}.gb"
            fasta_output_file = f"results/{sample}/genbanks/getorganelle/{annotation}.fasta"

            if os.path.exists(
                f"{annotations_path}/{annotation}/mitos2.fasta"
            ) and os.path.exists(f"{annotations_path}/{annotation}/result.gff"):

                create_genbank(
                    fasta_file=f"{annotations_path}/{annotation}/mitos2.fasta",
                    gff_file=f"{annotations_path}/{annotation}/result.gff",
                    output_file=gb_output_file,
                    code=code,
                )
