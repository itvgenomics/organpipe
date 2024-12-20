import os
import argparse


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Create a configuration file for organelle genome assembly"
    )

    parser.add_argument(
        "--organelle",
        help="Type of organelle (e.g., mito or chloro)",
    )
    parser.add_argument(
        "--sample",
        help="Sample name for the project",
    )
    parser.add_argument(
        "--genome_range",
        help="Genome range (e.g., 12000-22000)",
    )
    parser.add_argument(
        "--reads_length",
        help="Read length (e.g., 150)",
    )
    parser.add_argument(
        "--insert_size",
        help="Insert size (e.g., 300)",
    )
    parser.add_argument(
        "--forward",
        help="Path to forward reads file",
    )
    parser.add_argument(
        "--reverse",
        help="Path to reverse reads file",
    )

    parser.add_argument(
        "--max_memory",
        help="Maximum memory to be used (e.g., 8G)",
    )
    parser.add_argument(
        "--kmer",
        help="kmer size to build the hashtable (e.g., 33)",
    )

    args = parser.parse_args()
    return args


def create_config_file(
    organelle,
    sample,
    genome_range,
    reads_length,
    insert_size,
    forward,
    reverse,
    max_memory,
    kmer,
):

    config_file = f"""
    Project:
    -----------------------
    Project name          = {sample}
    Type                  = {organelle}
    Genome Range          = {genome_range}
    K-mer                 = {kmer}
    Max memory            = {max_memory}
    Extended log          =
    Save assembled reads  = yes
    Seed Input            = resources/initial_seed.fasta
    Reference sequence    = 
    Variance detection    =
    Chloroplast sequence  =

    Dataset 1:
    -----------------------
    Read Length           = {reads_length}
    Insert size           = {insert_size}
    Platform              = illumina
    Single/Paired         = PE
    Combined reads        =
    Forward reads         = {forward}
    Reverse reads         = {reverse}
    Store Hash            = yes

    Heteroplasmy:
    -----------------------
    MAF                   =
    HP exclude list       =
    PCR-free              =

    Optional:
    -----------------------
    Insert size auto      = yes
    Insert Range          =
    Insert Range strict   =
    Use Quality Scores    = no
    Extented log          = 1
    Output path           = results/{sample}/hashtable/kmer{kmer}/
    """

    with open(f"results/{sample}/hashtable/kmer{kmer}/hash_config.txt", "w+") as file:
        file.write(config_file)


if __name__ == "__main__":
    args = parse_arguments()
    create_config_file(
        args.organelle,
        args.sample,
        args.genome_range,
        args.reads_length,
        args.insert_size,
        args.forward,
        args.reverse,
        args.max_memory,
        args.kmer,
    )
