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
        "--max_memory",
        help="Maximum memory to be used (e.g., 8G)",
    )
    parser.add_argument(
        "--kmer",
        help="kmer to be used (e.g., 33)",
    )
    parser.add_argument(
        "--seed",
        help="Seed name to be used",
    )
    parser.add_argument("--reference", help="Path to the reference sequence", nargs="?")
    parser.add_argument(
        "--hashtable",
        help="Path to the hash table txt file",
    )
    parser.add_argument(
        "--hash2b",
        help="Path to the hash2b table txt file",
    )
    parser.add_argument(
        "--hash2c",
        help="Path to the hash2c table txt file",
    )
    args = parser.parse_args()
    return args


def create_config_file(
    organelle,
    sample,
    genome_range,
    reads_length,
    insert_size,
    max_memory,
    kmer,
    seed,
    reference,
    hashtable,
    hash2b,
    hash2c,
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
    Seed Input            = resources/{sample}/seeds/{seed}.fasta
    Reference sequence    = {reference}
    Variance detection    =
    Chloroplast sequence  =

    Dataset 1:
    -----------------------
    Read Length           = {reads_length}
    Insert size           = {insert_size}
    Platform              = illumina
    Single/Paired         = PE
    Combined reads        =
    Forward reads         = {hash2b}
    Reverse reads         = {hash2c}
    Store Hash            = {hashtable}

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
    Output path           = results/{sample}/novoplasty/{seed}/kmer{kmer}/
    """

    os.makedirs(f"results/{sample}/novoplasty/{seed}/kmer{kmer}", exist_ok=True)

    config_file = config_file.replace("= None", "= ")

    with open(
        f"results/{sample}/novoplasty/{seed}/kmer{kmer}/config.txt", "w+"
    ) as file:
        file.write(config_file)


if __name__ == "__main__":
    args = parse_arguments()
    create_config_file(
        organelle=args.organelle,
        sample=args.sample,
        genome_range=args.genome_range,
        reads_length=args.reads_length,
        insert_size=args.insert_size,
        max_memory=args.max_memory,
        kmer=args.kmer,
        seed=args.seed,
        reference=args.reference,
        hashtable=args.hashtable,
        hash2b=args.hash2b,
        hash2c=args.hash2c,
    )
