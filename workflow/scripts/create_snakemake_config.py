import yaml
import argparse
import os
import shutil
import pandas as pd
import subprocess
import csv
import logging
import sys

# Custom scripts
import get_seeds
import download_seeds
import download_references


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Run and parse mitochondrial genome assembly data"
    )

    parser.add_argument("--configfile", help="Path to the config file")

    args = parser.parse_args()

    return args

def detect_delimiter(file_path):
    with open(file_path, 'r') as file:
        first_line = file.readline()
        # Check for semicolon
        if ';' in first_line:
            return ';'
        else:
            return ','  # Default to comma

def create_dirs(sample):
    dirs = [
        f"resources/{sample}",
        f"resources/{sample}/rawreads",
        f"resources/{sample}/seeds",
    ]

    for dir in dirs:
        os.makedirs(dir, exist_ok=True)


def prepare_seeds(
    seed_format,
    seed_file,
    sample,
    search_ncbi,
    search_term,
    max_references,
    genes,
    organelle,
    feature,
    sequencing_type,
):
    if str(seed_format).lower() == "genbank":

        # Extract the sequences from the gb file
        genbank_file = seed_file

        if not os.path.exists(f"resources/{sample}/seeds/seeds.fasta"):

            valid_features = ["CDS", "rRNA", "tRNA"]
            if feature not in valid_features:
                raise ValueError(
                    f"Invalid feature: {feature}. Must be one of {valid_features}."
                )

            get_seeds.main(
                featuretype=feature,
                genbank_file=genbank_file,
                output_file=f"resources/{sample}/seeds/seeds.fasta",
            )

            # Split the seeds file into multiple sequences
            get_seeds.split_fasta(
                fasta_file=f"resources/{sample}/seeds/seeds.fasta",
                output_dir=f"resources/{sample}/seeds/",
            )

        # Extract the headers to use as seeds into snakemake config file
        seeds = get_seeds.get_fasta_headers(
            fasta_file=f"resources/{sample}/seeds/seeds.fasta"
        )

        return seeds

    elif str(seed_format).lower() == "fasta":

        if not os.path.exists(f"resources/{sample}/seeds/seeds.fasta"):

            shutil.copy(seed_file, f"resources/{sample}/seeds/seeds.fasta")
            # Split the seeds file into multiple sequences
            get_seeds.split_fasta(
                fasta_file=seed_file,
                output_dir=f"resources/{sample}/seeds/",
            )
        # Extract the headers to use as seeds into snakemake config file
        seeds = get_seeds.get_fasta_headers(
            fasta_file=f"resources/{sample}/seeds/seeds.fasta"
        )

        return seeds

    elif str(search_ncbi).lower() == "yes":
        if str(sequencing_type).lower() == "short":
            if not os.path.exists(f"resources/{sample}/seeds/seeds.fasta"):
                outpath = f"resources/{sample}/seeds/"
                download_seeds.main(
                    taxon=search_term,
                    outpath=outpath,
                    maxcount=max_references,
                    genes=genes,
                    organelle=organelle,
                )

                # Split the seeds file into multiple sequences
                get_seeds.split_fasta(
                    fasta_file=f"resources/{sample}/seeds/seeds.fasta",
                    output_dir=f"resources/{sample}/seeds/",
                )

            # Extract the headers to use as seeds into snakemake config file
            seeds = get_seeds.get_fasta_headers(
                fasta_file=f"resources/{sample}/seeds/seeds.fasta"
            )
            return seeds
        else:
            seeds = [
                f.split(".fasta")[0]
                for f in os.listdir(f"resources/{sample}/seeds/")
                if f.endswith(".fasta")
            ]

            if len(seeds) == 0:
                download_references.main(
                    species=search_term,
                    email="ncbiapirunner@gmail.com",
                    outfolder=f"resources/{sample}/seeds/",
                    min_length=14000,
                    n=int(max_references),
                )

            return seeds


def get_samples_ids(reads_path):
    # List all files in the directory
    file_list = os.listdir(reads_path)

    # Extract sample IDs
    sample_ids = set(
        filename.split("_")[0]
        for filename in file_list
        if filename.endswith(".fasta.gz")
    )

    return sample_ids


def get_samples_ids(reads_path):
    # List all files in the directory
    file_list = os.listdir(reads_path)
    sample_ids = []

    for filename in file_list:
        if str(filename).endswith(".fasta.gz"):
            sample_id = filename.split(".fasta.gz")[0]
            sample_ids.append(
                str(sample_id)
                .replace("_R1", "")
                .replace("_R2", "")
                .replace("_pair1", "")
                .replace("_pair2", "")
            )
        elif str(filename).endswith(".fastq.gz"):
            sample_id = filename.split(".fastq.gz")[0]
            sample_ids.append(
                str(sample_id)
                .replace("_R1", "")
                .replace("_R2", "")
                .replace("_pair1", "")
                .replace("_pair2", "")
            )

    return set(sample_ids)


def replace_none_with_empty_string(d):
    for key, value in d.items():
        if isinstance(value, dict):
            replace_none_with_empty_string(
                value
            )  # Recursively handle nested dictionaries
        elif value is None:
            d[key] = ""  # Replace None with an empty string


if __name__ == "__main__":
    args = parse_arguments()

    FORMAT = "%(asctime)s [%(levelname)s] %(message)s"

    logging.basicConfig(
        level=logging.INFO,
        stream=sys.stdout,
        format=FORMAT,
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    ## DEBUG
    # Print all arguments
    print("\nParsed Arguments:")
    for arg in vars(args):
        print(f"{arg}: {getattr(args, arg)}")

    if str(args.configfile).endswith(".yaml") or str(args.configfile).endswith(".yml"):
        with open(args.configfile, "r") as file:
            # Load the YAML content
            data = yaml.safe_load(file)

        if "all" in data["sample"]:
            sample_ids = get_samples_ids(data["reads_path"])
        else:
            sample_ids = [
                sample_id.strip() for sample_id in str(data["sample"]).split(",")
            ]


        samples_dict = {"samples": {}}
        for sample_id in sample_ids:
            logging.info(f"Preparing Sample: {sample_id}")

            sample = sample_id
            kmers = data["kmers"].split(",")

            create_dirs(sample)

            if data["search_genes"] == None:
                genes = ""
            else:
                genes = data["search_genes"].split(",")

            seeds = prepare_seeds(
                seed_format=data["seed_format"],
                seed_file=data["seed_file"],
                sample=sample_id,
                search_ncbi=data["search_ncbi"],
                search_term=data["search_term"],
                genes=genes,
                organelle=data["organelle"],
                max_references=data["max_references"],
                feature=data["feature"],
                sequencing_type=data["sequencing_type"],
            )

            samples_dict["samples"][sample] = {
                "seeds": seeds,
                "kmers": kmers,
                "reads_path": data["reads_path"],
                "organelle": data["organelle"],
                "genetic_code": (
                    int(data["genetic_code"])
                    if data["genetic_code"] != None
                    else str(data["genetic_code"])
                ),
                "sequencing_type": data["sequencing_type"],
                "genome_range": data["genome_range"],
                "annotation": data["annotation"],
                "adapterremoval": data["adapterremoval"],
                "seed_format": data["seed_format"],
                "seed_file": data["seed_file"],
                "feature": data["feature"],
                "search_ncbi": data["search_ncbi"],
                "search_genes": data["search_genes"],
                "search_term": data["search_term"],
                "max_memory": (
                    int(data["max_memory"]) if data["max_memory"] != None else ""
                ),
                "reads_length": (
                    int(data["reads_length"]) if data["reads_length"] != None else ""
                ),
                "insert_size": (
                    int(data["insert_size"]) if data["insert_size"] != None else ""
                ),
                "reads_path": data["reads_path"],
                "adapters": data["adapters"],
                "max_references": (
                    int(data["max_references"])
                    if data["max_references"] != None
                    else ""
                ),
                "minlength": (
                    int(data["minlength"]) if data["minlength"] != None else ""
                ),
                "minquality": (
                    int(data["minquality"]) if data["minquality"] != None else ""
                ),
                "reference": data["reference"],
                "run_nhmmer": data["run_nhmmer"],
                "nhmmer_db": data["nhmmer_db"],
                "run_images": data["run_images"],
            }

        replace_none_with_empty_string(samples_dict)

        with open("config/snakemake_config.yaml", "w") as file:
            yaml.dump(samples_dict, file)

    elif str(args.configfile).endswith(".csv"):
        delimiter = detect_delimiter(args.configfile)
        df = pd.read_csv(args.configfile, delimiter=delimiter)

        df.fillna("", inplace=True)

        # Initiate the samples dict
        samples_dict = {"samples": {}}

        for idx, row in df.iterrows():
            sample = row["sample"]

            logging.info(f"Preparing Sample: {sample}")

            kmers = row["kmers"].split(",")
            genes = row["search_genes"].split(",")

            create_dirs(sample)

            seeds = prepare_seeds(
                seed_format=row["seed_format"],
                seed_file=row["seed_file"],
                sample=sample,
                search_ncbi=row["search_ncbi"],
                search_term=row["search_term"],
                genes=genes,
                organelle=row["organelle"],
                max_references=row["max_references"],
                feature=row["feature"],
                sequencing_type=row["sequencing_type"],
            )

            samples_dict["samples"][sample] = {
                "seeds": seeds,
                "kmers": kmers,
                "reads_path": row["reads_path"],
                "organelle": row["organelle"],
                "genetic_code": (
                    int(row["genetic_code"])
                    if row["genetic_code"] != ""
                    else str(row["genetic_code"])
                ),
                "sequencing_type": row["sequencing_type"],
                "genome_range": row["genome_range"],
                "annotation": row["annotation"],
                "adapterremoval": row["adapterremoval"],
                "seed_format": row["seed_format"],
                "seed_file": row["seed_file"],
                "feature": row["feature"],
                "search_ncbi": row["search_ncbi"],
                "search_genes": row["search_genes"],
                "search_term": row["search_term"],
                "max_memory": (
                    int(row["max_memory"])
                    if row["max_memory"] != ""
                    else str(row["max_memory"])
                ),
                "reads_length": (
                    int(row["reads_length"])
                    if row["reads_length"] != ""
                    else str(row["reads_length"])
                ),
                "insert_size": (
                    int(row["insert_size"])
                    if row["insert_size"] != ""
                    else str(row["insert_size"])
                ),
                "reads_path": row["reads_path"],
                "adapters": row["adapters"],
                "max_references": (
                    int(row["max_references"])
                    if row["max_references"] != ""
                    else str(row["max_references"])
                ),
                "minlength": (
                    int(row["minlength"])
                    if row["minlength"] != ""
                    else str(row["minlength"])
                ),
                "minquality": (
                    int(row["minquality"])
                    if row["minquality"] != ""
                    else str(row["minquality"])
                ),
                "reference": row["reference"],
                "run_nhmmer": row["run_nhmmer"],
                "nhmmer_db": row["nhmmer_db"],
                "run_images": row["run_images"],
            }

        with open("config/snakemake_config.yaml", "w") as file:
            yaml.dump(samples_dict, file)
