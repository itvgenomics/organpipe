import yaml
import argparse
import os
import shutil
import pandas as pd
import logging
import sys

# Custom scripts
import get_seeds
import download_seeds
import download_references
from naming import SHORT_READ_SUFFIXES, LONG_READ_SUFFIXES, mate_marker, read_stem

# Fields written in lowercase, so the rules can compare them directly.
LOWERCASE_FIELDS = [
    "organelle",
    "sequencing_type",
    "run_trimming",
    "annotation",
    "run_nhmmer",
    "run_images",
    "run_novoplasty",
    "run_getorganelle",
    "search_ncbi",
]

INTEGER_FIELDS = [
    "genetic_code",
    "max_memory",
    "reads_length",
    "insert_size",
    "max_references",
    "minlength",
    "minquality",
    "n_references",
    "n_rounds",
]

TEXT_FIELDS = [
    "reads_path",
    "genome_range",
    "seed_format",
    "seed_file",
    "feature",
    "search_genes",
    "search_term",
    "adapters",
    "pacbio_adapters",
    "reference",
    "search_species",
    "database",
    "target_size",
    "extra_flags",
    "spades_kmers",
]


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Run and parse mitochondrial genome assembly data"
    )

    parser.add_argument("--configfile", help="Path to the config file")

    args = parser.parse_args()

    return args


def detect_delimiter(file_path):
    with open(file_path, "r") as file:
        first_line = file.readline()
        # Check for semicolon
        if ";" in first_line:
            return ";"
        else:
            return ","  # Default to comma


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
    search_species,
    n_references
):

    if sequencing_type == "short":
        if seed_format.lower() == "genbank":

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

            # Extract the headers to use them as seeds in the Snakemake config file
            seeds = get_seeds.get_fasta_headers(
                fasta_file=f"resources/{sample}/seeds/seeds.fasta"
            )

            return seeds

        elif seed_format.lower() == "fasta":

            if not os.path.exists(f"resources/{sample}/seeds/seeds.fasta"):

                shutil.copy(seed_file, f"resources/{sample}/seeds/seeds.fasta")
                # Split the seeds file into multiple sequences
                get_seeds.split_fasta(
                    fasta_file=seed_file,
                    output_dir=f"resources/{sample}/seeds/",
                )
            # Extract the headers to use them as seeds in the Snakemake config file
            seeds = get_seeds.get_fasta_headers(
                fasta_file=f"resources/{sample}/seeds/seeds.fasta"
            )

            return seeds

        elif search_ncbi == "yes":
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

        return []

    seeds = [
        f.split(".fasta")[0]
        for f in os.listdir(f"resources/{sample}/seeds/")
        if f.endswith(".fasta")
    ]

    if len(seeds) == 0:
        download_references.main(
            species=search_species,
            email="ncbiapirunner@gmail.com",
            outfolder=f"resources/{sample}/seeds/",
            min_length=10000,
            n=int(n_references),
        )

    seeds = [
        f.split(".fasta")[0]
        for f in os.listdir(f"resources/{sample}/seeds/")
        if f.endswith(".fasta")
    ]

    if len(seeds) == 0:
        print(
            f"\nERROR: No reference mitogenome was found on NCBI for sample {sample} "
            f"and search_species '{search_species}'.\n"
            "\n"
            "Please try one of the following:\n"
            "  - Use a broader taxon in search_species.\n"
            f"  - Copy a reference (<name>.fasta and <name>.gb) to resources/{sample}/seeds/.\n"
        )
        raise SystemExit(1)

    return seeds


def get_samples_ids(reads_path, sequencing_type):
    """Return the samples of `sample: all`: every read file or read pair in reads_path."""
    suffixes = SHORT_READ_SUFFIXES if sequencing_type == "short" else LONG_READ_SUFFIXES
    sample_ids = set()

    for filename in os.listdir(reads_path):
        if not filename.endswith(suffixes) or "Undetermined" in filename:
            continue

        marker = mate_marker(filename)
        if sequencing_type == "short" and marker and marker.group(1)[-1] == "1":
            sample_ids.add(filename[: marker.start()])
        elif sequencing_type == "long" and not marker:
            sample_ids.add(read_stem(filename))

    if not sample_ids:
        raise SystemExit(f"\nERROR: `sample: all`, but no read files were found in {reads_path}.\n")

    return sorted(sample_ids)


def read_config_rows(configfile):
    """Return one dict per sample, with every value as text and blank values as ""."""
    if configfile.endswith((".yaml", ".yml")):
        with open(configfile, "r") as file:
            data = yaml.safe_load(file)

        data = {key: "" if value is None else str(value) for key, value in data.items()}
        sequencing_type = data["sequencing_type"].strip().lower()

        if data["sample"].strip() == "all":
            sample_ids = get_samples_ids(data["reads_path"], sequencing_type)
        else:
            sample_ids = [sample.strip() for sample in data["sample"].split(",")]

        return [{**data, "sample": sample} for sample in sample_ids]

    delimiter = detect_delimiter(configfile)
    df = pd.read_csv(configfile, delimiter=delimiter, dtype=str, keep_default_na=False)
    return df.to_dict("records")


def to_int(value):
    return int(float(value)) if value != "" else ""


def split_list(value):
    return [item.strip() for item in value.split(",") if item.strip()]


def sample_config(row):
    """Build the Snakemake config entry of one sample, preparing its seeds on the way."""
    sample = row["sample"]

    logging.info(f"Preparing Sample: {sample}")

    create_dirs(sample)

    entry = {field: row[field].strip().lower() for field in LOWERCASE_FIELDS}
    entry.update({field: to_int(row[field]) for field in INTEGER_FIELDS})
    entry.update({field: row[field] for field in TEXT_FIELDS})

    entry["kmers"] = split_list(row["kmers"])
    entry["seeds"] = prepare_seeds(
        seed_format=row["seed_format"],
        seed_file=row["seed_file"],
        sample=sample,
        search_ncbi=entry["search_ncbi"],
        search_term=row["search_term"],
        genes=split_list(row["search_genes"]),
        organelle=entry["organelle"],
        max_references=entry["max_references"],
        feature=row["feature"],
        sequencing_type=entry["sequencing_type"],
        search_species=row["search_species"],
        n_references=row["n_references"],
    )

    return entry


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

    samples_dict = {
        "samples": {row["sample"]: sample_config(row) for row in read_config_rows(args.configfile)}
    }

    with open("config/snakemake_config.yaml", "w") as file:
        yaml.dump(samples_dict, file)
