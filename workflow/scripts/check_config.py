import yaml
import argparse
import os
import pandas as pd
import re
import logging
import sys

# Expected fields and their descriptions
REQUIRED_FIELDS = {
    "sample": "Specify the sample name or 'all' to process all samples in reads_path.",
    "reads_path": "Provide the directory containing the sequencing read files.",
    "organelle": "Set to 'mito' for mitochondrial or 'chloro' for chloroplast genomes.",
    "genetic_code": "Define the genetic code number used for annotation (e.g., 1 = Standard Code).",
    "reference": "Path to a reference genome FASTA file (optional, leave blank if none).",
    "sequencing_type": "Type of sequencing data: 'Short' (Illumina) or 'Long' (PacBio/ONT).",
    "genome_range": "Expected genome size range in kb (e.g., 15000-18000).",
    "run_trimming": "Set to 'Yes' to run Fastp for trimming, 'No' to skip.",
    "adapters": "Path to the text file listing adapter sequences (used if run_trimming is 'Yes').",
    "minlength": "Minimum read length to retain after trimming.",
    "minquality": "Minimum base quality score threshold for trimming.",
    "seed_format": "Format of seed file: 'fasta' or 'genbank'. Leave blank for NCBI search.",
    "seed_file": "Path to the seed file used to start the assembly (optional).",
    "feature": "If using GenBank seed, specify feature type: 'CDS', 'rRNA', or 'tRNA'.",
    "search_ncbi": "Set to 'Yes' to search NCBI for seed sequences, 'No' otherwise.",
    "search_genes": "List of genes to search on NCBI (e.g., COI,16S). No spaces allowed.",
    "search_term": "Taxon name used when querying NCBI (e.g., 'Amphisbaena').",
    "max_references": "Maximum number of reference sequences to download per gene.",
    "kmers": "List of k-mer sizes for assembly (e.g., 21,33,55). No spaces allowed.",
    "max_memory": "Maximum memory in GB allowed during assembly.",
    "reads_length": "Read length (bp) of input sequencing data.",
    "insert_size": "Average insert size (bp) for paired-end reads.",
    "annotation": "Set to 'Yes' to run genome annotation after assembly.",
    "run_nhmmer": "Set to 'Yes' to identify ncRNAs with nhmmer, 'No' to skip.",
    "nhmmer_db": "Path to HMM database file used with nhmmer.",
    "run_images": "Set to 'Yes' to generate visualizations like genome maps and depth plots.",
}


def detect_delimiter(file_path):
    with open(file_path, "r") as file:
        first_line = file.readline()
        # Check for semicolon
        if ";" in first_line:
            return ";"
        else:
            return ","


def parse_samples(sample_field):
    if not sample_field:
        return []
    return [s.strip() for s in str(sample_field).split(",") if s.strip()]


def check_required_fields(config, required_fields):
    return [field for field in required_fields if field not in config]


def check_spacing_issues(config):
    errors = []

    sample = config.get("sample")
    if sample and " " in str(sample):
        errors.append(
            "Field `sample` must not contain spaces (e.g., use `sample1` not `sample 1`)."
        )

    kmers = config.get("kmers")
    if isinstance(kmers, str) and " " in kmers:
        errors.append(
            "Field `kmers` must be a comma-separated list **without spaces** (e.g., `19,23,33`)."
        )

    genes = config.get("search_genes")
    if isinstance(genes, str) and " " in genes:
        errors.append(
            "Field `search_genes` must be a comma-separated list **without spaces** (e.g., `COI,16S,NAD4`)."
        )

    # Generic: leading/trailing whitespace
    for key, value in config.items():
        if isinstance(value, str) and value != value.strip():
            errors.append(
                f"Field `{key}` contains leading/trailing whitespace. Please remove extra spaces."
            )

    return errors


def print_missing_fields(missing_fields):
    logging.info("\n Missing Required Fields:")
    for field in missing_fields:
        raise ValueError(
            f" - `{field}`: {REQUIRED_FIELDS[field]} \n Tip: Add these fields to your config.yaml with appropriate values.\n"
        )


def print_spacing_errors(spacing_errors):
    logging.info("\n Spacing Issues Detected:")
    for msg in spacing_errors:
        raise ValueError(f" - {msg} \n Tip: Remove any spaces where not allowed.\n")


def check_reads_path(config):
    errors = []
    path = config.get("reads_path")
    sample_field = config.get("sample", "")
    samples = parse_samples(sample_field)

    if not path:
        return ["Field `reads_path` is not set."]
    if not os.path.isdir(path):
        return [f"`reads_path` '{path}' is not a valid directory."]

    files = os.listdir(path)
    fastq_fasta_files = [
        f for f in files if f.endswith(".fastq.gz") or f.endswith(".fasta.gz")
    ]

    if not fastq_fasta_files:
        errors.append(
            f"`reads_path` '{path}' does not contain any .fastq.gz or .fasta.gz files."
        )

    r1_files = [f for f in fastq_fasta_files if "_R1" in f or "_pair1" in f]
    r2_files = [f for f in fastq_fasta_files if "_R2" in f or "_pair2" in f]

    if not r1_files:
        errors.append("No read files found with `_R1` or `_pair1` in the filename.")
    if not r2_files:
        errors.append("No read files found with `_R2` or `_pair2` in the filename.")

    # Sample validation
    if samples and samples != ["all"]:
        for sample in samples:
            sample_matches = [f for f in fastq_fasta_files if sample in f]
            if not sample_matches:
                errors.append(
                    f"No read filenames in '{path}' contain the sample name '{sample}'. "
                    "Filenames must include each sample name."
                )

    return errors


def check_organelle(config):
    errors = []
    organelle = config.get("organelle", "").strip().lower()

    if organelle not in {"mito", "chloro", "Mito", "Chloro"}:
        errors.append(
            f"Invalid value for `organelle`: '{config.get('organelle')}'. "
            "Accepted values are: 'mito' or 'chloro'."
        )

    return errors


def check_genetic_code(config):
    errors = []
    valid_codes = {
        1,
        2,
        3,
        4,
        5,
        6,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
        21,
        22,
        23,
        24,
        25,
        26,
        27,
        28,
        29,
        30,
        31,
        32,
        33,
    }

    code = config.get("genetic_code")

    valid_codes_desc = """
        1. The Standard Code
        2. The Vertebrate Mitochondrial Code
        3. The Yeast Mitochondrial Code
        4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
        5. The Invertebrate Mitochondrial Code
        6. The Ciliate, Dasycladacean and Hexamita Nuclear Code
        9. The Echinoderm and Flatworm Mitochondrial Code
        10. The Euplotid Nuclear Code
        11. The Bacterial, Archaeal and Plant Plastid Code
        12. The Alternative Yeast Nuclear Code
        13. The Ascidian Mitochondrial Code
        14. The Alternative Flatworm Mitochondrial Code
        15. Blepharisma Nuclear Code
        16. Chlorophycean Mitochondrial Code
        21. Trematode Mitochondrial Code
        22. Scenedesmus obliquus Mitochondrial Code
        23. Thraustochytrium Mitochondrial Code
        24. Rhabdopleuridae Mitochondrial Code
        25. Candidate Division SR1 and Gracilibacteria Code
        26. Pachysolen tannophilus Nuclear Code
        27. Karyorelict Nuclear Code
        28. Condylostoma Nuclear Code
        29. Mesodinium Nuclear Code
        30. Peritrich Nuclear Code
        31. Blastocrithidia Nuclear Code
        32. Balanophoraceae Plastid Code
        33. Cephalodiscidae Mitochondrial UAA-Tyr Code
        """

    try:
        code_int = int(code)
        if code_int not in valid_codes:
            errors.append(
                f"Invalid `genetic_code`: '{code}'. Must be one of the valid numeric codes {valid_codes_desc}"
            )
    except (ValueError, TypeError):
        errors.append(
            f"`genetic_code` must be a number (e.g., 1, 5, 11), but got: '{code}'."
        )

    return errors


def check_reference_file(config):
    errors = []
    reference = str(config.get("reference", "")).strip()

    if reference not in ["None", "nan"]:
        if not os.path.isfile(reference):
            errors.append(
                f"`reference` is set to '{reference}', but this file does not exist."
            )
        else:
            valid_exts = (".fasta", ".fa", ".fna", ".fasta.gz", ".fa.gz")
            if not reference.endswith(valid_exts):
                errors.append(
                    f"`reference` file '{reference}' does not have a valid FASTA extension. "
                    f"Expected one of: {', '.join(valid_exts)}"
                )

    return errors


def check_sequencing_type(config):
    errors = []
    seq_type = config.get("sequencing_type", "").strip().lower()

    if seq_type not in {"short", "long", "Short", "Long"}:
        errors.append(
            f"Invalid value for `sequencing_type`: '{config.get('sequencing_type')}'. "
            "Accepted values are: 'short' or 'long'."
        )

    return errors


def check_genome_range(config):
    errors = []
    genome_range = str(config.get("genome_range", "")).strip()
    seq_type = config.get("sequencing_type", "").strip().lower()

    pattern = r"^\d+-\d+$"

    if seq_type in {"short", "Short"}:
        if genome_range in ["None", "nan"]:
            errors.append("`genome_range` field is empty.")
        elif not re.match(pattern, genome_range):
            errors.append(
                f"Invalid format for `genome_range`: '{genome_range}'. "
                "Expected format is 'min-max' (e.g., '14000-22000')."
            )
        else:
            min_val, max_val = map(int, genome_range.split("-"))
            if min_val > max_val:
                errors.append(
                    f"`genome_range` values are invalid: min '{min_val}' is greater than max '{max_val}'."
                )

        return errors

    else:
        return None


def check_yes_no_fields(config):
    errors = []
    fields = ["annotation", "run_trimming", "run_nhmmer", "run_images"]
    allowed = {"yes", "no", "Yes", "No"}

    for field in fields:
        value = str(config.get(field, "")).strip().lower()
        if value not in allowed:
            errors.append(
                f"Invalid value for `{field}`: '{config.get(field)}'. "
                "Accepted values are: 'yes' or 'no'."
            )
    return errors


def check_adapters_file(config):
    errors = []
    adapters = str(config.get("adapters", "")).strip()

    if adapters not in ["None", "nan"]:
        if not os.path.isfile(adapters):
            errors.append(
                f"`adapters` is set to '{adapters}', but this file does not exist."
            )
        else:
            valid_exts = (".fasta", ".fa", ".fna", ".txt")
            if not adapters.endswith(valid_exts):
                errors.append(
                    f"`adapters` file '{adapters}' does not have a valid FASTA extension. "
                    f"Expected one of: {', '.join(valid_exts)}"
                )

    return errors


def check_seed_format_and_file(config):
    errors = []
    seed_format = str(config.get("seed_format", "")).strip().lower()
    seed_file = str(config.get("seed_file", "")).strip()

    if seed_format == "fasta":
        if seed_file in ["None", "nan"]:
            errors.append(
                f"`seed_format` is set to '{seed_format}', but no `seed_file` is provided."
            )
        elif not os.path.isfile(seed_file):
            errors.append(
                f"`seed_file` is set to '{seed_file}', but this file does not exist."
            )

        valid_exts = (".fasta", ".fa", ".fna")
        if not seed_file.endswith(valid_exts):
            errors.append(
                f"`seed_file` is set to '{seed_file}', but it doesn't end with a valid FASTA extension "
                f"(expected one of: {', '.join(valid_exts)})."
            )

    elif seed_format == "genbank":
        feature = str(config.get("feature", "")).strip()
        if not feature or feature not in ["CDS", "rRNA", "tRNA"]:
            errors.append(
                f"`feature` is not set appropriate. Must be 'CDS', 'rRNA' or 'tRNA' (case sensitive)."
            )

        if seed_file in ["None", "nan"]:
            errors.append(
                f"`seed_format` is set to '{seed_format}', but no `seed_file` is provided."
            )
        elif not os.path.isfile(seed_file):
            errors.append(
                f"`seed_file` is set to '{seed_file}', but this file does not exist."
            )

        valid_exts = (".gb", ".gbk", ".genbank")
        if not seed_file.endswith(valid_exts):
            errors.append(
                f"`seed_file` is set to '{seed_file}', but it doesn't end with a valid GENBANK extension "
                f"(expected one of: {', '.join(valid_exts)})."
            )

    elif seed_format in ["None", "nan"]:
        errors.append(
            f"Invalid value for `seed_format`: '{config.get('seed_format')}'. "
            "Accepted values are: 'fasta' or 'genbank' (case sensitive)."
        )

    return errors


def check_run_trimming_requirements(config):
    errors = []
    run_trimming = str(config.get("run_trimming", "")).strip().lower()

    if run_trimming in ["Yes", "yes"]:
        minlength = str(config.get("minlength", ""))
        minquality = str(config.get("minquality", ""))

        if minlength in ["None", "nan"]:
            errors.append("`run_trimming` is 'Yes', but `minlength` is not set.")
        if minquality in ["None", "nan"]:
            errors.append("`run_trimming` is 'Yes', but `minquality` is not set.")

    return errors


def check_search_ncbi_requirements(config):
    errors = []
    search_ncbi = str(config.get("search_ncbi", "")).strip().lower()
    seq_type = config.get("sequencing_type", "").strip().lower()

    if search_ncbi not in {"yes", "no", "Yes", "No"}:
        errors.append(
            f"Invalid value for `search_ncbi`: '{config.get('search_ncbi')}'. "
            "Accepted values are: 'yes' or 'no'."
        )

    if search_ncbi == "yes":
        if str(config.get("search_genes", "")).strip() in ["None", "nan"]:
            errors.append("`search_ncbi` is 'Yes', but `search_genes` is not set.")
        if str(config.get("search_term", "")).strip() in ["None", "nan"]:
            errors.append("`search_ncbi` is 'Yes', but `search_term` is not set.")
        if str(config.get("max_references", "")).strip() in ["None", "nan"]:
            errors.append("`search_ncbi` is 'Yes', but `max_references` is not set.")

    if seq_type in ["Long", "long"]:
        if search_ncbi == "yes":
            errors.append(
                "`seq_type` is 'long', but `search_ncbi` is set to 'yes'. Please change it to 'no'"
            )

    return errors


def main():
    parser = argparse.ArgumentParser(
        description="Check required fields and formatting in a config.yaml file."
    )
    parser.add_argument(
        "-c", "--config", required=True, help="Path to the YAML configuration file."
    )
    args = parser.parse_args()

    FORMAT = "%(asctime)s [%(levelname)s] %(message)s"

    logging.basicConfig(
        level=logging.INFO,
        stream=sys.stdout,
        format=FORMAT,
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    config_path = args.config

    logging.info(f"\n Checking configuration file: {config_path}")

    all_configs = []

    if config_path.endswith(".csv"):
        delimiter = detect_delimiter(config_path)
        df = pd.read_csv(config_path, delimiter=delimiter)
        for idx, row in df.iterrows():
            config = row.to_dict()
            all_configs.append(config)

    elif config_path.endswith((".yaml", ".yml")):

        try:
            with open(config_path, "r") as f:
                config = yaml.safe_load(f)
                all_configs.append(config)

        except FileNotFoundError:
            logging.info(f"\n Error: File '{config_path}' not found.")
            return

        except yaml.YAMLError as e:
            logging.info(f"\n Error: Failed to parse YAML file. Details:\n{e}")
            return

        if not isinstance(config, dict):
            logging.info(
                "\n Error: The config file format is invalid (expected a top-level dictionary)."
            )
            return

    for config in all_configs:
        sample = config.get("sample")
        logging.info(f"INFO: Running Check for Sample {sample}")

        missing_fields = check_required_fields(config, REQUIRED_FIELDS)
        spacing_errors = check_spacing_issues(config)

        if missing_fields:
            print_missing_fields(missing_fields)

        if spacing_errors:
            print_spacing_errors(spacing_errors)

        read_path_errors = check_reads_path(config)

        if read_path_errors:
            logging.info("\n Problem with `reads_path` content:")
            for msg in read_path_errors:
                raise ValueError(f" - {msg}")

        organelle_errors = check_organelle(config)

        if organelle_errors:
            logging.info("\n Problem with `organelle` field:")
            for msg in organelle_errors:
                raise ValueError(f" - {msg}")

        genetic_code_errors = check_genetic_code(config)

        if genetic_code_errors:
            logging.info("\n Problem with `genetic_code` field:")
            for msg in genetic_code_errors:
                raise ValueError(f" - {msg}")

        reference_errors = check_reference_file(config)

        if reference_errors:
            logging.info("\n Problem with `reference` field:")
            for msg in reference_errors:
                raise ValueError(f" - {msg}")

        sequencing_errors = check_sequencing_type(config)

        if sequencing_errors:
            logging.info("\n Problem with `sequencing_type` field:")
            for msg in sequencing_errors:
                raise ValueError(f" - {msg}")

        genome_range_errors = check_genome_range(config)

        if genome_range_errors:
            logging.info("\n Problem with `genome_range` field:")
            for msg in genome_range_errors:
                raise ValueError(f" - {msg}")

        yes_no_errors = check_yes_no_fields(config)

        if yes_no_errors:
            logging.info("\n Problem with Yes/No fields:")
            for msg in yes_no_errors:
                raise ValueError(f" - {msg}")

        adapters_errors = check_adapters_file(config)

        if adapters_errors:
            logging.info("\n Problem with `reference` field:")
            for msg in adapters_errors:
                raise ValueError(f" - {msg}")

        adapter_req_errors = check_run_trimming_requirements(config)
        if adapter_req_errors:
            logging.info("\n Problem with Trimming requirements:")
            for msg in adapter_req_errors:
                raise ValueError(f" - {msg}")

        seed_errors = check_seed_format_and_file(config)

        if seed_errors:
            logging.info("\n Problem with `seed_format` / `seed_file` fields:")
            for msg in seed_errors:
                raise ValueError(f" - {msg}")

        search_errors = check_search_ncbi_requirements(config)

        if search_errors:
            logging.info("\n Problem with NCBI Search Configuration:")
            for msg in search_errors:
                raise ValueError(f" - {msg}")


if __name__ == "__main__":
    main()
