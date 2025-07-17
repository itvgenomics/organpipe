import yaml
import argparse
import os
import pandas as pd
import re
import logging
import sys

# Expected fields and their descriptions
REQUIRED_FIELDS = {
    "sample": "Specify the sample name to be used. If you want to process all samples in the reads_path, use 'all'. The reads filenames should include delimiters like _R1/_R2 or _pair1/_pair2 to distinguish paired-end reads. For short reads, ensure that only the sequencing read files (in fastq.gz format) are present in reads_path. For long reads, ensure that the reads are in fasta format or in fastq.gz format.",
    "reads_path": "Provide the directory path where the sequencing read files (_R1/_R2 or _pair1/_pair2) are located.",
    "organelle": "Specify the type of organelle genome you are assembling: 'mito' for mitochondrial or 'chloro' for chloroplast.",
    "genetic_code": "Define the genetic code to be used for mitochondrial genome annotation (e.g., 1 for Standard Code).",
    "reference": "If a reference genome is available, specify its path as a fasta file. Leave blank if no reference is available.",
    "sequencing_type": "Indicate the type of sequencing data: 'Short' for Illumina short reads or 'Long' for PacBio/ONT reads.",
    "genome_range": "Specify the expected genome size range in kilobases (e.g., '15000-18000' for mitochondria).",
    "run_trimming": "Choose whether to run the Fastp software for adapter trimming and quality control. For long reads, only .fastq.gz files are supported. Options: 'Yes' to perform trimming, 'No' to skip this step.",
    "adapters": "Provide the path to a text file listing adapter sequences to be removed by Fastp. You can leave this blank if Fastp can identify your adapters automatically.",
    "minlength": "Set the minimum read length to retain after trimming. Reads shorter than this length will be discarded.",
    "minquality": "Set the minimum base quality score threshold for trimming low-quality bases from reads.",
    "pacbio_adapters": "If you are using PacBio reads and want to trim them, make sure your files are in fastq.gz format. Change the following variables you want to customize cutadapt settings. Example: '-b ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT -b ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT'",
    "seed_format": "Specify the path to a seed file to initialize assembly. Supported formats: 'fasta' or 'genbank'. For long-read sequencing, leave seed_format and feature blank.",
    "seed_file": "Path to the seed file to initialize assembly.",
    "feature": "If a GenBank file is provided as the seed file, specify the feature type to use for assembly: 'CDS', 'rRNA', or 'tRNA'.",
    "search_ncbi": "If no seed file is provided, indicate whether to search NCBI for reference sequences. Options: 'Yes' to search NCBI, 'No' to skip NCBI search. Ensure seed_format is blank if 'Yes' is selected.",
    "search_genes": "List the specific genes (e.g., 'COI,16S,ATP6') to search for on NCBI if NCBI search is enabled.",
    "search_term": "Specify the taxon name to use as a query when searching NCBI (e.g., 'Amphisbaena').",
    "max_references": "Set the maximum number of reference sequences to download per gene during the NCBI search.",
    "kmers": "Specify the list of k-mer sizes to be used in genome assembly (e.g., 19,23,28,33,39).",
    "max_memory": "Set the maximum amount of memory (in GB) that can be used during assembly.",
    "reads_length": "Provide the read length (in base pairs) and average insert size (in base pairs) of the sequencing data.",
    "insert_size": "Provide the read length (in base pairs) and average insert size (in base pairs) of the sequencing data.",
    "annotation": "Indicate whether to run the annotation pipeline for the assembled genome. Options: 'Yes' or 'No'.",
    "run_nhmmer": "Specify whether to run nhmmer for identifying non-coding RNA (ncRNA) and intergenic regions. Options: 'Yes' or 'No'.",
    "nhmmer_db": "Provide the path to the HMM database file to be used with nhmmer. We are currently working on a better database for improved accuracy. A parsing step will be implemented to enhance this process.",
    "run_images": "Indicate whether to generate visualizations, such as OGDraw diagrams, depth plots, and recruitment plots. Options: 'Yes' or 'No'.",
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
    sequencing_type = config.get("sequencing_type", "").strip().lower()

    if not path:
        return ["Field `reads_path` is not set."]
    if not os.path.isdir(path):
        return [f"`reads_path` '{path}' is not a valid directory."]

    files = os.listdir(path)

    if sequencing_type == "short":
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

    elif sequencing_type == "long":
        fasta_files = [
            f for f in files if f.endswith(".fasta") or f.endswith(".fastq.gz")
        ]

        if not fasta_files:
            errors.append(
                f"`reads_path` '{path}' does not contain any .fasta or .fastq.gz files for long reads."
            )

        # Sample validation
        if samples and samples != ["all"]:
            for sample in samples:
                sample_matches = [f for f in fasta_files if sample in f]
                if not sample_matches:
                    errors.append(
                        f"No .fasta.gz filenames in '{path}' contain the sample name '{sample}'. "
                        "Filenames must include each sample name."
                    )
    else:
        errors.append(
            f"`sequencing_type` must be 'short' or 'long', but found: '{config.get('sequencing_type')}'"
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

    sequencing_type = str(config.get("sequencing_type", "")).strip().lower()

    for field in fields:
        if field == "run_trimming":
            if sequencing_type == "long":
                allowed = {"None", "nan", "no", "No", "none", "Yes", "yes"}
                value = str(config.get(field, "")).strip().lower()
                if value not in allowed:
                    errors.append(
                        f"Invalid value for `{field}`: '{config.get(field)}'. "
                        "Accepted values are: empty, yes or 'no'."
                    )

            elif sequencing_type == "short":
                value = str(config.get(field, "")).strip().lower()
                if value not in allowed:
                    errors.append(
                        f"Invalid value for `{field}`: '{config.get(field)}'. "
                        "Accepted values are: 'yes' or 'no'."
                    )

        else:
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
    sequencing_type = config.get("sequencing_type", "").strip().lower()

    if sequencing_type == "short":
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

    if sequencing_type == "long":
        if seed_format not in ["None", "nan"] and seed_file not in ["None", "nan"]:

            errors.append(
                f"`sequencing_type` is set to '{sequencing_type}', so `seed_format` and `seed_file` must be empty. Set `search_ncbi` to 'yes' to search NCBI for seed sequences."
            )

            return errors


def check_run_trimming_requirements(config):
    errors = []
    run_trimming = str(config.get("run_trimming", "")).strip().lower()
    sequencing_type = str(config.get("sequencing_type", "")).strip().lower()
    pacbio_adapters = config.get("pacbio_adapters", "").strip()

    if run_trimming in ["Yes", "yes"]:
        if sequencing_type == "short":
            minlength = str(config.get("minlength", ""))
            minquality = str(config.get("minquality", ""))

            if minlength in ["None", "nan"]:
                errors.append("`run_trimming` is 'Yes', but `minlength` is not set.")
            if minquality in ["None", "nan"]:
                errors.append("`run_trimming` is 'Yes', but `minquality` is not set.")

        elif sequencing_type == "long":
            if pacbio_adapters in ["None", "nan"]:
                errors.append(
                    "`run_trimming` is 'Yes', but `pacbio_adapters` is not set."
                )

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

    if seq_type == "long":
        if search_ncbi == "yes":

            if str(config.get("search_term", "")).strip() in ["None", "nan"]:
                errors.append("`search_ncbi` is 'Yes', but `search_term` is not set.")

            if str(config.get("max_references", "")).strip() in ["None", "nan"]:
                errors.append(
                    "`search_ncbi` is 'Yes', but `max_references` is not set."
                )

        else:
            errors.append(
                f"Invalid value for `search_ncbi`: '{config.get('search_ncbi')}'."
                f"If your sequencing_type is set to `long`, you need to run with search_ncbi: `yes`."
            )

    elif seq_type == "short":
        if search_ncbi == "yes":
            if str(config.get("search_genes", "")).strip() in ["None", "nan"]:
                errors.append("`search_ncbi` is 'Yes', but `search_genes` is not set.")

            if str(config.get("search_term", "")).strip() in ["None", "nan"]:
                errors.append("`search_ncbi` is 'Yes', but `search_term` is not set.")

            if str(config.get("max_references", "")).strip() in ["None", "nan"]:
                errors.append(
                    "`search_ncbi` is 'Yes', but `max_references` is not set."
                )

    return errors


def check_pacbio_adapters(config):
    errors = []
    seq_type = config.get("sequencing_type", "").strip().lower()
    adapter_value = config.get("pacbio_adapters", "").strip()

    if seq_type == "long":
        if not adapter_value:
            errors.append(
                "`sequencing_type` is 'long' but `pacbio_adapters` is not set."
            )
        else:
            # Accept multiple -b flags, each followed by a DNA sequence
            pattern = r"^(?:-b\s+[ACGT]+(?:\s+|$))+$"
            if not re.match(pattern, adapter_value):
                errors.append(
                    f"`pacbio_adapters` value '{adapter_value}' is not valid.\n"
                    "It must contain one or more adapter sequences in this format:\n"
                    "  -b ACTG... -b ACTG... (valid characters: A, C, G, T)"
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

        pacbio_adapters = check_pacbio_adapters(config)

        if pacbio_adapters:
            logging.info("\n Problem with `pacbio_adapters` setting:")
            for msg in pacbio_adapters:
                logging.info(f" - {msg}")


if __name__ == "__main__":
    main()
