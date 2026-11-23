import yaml
import argparse
import os
import pandas as pd
import re
import logging
import sys
from config_values import as_text, read_layout
from naming import find_long_reads, find_short_reads, find_single_reads

# Samples must match the `sample` wildcard in workflow/Snakefile.
VALID_SAMPLE_NAME = r"^[A-Za-z0-9_]+$"

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
    "run_images": "Indicate whether to generate visualizations, such as OGDraw diagrams, depth plots, and recruitment plots. Options: 'Yes' or 'No'.",
    "search_species": "Specify the taxon name to use as a query when searching NCBI (e.g., 'Amphisbaena') for complete mitogenome references.",
    "n_references": "Set the maximum number of reference sequences to download during the NCBI search (e.g., 5).",
    "run_novoplasty": "Indicate whether to run NOVOPlasty for genome assembly. Options: 'Yes' to run NOVOPlasty, 'No' to skip this step.",
    "run_getorganelle": "Indicate whether to run GetOrganelle for genome assembly. Options: 'Yes' to run GetOrganelle, 'No' to skip this step.",
    "database": "GetOrganelle organelle type (e.g., 'embplant_pt', 'animal_mt'). Required if run_getorganelle is 'Yes'.",
    "n_rounds": "Maximum number of extending rounds (suggested: >=2).",
    "target_size": "Hypothetical value(s) of target genome size.",
    "spades_kmers": "SPAdes kmer settings. Use the same format as in SPAdes. e.g: 21,55,85,115",
    "extra_flags": "Add any additional flags to be passed to the GetOrganelle command line.",
}


def detect_delimiter(file_path):
    with open(file_path, "r") as file:
        first_line = file.readline()
        # Check for semicolon
        if ";" in first_line:
            return ";"
        else:
            return ","


def get_text(config, field):
    """Return the value of a field as text, with blank values as ""."""
    return as_text(config.get(field))


def runs_novoplasty(config):
    """NOVOPlasty only runs on short reads. Its fields are ignored when it does not run."""
    return (
        get_text(config, "sequencing_type").lower() == "short"
        and get_text(config, "run_novoplasty").lower() == "yes"
    )


def parse_samples(sample_field):
    if not sample_field:
        return []
    return [s.strip() for s in str(sample_field).split(",") if s.strip()]


def check_required_fields(config, required_fields):
    return [field for field in required_fields if field not in config]


def check_sample_names(config):
    errors = []
    samples = parse_samples(config.get("sample"))

    if samples != ["all"]:
        for sample in samples:
            if not re.match(VALID_SAMPLE_NAME, sample):
                errors.append(
                    f"Sample name '{sample}' is not valid. Use only letters, numbers and underscores."
                )

    return errors


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

    spades_kmers = config.get("spades_kmers")
    if isinstance(spades_kmers, str) and " " in spades_kmers:
        errors.append(
            "Field `spades_kmers` must be a comma-separated list **without spaces** (e.g., `21,55,85,115`)."
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

        if read_layout(config) == "single":
            if samples and samples != ["all"]:
                for sample in samples:
                    try:
                        find_single_reads(sample, path)
                    except ValueError as e:
                        errors.append(str(e))
            return errors

        r1_files = [f for f in fastq_fasta_files if "_R1" in f or "_pair1" in f]
        r2_files = [f for f in fastq_fasta_files if "_R2" in f or "_pair2" in f]

        if not r1_files:
            errors.append("No read files found with `_R1` or `_pair1` in the filename.")
        if not r2_files:
            errors.append("No read files found with `_R2` or `_pair2` in the filename.")

        # Sample validation
        if samples and samples != ["all"]:
            for sample in samples:
                try:
                    find_short_reads(sample, path)
                except ValueError as e:
                    errors.append(str(e))

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
            run_trimming = get_text(config, "run_trimming").lower()
            for sample in samples:
                try:
                    reads = find_long_reads(sample, path)
                except ValueError as e:
                    errors.append(str(e))
                    continue

                if run_trimming == "yes" and not reads.endswith(".fastq.gz"):
                    errors.append(
                        f"`run_trimming` is 'Yes', but the reads of sample '{sample}' are not in "
                        f".fastq.gz format ({reads}). Trimming is only supported for .fastq.gz files."
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

    code = get_text(config, "genetic_code")
    sequencing_type = get_text(config, "sequencing_type").lower()
    organelle = get_text(config, "organelle").lower()
    annotation = get_text(config, "annotation").lower()

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

    if not code:
        # Used by MitoHiFi and by the MITOS2 annotation of mitochondria
        if sequencing_type == "long" or (organelle == "mito" and annotation == "yes"):
            errors.append(
                "`genetic_code` is not set. It is required for long reads and for the annotation of mitochondria."
            )
    elif not re.fullmatch(r"\d+(\.0)?", code):
        errors.append(
            f"`genetic_code` must be a number (e.g., 1, 5, 11), but got: '{code}'."
        )
    elif int(float(code)) not in valid_codes:
        errors.append(
            f"Invalid `genetic_code`: '{code}'. Must be one of the valid numeric codes {valid_codes_desc}"
        )

    return errors


def check_reference_file(config):
    errors = []
    reference = get_text(config, "reference")

    # Only NOVOPlasty uses the reference
    if reference and runs_novoplasty(config):
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


def check_read_layout(config):
    """Short reads come in _R1/_R2 files (paired) or in one file of single-end or merged reads (single)."""
    errors = []
    layout = read_layout(config)

    if layout not in {"paired", "single"}:
        errors.append(
            f"Invalid value for `read_layout`: '{config.get('read_layout')}'. "
            "Accepted values are: 'paired', 'single' or empty (paired)."
        )
    elif layout == "single" and runs_novoplasty(config):
        errors.append(
            "`read_layout` is 'single', but `run_novoplasty` is 'Yes'. NOVOPlasty needs paired reads: "
            "set `run_novoplasty` to 'No' and assemble the single-end or merged reads with GetOrganelle."
        )

    return errors


def check_genome_range(config):
    errors = []
    genome_range = get_text(config, "genome_range")

    pattern = r"^\d+-\d+$"

    if runs_novoplasty(config):
        if not genome_range:
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


def check_yes_no_fields(config):
    errors = []
    fields = ["annotation", "run_trimming", "run_nhmmer", "run_images", "run_novoplasty", "run_getorganelle"]

    sequencing_type = get_text(config, "sequencing_type").lower()

    for field in fields:
        value = get_text(config, field).lower()

        # For long reads only `annotation` is required; the steps of short reads can be left blank
        if sequencing_type == "long" and field != "annotation":
            if value not in {"", "yes", "no"}:
                errors.append(
                    f"Invalid value for `{field}`: '{config.get(field)}'. "
                    "Accepted values are: empty, 'yes' or 'no'."
                )
        elif value not in {"yes", "no"}:
            errors.append(
                f"Invalid value for `{field}`: '{config.get(field)}'. "
                "Accepted values are: 'yes' or 'no'."
            )

    return errors


def check_adapters_file(config):
    errors = []
    adapters = get_text(config, "adapters")

    if adapters:
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
    seed_format = get_text(config, "seed_format").lower()
    seed_file = get_text(config, "seed_file")
    sequencing_type = get_text(config, "sequencing_type").lower()
    search_ncbi = get_text(config, "search_ncbi").lower()

    if sequencing_type == "short":
        # Short-read seeds are only used by NOVOPlasty
        if not runs_novoplasty(config):
            return errors

        if seed_format == "fasta":
            if not seed_file:
                errors.append(
                    f"`seed_format` is set to '{seed_format}', but no `seed_file` is provided."
                )
            elif not os.path.isfile(seed_file):
                errors.append(
                    f"`seed_file` is set to '{seed_file}', but this file does not exist."
                )

            valid_exts = (".fasta", ".fa", ".fna")
            if seed_file and not seed_file.endswith(valid_exts):
                errors.append(
                    f"`seed_file` is set to '{seed_file}', but it doesn't end with a valid FASTA extension "
                    f"(expected one of: {', '.join(valid_exts)})."
                )

        elif seed_format == "genbank":
            feature = get_text(config, "feature")
            if feature not in ["CDS", "rRNA", "tRNA"]:
                errors.append(
                    f"`feature` is not set appropriate. Must be 'CDS', 'rRNA' or 'tRNA' (case sensitive)."
                )

            if not seed_file:
                errors.append(
                    f"`seed_format` is set to '{seed_format}', but no `seed_file` is provided."
                )
            elif not os.path.isfile(seed_file):
                errors.append(
                    f"`seed_file` is set to '{seed_file}', but this file does not exist."
                )

            valid_exts = (".gb", ".gbk", ".genbank")
            if seed_file and not seed_file.endswith(valid_exts):
                errors.append(
                    f"`seed_file` is set to '{seed_file}', but it doesn't end with a valid GENBANK extension "
                    f"(expected one of: {', '.join(valid_exts)})."
                )

        elif seed_format:
            errors.append(
                f"Invalid value for `seed_format`: '{config.get('seed_format')}'. "
                "Accepted values are: 'fasta', 'genbank' or empty."
            )

        elif search_ncbi != "yes":
            errors.append(
                "`seed_format` is empty, but `search_ncbi` is not 'yes'. "
                "Either provide a seed file or set `search_ncbi` to 'yes' to search NCBI for seed sequences."
            )

    elif sequencing_type == "long":
        if seed_format and seed_file:
            errors.append(
                f"`sequencing_type` is set to '{sequencing_type}', so `seed_format` and `seed_file` must be empty. Set `search_ncbi` to 'yes' to search NCBI for seed sequences."
            )

    return errors


def check_run_trimming_requirements(config):
    errors = []
    run_trimming = get_text(config, "run_trimming").lower()
    sequencing_type = get_text(config, "sequencing_type").lower()
    pacbio_adapters = get_text(config, "pacbio_adapters")

    if run_trimming == "yes":
        if sequencing_type == "short":
            if not get_text(config, "minlength"):
                errors.append("`run_trimming` is 'Yes', but `minlength` is not set.")
            if not get_text(config, "minquality"):
                errors.append("`run_trimming` is 'Yes', but `minquality` is not set.")

        elif sequencing_type == "long":
            if not pacbio_adapters:
                errors.append(
                    "`run_trimming` is 'Yes', but `pacbio_adapters` is not set."
                )

    return errors


def check_search_ncbi_requirements(config):
    errors = []
    search_ncbi = get_text(config, "search_ncbi").lower()

    # The NCBI search of short reads only provides seeds to NOVOPlasty. A blank value is "no".
    if runs_novoplasty(config):

        if search_ncbi not in {"", "yes", "no"}:
            errors.append(
                f"Invalid value for `search_ncbi`: '{config.get('search_ncbi')}'. "
                "Accepted values are: 'yes', 'no' or empty."
            )

        if search_ncbi == "yes":
            if not get_text(config, "search_genes"):
                errors.append("`search_ncbi` is 'Yes', but `search_genes` is not set.")

            if not get_text(config, "search_term"):
                errors.append("`search_ncbi` is 'Yes', but `search_term` is not set.")

            if not get_text(config, "max_references"):
                errors.append(
                    "`search_ncbi` is 'Yes', but `max_references` is not set."
                )

    return errors

def check_pacbio_adapters(config):
    errors = []
    seq_type = get_text(config, "sequencing_type").lower()
    run_trimming = get_text(config, "run_trimming").lower()
    adapter_value = get_text(config, "pacbio_adapters").strip()

    if seq_type == "long" and run_trimming == "yes":
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


def check_assembler_requirements(config):
    errors = []
    seq_type = get_text(config, "sequencing_type").lower()

    if seq_type == "short":
        run_novoplasty = get_text(config, "run_novoplasty").lower()
        run_getorganelle = get_text(config, "run_getorganelle").lower()

        if run_novoplasty != "yes" and run_getorganelle != "yes":
            errors.append("For short reads, at least one assembler (`run_novoplasty` or `run_getorganelle`) must be set to 'Yes'.")

        if run_novoplasty == "yes":
            max_memory = get_text(config, "max_memory")
            kmers = get_text(config, "kmers")
            if not get_text(config, "reads_length"):
                errors.append("`run_novoplasty` is 'Yes', but `reads_length` is not set.")
            if not get_text(config, "insert_size"):
                errors.append("`run_novoplasty` is 'Yes', but `insert_size` is not set.")
            if not re.match(r"^[1-9]\d*(\.0)?$", max_memory):
                errors.append(
                    f"`run_novoplasty` is 'Yes', but `max_memory` is '{max_memory}'. "
                    "It must be a whole number of GB (e.g., 40)."
                )
            # Without k-mers NOVOPlasty has no job to run, and the run would finish without assemblies
            if not re.fullmatch(r"[1-9]\d*(,[1-9]\d*)*", kmers):
                errors.append(
                    f"`run_novoplasty` is 'Yes', but `kmers` is '{kmers}'. "
                    "It must be a comma-separated list of k-mer sizes (e.g., 19,23,33)."
                )

        if run_getorganelle == "yes":
            database = get_text(config, "database")
            if not database:
                errors.append("`run_getorganelle` is 'Yes', but `database` is not set.")
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
        # Read as text, like create_snakemake_config.py, so empty cells are "" and numbers keep their form
        df = pd.read_csv(config_path, delimiter=delimiter, dtype=str, keep_default_na=False)
        all_configs = df.to_dict("records")

    elif config_path.endswith((".yaml", ".yml")):

        try:
            with open(config_path, "r") as f:
                config = yaml.safe_load(f)
                all_configs.append(config)

        except FileNotFoundError:
            logging.info(f"\n Error: File '{config_path}' not found.")
            sys.exit(1)

        except yaml.YAMLError as e:
            logging.info(f"\n Error: Failed to parse YAML file. Details:\n{e}")
            sys.exit(1)

        if not isinstance(config, dict):
            logging.info(
                "\n Error: The config file format is invalid (expected a top-level dictionary)."
            )
            sys.exit(1)

    # Blank fields can be written as `field:`, `field: ""` or an empty CSV cell, and all mean the same
    all_configs = [{key: as_text(value) for key, value in config.items()} for config in all_configs]

    for config in all_configs:
        sample = config.get("sample")
        logging.info(f"INFO: Running Check for Sample {sample}")

        missing_fields = check_required_fields(config, REQUIRED_FIELDS)
        spacing_errors = check_spacing_issues(config)

        if missing_fields:
            print_missing_fields(missing_fields)

        if spacing_errors:
            print_spacing_errors(spacing_errors)

        sample_name_errors = check_sample_names(config)

        if sample_name_errors:
            logging.info("\n Problem with `sample` field:")
            for msg in sample_name_errors:
                raise ValueError(f" - {msg}")

        read_layout_errors = check_read_layout(config)

        if read_layout_errors:
            logging.info("\n Problem with `read_layout` field:")
            for msg in read_layout_errors:
                raise ValueError(f" - {msg}")

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
            logging.info("\n Problem with `adapters` field:")
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
                raise ValueError(f" - {msg}")

        assembler_errors = check_assembler_requirements(config)
        if assembler_errors:
            logging.info("\n Problem with Assembler requirements:")
            for msg in assembler_errors:
                raise ValueError(f" - {msg}")


if __name__ == "__main__":
    main()
