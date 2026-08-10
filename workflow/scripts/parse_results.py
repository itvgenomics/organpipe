import pandas as pd
from Bio import SeqIO
import argparse
import os
import warnings
import glob
import openpyxl
import sqlite3
import logging
import sys
import yaml
import shutil
from Bio.SeqRecord import SeqRecord
from pathlib import Path


# ── Constants ──────────────────────────────────────────────────────────────────

HMMER_COLUMNS = [
    "target name",
    "accession_1",
    "query name",
    "accession_2",
    "E-value_full",
    "score_full",
    "bias_full",
    "E-value_best",
    "score_best",
    "bias_best",
    "exp",
    "reg",
    "clu",
    "ov",
    "env",
    "dom",
    "rep",
    "inc",
    "description of target",
]

# (relative-src-template, relative-dst-template) pairs for MITOS2 file copying.
# Placeholders are formatted with assembly=<assembly>.
MITOS2_FILES = [
    ("ignored.mitos", "{assembly}.ignored.mitos"),
    ("result.bed", "{assembly}.result.bed"),
    ("result.faa", "{assembly}.result.faa"),
    ("result.fas", "{assembly}.result.fas"),
    ("result.geneorder", "{assembly}.result.geneorder"),
    ("result.gff", "{assembly}.result.gff"),
    ("result.mitos", "{assembly}.result.mitos"),
    ("result.seq", "{assembly}.result.seq"),
    ("stst.dat", "{assembly}.stst.dat"),
    (
        "mitfi-global/sequence.fas-0_tRNAout.nc",
        "{assembly}_tRNAout.nc",
    ),
    (
        "mitfi-global/sequence.fas-0_rRNAout.nc",
        "{assembly}_rRNAout.nc",
    ),
]


# ── CLI ────────────────────────────────────────────────────────────────────────


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--sample",
        help="Sample ID to parse the data.",
        required=True,
    )
    args = parser.parse_args()
    return args


def get_sample_keys(yaml_file):
    with open(yaml_file, "r") as file:
        config = yaml.safe_load(file)  # Load the YAML file
        if "samples" in config:
            return config["samples"].keys()  # Get keys from the 'samples' section
        else:
            return None  # or raise an error if 'samples' is not found


# ── Path context helpers ───────────────────────────────────────────────────────


def get_path_context(logfile, config=None):
    """Extract sample, assembler, seed, kmer, and assembly from a
    ``results/{sample}/{tool}/{assembler}/[{folder}/]{assembly}/...`` path.

    Folder naming convention (novoplasty): ``{seed}_kmer{kmer}``.

    Parameters
    ----------
    logfile : str | Path
    config : dict | None
        Snakemake config.  Required for getorganelle seed/kmer lookup.

    Returns
    -------
    tuple[str, str, str, str, str]
        ``(sample, assembler, seed, kmer, assembly)``
    """
    parts = Path(logfile).parts
    sample = parts[1]
    assembler = parts[3]
    assembly = Path(logfile).parent.name

    if assembler == "novoplasty":
        folder = parts[4]
        seed = folder.split("_kmer")[0]
        kmer = folder.split("_kmer")[1]
    elif assembler == "getorganelle" and config is not None:
        seed = config["samples"][sample]["database"].lower()
        kmer = config["samples"][sample]["spades_kmers"].lower()
    else:
        seed = kmer = ""

    return sample, assembler, seed, kmer, assembly


def get_novoplasty_path_context(logfile):
    """Extract sample, seed, kmer from a
    ``results/{sample}/novoplasty/{seed}/kmer{kmer}/...`` path.

    Returns
    -------
    tuple[str, str, str]
        ``(sample, seed, kmer)``
    """
    parts = Path(logfile).parts
    sample = parts[1]
    seed = parts[3]
    kmer = parts[4].split("kmer")[1]
    return sample, seed, kmer


def get_pilon_path_context(logfile, config):
    """Extract sample, assembler, seed, kmer, and assembly from a pilon log path.

    Pilon's novoplasty assembly folder is named ``{sample}_{seed}_{kmer}_...``,
    which is a different convention from other tools.

    Returns
    -------
    tuple[str, str, str, str, str]
        ``(sample, assembler, seed, kmer, assembly)``
    """
    parts = Path(logfile).parts
    sample = parts[1]
    assembler = parts[3]
    assembly = parts[4]  # directory containing pilon.log

    if assembler == "novoplasty":
        seed = "_".join(parts[4].split(sample)[1].split("_")[1:-2])
        kmer = parts[4].split("_")[-2]
    elif assembler == "getorganelle":
        seed = config["samples"][sample]["database"].lower()
        kmer = config["samples"][sample]["spades_kmers"].lower()
    else:
        seed = kmer = ""

    return sample, assembler, seed, kmer, assembly


# ── Utilities ──────────────────────────────────────────────────────────────────


def concat_csv_files(directory):
    """Read all CSV files in *directory* and return a concatenated DataFrame.

    Returns ``None`` if no CSV files are found or if concatenation fails.
    """
    csv_files = list(Path(directory).glob("*.csv"))
    if not csv_files:
        logging.warning(f"No CSV files found in: {directory}")
        return None
    try:
        return pd.concat(
            [pd.read_csv(f) for f in csv_files], ignore_index=True
        )
    except Exception as e:
        logging.error(f"Error concatenating CSVs in {directory}: {e}")
        return None


def convert_genbank_to_fasta(genbank_file, fasta_file):
    try:
        with open(fasta_file, "w") as output_handle:
            # Parse the GenBank file and write to FASTA
            for record in SeqIO.parse(genbank_file, "genbank"):
                SeqIO.write(record, output_handle, "fasta")
        logging.info(f"Successfully converted {genbank_file} to {fasta_file}")
    except Exception as e:
        logging.error(f"An error occurred: {e}")


def ignore_symlinks(src, names):
    return [name for name in names if os.path.islink(os.path.join(src, name))]


def extract_features(genbank_file):
    all_data = []
    seen_sequences = set()  # Set to store sequences that have already been added

    for record in SeqIO.parse(genbank_file, "genbank"):
        for feature in record.features:
            if feature.type != "source":
                data = {}
                acronym = feature.qualifiers.get("gene", [""])[0]

                if not acronym:
                    acronym = feature.qualifiers.get("product", [""])[0]

                sequence = str(
                    feature.extract(record.seq)
                )  # Extract sequence as string

                # Add the data only if the sequence is unique
                if sequence not in seen_sequences:
                    data["feature"] = str(acronym).strip()
                    data["seq"] = sequence
                    all_data.append(data)

                    # Mark this sequence as seen
                    seen_sequences.add(sequence)

    return all_data


# ── NOVOPlasty ─────────────────────────────────────────────────────────────────


def parse_novoplasty(logfile):
    logging.info(f"Parsing NOVOPlasty log file: {logfile}")

    sample, seed, kmer = get_novoplasty_path_context(logfile)

    # Create an empty dictionary to store the parsed key-value pairs
    parsed_data = {}
    contigs = {}

    with open(logfile, "r") as file:
        for line in file:
            # Skip lines that are dashes or empty
            if line.strip() == "" or line.startswith("-"):
                continue

            # Split the line into key and value
            if ":" in line:
                key, value = line.split(":", 1)
                key = key.strip()
                value = value.strip()

                # Process the value if it's a percentage or contains other symbols
                if "%" in value:
                    value = str(value)
                elif value.isdigit():
                    value = int(value)
                elif "bp" in value:  # Keep 'bp' values as-is
                    pass
                else:
                    try:
                        value = float(value)
                    except ValueError:
                        pass  # Keep as string if it can't be converted

                # Handle contig keys
                if key.startswith("Contig"):
                    contigs[key] = value
                else:
                    parsed_data[key] = value

    # Combine contigs into a single key if any exist
    if contigs:
        parsed_data["Contigs"] = contigs

    parsed_data["Circularized"] = "No"
    for file in os.listdir(os.path.dirname(logfile)):
        if "Circularized_" in file or "Option_" in file:
            parsed_data["Circularized"] = "Yes"
            break

    parsed_data["Sample"] = sample
    parsed_data["Seed"] = seed
    parsed_data["kmer"] = kmer

    df = pd.DataFrame.from_dict(parsed_data, orient="index").T

    columns = [
        "Sample",
        "Seed",
        "kmer",
        "Total contigs",
        "Largest contig",
        "Smallest contig",
        "Average insert size",
        "Total reads",
        "Aligned reads",
        "Assembled reads",
        "Organelle genome %",
        "Average organelle coverage",
        "Contigs",
        "Circularized",
    ]

    df = df.reindex(columns=columns, fill_value="")
    for char in ["{", "}", "'"]:
        df["Contigs"] = df["Contigs"].apply(lambda x: str(x).replace(char, ""))

    df = df.fillna("")

    out_dir = Path(f"workflow/reports/{sample}/novoplasty")
    out_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_dir / f"{seed}_kmer{kmer}.csv", index=False)


def get_novoplasty_files(logfile):
    sample, seed, kmer = get_novoplasty_path_context(logfile)

    logging.info(
        f"Getting NOVOPlasty files for sample: {sample}, seed: {seed}, kmer: {kmer}"
    )

    novoplasty_outdir = Path(f"workflow/reports/{sample}/files/{seed}/kmer{kmer}")
    novoplasty_outdir.mkdir(parents=True, exist_ok=True)

    base = Path(f"results/{sample}/novoplasty/{seed}/kmer{kmer}")
    shutil.copy(
        base / f"Assembled_reads_{sample}_R1.fasta",
        novoplasty_outdir / "Assembled_reads_R1.fasta",
    )
    shutil.copy(
        base / f"Assembled_reads_{sample}_R2.fasta",
        novoplasty_outdir / "Assembled_reads_R2.fasta",
    )
    shutil.copy(
        base / f"log_{sample}.txt",
        novoplasty_outdir / "log_novoplasty.txt",
    )


# ── MITOS2 ─────────────────────────────────────────────────────────────────────


def parse_mitos2(logfile, config):
    logging.info(f"Parsing MITOS2 files: {logfile}")

    sample, assembler, seed, kmer, assembly = get_path_context(logfile, config)

    out_dir = Path(f"workflow/reports/{sample}/mitos2")
    out_dir.mkdir(parents=True, exist_ok=True)

    mitos_df = pd.DataFrame(
        columns=[
            "sample",
            "seed",
            "kmer",
            "assembly",
            "gene_order",
            "missing",
            "duplicated",
            "forward_genes",
            "forward_count",
            "reverse_genes",
            "reverse_count",
            "transporter_genes",
            "transporter_count",
            "ribosomal_genes",
            "ribosomal_count",
            "origins",
            "origins_count",
            "coding_genes",
            "coding_count",
        ]
    )

    mitos_df["sample"] = [sample]
    mitos_df["seed"] = [seed]
    mitos_df["kmer"] = [kmer]
    mitos_df["assembly"] = [assembly]

    # Parse missing / duplicated entries from the log file
    targets = {"missing": "", "duplicated": ""}
    with open(logfile, "r") as file:
        for line in file:
            for key in targets:
                if key in line and ":" in line:
                    _, value = line.split(":", 1)
                    targets[key] = value.strip()

    for key, value in targets.items():
        mitos_df[key] = value

    # Get gene order file
    assembly_dir = Path(logfile).parent
    geneorder_file = assembly_dir / "result.geneorder"

    if geneorder_file.exists():
        with open(geneorder_file, "r") as file:
            transporters = []
            ribosomais = []
            codings = []
            forward = []
            reverse = []
            origins = []

            data = file.read()
            gene_line = data.split("\n")[1]
            mitos_df["gene_order"] = gene_line

            # Get forward and reverse genes
            for cds in gene_line.split():
                if cds.startswith("-"):
                    reverse.append(cds)
                else:
                    forward.append(cds)

            # Get type of gene
            for cds in gene_line.split():
                if cds.startswith(("t", "-t")):
                    transporters.append(cds)
                elif cds.startswith(("r", "-r")):
                    ribosomais.append(cds)
                elif cds.startswith(("OH", "-OH", "OL", "-OL")):
                    origins.append(cds)
                else:
                    codings.append(cds)

        mitos_df["forward_genes"] = [forward]
        mitos_df["forward_count"] = len(forward)
        mitos_df["reverse_genes"] = [reverse]
        mitos_df["reverse_count"] = len(reverse)
        mitos_df["transporter_genes"] = [transporters]
        mitos_df["transporter_count"] = len(transporters)
        mitos_df["ribosomal_genes"] = [ribosomais]

        unique_rRNA = {"rrnL": 0, "rrnS": 0}  # Initialize counts for rrnL and rrnS
        for item in ribosomais:
            if "rrnL" in item:
                unique_rRNA["rrnL"] = 1
            if "rrnS" in item:
                unique_rRNA["rrnS"] = 1
        mitos_df["ribosomal_count"] = sum(unique_rRNA.values())

        mitos_df["origins"] = [origins]

        unique_origins = {"OH": 0, "OL": 0}  # Initialize counts for OH and OL
        for item in origins:
            if "OH" in item:
                unique_origins["OH"] = 1
            if "OL" in item:
                unique_origins["OL"] = 1
        mitos_df["origins_count"] = sum(unique_origins.values())

        mitos_df["coding_genes"] = [codings]
        mitos_df["coding_count"] = len(codings)

        # Convert the lists to strings and remove unwanted characters
        columns_to_correct = [
            "ribosomal_genes",
            "origins",
            "coding_genes",
            "transporter_genes",
            "reverse_genes",
            "forward_genes",
        ]

        for column in columns_to_correct:
            mitos_df[column] = (
                mitos_df[column]
                .apply(str)
                .str.replace("[", "", regex=False)
                .str.replace("]", "", regex=False)
                .str.replace("'", "", regex=False)
            )

    mitos_df.to_csv(out_dir / f"{assembly}.csv", index=False)


def get_mitos2_files(logfile, config):
    sample, assembler, seed, kmer, assembly = get_path_context(logfile, config)

    if assembler == "novoplasty":
        logging.info(
            f"Getting MITOS2 files for sample: {sample}, seed: {seed}, kmer: {kmer}, assembly: {assembly}"
        )
        src_base = Path(
            f"results/{sample}/mitos2/novoplasty/{seed}_kmer{kmer}/{assembly}"
        )
        dst_base = Path(f"workflow/reports/{sample}/files/{seed}/kmer{kmer}")

    elif assembler == "getorganelle":
        logging.info(
            f"Getting MITOS2 files for sample: {sample}, assembly: {assembly}"
        )
        src_base = Path(f"results/{sample}/mitos2/getorganelle/{assembly}")
        dst_base = Path(f"workflow/reports/{sample}/files/{assembly}")

    else:
        return

    dst_base.mkdir(parents=True, exist_ok=True)

    for src_template, dst_template in MITOS2_FILES:
        src = src_base / src_template.format(assembly=assembly)
        dst = dst_base / dst_template.format(assembly=assembly)
        shutil.copy(src, dst)


# ── Pilon ──────────────────────────────────────────────────────────────────────


def parse_pilon(logfile, config):
    logging.info(f"Parsing Pilon data for: {logfile}")

    sample, assembler, seed, kmer, assembly = get_pilon_path_context(logfile, config)

    pilon_df = pd.DataFrame(
        columns=[
            "sample",
            "seed",
            "kmer",
            "assembly",
            "genome_size",
            "reads",
            "filtered",
            "mapped",
            "propper",
            "stray",
            "fr",
            "insert_max",
            "coverage",
            "minDepth",
            "confirmed",
            "corrected",
            "mean_frags_coverage",
            "mean_total_coverage",
            "changes_number",
            "changes",
        ]
    )

    pilon_df["sample"] = [sample]
    pilon_df["seed"] = [seed]
    pilon_df["kmer"] = [kmer]
    pilon_df["assembly"] = [assembly]

    # Initialize all fields to avoid referencing unassigned variables
    genome_size = reads = filtered = mapped = propper = stray = ""
    fr = insert_max = coverage = minDepth = confirmed = corrected = ""
    mean_frags_coverage = mean_total_coverage = ""

    with open(logfile, "r") as file:
        for line in file:
            if "Input genome size" in line:
                genome_size = line.split(":")[1].replace("\n", "").strip()
            elif line.startswith(f"results/{sample}/pilon/{assembler}"):
                data = line.split(".bam:")
                if len(data) == 7:
                    reads = str(data).split(",")[1].strip().split(" ")[1]
                    filtered = str(data).split(",")[2].strip().split(" ")[0]
                    mapped = str(data).split(",")[3].strip().split(" ")[0]
                    propper = str(data).split(",")[4].strip().split(" ")[0]
                    stray = str(data).split(",")[5].strip().split(" ")[0]
                    fr = str(data).split(",")[6].strip().split(" ")[1:]
                    insert_max = str(data).split(",")[7].strip().split(" ")[1]
                else:
                    reads = str(data).split(",")[1].strip().split(" ")[1]
                    filtered = str(data).split(",")[2].strip().split(" ")[0]
                    mapped = str(data).split(",")[3].strip().split(" ")[0]
                    propper = str(data).split(",")[4].strip().split(" ")[0]
                    stray = str(data).split(",")[5].strip().split(" ")[0]
                    insert_max = str(data).split(",")[6].strip().split(" ")[1]
            elif line.startswith("Total Reads:"):
                coverage = line.split(",")[1].split(":")[1].strip()
                minDepth = line.split(",")[2].split(":")[1].strip()
            elif line.startswith("Confirmed"):
                confirmed = line.split(" ", maxsplit=1)[1].replace("\n", "")
            elif line.startswith("Corrected"):
                corrected = line.split(" ", maxsplit=1)[1].replace("\n", "")
            elif line.startswith("Mean frags coverage"):
                mean_frags_coverage = line.split(": ")[1].replace("\n", "")
            elif line.startswith("Mean total coverage"):
                mean_total_coverage = line.split(": ")[1].replace("\n", "")

    pilon_df["genome_size"] = genome_size
    pilon_df["reads"] = reads
    pilon_df["filtered"] = filtered
    pilon_df["mapped"] = mapped
    pilon_df["propper"] = propper
    pilon_df["stray"] = stray
    pilon_df["fr"] = [fr]
    pilon_df["insert_max"] = insert_max
    pilon_df["coverage"] = coverage
    pilon_df["minDepth"] = minDepth
    pilon_df["confirmed"] = confirmed
    pilon_df["corrected"] = corrected
    pilon_df["mean_frags_coverage"] = mean_frags_coverage
    pilon_df["mean_total_coverage"] = mean_total_coverage

    # Get number of changes
    changes = []
    changes_file = Path(
        f"results/{sample}/pilon/{assembler}/{assembly}/{assembly}.changes"
    )
    with open(changes_file, "r") as file:
        lines = file.readlines()  # Read all lines at once
        line_count = len(lines)  # Count the lines
        if line_count > 0:
            for line in lines:  # Iterate over the list of lines
                change_parts = line.strip().split()
                original_coord = change_parts[0]
                new_coord = change_parts[1]
                original_seq = change_parts[2]
                new_seq = change_parts[3]

                string = f"original_coord: {original_coord} new_coord: {new_coord} original_seq: {original_seq} new_seq: {new_seq}"
                changes.append(string)

    pilon_df["changes_number"] = line_count
    pilon_df["changes"] = str(changes)

    out_dir = Path(f"workflow/reports/{sample}/pilon")
    out_dir.mkdir(parents=True, exist_ok=True)
    pilon_df.to_csv(out_dir / f"{assembly}.csv", index=False)


def get_pilon_files(logfile, config):
    sample, assembler, seed, kmer, assembly = get_pilon_path_context(logfile, config)

    src_base = Path(f"results/{sample}/pilon/{assembler}/{assembly}")

    if assembler == "novoplasty":
        logging.info(
            f"Getting Pilon files for sample: {sample}, seed: {seed}, kmer: {kmer}, assembly: {assembly}"
        )
        dst_base = Path(f"workflow/reports/{sample}/files/{seed}/kmer{kmer}")

    elif assembler == "getorganelle":
        logging.info(
            f"Getting Pilon files for sample: {sample}, assembly: {assembly}"
        )
        dst_base = Path(f"workflow/reports/{sample}/files/{assembly}")

    else:
        return

    dst_base.mkdir(parents=True, exist_ok=True)
    shutil.copy(src_base / f"{assembly}.changes", dst_base / f"{assembly}_pilon.changes")
    shutil.copy(src_base / "pilon.log", dst_base / f"{assembly}_pilon.log")


# ── CPGAVAS2 ───────────────────────────────────────────────────────────────────


def parse_cpgavas2_report_table(logfile, config):
    logging.info(f"Parsing CPGAVAS2 Report Table: {logfile}")

    sample, assembler, seed, kmer, assembly = get_path_context(logfile, config)

    output_dirs = [
        Path(f"workflow/reports/{sample}/cpgavas2/{assembler}/report_table/gene_composition"),
        Path(f"workflow/reports/{sample}/cpgavas2/{assembler}/report_table/intron_exon"),
        Path(f"workflow/reports/{sample}/cpgavas2/{assembler}/report_table/codon_usage"),
    ]
    for d in output_dirs:
        d.mkdir(parents=True, exist_ok=True)

    # Read the text and split it into sections
    with open(logfile, "r") as file:
        text = file.read()

    sections = text.split(
        "__________________________________________________________________________________"
    )

    # Parse Table 1: Gene Composition
    try:
        table1_data = [line.split("\t") for line in sections[2].strip().split("\n")]
        df1 = pd.DataFrame(
            table1_data,
            columns=["Category of genes", "Group of genes", "Name of genes"],
        )
        df1.dropna(inplace=True)
        df1.insert(0, "Kmer", kmer)
        df1.insert(0, "Seed", seed)
        df1.insert(0, "Assembly", assembly)
        df1.insert(0, "Sample", sample)
        df1.to_csv(output_dirs[0] / f"{assembly}_gene_composition.csv", index=False)
    except Exception as e:
        logging.error(f"Warning while parsing Table: {e}")

    # Parse Table 2: Lengths of Introns and Exons
    try:
        table2_data = [line.split("\t") for line in sections[5].strip().split("\n")]
        df2 = pd.DataFrame(
            table2_data,
            columns=[
                "Gene", "Strand", "Start", "End",
                "ExonI", "IntronI", "ExonII", "IntronII", "ExonIII",
            ],
        )
        df2.insert(0, "Kmer", kmer)
        df2.insert(0, "Seed", seed)
        df2.insert(0, "Assembly", assembly)
        df2.insert(0, "Sample", sample)
        df2.to_csv(output_dirs[1] / f"{assembly}_intron_exon.csv", index=False)
    except Exception as e:
        logging.error(f"Warning while parsing Table: {e}")

    # Parse Table 3: Codon Usage
    try:
        table3_data = [line.split("\t") for line in sections[7].strip().split("\n")[1:]]
        df3 = pd.DataFrame(
            table3_data, columns=["Codon", "Amino acid", "Frequency", "Number"]
        )
        df3.insert(0, "Kmer", kmer)
        df3.insert(0, "Seed", seed)
        df3.insert(0, "Assembly", assembly)
        df3.insert(0, "Sample", sample)
        df3.to_csv(output_dirs[2] / f"{assembly}_codon_usage.csv", index=False)
    except Exception as e:
        logging.error(f"Warning while parsing Table: {e}")


def parse_cpgavas2_problems(logfile, config):
    logging.info(f"Parsing Problems: {logfile}")

    sample, assembler, seed, kmer, assembly = get_path_context(logfile, config)

    out_dir = Path(f"workflow/reports/{sample}/cpgavas2/{assembler}/problems")
    out_dir.mkdir(parents=True, exist_ok=True)

    with open(logfile, "r") as file:
        lines = file.readlines()
        text = [
            line
            for line in lines
            if not line.startswith("#") and not line.startswith("Possible")
        ]

        lines = "".join(text).split("\n")[:-1]
        data = [line.split("\t", maxsplit=1) for line in lines]

        df = pd.DataFrame(data, columns=["identifier", "problem"])
        df["problem"] = df["problem"].apply(lambda x: str(x).replace("\t", " "))

        df.insert(0, "Kmer", kmer)
        df.insert(0, "Seed", seed)
        df.insert(0, "Assembly", assembly)
        df.insert(0, "Sample", sample)

        df.to_csv(out_dir / f"{assembly}_problems.csv", index=False)


def get_cpgavas2_files(logfile, config):
    sample, assembler, seed, kmer, assembly = get_path_context(logfile, config)

    if assembler == "novoplasty":
        logging.info(
            f"Copying CPGAVAS2 files for sample: {sample}, seed: {seed}, kmer: {kmer}, assembly: {assembly}"
        )
        source_path = Path(f"results/{sample}/cpgavas2/novoplasty/{seed}_kmer{kmer}/{assembly}")
        destination_path = Path(f"workflow/reports/{sample}/files/cpgavas2/novoplasty/{assembly}")

    elif assembler == "getorganelle":
        logging.info(
            f"Copying CPGAVAS2 files for sample: {sample}, assembly: {assembly}"
        )
        source_path = Path(f"results/{sample}/cpgavas2/getorganelle/{assembly}")
        destination_path = Path(f"workflow/reports/{sample}/files/cpgavas2/getorganelle/{assembly}")

    else:
        return

    if not destination_path.exists():
        shutil.copytree(source_path, destination_path, ignore=ignore_symlinks)


# ── nhmmer ─────────────────────────────────────────────────────────────────────


def parse_nhmmer(logfile, subdir, config=None):
    """Parse an nhmmer tblout file and write a CSV to the reports tree.

    This single function replaces the previous four separate functions:
    ``parse_ncRNA_nhmmer``, ``parse_intergenes_nhmmer``,
    ``parse_ncRNA_nhmmer_long``, and ``parse_intergenes_nhmmer_long``.

    Parameters
    ----------
    logfile : str | Path
        Path to the tblout file.
    subdir : str
        Output sub-directory under ``workflow/reports/{sample}/nhmmer/``.
        Typically ``'ncRNA'`` or ``'intergenes'``.
    config : dict | None
        Snakemake config dict.  When provided the function uses the
        short-read path structure (assembler at ``parts[3]``).
        When ``None`` it uses the long-read structure (seed at ``parts[3]``).
    """
    parts = Path(logfile).parts
    sample = parts[1]

    out_dir = Path(f"workflow/reports/{sample}/nhmmer/{subdir}")
    out_dir.mkdir(parents=True, exist_ok=True)

    if config is not None:
        # Short reads: results/{sample}/nhmmer/{assembler}/{folder}/{assembly}/file
        assembler = parts[3]
        assembly = Path(logfile).parent.name
        if assembler == "novoplasty":
            folder = parts[4]
            seed = folder.split("_kmer")[0]
            kmer = folder.split("_kmer")[1]
        elif assembler == "getorganelle":
            seed = config["samples"][sample]["database"].lower()
            kmer = config["samples"][sample]["spades_kmers"].lower()
        else:
            seed = kmer = ""
        out_name = assembly
    else:
        # Long reads: results/{sample}/nhmmer/{seed}/file
        seed = parts[3]
        kmer = assembly = None
        out_name = seed

    logging.info(f"Parsing {subdir} NHMMER data: {logfile}")

    try:
        df = pd.read_csv(logfile, sep=r"\s+", comment="#", header=None)
        df.columns = HMMER_COLUMNS

        if config is not None:
            df.insert(0, "kmer", kmer)
            df.insert(0, "seed", seed)
            df.insert(0, "assembly", assembly)
            df.insert(0, "sample", sample)
        else:
            df.insert(0, "seed", seed)
            df.insert(0, "sample", sample)

        df.to_csv(out_dir / f"{out_name}.csv", index=False)

    except Exception as e:
        logging.error(f"Warning while parsing Table: {e}")


# ── MitoHiFi ───────────────────────────────────────────────────────────────────


def parse_mitohifi(statsfile):
    logging.info(f"Parsing MitoHifi for file: {statsfile}")

    parts = Path(statsfile).parts
    sample = parts[1]
    seed = parts[3]

    df = pd.read_csv(statsfile, sep="\t", comment="#")
    df.insert(0, "seed", seed)
    df.insert(0, "sample", sample)

    out_dir = Path(f"workflow/reports/{sample}/mitohifi")
    out_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_dir / f"{seed}.csv", index=False)


# ── Pipeline section helpers ───────────────────────────────────────────────────


def copy_images(sample, config):
    """Copy PNG images (depth plots, recruitment plots, OGDraw) into the reports tree."""
    logging.info(f"Getting all images from sample: {sample}")

    assemblers = []
    if config["samples"][sample]["run_novoplasty"].lower() == "yes":
        assemblers.append("novoplasty")
    if config["samples"][sample]["run_getorganelle"].lower() == "yes":
        assemblers.append("getorganelle")

    for assembler in assemblers:
        for root, dirs, files in os.walk(f"results/{sample}/images/{assembler}"):
            for file in files:
                if not file.endswith(".png"):
                    continue

                if ".depth." in file:
                    dest_dir = Path(f"workflow/reports/{sample}/images/depth")
                elif "_recruitment_plot." in file:
                    dest_dir = Path(f"workflow/reports/{sample}/images/recruitment_plot")
                else:
                    dest_dir = Path(f"workflow/reports/{sample}/images/ogdraw")

                dest_dir.mkdir(parents=True, exist_ok=True)
                logging.info(f"Copying file: {file}")
                shutil.copy(os.path.join(root, file), dest_dir)


def copy_genbanks_novoplasty(sample, config):
    """Copy GenBank files produced by the NOVOPlasty assembly path."""
    logging.info(f"Getting GenBank file for sample: {sample}")

    seeds = list(config["samples"][sample]["seeds"])

    gb_files = []
    for root, dirs, files in os.walk(f"results/{sample}/genbanks/novoplasty"):
        for file in files:
            if file.endswith(".gb"):
                gb_files.append(file)

    for seed in seeds:
        for file in gb_files:
            if f"_{seed}_" in file:
                logging.info(f"Copying GenBank file: {file}")
                dest_dir = Path(f"workflow/reports/{sample}/genbanks/{seed}")
                dest_dir.mkdir(parents=True, exist_ok=True)
                shutil.copy(f"results/{sample}/genbanks/novoplasty/{file}", dest_dir)

        Path(f"workflow/reports/{sample}/fastas/{seed}").mkdir(parents=True, exist_ok=True)


def copy_genbanks_getorganelle(sample, config):
    """Copy GenBank files produced by the GetOrganelle assembly path."""
    gb_files = []
    for root, dirs, files in os.walk(f"results/{sample}/genbanks/getorganelle"):
        for file in files:
            if file.endswith(".gb"):
                gb_files.append(file)

    organelle = config["samples"][sample]["organelle"].lower()
    if organelle == "mito":
        base_dir = Path(f"results/{sample}/mitos2/getorganelle")
    elif organelle == "chloro":
        base_dir = Path(f"results/{sample}/cpgavas2/getorganelle")

    assemblies = sorted(d.name for d in base_dir.iterdir() if d.is_dir())

    logging.info(f"Base dir: {base_dir}, assemblies: {assemblies}")
    for assembly in assemblies:
        logging.info(f"Processing assembly: {assembly}")
        for file in gb_files:
            logging.info(f"Checking file: {file} against assembly: {assembly}")
            if assembly in str(file):
                logging.info(f"Copying GenBank file: {file}")
                dest_dir = Path(f"workflow/reports/{sample}/genbanks/{assembly}")
                dest_dir.mkdir(parents=True, exist_ok=True)
                shutil.copy(f"results/{sample}/genbanks/getorganelle/{file}", dest_dir)

        Path(f"workflow/reports/{sample}/fastas/{assembly}").mkdir(parents=True, exist_ok=True)


def convert_genbanks_to_fastas(sample):
    """Convert all GenBank files in the reports tree to FASTA."""
    logging.info(f"Getting FASTA file for sample: {sample}")

    genbanks_dir = Path(f"workflow/reports/{sample}/genbanks")
    if not genbanks_dir.exists():
        return

    for genbank_subdir in genbanks_dir.iterdir():
        if not genbank_subdir.is_dir():
            continue
        for gb_file in genbank_subdir.iterdir():
            logging.info(f"Getting FASTA file for sample: {gb_file.name}")
            fasta_file = (
                str(gb_file)
                .replace("/genbanks/", "/fastas/")
                .replace(".gb", ".fasta")
            )
            convert_genbank_to_fasta(gb_file, fasta_file)


def copy_assembly_fastas(sample, config):
    """Copy raw assembly FASTA files when annotation is disabled.

    Used as a fallback when ``annotation`` is set to ``"No"`` and no GenBank
    files are produced.  Copies from:

    * NOVOPlasty: ``results/{sample}/assemblies/novoplasty/``
    * GetOrganelle: ``results/{sample}/assemblies/getorganelle/*.fasta``

    All files land in ``workflow/reports/{sample}/fastas/``.
    """
    sample_cfg = config["samples"][sample]
    dest_dir = Path(f"workflow/reports/{sample}/fastas")
    dest_dir.mkdir(parents=True, exist_ok=True)

    if sample_cfg["run_novoplasty"].lower() == "yes":
        src_dir = Path(f"results/{sample}/assemblies/novoplasty")
        if src_dir.exists():
            for fasta_file in src_dir.iterdir():
                if fasta_file.is_file():
                    logging.info(f"Copying NOVOPlasty assembly: {fasta_file.name}")
                    shutil.copy(fasta_file, dest_dir / fasta_file.name)
        else:
            logging.warning(f"NOVOPlasty assembly directory not found: {src_dir}")

    if sample_cfg["run_getorganelle"].lower() == "yes":
        src_dir = Path(f"results/{sample}/assemblies/getorganelle")
        if src_dir.exists():
            for fasta_file in src_dir.glob("*.fasta"):
                logging.info(f"Copying GetOrganelle assembly: {fasta_file.name}")
                shutil.copy(fasta_file, dest_dir / fasta_file.name)
        else:
            logging.warning(f"GetOrganelle assembly directory not found: {src_dir}")


def extract_genes_from_genbanks(sample, config):
    """Extract per-gene FASTA files from all GenBank files in the reports tree.

    Applies to both short-read (mito only, skips rotated) and long-read samples.
    Removes a stale ``genes/`` directory before writing fresh output.
    """
    genbanks_dir = Path(f"workflow/reports/{sample}/genbanks")
    if not genbanks_dir.exists():
        return

    seq_type = config["samples"][sample]["sequencing_type"].lower()
    organelle = config["samples"][sample]["organelle"].lower()
    genes_dir = Path(f"workflow/reports/{sample}/genes")

    # Start fresh – remove any stale gene files from a previous run
    if genes_dir.exists():
        shutil.rmtree(genes_dir)

    for root, dirs, files in os.walk(genbanks_dir):
        for file in files:
            # For short reads, only mito genbanks are processed; rotated files are skipped
            if seq_type == "short":
                if organelle != "mito" or ".rotated.gb" in file:
                    continue

            filepath = os.path.join(root, file)
            logging.info(f"Getting genes from GenBank file: {filepath}")

            all_data = extract_features(filepath)
            assembly = file.replace(".gb", "")

            for data in all_data:
                if seq_type == "short":
                    feature = str(data["feature"]).split("(")[0].strip()
                else:
                    feature = data["feature"]
                seq = data["seq"]

                genes_dir.mkdir(parents=True, exist_ok=True)
                gene_fasta = genes_dir / f"{feature}.fasta"
                with open(gene_fasta, "a" if gene_fasta.exists() else "w") as output_fasta:
                    output_fasta.write(f">{feature}_{assembly}\n{seq}\n")


def generate_summary_mito(sample, config):
    """Write ``summary.csv`` for mitochondrial short-read samples."""
    pilon_csv = Path(f"workflow/reports/{sample}/pilon.csv")
    mitos2_csv = Path(f"workflow/reports/{sample}/mitos2.csv")

    if not (pilon_csv.exists() and mitos2_csv.exists()):
        return

    logging.info(f"Writing summary.csv for sample: {sample}")

    df_mitos2 = pd.read_csv(mitos2_csv)[
        ["sample", "seed", "kmer", "assembly",
         "transporter_count", "ribosomal_count", "origins_count", "coding_count"]
    ]
    df_pilon = pd.read_csv(pilon_csv)[["assembly", "genome_size", "changes_number"]]

    df_summary = pd.merge(df_mitos2, df_pilon, on="assembly")
    df_summary.to_csv(f"workflow/reports/{sample}/summary.csv", index=False)


def generate_summary_chloro(sample, config):
    """Write ``summary.csv`` for chloroplast short-read samples."""
    base_dir = Path(f"workflow/reports/{sample}/genbanks/")
    pilon_csv = Path(f"workflow/reports/{sample}/pilon.csv")

    if not (pilon_csv.exists() and base_dir.exists()):
        return

    logging.info(f"Writing summary.csv for sample: {sample}")

    rows = []
    for root, dirs, files in os.walk(base_dir):
        seed = os.path.basename(root)
        for file in files:
            if not file.endswith(("chloe.gb", "cpgavas2.gb")):
                continue

            filepath = os.path.join(root, file)
            assembly = os.path.splitext(file)[0]

            try:
                kmer = int(assembly.split("_")[-2])
            except Exception:
                seed = config["samples"][sample]["database"]
                kmer = config["samples"][sample]["spades_kmers"]

            transporter_count = ribosomal_count = coding_count = genome_size = 0
            for record in SeqIO.parse(filepath, "genbank"):
                genome_size += len(record.seq)
                for feature in record.features:
                    if feature.type == "CDS":
                        coding_count += 1
                    elif feature.type == "rRNA":
                        ribosomal_count += 1
                    elif feature.type == "tRNA":
                        transporter_count += 1

            rows.append({
                "Sample": sample,
                "Seed": seed,
                "kmer": kmer,
                "assembly": assembly,
                "annotation_tool": "chloe" if "chloe" in file else "cpgavas2",
                "transporter_count": transporter_count,
                "ribosomal_count": ribosomal_count,
                "coding_count": coding_count,
                "genome_size": genome_size,
            })

    if rows:
        pd.DataFrame(rows).to_csv(f"workflow/reports/{sample}/summary.csv", index=False)
    else:
        logging.warning(
            f"No chloe.gb or cpgavas2.gb files found for sample {sample} to generate summary."
        )


# ── Top-level pipeline orchestrators ───────────────────────────────────────────


def process_short_reads(sample, config):
    """Run all parsing and file-gathering steps for short-read samples."""
    sample_cfg = config["samples"][sample]
    organelle = sample_cfg["organelle"].lower()
    annotation = sample_cfg["annotation"].lower()

    if sample_cfg["run_trimming"].lower() == "yes":
        logging.info(f"Getting Fastp files for sample: {sample}")
        os.makedirs(f"workflow/reports/{sample}/files/", exist_ok=True)
        shutil.copy(
            f"resources/{sample}/rawreads/fastp.html",
            f"workflow/reports/{sample}/files/fastp.html",
        )

    # NOVOPlasty
    if sample_cfg["run_novoplasty"].lower() == "yes":
        logging.info(f"Parsing NOVOPlasty sample: {sample}")
        for root, dirs, files in os.walk(f"results/{sample}/novoplasty"):
            for file in files:
                if "log_" in file and "_extended" not in file:
                    parse_novoplasty(os.path.join(root, file))
                    get_novoplasty_files(os.path.join(root, file))

        logging.info(f"Joining NOVOPlasty .csv files for sample: {sample}")
        combined_df = concat_csv_files(f"workflow/reports/{sample}/novoplasty")
        if combined_df is not None:
            combined_df = combined_df.sort_values(by=["Seed", "kmer"])
            combined_df.reset_index(drop=True, inplace=True)
            combined_df.to_csv(f"workflow/reports/{sample}/novoplasty.csv", index=False)

    # Annotation: MITOS2 or CPGAVAS2
    if organelle == "mito" and annotation == "yes":
        logging.info(f"Parsing MITOS2 sample: {sample}")
        for root, dirs, files in os.walk(f"results/{sample}/mitos2"):
            for file in files:
                if "mitos.log" in file:
                    parse_mitos2(os.path.join(root, file), config)
                    get_mitos2_files(os.path.join(root, file), config)

        logging.info(f"Joining MITOS2 .csv files for sample: {sample}")
        mitos2_dir = Path(f"workflow/reports/{sample}/mitos2")
        if mitos2_dir.exists():
            combined_df = concat_csv_files(mitos2_dir)
            if combined_df is not None:
                combined_df = combined_df.sort_values(by=["seed", "kmer"])
                combined_df.reset_index(drop=True, inplace=True)
                combined_df.to_csv(f"workflow/reports/{sample}/mitos2.csv", index=False)

    elif organelle == "chloro" and annotation == "yes":
        logging.info(f"Parsing CPGAVAS2 sample: {sample}")
        for root, dirs, files in os.walk(f"results/{sample}/cpgavas2"):
            for file in files:
                if "_reportTable.txt" in file:
                    parse_cpgavas2_report_table(os.path.join(root, file), config)
                    get_cpgavas2_files(os.path.join(root, file), config)
                elif ".annotation_with_problems.txt" in file:
                    parse_cpgavas2_problems(os.path.join(root, file), config)

        for report_type, out_name in [
            ("codon_usage", "cpgavas2_codon_usage"),
            ("intron_exon", "cpgavas2_intron_exon"),
            ("gene_composition", "cpgavas2_gene_composition"),
        ]:
            logging.info(f"Joining CPGAVAS2 Report Table {report_type} .csv files for sample: {sample}")
            combined_df = concat_csv_files(
                f"workflow/reports/{sample}/cpgavas2/report_table/{report_type}"
            )
            if combined_df is not None and not combined_df.empty:
                combined_df = combined_df.sort_values(by=["Seed", "Kmer"])
                combined_df.reset_index(drop=True, inplace=True)
                combined_df.to_csv(f"workflow/reports/{sample}/{out_name}.csv", index=False)

        logging.info(f"Joining CPGAVAS2 Problems .csv files for sample: {sample}")
        combined_df = concat_csv_files(f"workflow/reports/{sample}/cpgavas2/problems")
        if combined_df is not None and not combined_df.empty:
            combined_df = combined_df.sort_values(by=["Seed", "Kmer"])
            combined_df.reset_index(drop=True, inplace=True)
            combined_df.to_csv(f"workflow/reports/{sample}/cpgavas2_problems.csv", index=False)

    # Pilon
    logging.info(f"Parsing PILON sample: {sample}")
    for root, dirs, files in os.walk(f"results/{sample}/pilon"):
        for file in files:
            if file == "pilon.log":
                parse_pilon(os.path.join(root, file), config)
                get_pilon_files(os.path.join(root, file), config)

    logging.info(f"Joining PILON .csv files for sample: {sample}")
    pilon_dir = Path(f"workflow/reports/{sample}/pilon")
    if pilon_dir.exists():
        combined_df = concat_csv_files(pilon_dir)
        if combined_df is not None:
            combined_df = combined_df.sort_values(by=["seed", "kmer"])
            combined_df.reset_index(drop=True, inplace=True)
            combined_df.to_csv(f"workflow/reports/{sample}/pilon.csv", index=False)

    # nhmmer
    if sample_cfg["run_nhmmer"].lower() == "yes":
        for nhmmer_file, subdir in [
            ("rRNA-tRNA.tblout.out", "ncRNA"),
            ("intergenes_filter.tblout.out", "intergenes"),
        ]:
            logging.info(f"Parsing {subdir} NHMMER results for sample: {sample}")
            for root, dirs, files in os.walk(f"results/{sample}/nhmmer"):
                for file in files:
                    if file == nhmmer_file:
                        parse_nhmmer(os.path.join(root, file), subdir, config)

            logging.info(f"Joining {subdir} NHMMER .csv files for sample: {sample}")
            try:
                combined_df = concat_csv_files(f"workflow/reports/{sample}/nhmmer/{subdir}")
                if combined_df is not None:
                    combined_df = combined_df.sort_values(by=["seed", "kmer"])
                    combined_df.reset_index(drop=True, inplace=True)
                    combined_df.to_csv(
                        f"workflow/reports/{sample}/nhmmer_{subdir}.csv", index=False
                    )
            except Exception as e:
                logging.error(f"Error processing {subdir} NHMMER results for sample {sample}: {e}")

    # Images
    if sample_cfg["run_images"].lower() == "yes":
        copy_images(sample, config)

    # GenBank / FASTA / gene extraction
    if annotation == "yes":
        logging.info(f"Getting GenBank file for sample: {sample}")
        if sample_cfg["run_novoplasty"].lower() == "yes":
            copy_genbanks_novoplasty(sample, config)

        if sample_cfg["run_getorganelle"].lower() == "yes":
            copy_genbanks_getorganelle(sample, config)

        convert_genbanks_to_fastas(sample)
        extract_genes_from_genbanks(sample, config)
    else:
        # No annotation → no GenBank files; copy raw assembly FASTAs instead
        logging.info(f"Annotation disabled — copying raw assembly FASTAs for sample: {sample}")
        copy_assembly_fastas(sample, config)

    # Summary
    if organelle == "mito":
        generate_summary_mito(sample, config)
    elif organelle == "chloro":
        generate_summary_chloro(sample, config)


def process_long_reads(sample, config):
    """Run all parsing and file-gathering steps for long-read samples."""
    logging.info(f"Parsing MitoHifi from sample: {sample}")

    for annotation in os.listdir(f"results/{sample}/mitohifi"):
        mitohifi_dir = Path(f"results/{sample}/mitohifi/{annotation}")
        stats_file = mitohifi_dir / "contigs_stats.tsv"

        if stats_file.exists():
            parse_mitohifi(stats_file)

        logging.info(f"Copying files from sample: {sample}")
        for file in os.listdir(mitohifi_dir):
            source_file = mitohifi_dir / file
            new_filename = f"{annotation}_{file}"

            if file.endswith(".png"):
                dest_dir = Path(f"workflow/reports/{sample}/images")
            elif file.endswith(".gb"):
                dest_dir = Path(f"workflow/reports/{sample}/genbanks")
            elif file.endswith(".fasta"):
                dest_dir = Path(f"workflow/reports/{sample}/fastas")
            else:
                continue

            dest_dir.mkdir(parents=True, exist_ok=True)
            logging.info(f"Copying file: {file}")
            shutil.copy2(source_file, dest_dir / new_filename)

    logging.info(f"Joining MitoHifi .csv files for sample: {sample}")
    combined_df = concat_csv_files(f"workflow/reports/{sample}/mitohifi")
    if combined_df is not None:
        combined_df = combined_df.sort_values(by=["seed"])
        combined_df.reset_index(drop=True, inplace=True)
        combined_df.to_csv(f"workflow/reports/{sample}/mitohifi.csv", index=False)

    # nhmmer (long-read mode: no config argument -> seed at parts[3])
    if config["samples"][sample]["run_nhmmer"].lower() == "yes":
        for nhmmer_file, subdir in [
            ("rRNA-tRNA.tblout.out", "ncRNA"),
            ("intergenes_filter.tblout.out", "intergenes"),
        ]:
            logging.info(f"Parsing {subdir} NHMMER results for sample: {sample}")
            for root, dirs, files in os.walk(f"results/{sample}/nhmmer"):
                for file in files:
                    if file == nhmmer_file:
                        parse_nhmmer(os.path.join(root, file), subdir)  # no config = long-read mode

            logging.info(f"Joining {subdir} NHMMER .csv files for sample: {sample}")
            try:
                combined_df = concat_csv_files(f"workflow/reports/{sample}/nhmmer/{subdir}")
                if combined_df is not None:
                    combined_df = combined_df.sort_values(by=["seed"])
                    combined_df.reset_index(drop=True, inplace=True)
                    combined_df.to_csv(
                        f"workflow/reports/{sample}/nhmmer_{subdir}.csv", index=False
                    )
            except Exception as e:
                logging.error(f"Error processing {subdir} NHMMER results for sample {sample}: {e}")

    extract_genes_from_genbanks(sample, config)


# ── Entry point ────────────────────────────────────────────────────────────────


if __name__ == "__main__":
    FORMAT = "%(asctime)s [%(levelname)s] %(message)s"

    logging.basicConfig(
        level=logging.INFO,
        stream=sys.stdout,
        format=FORMAT,
        datefmt="%d-%m-%Y %H:%M:%S",
    )

    args = parse_arguments()
    sample_id = args.sample

    with open("config/snakemake_config.yaml", "r") as config_file:
        config = yaml.safe_load(config_file)

    for sample in config["samples"].keys():
        if sample != sample_id:
            continue

        seq_type = config["samples"][sample]["sequencing_type"].lower()
        if seq_type == "short":
            process_short_reads(sample, config)
        elif seq_type == "long":
            process_long_reads(sample, config)
