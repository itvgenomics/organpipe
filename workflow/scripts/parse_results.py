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


def get_sample_keys(yaml_file):
    with open(yaml_file, "r") as file:
        config = yaml.safe_load(file)  # Load the YAML file
        if "samples" in config:
            return config["samples"].keys()  # Get keys from the 'samples' section
        else:
            return None  # or raise an error if 'samples' is not found


def parse_novoplasty(logfile):
    logging.info(f"Parsing NOVOPlasty log file: {logfile}")

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
                    value = value
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

    sample = str(logfile).split("/")[1]
    seed = str(logfile).split("/")[3]
    kmer = str(logfile).split("/")[4].split("kmer")[1]

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

    if not os.path.exists(f"workflow/reports/{sample}/novoplasty"):
        os.makedirs(f"workflow/reports/{sample}/novoplasty", exist_ok=True)

    df.to_csv(
        f"workflow/reports/{sample}/novoplasty/{seed}_kmer{kmer}.csv", index=False
    )


def concat_csv_files(directory):
    # Create an empty list to hold the DataFrames
    dataframes = []

    try:
        # Loop through all files in the directory
        for filename in os.listdir(directory):
            if filename.endswith(".csv"):
                # Construct full file path
                file_path = os.path.join(directory, filename)
                # Read the CSV file and append the DataFrame to the list
                df = pd.read_csv(file_path)
                dataframes.append(df)

        try:
            # Concatenate all DataFrames in the list into a single DataFrame
            combined_df = pd.concat(dataframes, ignore_index=True)
            return combined_df

        except Exception as e:
            logging.error(f"Warning while parsing Table: {e}")

    except Exception as e:
        logging.error(
            f"Warning: {e}. This might be due to an empty directory or no CSV files found."
        )


def parse_mitos2(logfile):
    logging.info(f"Parsing MITOS2 files: {logfile}")

    sample = str(logfile).split("/")[1]
    seed = str(logfile).split("/")[3].split("_kmer")[0]
    kmer = str(logfile).split("/")[3].split("_kmer")[1]
    assembly_dir = os.path.dirname(logfile)
    assembly = str(assembly_dir).split("/")[-1]

    if not os.path.exists(f"workflow/reports/{sample}/mitos2"):
        os.makedirs(f"workflow/reports/{sample}/mitos2", exist_ok=True)

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
            "ribossomal_genes",
            "ribossomal_count",
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

    # Read the file contents into a variable
    with open(logfile, "r") as file:
        for line in file:
            if "missing" in line:
                key, value = line.split(":")
            else:
                value = ""
                key = "missing"
                mitos_df["missing"] = ""
            mitos_df[key] = value
            if "duplicated" in line:
                key, value = line.split(":")
            else:
                value = ""
                key = "duplicated"
                mitos_df["duplicated"] = ""
            mitos_df[key] = value

    # get gene order file
    with open(f"{assembly_dir}/result.geneorder", "r") as file:
        data = file.read()
        mitos_df["gene_order"] = data.split("\n")[1]

    # get gene order file
    with open(f"{assembly_dir}/result.geneorder", "r") as file:
        transporters = []
        ribossomals = []
        codings = []
        forward = []
        reverse = []
        origins = []

        data = file.read()
        mitos_df["gene_order"] = data.split("\n")[1]

        # Get forward and reverse genes
        for cds in data.split("\n")[1].split():
            if cds.startswith("-"):
                reverse.append(cds)
            else:
                forward.append(cds)

        # Get type of gene
        for cds in data.split("\n")[1].split():
            if cds.startswith("t") | cds.startswith("-t"):
                transporters.append(cds)
            elif cds.startswith("r") | cds.startswith("-r"):
                ribossomals.append(cds)
            elif (
                cds.startswith("OH")
                | cds.startswith("OL")
                | cds.startswith("-OH")
                | cds.startswith("-OL")
            ):
                origins.append(cds)
            else:
                codings.append(cds)

    mitos_df["forward_genes"] = [forward]
    mitos_df["forward_count"] = len(forward)
    mitos_df["reverse_genes"] = [reverse]
    mitos_df["reverse_count"] = len(reverse)
    mitos_df["transporter_genes"] = [transporters]
    mitos_df["transporter_count"] = len(transporters)
    mitos_df["ribossomal_genes"] = [ribossomals]

    unique_counts = {
        "rrnL": 0,
        "rrnS": 0,
    }  # Initialize counts for rrnL and rrnS
    for item in ribossomals:
        if "rrnL" in item:
            unique_counts["rrnL"] = 1
        if "rrnS" in item:
            unique_counts["rrnS"] = 1

    count = sum(unique_counts.values())
    mitos_df["ribossomal_count"] = count

    mitos_df["origins"] = [origins]

    unique_counts = {"OH": 0, "OL": 0}  # Initialize counts for OH and OL
    for item in origins:
        if "OH" in item:
            unique_counts["OH"] = 1
        if "OL" in item:
            unique_counts["OL"] = 1

    count = sum(unique_counts.values())
    mitos_df["origins_count"] = count

    mitos_df["coding_genes"] = [codings]
    mitos_df["coding_count"] = len(codings)

    # Convert the lists to strings and remove not wanted character
    columns_to_correct = [
        "ribossomal_genes",
        "origins",
        "coding_genes",
        "transporter_genes",
        "reverse_genes",
        "forward_genes",
        "coding_genes",
    ]

    for column in columns_to_correct:
        mitos_df[column] = mitos_df[column].apply(lambda x: str(x))
        mitos_df[column] = (
            mitos_df[column]
            .str.replace("[", "")
            .str.replace("]", "")
            .str.replace("'", "")
        )

    mitos_df.to_csv(f"workflow/reports/{sample}/mitos2/{assembly}.csv", index=False)


def parse_pilon(logfile):
    logging.info(f"Parsing Pilon data for: {logfile}")
    sample = str(logfile).split("/")[1]
    seed = "_".join(str(logfile).split("/")[3].split(sample)[1].split("_")[1:-2])
    kmer = str(logfile).split("/")[3].split("_")[-2]
    assembly_dir = os.path.dirname(logfile)
    assembly = str(assembly_dir).split("/")[-1]

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
            "max",
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
    with open(logfile, "r") as file:
        for line in file:
            if "Input genome size" in line:
                genome_size = line.split(":")[1].replace("\n", "").strip()
            elif line.startswith(f"results/{sample}/pilon"):
                data = line.split(".bam:")
                if len(data) == 7:
                    reads = str(data).split(",")[1].strip().split(" ")[1]
                    filtered = str(data).split(",")[2].strip().split(" ")[0]
                    mapped = str(data).split(",")[3].strip().split(" ")[0]
                    propper = str(data).split(",")[4].strip().split(" ")[0]
                    stray = str(data).split(",")[5].strip().split(" ")[0]
                    fr = str(data).split(",")[6].strip().split(" ")[1:]
                    max = str(data).split(",")[7].strip().split(" ")[1]
                else:
                    reads = str(data).split(",")[1].strip().split(" ")[1]
                    filtered = str(data).split(",")[2].strip().split(" ")[0]
                    mapped = str(data).split(",")[3].strip().split(" ")[0]
                    propper = str(data).split(",")[4].strip().split(" ")[0]
                    stray = str(data).split(",")[5].strip().split(" ")[0]
                    max = str(data).split(",")[6].strip().split(" ")[1]

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

            try:
                genome_size
            except NameError:
                genome_size = ""
            try:
                data
            except NameError:
                data = ""
            try:
                reads
            except NameError:
                reads = ""
            try:
                filtered
            except NameError:
                filtered = ""
            try:
                mapped
            except NameError:
                mapped = ""
            try:
                propper
            except NameError:
                propper = ""
            try:
                stray
            except NameError:
                stray = ""
            try:
                fr
            except NameError:
                fr = ""
            try:
                max
            except NameError:
                max = ""
            try:
                coverage
            except NameError:
                coverage = ""
            try:
                minDepth
            except NameError:
                minDepth = ""
            try:
                confirmed
            except NameError:
                confirmed = ""
            try:
                corrected
            except NameError:
                corrected = ""
            try:
                mean_frags_coverage
            except NameError:
                mean_frags_coverage = ""
            try:
                mean_total_coverage
            except NameError:
                mean_total_coverage = ""

        pilon_df["genome_size"] = genome_size
        pilon_df["reads"] = reads
        pilon_df["filtered"] = filtered
        pilon_df["mapped"] = mapped
        pilon_df["propper"] = propper
        pilon_df["stray"] = stray
        pilon_df["fr"] = [fr]
        pilon_df["max"] = max
        pilon_df["coverage"] = coverage
        pilon_df["minDepth"] = minDepth
        pilon_df["confirmed"] = confirmed
        pilon_df["corrected"] = corrected
        pilon_df["mean_frags_coverage"] = mean_frags_coverage
        pilon_df["mean_total_coverage"] = mean_total_coverage

        # Get number of changes
        changes = []
        with open(f"results/{sample}/pilon/{assembly}/{assembly}.changes", "r") as file:
            lines = file.readlines()  # Read all lines at once
            line_count = len(lines)  # Count the lines
            if line_count > 0:
                for line in lines:  # Iterate over the list of lines
                    parts = line.strip().split()
                    original_coord = parts[0]
                    new_coord = parts[1]
                    original_seq = parts[2]
                    new_seq = parts[3]

                    string = f"original_coord: {original_coord} new_coord: {new_coord} original_seq: {original_seq} new_seq: {new_seq}"
                    changes.append(string)

        pilon_df["changes_number"] = line_count
        pilon_df["changes"] = str(changes)

        if not os.path.exists(f"workflow/reports/{sample}/pilon"):
            os.makedirs(f"workflow/reports/{sample}/pilon", exist_ok=True)

        pilon_df.to_csv(f"workflow/reports/{sample}/pilon/{assembly}.csv", index=False)


def parse_cpgavas2_report_table(logfile):
    logging.info(f"Parsing CPGAVAS2 Report Table: {logfile}")

    sample = str(logfile).split("/")[1]
    seed = str(logfile).split("/")[3].split("_kmer")[0]
    kmer = str(logfile).split("/")[3].split("_kmer")[1]
    assembly_dir = os.path.dirname(logfile)
    assembly = str(assembly_dir).split("/")[-1]

    output_dirs = [
        f"workflow/reports/{sample}/cpgavas2/report_table/gene_composition",
        f"workflow/reports/{sample}/cpgavas2/report_table/intron_exon",
        f"workflow/reports/{sample}/cpgavas2/report_table/codon_usage",
    ]

    for output_dir in output_dirs:
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

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

            df1.to_csv(
                f"workflow/reports/{sample}/cpgavas2/report_table/gene_composition/{assembly}_gene_composition.csv",
                index=False,
            )
        except Exception as e:
            logging.error(f"Warning while parsing Table: {e}")

        # Parse Table 2: Lengths of Introns and Exons
        try:
            table2_data = [line.split("\t") for line in sections[5].strip().split("\n")]
            df2 = pd.DataFrame(
                table2_data,
                columns=[
                    "Gene",
                    "Strand",
                    "Start",
                    "End",
                    "ExonI",
                    "IntronI",
                    "ExonII",
                    "IntronII",
                    "ExonIII",
                ],
            )
            df2.insert(0, "Kmer", kmer)
            df2.insert(0, "Seed", seed)
            df2.insert(0, "Assembly", assembly)
            df2.insert(0, "Sample", sample)

            df2.to_csv(
                f"workflow/reports/{sample}/cpgavas2/report_table/intron_exon/{assembly}_intron_exon.csv",
                index=False,
            )
        except Exception as e:
            logging.error(f"Warning while parsing Table: {e}")

        # Parse Table 3: Codon Usage
        try:
            table3_data = [
                line.split("\t") for line in sections[7].strip().split("\n")[1:]
            ]
            df3 = pd.DataFrame(
                table3_data, columns=["Codon", "Amino acid", "Frequency", "Number"]
            )
            df3.insert(0, "Kmer", kmer)
            df3.insert(0, "Seed", seed)
            df3.insert(0, "Assembly", assembly)
            df3.insert(0, "Sample", sample)

            df3.to_csv(
                f"workflow/reports/{sample}/cpgavas2/report_table/codon_usage/{assembly}_codon_usage.csv",
                index=False,
            )
        except Exception as e:
            logging.error(f"Warning while parsing Table: {e}")


def parse_cpgavas2_problems(logfile):
    sample = str(logfile).split("/")[1]
    seed = str(logfile).split("/")[3].split("_kmer")[0]
    kmer = str(logfile).split("/")[3].split("_kmer")[1]
    assembly_dir = os.path.dirname(logfile)
    assembly = str(assembly_dir).split("/")[-1]

    logging.info(f"Parsing Problems: {logfile}")

    if not os.path.exists(f"workflow/reports/{sample}/cpgavas2/problems"):
        os.makedirs(f"workflow/reports/{sample}/cpgavas2/problems")

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

        df.to_csv(
            f"workflow/reports/{sample}/cpgavas2/problems/{assembly}_problems.csv",
            index=False,
        )


def parse_ncRNA_nhmmer(logfile):
    logging.info(f"Parsing rRNA-tRNA NHMMER data: {logfile}")

    sample = str(logfile).split("/")[1]
    seed = str(logfile).split("/")[3].split("_kmer")[0]
    kmer = str(logfile).split("/")[3].split("_kmer")[1]
    assembly_dir = os.path.dirname(logfile)
    assembly = str(assembly_dir).split("/")[-1]

    if not os.path.exists(f"workflow/reports/{sample}/nhmmer/ncRNA/"):
        os.makedirs(f"workflow/reports/{sample}/nhmmer/ncRNA/", exist_ok=True)

    try:
        df = pd.read_csv(
            logfile,
            delim_whitespace=True,
            comment="#",
        )

        hmmercolumns = [
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

        df.columns = hmmercolumns

        df.insert(0, "kmer", kmer)
        df.insert(0, "seed", seed)
        df.insert(0, "assembly", assembly)
        df.insert(0, "sample", sample)

        df.to_csv(f"workflow/reports/{sample}/nhmmer/ncRNA/{assembly}.csv", index=False)

    except Exception as e:
        logging.error(f"Warning while parsing Table: {e}")


def parse_intergenes_nhmmer(logfile):
    logging.info(f"Parsing intergenes NHMMER data: {logfile}")

    sample = str(logfile).split("/")[1]
    seed = str(logfile).split("/")[3].split("_kmer")[0]
    kmer = str(logfile).split("/")[3].split("_kmer")[1]
    assembly_dir = os.path.dirname(logfile)
    assembly = str(assembly_dir).split("/")[-1]

    if not os.path.exists(f"workflow/reports/{sample}/nhmmer/intergenes/"):
        os.makedirs(f"workflow/reports/{sample}/nhmmer/intergenes/", exist_ok=True)

    try:
        df = pd.read_csv(
            logfile,
            delim_whitespace=True,
            comment="#",
        )

        hmmercolumns = [
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

        df.columns = hmmercolumns

        df.insert(0, "kmer", kmer)
        df.insert(0, "seed", seed)
        df.insert(0, "assembly", assembly)
        df.insert(0, "sample", sample)

        df.to_csv(
            f"workflow/reports/{sample}/nhmmer/intergenes/{assembly}.csv", index=False
        )

    except Exception as e:
        logging.error(f"Warning while parsing Table: {e}")


def parse_mitohifi(statsfile):
    logging.info(f"Parsing MitoHifi for file: {statsfile}")

    df = pd.read_csv(statsfile, sep="\t", comment="#")

    sample = str(statsfile).split("/")[1]
    seed = str(statsfile).split("/")[3]
    df.insert(0, "seed", seed)
    df.insert(0, "sample", sample)

    if not os.path.exists(f"workflow/reports/{sample}/mitohifi"):
        os.makedirs(f"workflow/reports/{sample}/mitohifi")

    df.to_csv(f"workflow/reports/{sample}/mitohifi/{seed}.csv", index=False)


def parse_ncRNA_nhmmer_long(logfile):
    logging.info(f"Parsing rRNA-tRNA NHMMER data: {logfile}")

    sample = str(logfile).split("/")[1]
    seed = str(logfile).split("/")[3]

    if not os.path.exists(f"workflow/reports/{sample}/nhmmer/ncRNA/"):
        os.makedirs(f"workflow/reports/{sample}/nhmmer/ncRNA/", exist_ok=True)

    try:
        df = pd.read_csv(
            logfile,
            delim_whitespace=True,
            comment="#",
        )

        hmmercolumns = [
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

        df.columns = hmmercolumns

        df.insert(0, "seed", seed)
        df.insert(0, "sample", sample)

        df.to_csv(f"workflow/reports/{sample}/nhmmer/ncRNA/{seed}.csv", index=False)

    except Exception as e:
        logging.error(f"Warning while parsing Table: {e}")


def parse_intergenes_nhmmer_long(logfile):
    logging.info(f"Parsing intergenes NHMMER data: {logfile}")

    sample = str(logfile).split("/")[1]
    seed = str(logfile).split("/")[3]

    if not os.path.exists(f"workflow/reports/{sample}/nhmmer/intergenes/"):
        os.makedirs(f"workflow/reports/{sample}/nhmmer/intergenes/", exist_ok=True)

    try:
        df = pd.read_csv(
            logfile,
            delim_whitespace=True,
            comment="#",
        )

        hmmercolumns = [
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

        df.columns = hmmercolumns

        df.insert(0, "seed", seed)
        df.insert(0, "sample", sample)

        df.to_csv(
            f"workflow/reports/{sample}/nhmmer/intergenes/{seed}.csv", index=False
        )

    except Exception as e:
        logging.error(f"Warning while parsing Table: {e}")


def convert_genbank_to_fasta(genbank_file, fasta_file):
    try:
        with open(fasta_file, "w") as output_handle:
            # Parse the GenBank file and write to FASTA
            for record in SeqIO.parse(genbank_file, "genbank"):
                SeqIO.write(record, output_handle, "fasta")
        logging.info(f"Successfully converted {genbank_file} to {fasta_file}")
    except Exception as e:
        logging.error(f"An error occurred: {e}")


def get_novoplasty_files(logfile):
    sample = str(logfile).split("/")[1]
    seed = str(logfile).split("/")[3]
    kmer = str(logfile).split("/")[4].split("kmer")[1]

    logging.info(
        f"Getting NOVOPlasty files for sample: {sample}, seed: {seed}, kmer: {kmer}"
    )

    novoplasty_outdir = f"workflow/reports/{sample}/files/{seed}/kmer{kmer}"

    if not os.path.exists(novoplasty_outdir):
        os.makedirs(novoplasty_outdir)

    shutil.copy(
        f"results/{sample}/novoplasty/{seed}/kmer{kmer}/Assembled_reads_{sample}_R1.fasta",
        f"{novoplasty_outdir}/Assembled_reads_R1.fasta",
    )
    shutil.copy(
        f"results/{sample}/novoplasty/{seed}/kmer{kmer}/Assembled_reads_{sample}_R2.fasta",
        f"{novoplasty_outdir}/Assembled_reads_R2.fasta",
    )
    shutil.copy(
        f"results/{sample}/novoplasty/{seed}/kmer{kmer}/log_{sample}.txt",
        f"{novoplasty_outdir}/log_novoplasty.txt",
    )


def get_mitos2_files(logfile):
    sample = str(logfile).split("/")[1]
    seed = str(logfile).split("/")[3].split("_kmer")[0]
    kmer = str(logfile).split("/")[3].split("_kmer")[1]
    assembly_dir = os.path.dirname(logfile)
    assembly = str(assembly_dir).split("/")[-1]

    logging.info(
        f"Getting NOVOPlasty files for sample: {sample}, seed: {seed}, kmer: {kmer}, assembly: {assembly}"
    )

    mitos2_outdir = f"workflow/reports/{sample}/files/{seed}/kmer{kmer}"

    if not os.path.exists(mitos2_outdir):
        os.makedirs(mitos2_outdir)

    shutil.copy(
        f"results/{sample}/mitos2/{seed}_kmer{kmer}/{assembly}/ignored.mitos",
        f"{mitos2_outdir}/{assembly}.ignored.mitos",
    )
    shutil.copy(
        f"results/{sample}/mitos2/{seed}_kmer{kmer}/{assembly}/result.bed",
        f"{mitos2_outdir}/{assembly}.result.bed",
    )
    shutil.copy(
        f"results/{sample}/mitos2/{seed}_kmer{kmer}/{assembly}/result.faa",
        f"{mitos2_outdir}/{assembly}.result.faa",
    )
    shutil.copy(
        f"results/{sample}/mitos2/{seed}_kmer{kmer}/{assembly}/result.fas",
        f"{mitos2_outdir}/{assembly}.result.fas",
    )
    shutil.copy(
        f"results/{sample}/mitos2/{seed}_kmer{kmer}/{assembly}/result.geneorder",
        f"{mitos2_outdir}/{assembly}.result.geneorder",
    )
    shutil.copy(
        f"results/{sample}/mitos2/{seed}_kmer{kmer}/{assembly}/result.gff",
        f"{mitos2_outdir}/{assembly}.result.gff",
    )
    shutil.copy(
        f"results/{sample}/mitos2/{seed}_kmer{kmer}/{assembly}/result.mitos",
        f"{mitos2_outdir}/{assembly}.result.mitos",
    )
    shutil.copy(
        f"results/{sample}/mitos2/{seed}_kmer{kmer}/{assembly}/result.seq",
        f"{mitos2_outdir}/{assembly}.result.seq",
    )
    shutil.copy(
        f"results/{sample}/mitos2/{seed}_kmer{kmer}/{assembly}/stst.dat",
        f"{mitos2_outdir}/{assembly}.stst.dat",
    )
    shutil.copy(
        f"results/{sample}/mitos2/{seed}_kmer{kmer}/{assembly}/mitfi-global/sequence.fas-0_tRNAout.nc",
        f"{mitos2_outdir}/{assembly}_tRNAout.nc",
    )
    shutil.copy(
        f"results/{sample}/mitos2/{seed}_kmer{kmer}/{assembly}/mitfi-global/sequence.fas-0_rRNAout.nc",
        f"{mitos2_outdir}/{assembly}_rRNAout.nc",
    )


def get_pilon_files(logfile):

    sample = str(logfile).split("/")[1]
    seed = "_".join(str(logfile).split("/")[3].split(sample)[1].split("_")[1:-2])
    kmer = str(logfile).split("/")[3].split("_")[-2]
    assembly_dir = os.path.dirname(logfile)
    assembly = str(assembly_dir).split("/")[-1]

    logging.info(
        f"Getting NOVOPlasty files for sample: {sample}, seed: {seed}, kmer: {kmer}, assembly: {assembly}"
    )

    pilon_outdir = f"workflow/reports/{sample}/files/{seed}/kmer{kmer}"

    if not os.path.exists(pilon_outdir):
        os.makedirs(pilon_outdir)

    shutil.copy(
        f"results/{sample}/pilon/{assembly}/{assembly}.changes",
        f"{pilon_outdir}/{assembly}_pilon.changes",
    )
    shutil.copy(
        f"results/{sample}/pilon/{assembly}/pilon.log",
        f"{pilon_outdir}/{assembly}_pilon.log",
    )


def ignore_symlinks(src, names):
    return [name for name in names if os.path.islink(os.path.join(src, name))]


def get_cpgavas2_files(logfile):
    sample = str(logfile).split("/")[1]
    seed = str(logfile).split("/")[3].split("_kmer")[0]
    kmer = str(logfile).split("/")[3].split("_kmer")[1]
    assembly_dir = os.path.dirname(logfile)
    assembly = str(assembly_dir).split("/")[-1]

    logging.info(
        f"Copying CPGAVAS2 files for sample: {sample}, seed: {seed}, kmer: {kmer}, assembly: {assembly}"
    )

    source_path = f"results/{sample}/cpgavas2/{seed}_kmer{kmer}/{assembly}"
    destination_path = f"workflow/reports/{sample}/files/cpgavas2/{assembly}"

    if not os.path.exists(destination_path):
        shutil.copytree(source_path, destination_path, ignore=ignore_symlinks)


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


if __name__ == "__main__":
    FORMAT = "%(asctime)s [%(levelname)s] %(message)s"

    logging.basicConfig(
        level=logging.INFO,
        stream=sys.stdout,
        format=FORMAT,
        datefmt="%d-%m-%Y %H:%M:%S",
    )

    with open("config/snakemake_config.yaml", "r") as config_file:
        config = yaml.safe_load(config_file)

        sample_keys = config["samples"].keys()

        for sample in sample_keys:
            if config["samples"][sample]["sequencing_type"].lower() == "short":
                logging.info(f"Parsing NOVOPlasty sample: {sample}")
                for root, dirs, files in os.walk(f"results/{sample}/novoplasty"):
                    for file in files:
                        if "log_" in file:
                            parse_novoplasty(os.path.join(root, file))
                            get_novoplasty_files(os.path.join(root, file))

                logging.info(f"Joining NOVOPlasty .csv files for sample: {sample}")
                combined_df = concat_csv_files(f"workflow/reports/{sample}/novoplasty")

                combined_df = combined_df.sort_values(by=["Seed", "kmer"])
                combined_df.reset_index(drop=True, inplace=True)

                combined_df.to_csv(
                    f"workflow/reports/{sample}/novoplasty.csv", index=False
                )

                if config["samples"][sample]["organelle"].lower() == "mito":
                    logging.info(f"Parsing MITOS2 sample: {sample}")
                    for root, dirs, files in os.walk(f"results/{sample}/mitos2"):
                        for file in files:
                            if "mitos.log" in file:
                                parse_mitos2(os.path.join(root, file))
                                get_mitos2_files(os.path.join(root, file))

                    logging.info(f"Joining MITOS2 .csv files for sample: {sample}")

                    if os.path.exists(f"workflow/reports/{sample}/mitos2"):

                        combined_df = concat_csv_files(
                            f"workflow/reports/{sample}/mitos2"
                        )

                        combined_df = combined_df.sort_values(by=["seed", "kmer"])
                        combined_df.reset_index(drop=True, inplace=True)

                        combined_df.to_csv(
                            f"workflow/reports/{sample}/mitos2.csv", index=False
                        )

                elif config["samples"][sample]["organelle"].lower() == "chloro":
                    logging.info(f"Parsing CPGAVAS2 sample: {sample}")
                    for root, dirs, files in os.walk(f"results/{sample}/cpgavas2"):
                        for file in files:
                            if "_reportTable.txt" in file:
                                parse_cpgavas2_report_table(os.path.join(root, file))
                                get_cpgavas2_files(os.path.join(root, file))

                            elif ".annotation_with_problems.txt" in file:
                                parse_cpgavas2_problems(os.path.join(root, file))

                    logging.info(
                        f"Joining CPGAVAS2 Report Table Condon Usage .csv files for sample: {sample}"
                    )

                    combined_df = concat_csv_files(
                        f"workflow/reports/{sample}/cpgavas2/report_table/codon_usage"
                    )

                    if combined_df is not None and not combined_df.empty:
                        combined_df = combined_df.sort_values(by=["Seed", "Kmer"])
                        combined_df.reset_index(drop=True, inplace=True)

                        combined_df.to_csv(
                            f"workflow/reports/{sample}/cpgavas2_codon_usage.csv",
                            index=False,
                        )

                    logging.info(
                        f"Joining CPGAVAS2 Report Table Intron Exon .csv files for sample: {sample}"
                    )

                    combined_df = concat_csv_files(
                        f"workflow/reports/{sample}/cpgavas2/report_table/intron_exon"
                    )

                    if combined_df is not None and not combined_df.empty:
                        combined_df = combined_df.sort_values(by=["Seed", "Kmer"])
                        combined_df.reset_index(drop=True, inplace=True)

                        combined_df.to_csv(
                            f"workflow/reports/{sample}/cpgavas2_intron_exon.csv",
                            index=False,
                        )

                    logging.info(
                        f"Joining CPGAVAS2 Report Table Gene Composition .csv files for sample: {sample}"
                    )

                    combined_df = concat_csv_files(
                        f"workflow/reports/{sample}/cpgavas2/report_table/gene_composition"
                    )

                    if combined_df is not None and not combined_df.empty:
                        combined_df = combined_df.sort_values(by=["Seed", "Kmer"])
                        combined_df.reset_index(drop=True, inplace=True)

                        combined_df.to_csv(
                            f"workflow/reports/{sample}/cpgavas2_gene_composition.csv",
                            index=False,
                        )

                    logging.info(
                        f"Joining CPGAVAS2 Problems .csv files for sample: {sample}"
                    )

                    combined_df = concat_csv_files(
                        f"workflow/reports/{sample}/cpgavas2/problems"
                    )

                    if combined_df is not None and not combined_df.empty:
                        combined_df = combined_df.sort_values(by=["Seed", "Kmer"])
                        combined_df.reset_index(drop=True, inplace=True)

                        combined_df.to_csv(
                            f"workflow/reports/{sample}/cpgavas2_problems.csv",
                            index=False,
                        )

                logging.info(f"Parsing PILON sample: {sample}")
                for root, dirs, files in os.walk(f"results/{sample}/pilon"):
                    for file in files:
                        if file == "pilon.log":
                            parse_pilon(os.path.join(root, file))
                            get_pilon_files(os.path.join(root, file))

                logging.info(f"Joining PILON .csv files for sample: {sample}")

                if os.path.exists(f"workflow/reports/{sample}/pilon"):
                    combined_df = concat_csv_files(f"workflow/reports/{sample}/pilon")

                    combined_df = combined_df.sort_values(by=["seed", "kmer"])
                    combined_df.reset_index(drop=True, inplace=True)

                    combined_df.to_csv(
                        f"workflow/reports/{sample}/pilon.csv", index=False
                    )

                if config["samples"][sample]["run_nhmmer"].lower() == "yes":
                    logging.info(f"Parsing ncRNA NHMMER results for sample: {sample}")
                    for root, dirs, files in os.walk(f"results/{sample}/nhmmer"):
                        for file in files:
                            if file == "rRNA-tRNA.tblout.out":
                                parse_ncRNA_nhmmer(os.path.join(root, file))

                    logging.info(
                        f"Joining ncRNA NHMMER .csv files for sample: {sample}"
                    )
                    try:
                        combined_df = concat_csv_files(
                            f"workflow/reports/{sample}/nhmmer/ncRNA"
                        )

                        combined_df = combined_df.sort_values(by=["seed", "kmer"])
                        combined_df.reset_index(drop=True, inplace=True)

                        combined_df.to_csv(
                            f"workflow/reports/{sample}/nhmmer_ncRNA.csv",
                            index=False,
                        )
                    except Exception as e:
                        logging.error(
                            f"Error processing ncRNA NHMMER results for sample {sample}: {e}"
                        )

                    logging.info(
                        f"Parsing intergenes NHMMER results for sample: {sample}"
                    )
                    for root, dirs, files in os.walk(f"results/{sample}/nhmmer"):
                        for file in files:
                            if file == "intergenes_filter.tblout.out":
                                parse_intergenes_nhmmer(os.path.join(root, file))

                    logging.info(
                        f"Joining intergenes NHMMER .csv files for sample: {sample}"
                    )
                    try:
                        combined_df = concat_csv_files(
                            f"workflow/reports/{sample}/nhmmer/intergenes"
                        )

                        combined_df = combined_df.sort_values(by=["seed", "kmer"])
                        combined_df.reset_index(drop=True, inplace=True)

                        combined_df.to_csv(
                            f"workflow/reports/{sample}/nhmmer_intergenes.csv",
                            index=False,
                        )
                    except Exception as e:
                        logging.error(
                            f"Error processing intergenes NHMMER results for sample {sample}: {e}"
                        )

                if config["samples"][sample]["run_images"].lower() == "yes":
                    logging.info(f"Getting all images from sample: {sample}")
                    for root, dirs, files in os.walk(f"results/{sample}/images"):
                        for file in files:
                            if file.endswith(".png") and ".depth." in file:
                                depth_img_dir = (
                                    f"workflow/reports/{sample}/images/depth"
                                )

                                if not os.path.exists(depth_img_dir):
                                    os.makedirs(depth_img_dir, exist_ok=True)

                                logging.info(f"Copying file: {file}")
                                shutil.copy(os.path.join(root, file), depth_img_dir)

                            elif file.endswith(".png") and "_recruitment_plot." in file:
                                rec_plot_dir = (
                                    f"workflow/reports/{sample}/images/recruitment_plot"
                                )

                                if not os.path.exists(rec_plot_dir):
                                    os.makedirs(rec_plot_dir, exist_ok=True)

                                logging.info(f"Copying file: {file}")
                                shutil.copy(os.path.join(root, file), rec_plot_dir)

                            elif file.endswith(".png"):
                                ogdraw_dir = f"workflow/reports/{sample}/images/ogdraw"

                                if not os.path.exists(ogdraw_dir):
                                    os.makedirs(ogdraw_dir, exist_ok=True)

                                logging.info(f"Copying file: {file}")
                                shutil.copy(os.path.join(root, file), ogdraw_dir)

                logging.info(f"Getting GenBank file for sample: {sample}")

                kmers = [kmer for kmer in config["samples"][sample]["kmers"]]
                seeds = [seed for seed in config["samples"][sample]["seeds"]]

                gb_files = []
                for root, dirs, files in os.walk(f"results/{sample}/genbanks"):
                    for file in files:
                        if file.endswith(".gb"):
                            gb_files.append(file)

                for seed in seeds:
                    for file in gb_files:
                        if f"_{seed}_" in file:
                            logging.info(f"Copying GenBank file: {file}")

                            if not os.path.exists(
                                f"workflow/reports/{sample}/genbanks/{seed}"
                            ):
                                os.makedirs(
                                    f"workflow/reports/{sample}/genbanks/{seed}"
                                )

                            shutil.copy(
                                f"results/{sample}/genbanks/{file}",
                                f"workflow/reports/{sample}/genbanks/{seed}/",
                            )

                    if not os.path.exists(f"workflow/reports/{sample}/fastas/{seed}"):
                        os.makedirs(
                            f"workflow/reports/{sample}/fastas/{seed}", exist_ok=True
                        )

                logging.info(f"Getting FASTA file for sample: {sample}")

                if os.path.exists(f"workflow/reports/{sample}/genbanks"):
                    genbanks_dir = f"workflow/reports/{sample}/genbanks"

                    for genbanks in os.listdir(genbanks_dir):
                        for file in os.listdir(os.path.join(genbanks_dir, genbanks)):
                            logging.info(f"Getting FASTA file for sample: {file}")

                            genbank = os.path.join(genbanks_dir, genbanks, file)
                            fasta_file = genbank.replace(
                                "/genbanks/", "/fastas/"
                            ).replace(".gb", ".fasta")

                            convert_genbank_to_fasta(genbank, fasta_file)

            elif config["samples"][sample]["sequencing_type"].lower() == "long":
                logging.info(f"Parsing MitoHifi from sample: {sample}")
                for annotation in os.listdir(f"results/{sample}/mitohifi"):
                    mitohifi_dir = f"results/{sample}/mitohifi/{annotation}"

                    if os.path.exists(
                        f"results/{sample}/mitohifi/{annotation}/contigs_stats.tsv"
                    ):
                        parse_mitohifi(
                            f"results/{sample}/mitohifi/{annotation}/contigs_stats.tsv"
                        )

                    logging.info(f"Copying files from sample: {sample}")
                    for file in os.listdir(mitohifi_dir):
                        if file.endswith(".png"):
                            logging.info(f"Copying file: {file}")
                            source_file = os.path.join(mitohifi_dir, file)
                            new_filename = f"{annotation}_{file}"
                            dest_file = os.path.join(
                                f"workflow/reports/{sample}/images", new_filename
                            )

                            if not os.path.exists(f"workflow/reports/{sample}/images"):
                                os.makedirs(
                                    f"workflow/reports/{sample}/images", exist_ok=True
                                )

                            shutil.copy2(source_file, dest_file)

                        elif file.endswith(".gb"):
                            logging.info(f"Copying file: {file}")
                            source_file = os.path.join(mitohifi_dir, file)
                            new_filename = f"{annotation}_{file}"
                            dest_file = os.path.join(
                                f"workflow/reports/{sample}/genbanks", new_filename
                            )

                            if not os.path.exists(
                                f"workflow/reports/{sample}/genbanks"
                            ):
                                os.makedirs(
                                    f"workflow/reports/{sample}/genbanks", exist_ok=True
                                )

                            shutil.copy2(source_file, dest_file)

                        elif file.endswith(".fasta"):
                            logging.info(f"Copying file: {file}")
                            source_file = os.path.join(mitohifi_dir, file)
                            new_filename = f"{annotation}_{file}"
                            dest_file = os.path.join(
                                f"workflow/reports/{sample}/fastas", new_filename
                            )

                            if not os.path.exists(f"workflow/reports/{sample}/fastas"):
                                os.makedirs(
                                    f"workflow/reports/{sample}/fastas", exist_ok=True
                                )

                            shutil.copy2(source_file, dest_file)

                logging.info(
                    f"Joining intergenes MitoHifi .csv files for sample: {sample}"
                )

                combined_df = concat_csv_files(f"workflow/reports/{sample}/mitohifi")

                combined_df = combined_df.sort_values(by=["seed"])
                combined_df.reset_index(drop=True, inplace=True)

                combined_df.to_csv(
                    f"workflow/reports/{sample}/mitohifi.csv",
                    index=False,
                )

                if config["samples"][sample]["run_nhmmer"].lower() == "yes":
                    logging.info(f"Parsing ncRNA NHMMER results for sample: {sample}")
                    for root, dirs, files in os.walk(f"results/{sample}/nhmmer"):
                        for file in files:
                            if file == "rRNA-tRNA.tblout.out":
                                parse_ncRNA_nhmmer_long(os.path.join(root, file))

                    logging.info(
                        f"Joining ncRNA NHMMER .csv files for sample: {sample}"
                    )
                    try:
                        combined_df = concat_csv_files(
                            f"workflow/reports/{sample}/nhmmer/ncRNA"
                        )

                        combined_df = combined_df.sort_values(by=["seed"])
                        combined_df.reset_index(drop=True, inplace=True)

                        combined_df.to_csv(
                            f"workflow/reports/{sample}/nhmmer_ncRNA.csv",
                            index=False,
                        )
                    except Exception as e:
                        logging.error(
                            f"Error processing ncRNA NHMMER results for sample {sample}: {e}"
                        )

                    logging.info(
                        f"Parsing intergenes NHMMER results for sample: {sample}"
                    )
                    for root, dirs, files in os.walk(f"results/{sample}/nhmmer"):
                        for file in files:
                            if file == "intergenes_filter.tblout.out":
                                parse_intergenes_nhmmer_long(os.path.join(root, file))

                    logging.info(
                        f"Joining intergenes NHMMER .csv files for sample: {sample}"
                    )
                    try:
                        combined_df = concat_csv_files(
                            f"workflow/reports/{sample}/nhmmer/intergenes"
                        )

                        combined_df = combined_df.sort_values(by=["seed"])
                        combined_df.reset_index(drop=True, inplace=True)

                        combined_df.to_csv(
                            f"workflow/reports/{sample}/nhmmer_intergenes.csv",
                            index=False,
                        )
                    except Exception as e:
                        logging.error(
                            f"Error processing intergenes NHMMER results for sample {sample}: {e}"
                        )

            if os.path.exists(f"workflow/reports/{sample}/genes"):
                shutil.rmtree(
                    f"workflow/reports/{sample}/genes",
                )

            if os.path.exists(f"workflow/reports/{sample}/genbanks"):
                for root, dirs, files in os.walk(f"workflow/reports/{sample}/genbanks"):
                    if config["samples"][sample]["sequencing_type"].lower() == "short":

                        for file in files:
                            if (
                                config["samples"][sample]["organelle"].lower() == "mito"
                                and ".rotated.gb" not in file
                            ):

                                logging.info(
                                    f"Getting genes from GenBank file: {os.path.join(root, file)}"
                                )
                                all_data = extract_features(os.path.join(root, file))
                                for data in all_data:
                                    assembly = file.split("/")[-1].replace(".gb", "")
                                    feature = data["feature"]
                                    seq = data["seq"]

                                    if not os.path.exists(
                                        f"workflow/reports/{sample}/genes"
                                    ):
                                        os.makedirs(
                                            f"workflow/reports/{sample}/genes",
                                            exist_ok=True,
                                        )

                                    file_path = f"workflow/reports/{sample}/genes/{feature}.fasta"
                                    with open(
                                        file_path,
                                        "a" if os.path.exists(file_path) else "w",
                                    ) as output_fasta:
                                        output_fasta.write(
                                            f">{feature}_{assembly}\n{seq}\n"
                                        )

                    elif config["samples"][sample]["sequencing_type"].lower() == "long":
                        for root, dirs, files in os.walk(
                            f"workflow/reports/{sample}/genbanks"
                        ):
                            for file in files:

                                logging.info(
                                    f"Getting genes from GenBank file: {os.path.join(root, file)}"
                                )

                                all_data = extract_features(os.path.join(root, file))

                                for data in all_data:
                                    assembly = file.split("/")[-1].replace(".gb", "")
                                    feature = data["feature"]
                                    seq = data["seq"]

                                    if not os.path.exists(
                                        f"workflow/reports/{sample}/genes"
                                    ):
                                        os.makedirs(
                                            f"workflow/reports/{sample}/genes",
                                            exist_ok=True,
                                        )

                                    file_path = f"workflow/reports/{sample}/genes/{feature}.fasta"
                                    with open(
                                        file_path,
                                        "a" if os.path.exists(file_path) else "w",
                                    ) as output_fasta:
                                        output_fasta.write(
                                            f">{feature}_{assembly}\n{seq}\n"
                                        )

            if config["samples"][sample]["sequencing_type"].lower() == "short":
                if config["samples"][sample]["organelle"].lower() == "mito":
                    if (
                        os.path.exists(f"workflow/reports/{sample}/novoplasty.csv")
                        and os.path.exists(f"workflow/reports/{sample}/pilon.csv")
                        and os.path.exists(f"workflow/reports/{sample}/mitos2.csv")
                    ):

                        logging.info(f"Writting abstract.csv for sample: {sample}")

                        df_novoplasty = pd.read_csv(
                            f"workflow/reports/{sample}/novoplasty.csv"
                        )
                        df_mitos2 = pd.read_csv(f"workflow/reports/{sample}/mitos2.csv")

                        df_abstract = pd.merge(
                            df_novoplasty,
                            df_mitos2,
                            left_on=["Sample", "Seed", "kmer"],
                            right_on=["sample", "seed", "kmer"],
                        )

                        columns_to_keep = [
                            "Sample",
                            "Seed",
                            "kmer",
                            "assembly",
                            "transporter_count",
                            "ribossomal_count",
                            "origins_count",
                            "coding_count",
                        ]
                        df_abstract = df_abstract[columns_to_keep]

                        df_pilon = pd.read_csv(f"workflow/reports/{sample}/pilon.csv")

                        columns_to_keep = ["assembly", "genome_size", "changes_number"]
                        df_pilon = df_pilon[columns_to_keep]

                        df_abstract = pd.merge(
                            df_abstract,
                            df_pilon,
                            left_on="assembly",
                            right_on="assembly",
                        )
                        df_abstract.to_csv(
                            f"workflow/reports/{sample}/abstract.csv", index=False
                        )
