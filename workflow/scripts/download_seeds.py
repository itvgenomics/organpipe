import argparse
import os
import logging
import sys
import textwrap
import time
from io import StringIO
from Bio import Entrez
from Bio import SeqIO


# Create a Function to read the kmers list
def parse_csv_list(gene_list):
    try:
        return [str(x) for x in gene_list.split(",")]
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"{gene_list} is not a comma-separated list of integers"
        )


def parse_arguments():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(
            """
            This Python script is designed to search the NCBI database and download seed sequences to use in a PipePlasty run. It retrieves genetic information in the form of FASTA files based on specified taxon and gene list.

            Command-line Arguments:
                --taxon: Specify the taxon to search.
                --outpath: Set the directory to write output files.
                --maxcount: Number of FASTA files to download for each gene (default: 3).
                --genes: Comma-separated list of genes to search for.
                --organelle: Organelle type [mito or chloro]

            Exemple:

            python findSeed.py --taxon Amphisbaena --outpath /path/to/output --maxcount 5 --genes COI,ATP6,16S,CYTB --organelle mito

            Output:

            The script generates multiple FASTA files for each gene and organism combination. The files are stored in the specified output directory.

                - <gene>_<organism>.fasta: Individual FASTA files for each gene and organism combination.
                - seeds.fasta: Aggregated file containing all downloaded seed sequences.

            Ensure that the specified taxon and gene names are accurate for successful search results.
            """
        ),
    )

    parser.add_argument("--taxon", required=True, help="Taxon to search")
    parser.add_argument(
        "--outpath", required=True, help="The directory to write output files"
    )
    parser.add_argument(
        "--maxcount",
        type=int,
        default=3,
        help="Number of fasta files to download for each gene",
    )

    parser.add_argument(
        "--genes", type=parse_csv_list, help="Comma Separated Gene List", required=True
    )
    parser.add_argument(
        "--organelle",
        choices=['"mito"', '"chloro"'],
        help="Organelle Type",
        required=True,
    )

    args = parser.parse_args()

    return args


def fetch_gb(db, id_str, rettype, retmode, max_retries=5, delay=5):
    for attempt in range(max_retries):
        try:
            handle = Entrez.efetch(db=db, id=id_str, rettype=rettype, retmode=retmode)
            data = handle.read()
            handle.close()
            return data
        except Exception as e:
            logging.warning(f"Attempt {attempt + 1} failed with error: {e}")
            if attempt < max_retries - 1:
                time.sleep(delay)
            else:
                raise


def find_seed(organism, gene, max_n, outpath, organelle):
    organisms = []

    n = 0

    if organelle == "mito":
        term = (
            f"({organism}[Organism] AND {gene}[All Fields] AND mitochondrion[filter])"
        )
    elif organelle == "chloro":
        term = f"({organism}[Organism] AND {gene}[All Fields] AND chloroplast[filter])"

    logging.info(f"NCBI Nucleotide Query: {term}")
    logging.info("Getting IDs...")

    handle = Entrez.esearch(db="nucleotide", term=term, idtype="acc", retmax=250)
    id_record = Entrez.read(handle)
    handle.close()

    id_list = id_record["IdList"]

    if len(id_list) > 248:
        logging.info(f"Found 250+ IDs.")
    else:
        logging.info(f"Found {len(id_list)} IDs.")

    id_str = ",".join(id_list)

    logging.info(f"Getting GB files...")

    gb_record = fetch_gb(db="nucleotide", id_str=id_str, rettype="gb", retmode="text")

    if gb_record:
        logging.info(f"Writing FASTA files.")

        with open(f"{outpath}/{gene}.gb", "w") as gb_file:
            for record in gb_record:
                gb_file.write(record)

        for record in SeqIO.parse(
            f"{outpath}/{gene}.gb",
            "genbank",
        ):
            organism_name = record.annotations.get("organism", "")
            sequence = str(record.seq)

            if organism_name not in organisms:
                organisms.append(organism_name)

                organism_name = (
                    organism_name.replace(" ", "_")
                    .replace("/", ".")
                    .replace(":", ".")
                    .replace("'", "")
                    .replace("(", "")
                    .replace(")", "")
                    .replace("_aff.", "")
                )

                fasta_file = f"{outpath}/{gene}_{organism_name}.fasta"

                if len(sequence) > 3000:
                    for feature in record.features:
                        gene_name = str(feature.qualifiers.get("gene"))

                        if feature.type == "gene" and (
                            gene.lower() in gene_name or gene.upper() in gene_name
                        ):
                            gene_sequence = str(feature.location.extract(record).seq)
                            with open(fasta_file, "w") as output:
                                output.write(f">{gene}_{organism_name}\n")
                                output.write(gene_sequence)

                            break
                else:
                    with open(fasta_file, "w") as output:
                        output.write(f">{gene}_{organism_name}\n")
                        output.write(sequence)

                n += 1

            if n == max_n:
                return logging.info(
                    f"Got {len(organisms)} distinct organisms and {n} sequences for gene {gene}"
                )

    return logging.info(
        f"Got {len(organisms)} distinct organisms and {n} sequences for gene {gene}"
    )


def main(taxon, outpath, maxcount, genes, organelle):
    # args = parse_arguments()

    Entrez.email = "ncbi.itvds@pq.itv.org"
    Entrez.api_key = "f178b74aee7a97072b0dc4e10fc75b3b5508"

    FORMAT = "%(asctime)s [%(levelname)s] %(message)s"

    logging.basicConfig(
        level=logging.INFO,
        stream=sys.stdout,
        format=FORMAT,
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    for gene in genes:
        logging.info(f"Running search for gene: {gene}.")
        find_seed(f"{taxon}", gene, maxcount, outpath, organelle)

    logging.info("Creating the seeds.fasta file.")

    # Iterate over each file in the input folder
    for filename in os.listdir(outpath):
        if filename.endswith(".fasta"):
            filepath = os.path.join(outpath, filename)
            # Append the contents of the current file to the output file
            with open(filepath, "r") as infile, open(
                f"{outpath}/seeds.fasta", "a"
            ) as outfile:
                records = SeqIO.read(infile, "fasta")
                SeqIO.write(records, outfile, "fasta")

    logging.info("Done!")
