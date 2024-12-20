import seaborn as sns
from matplotlib import pyplot as plt, lines
import argparse
import os
from Bio import SeqIO
import pandas as pd
import warnings


def parse_arguments():
    # Set up argument parsing
    parser = argparse.ArgumentParser(
        description="Convert GenBank files to FASTA format."
    )
    parser.add_argument(
        "--parse_gb", action="store_true", help="Trigger GenBank to FASTA conversion."
    )
    parser.add_argument(
        "--recruitment_plot", action="store_true", help="Trigger Recruitment Plot."
    )
    parser.add_argument("--depth", action="store_true", help="Trigger Depth Plot.")
    parser.add_argument(
        "--sample",
        help="Sample dir",
        required=False,
        nargs="?",
    )
    parser.add_argument(
        "--kmer",
        help="kmer used",
        nargs="?",
    )
    parser.add_argument(
        "--seed",
        help="Seed used",
        nargs="?",
    )
    parser.add_argument(
        "--organelle",
        help="Type of organelle (e.g., mito or chloro)",
    )
    parser.add_argument(
        "--blastn_fasta",
        help="Fasta file used to perform makeblastdb",
        nargs="?",
    )
    parser.add_argument(
        "--blastn_out",
        help="Blastn output tabular outfmt6 file",
        nargs="?",
    )
    parser.add_argument(
        "--depth_bam",
        help="Depth output file from samtools depth",
        nargs="?",
    )

    args = parser.parse_args()

    return args


def count_nucleotides(fasta_file):
    total_count = 0

    with open(fasta_file, "r") as file:
        for record in SeqIO.parse(file, "fasta"):
            sequence = record.seq
            total_count += len(sequence)

    return total_count


def plot_recruitment_plot(
    blast_tabular_output,
    output_plot,
    plot_title,
    genome_size,
    minimum_identity_displayed=70,
    min_alignment_length=100,
    identity_species_cutoff=95,
):
    """Parse BLAST-like tabular output and create a recruitment plot for the alignments.

    Args:
        blast_tabular_output (str): Path to input BLAST-like alignment file.
        output_plot (str): Path to where the image output in PNG will be saved.
        plot_title (str): Plot title (e.g., Genome Name).
        genome_size (int): Genome size.
        minimum_identity_displayed (float): Y-axis lower bound (default % identity 70).
        min_alignment_length (int): Minimum alignment length allowed (default = 100 bp).
        identity_species_cutoff (float): Species line draw on plot (default = 95).
    """
    # Set Seaborn theme for aesthetics
    sns.set_theme(color_codes=True)

    # Define column names for the BLAST-like tabular input
    columns = [
        "query_id",
        "subject_id",
        "identity",
        "alignment_length",
        "mismatches",
        "gap_opens",
        "q_start",
        "q_end",
        "s_start",
        "s_end",
        "e_value",
        "bit_score",
    ]

    # Read the tabular file into a Pandas DataFrame
    data = pd.read_csv(blast_tabular_output, sep="\t", header=None, names=columns)

    # Filter the data based on alignment length and identity thresholds
    data = data[
        (data["alignment_length"] >= min_alignment_length)
        & (data["identity"] >= minimum_identity_displayed)
    ]

    # Sort the start and stop positions for consistent plotting
    data["start"] = data[["s_start", "s_end"]].min(axis=1)
    data["stop"] = data[["s_start", "s_end"]].max(axis=1)

    # Initialize the plot
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.set_xlim(0, genome_size + 1)
    ax.set_ylim(minimum_identity_displayed, 100)

    # Plot each alignment as a line
    for _, row in data.iterrows():
        ax.add_line(
            lines.Line2D(
                [row["start"], row["stop"]],
                [row["identity"], row["identity"]],
                color="#89896E",
            )
        )

    # Add species cutoff line
    ax.add_line(
        lines.Line2D(
            [0, genome_size + 1],
            [identity_species_cutoff, identity_species_cutoff],
            color="k",
        )
    )

    # Set axis labels and title
    ax.set_ylabel("% Identity")
    ax.set_xlabel("Genome Position (bp)")
    ax.set_title(plot_title)

    # Save the plot
    fig.savefig(output_plot, format="png", dpi=400)
    plt.close()


def parse_genbank_to_fasta(input_file, fasta_dir):
    # Parse the GenBank file
    with open(input_file, "r") as gb_file:
        for record in SeqIO.parse(gb_file, "genbank"):
            file_basename = os.path.basename(input_file)
            fasta_header = os.path.splitext(file_basename)[0]
            # Construct the output FASTA filename
            fasta_filename = file_basename.replace(".gb", ".fasta")

            record.id = fasta_header
            record.name = fasta_header
            record.description = ""  # Optional: leave the description blank

            # Write the sequence to a FASTA file
            with open(os.path.join(fasta_dir, fasta_filename), "w") as fasta_file:
                SeqIO.write(record, fasta_file, "fasta")


# The following code is found in https://onestopdataanalysis.com/depth-plot/
def parse_depth(depth_input, genome_size):
    """Parse depth file.

    Args:
        depth_input (str): Path to depth file.
        genome_size (int): Genome size.

    Returns:
        list: List with depth.

    """
    depth = [0] * genome_size
    references = set()

    with open(depth_input) as depth_object:
        for row in depth_object:
            genome_id, position, depth_count = row.split()

            references.add(genome_id)

            if len(references) > 1:
                raise Exception(" This script only handles one genome - contig.")

            depth[int(position)] = int(depth_count)

    return depth


def plot_depth(depth_report, output_name, plot_title, genome_size, normalize=False):
    """Plot genome Depth across genome.

    Args:
        depth_report (str): Path to samtool's depth file.
        output_name (str): Path to output PNG image.
        plot_title (str): Plot title.
        genome_size (int): Genome size.
        normalize (bool): If `True`, normalizes the depth by the largest depth (default = `False`).
        depth_cut_off (int): Plot a line to represent a targeted depth (default = 20).

    """
    data = parse_depth(depth_report, genome_size)

    y_label = "Normalized Depth" if normalize else "Depth"
    data = [xx / max(data) for xx in data] if normalize else data

    sns.set_theme(color_codes=True)
    plt.title(plot_title)
    ax = plt.subplot(111)

    data = data[1:]

    plt.plot(range(len(data)), data)
    plt.xlabel("Genome Position (bp)")
    plt.ylabel(y_label)

    plt.savefig(output_name, bbox_inches="tight", dpi=400)
    plt.close()


if __name__ == "__main__":
    args = parse_arguments()

    sample = args.sample
    seed = args.seed
    kmer = args.kmer
    organelle = args.organelle
    blastn_fasta = args.blastn_fasta
    blastn_out = args.blastn_out
    depth_bam = args.depth_bam

    if args.parse_gb:
        genbank_dir = f"results/{sample}/genbanks"
        output_dir = f"results/{sample}/images"

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        for genbank in os.listdir(genbank_dir):
            if organelle == "mito":
                if genbank.endswith(".gb") and seed in genbank and kmer in genbank:
                    fasta_dir = f"results/{sample}/images/{seed}_kmer{kmer}"

                    if not os.path.exists(fasta_dir):
                        os.makedirs(fasta_dir)

                    parse_genbank_to_fasta(
                        os.path.join(genbank_dir, genbank), fasta_dir
                    )
            elif organelle == "chloro":
                if (
                    genbank.endswith(".gb")
                    and seed in genbank
                    and kmer in genbank
                    and "cpgavas2" in genbank
                ):
                    fasta_dir = f"results/{sample}/images/{seed}_kmer{kmer}"

                    if not os.path.exists(fasta_dir):
                        os.makedirs(fasta_dir)

                    parse_genbank_to_fasta(
                        os.path.join(genbank_dir, genbank), fasta_dir
                    )

    if args.recruitment_plot:
        try:
            genomesize = count_nucleotides(blastn_fasta)

            plot_recruitment_plot(
                blast_tabular_output=blastn_out,
                output_plot=blastn_out.replace(".blastn.out", "_recruitment_plot.png"),
                plot_title="Recruitment Plot",
                genome_size=int(genomesize),
                minimum_identity_displayed=80,
                min_alignment_length=100,
                identity_species_cutoff=95,
            )
        except Exception as e:
            warnings.warn(
                f"An error occurred: {e}. The script will continue running.",
                RuntimeWarning,
            )

    if args.depth:
        if organelle == "mito":
            try:
                fasta_file = depth_bam.replace(".depth", ".fasta")
                genomesize = count_nucleotides(fasta_file)

                plot_depth(
                    depth_report=depth_bam,
                    output_name=f"{depth_bam}.png",
                    plot_title="Depth Plot",
                    genome_size=int(genomesize) + 1,
                    normalize=False,
                )
            except Exception as e:
                warnings.warn(
                    f"An error occurred: {e}. The script will continue running.",
                    RuntimeWarning,
                )

        elif organelle == "chloro":
            print("depth_bam", depth_bam)

            if ".cpgavas2.rotated." in depth_bam:
                fasta_file = depth_bam.replace(".depth", ".fasta")
            else:
                fasta_file = depth_bam.replace(".depth", ".cpgavas2.fasta")

            try:
                genomesize = count_nucleotides(fasta_file)

                print("fasta_file", fasta_file)

                plot_depth(
                    depth_report=depth_bam,
                    output_name=f"{depth_bam}.png",
                    plot_title="Depth Plot",
                    genome_size=int(genomesize) + 1,
                    normalize=False,
                )
            except Exception as e:
                warnings.warn(
                    f"An error occurred: {e}. The script will continue running.",
                    RuntimeWarning,
                )
