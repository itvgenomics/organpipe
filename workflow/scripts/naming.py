"""File name conventions and warnings shared by the rules and the scripts."""

import os
import re

SHORT_READ_SUFFIXES = (".fastq.gz", ".fasta.gz")
LONG_READ_SUFFIXES = (".fastq.gz", ".fasta")

# Read pair marker in a file name: sample_R1.fastq.gz, sample_S1_R2_001.fastq.gz...
MATE_MARKER = re.compile(r"_(R[12]|pair[12])(?=[._])")

# parse_results.py collects the log lines that start with this into warnings.txt.
WARNING_PREFIX = "OrganPipe WARNING:"

# A run is one execution of an assembler, and it can produce several assemblies.
# NOVOPlasty has one run per seed and k-mer ({seed}_kmer{kmer}); GetOrganelle has one per sample.
GETORGANELLE_RUN = "all"


def sif_name(name, image):
    """File name of a container image: <name>_<tag>.sif"""
    tag = image.rsplit(":", 1)[1] if ":" in image.split("/")[-1] else "latest"
    return f"{name}_{re.sub(r'[^A-Za-z0-9._-]', '_', tag)}.sif"


def warn(message):
    """Report a problem that does not stop the run but should reach the user."""
    print(f"{WARNING_PREFIX} {message}", flush=True)


def has_sample_name(sample, name):
    """True if sample is a whole name inside name: hs1 is in run_hs1_S1, not in hs10."""
    return re.search(rf"(?<![A-Za-z0-9]){re.escape(sample)}(?![A-Za-z0-9])", name) is not None


def mate_marker(filename):
    """Return the last read pair marker of a file name, or None."""
    markers = list(MATE_MARKER.finditer(filename))
    return markers[-1] if markers else None


def read_stem(filename):
    """Return a read file name without its extension."""
    return re.sub(r"\.(fastq\.gz|fasta\.gz|fasta)$", "", filename)


def pick_read_file(sample, description, reads_path, exact, partial):
    """Return the only matching file. A file named exactly after the sample wins."""
    files = exact or partial
    if len(files) == 1:
        return os.path.join(reads_path, files[0])
    if not files:
        raise ValueError(f"No {description} found for sample '{sample}' in '{reads_path}'.")
    raise ValueError(
        f"More than one {description} found for sample '{sample}' in '{reads_path}': "
        f"{', '.join(files)}. Keep one file per sample, or merge them first."
    )


def find_short_reads(sample, reads_path):
    """Return the forward and reverse read files of a sample."""
    exact = {"1": [], "2": []}
    partial = {"1": [], "2": []}

    for filename in sorted(os.listdir(reads_path)):
        marker = mate_marker(filename)
        if not marker or not filename.endswith(SHORT_READ_SUFFIXES):
            continue

        prefix = filename[: marker.start()]
        mate = marker.group(1)[-1]
        if prefix == sample:
            exact[mate].append(filename)
        elif has_sample_name(sample, prefix):
            partial[mate].append(filename)

    return {
        "forward": pick_read_file(sample, "forward read file", reads_path, exact["1"], partial["1"]),
        "reverse": pick_read_file(sample, "reverse read file", reads_path, exact["2"], partial["2"]),
    }


def find_long_reads(sample, reads_path):
    """Return the read file of a long-read sample."""
    exact = []
    partial = []

    for filename in sorted(os.listdir(reads_path)):
        if mate_marker(filename) or not filename.endswith(LONG_READ_SUFFIXES):
            continue

        stem = read_stem(filename)
        if stem == sample:
            exact.append(filename)
        elif has_sample_name(sample, stem):
            partial.append(filename)

    return pick_read_file(sample, "read file", reads_path, exact, partial)
