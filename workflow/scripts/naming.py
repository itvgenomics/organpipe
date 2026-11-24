"""File name conventions and warnings shared by the rules and the scripts."""

import os
import re
import warnings

SHORT_READ_SUFFIXES = (".fastq.gz", ".fasta.gz")
LONG_READ_SUFFIXES = (".fastq.gz", ".fasta")

# Read pair marker in a file name: sample_R1.fastq.gz, sample_S1_R2_001.fastq.gz...
MATE_MARKER = re.compile(r"_(R[12]|pair[12])(?=[._])")

# parse_results.py collects the log lines that start with this into warnings.txt.
WARNING_PREFIX = "OrganPipe WARNING:"

# A run is one execution of an assembler, and it can produce several assemblies.
# NOVOPlasty has one run per seed and k-mer ({seed}_kmer{kmer}); GetOrganelle has one per sample.
GETORGANELLE_RUN = "all"


# ── GenBank names ──────────────────────────────────────────────────────────────
# The LOCUS name is the document name in Geneious and similar tools, and the
# DEFINITION is its description. The GenBank format accepts LOCUS names longer than
# 16 characters since release 229 (the other fields move to the right). Above this
# length only the seed name is shortened, keeping its number, so names stay unique.
LOCUS_MAX_LENGTH = 50

ORGANELLE_NAME = {"mito": "mitochondrion", "chloro": "chloroplast"}
ASSEMBLER_NAME = {"novoplasty": "NOVOPlasty", "getorganelle": "GetOrganelle", "mitohifi": "MitoHiFi"}
ANNOTATOR_NAME = {"mitos2": "MITOS2", "chloe": "Chloe", "cpgavas2": "CPGAVAS2"}
# First gene of the rotated GenBank files (rotate_genbank.py)
START_GENE = {"mito": "rrnS", "chloro": "psbA"}

# NOVOPlasty assembly: {status}_{n}_{sample}_{seed}_{kmer}_{index} (README, "Assembly names")
NOVOPLASTY_RUN = re.compile(r"(?P<seed>.+)_kmer(?P<kmer>\d+)")


def locus_text(text):
    """Keep only the characters allowed in a LOCUS name."""
    return re.sub(r"[^A-Za-z0-9_.-]", "_", text)


def shorten_seed(seed, excess):
    """Remove excess characters from the end of a seed name, keeping its number ({number}-)."""
    keep = len(seed.split("-", 1)[0]) + 1 if "-" in seed else 1
    return seed[: max(len(seed) - excess, keep)].rstrip("_.")


def novoplasty_parts(sample, run, assembly):
    """Seed, k-mer and short assembly tag ({status}{n}_{index}) of a NOVOPlasty assembly."""
    match = NOVOPLASTY_RUN.fullmatch(run)
    if not match:
        return None
    seed, kmer = match.group("seed"), match.group("kmer")
    head, found, index = assembly.partition(f"_{sample}_{seed}_{kmer}_")
    status = re.fullmatch(r"([A-Za-z]+)_(\d+)", head) if found else None
    if not status:
        return None
    return seed, kmer, f"{status.group(1)}{status.group(2)}_{index}"


def genbank_identity(sample, organelle, assembler, run, assembly=None, annotator=None, rotated=False):
    """LOCUS name and DEFINITION of the GenBank file of an assembly.

    NOVOPlasty:   {sample}_{seed}_k{kmer}_{status}{n}_{index}, e.g. DNA1_1-rbcL_1858743_k33_OPT1_1
    GetOrganelle: {sample}_{assembly},                         e.g. DNA1_embplant_pt_complete_1
    MitoHiFi:     {sample}_{reference}_MitoHiFi
    The chloroplast annotation tool (_chloe, _cpgavas2) and _rot (rotated) are added at the end.
    """
    suffix = ""
    if annotator in ("chloe", "cpgavas2"):
        suffix += f"_{annotator}"
    if rotated:
        suffix += "_rot"

    definition = f"{sample} {ORGANELLE_NAME.get(organelle, organelle)}, {ASSEMBLER_NAME.get(assembler, assembler)} assembly"
    parts = novoplasty_parts(sample, run, assembly) if assembler == "novoplasty" else None

    if parts:
        seed, kmer, tag = parts
        rest = f"_k{kmer}_{tag}{suffix}"
        name = locus_text(f"{sample}_{seed}{rest}")
        excess = len(name) - LOCUS_MAX_LENGTH
        if excess > 0:
            name = locus_text(f"{sample}_{shorten_seed(seed, excess)}{rest}")
        definition += f" {assembly} (seed {seed}, k-mer {kmer})"
    elif assembler == "mitohifi":
        name = locus_text(f"{sample}_{run}_MitoHiFi{suffix}")
        definition += f" (reference {run})"
    else:
        name = locus_text(f"{sample}_{assembly}{suffix}")
        definition += f" {assembly}"

    if annotator in ANNOTATOR_NAME:
        definition += f", annotated with {ANNOTATOR_NAME[annotator]}"
    if rotated:
        definition += f", rotated to start at {START_GENE.get(organelle, 'the start gene')}"
    return name, definition


def set_genbank_identity(record, name, definition):
    """Write the LOCUS name and DEFINITION into a Biopython record."""
    record.name = name
    record.id = name
    record.description = definition
    return record


def write_genbank(record, path):
    """Write a GenBank file, without Biopython's warning about LOCUS names over 16 characters."""
    from Bio import SeqIO

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Increasing length of locus line")
        SeqIO.write(record, path, "genbank")


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


def find_read_file(sample, reads_path, suffixes):
    """Return the only read file of a sample whose reads are in one file (no _R1/_R2 marker)."""
    exact = []
    partial = []

    for filename in sorted(os.listdir(reads_path)):
        if mate_marker(filename) or not filename.endswith(suffixes):
            continue

        stem = read_stem(filename)
        if stem == sample:
            exact.append(filename)
        elif has_sample_name(sample, stem):
            partial.append(filename)

    return pick_read_file(sample, "read file", reads_path, exact, partial)


def find_long_reads(sample, reads_path):
    """Return the read file of a long-read sample."""
    return find_read_file(sample, reads_path, LONG_READ_SUFFIXES)


def find_single_reads(sample, reads_path):
    """Return the read file of a short-read sample with read_layout: single (single-end or merged reads)."""
    return find_read_file(sample, reads_path, SHORT_READ_SUFFIXES)
