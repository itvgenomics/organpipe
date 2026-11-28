import re
import sys
import yaml

sys.path.insert(0, os.path.join(workflow.basedir, "scripts"))
from naming import (
    GETORGANELLE_RUN, START_GENE, find_long_reads, find_short_reads, find_single_reads, genbank_identity, sif_name,
)

# Container images, pulled by workflow/scripts/singularity.py.
with open(os.path.join(workflow.basedir, "containers.yaml")) as handle:
    CONTAINERS = yaml.safe_load(handle)

# Name of a MitoHiFi result in the paths shared with the assemblies of the other assemblers.
MITOHIFI_ASSEMBLY = "final_mitogenome"

# Sequences searched with nhmmscan for each assembly.
NHMMER_REGIONS = ["rRNA-tRNA", "intergenes_filter"]

# Suffix of the GenBank files of each annotation tool.
ANNOTATOR_SUFFIX = {"mitos2": "", "chloe": ".chloe", "cpgavas2": ".cpgavas2"}

def sif(name):
    """Path of an image in workflow/containers.yaml."""
    return f"{config['sif_dir']}/{sif_name(name, CONTAINERS[name])}"

def sample_param(field):
    """Function for `params:` that returns a config field of the sample."""
    return lambda wildcards: config["samples"][wildcards.sample][field]

def is_yes(sample, field):
    return config["samples"][sample][field] == "yes"

def single_end(sample):
    """Short reads in one file of single-end or merged reads (read_layout: single), assembled with GetOrganelle -u."""
    return config["samples"][sample].get("read_layout", "paired") == "single"

def only_samples(samples):
    """Wildcard constraint that matches only these samples, and no sample when the list is empty."""
    return "|".join(re.escape(sample) for sample in samples) or "(?!)"

# Paired and single-end samples have their own trimming rule, and both write the same fastp report.
PAIRED_SAMPLES = only_samples([sample for sample in SAMPLES if not single_end(sample)])
SINGLE_END_SAMPLES = only_samples([sample for sample in SAMPLES if single_end(sample)])

def output_dir(wildcards, output):
    """Function for `params:` that returns the directory of the first output."""
    return os.path.dirname(output[0])

def annotators(sample):
    """Annotation tools of a sample: MITOS2 for mitochondria, and those of run_cpgavas2 and run_chloe for chloroplasts.

    The first one annotates the sequences of the plots and of the nhmmer search: CPGAVAS2 when it runs.
    """
    if config["samples"][sample]["organelle"] == "mito":
        return ["mitos2"]
    return [annotator for annotator in ["cpgavas2", "chloe"] if is_yes(sample, f"run_{annotator}")]

def genbank_path(wildcards, annotator, rotated=False):
    suffix = ANNOTATOR_SUFFIX[annotator] + (".rotated" if rotated else "")
    return (
        f"results/{wildcards.sample}/genbanks/{wildcards.assembler}/{wildcards.run}/"
        f"{wildcards.assembly}{suffix}.gb"
    )

def genbank_name(annotator=None, rotated=False):
    """Function for `params:` with the LOCUS name and DEFINITION of a GenBank file.

    Returns {"name": ..., "definition": ...}. Without annotator, the {annotator} wildcard is used.
    """
    def identity(wildcards):
        name, definition = genbank_identity(
            wildcards.sample,
            config["samples"][wildcards.sample]["organelle"],
            wildcards.assembler,
            wildcards.run,
            wildcards.assembly,
            annotator or wildcards.annotator,
            rotated,
        )
        return {"name": name, "definition": definition}
    return identity

def short_reads(wildcards):
    """Trimmed or raw reads of a sample, for the assemblers: r1 and r2, or single with read_layout: single."""
    prefix = f"resources/{wildcards.sample}/rawreads/{wildcards.sample}"
    trimmed = is_yes(wildcards.sample, "run_trimming")
    if single_end(wildcards.sample):
        return {"single": f"{prefix}.SE.trimmed.gz" if trimmed else f"{prefix}_SE.fastq.gz"}
    if trimmed:
        return {"r1": f"{prefix}.R1.trimmed.gz", "r2": f"{prefix}.R2.trimmed.gz"}
    return {"r1": f"{prefix}_R1.fastq.gz", "r2": f"{prefix}_R2.fastq.gz"}

def assembled_reads(wildcards):
    """Reads each assembler kept for a run, for error correction and the plots.

    `reads` is the pair of files, or one file for single-end reads: bwa-mem2 takes either.
    """
    if wildcards.assembler == "getorganelle":
        prefix = f"results/{wildcards.sample}/getorganelle/extended"
        if single_end(wildcards.sample):
            # With -u alone, GetOrganelle writes the reads it kept to extended_1.fq
            return {"reads": [f"{prefix}_1.fq"]}
        return {"reads": [f"{prefix}_1_paired.fq", f"{prefix}_2_paired.fq"]}

    seed, kmer = wildcards.run.rsplit("_kmer", 1)
    prefix = f"results/{wildcards.sample}/novoplasty/{seed}/kmer{kmer}/Assembled_reads_{wildcards.sample}"
    return {"reads": [f"{prefix}_R1.fasta", f"{prefix}_R2.fasta"]}
