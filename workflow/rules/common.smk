import re
import sys
import yaml

sys.path.insert(0, os.path.join(workflow.basedir, "scripts"))
from naming import GETORGANELLE_RUN, find_long_reads, find_short_reads, sif_name

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

def output_dir(wildcards, output):
    """Function for `params:` that returns the directory of the first output."""
    return os.path.dirname(output[0])

def annotators(sample):
    """Annotation tools of a sample."""
    return ["mitos2"] if config["samples"][sample]["organelle"] == "mito" else ["cpgavas2", "chloe"]

def genbank_path(wildcards, annotator, rotated=False):
    suffix = ANNOTATOR_SUFFIX[annotator] + (".rotated" if rotated else "")
    return (
        f"results/{wildcards.sample}/genbanks/{wildcards.assembler}/{wildcards.run}/"
        f"{wildcards.assembly}{suffix}.gb"
    )

def short_reads(wildcards):
    """Trimmed or raw reads of a sample, for the assemblers."""
    prefix = f"resources/{wildcards.sample}/rawreads/{wildcards.sample}"
    if is_yes(wildcards.sample, "run_trimming"):
        return {"r1": f"{prefix}.R1.trimmed.gz", "r2": f"{prefix}.R2.trimmed.gz"}
    return {"r1": f"{prefix}_R1.fastq.gz", "r2": f"{prefix}_R2.fastq.gz"}

def assembled_reads(wildcards):
    """Reads each assembler kept for a run, for error correction and the plots."""
    if wildcards.assembler == "getorganelle":
        prefix = f"results/{wildcards.sample}/getorganelle/extended"
        return {"r1": f"{prefix}_1_paired.fq", "r2": f"{prefix}_2_paired.fq"}

    seed, kmer = wildcards.run.rsplit("_kmer", 1)
    prefix = f"results/{wildcards.sample}/novoplasty/{seed}/kmer{kmer}/Assembled_reads_{wildcards.sample}"
    return {"r1": f"{prefix}_R1.fasta", "r2": f"{prefix}_R2.fasta"}
