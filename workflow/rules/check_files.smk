import os
import subprocess

def get_shortreads(sample):
    files_list = []
    for root, dirs, files in os.walk(config["samples"][sample]["reads_path"]):
        for file in files:
            file_path = os.path.join(root, file)
            files_list.append(file_path)

    reads = [file for file in files_list if sample in file]
    forward_reads = [file for file in reads if "_R1" in file or "pair1" in file]
    reverse_reads = [file for file in reads if "_R2" in file or "pair2" in file]
    return {"forward": forward_reads[0], "reverse": reverse_reads[0]}

def get_longreads(sample):
    files_list = []
    for root, dirs, files in os.walk(config["samples"][sample]["reads_path"]):
        for file in files:
            file_path = os.path.join(root, file)
            files_list.append(file_path)

    long_reads = [file for file in files_list if sample in file]
    return {"long_reads": long_reads[0]}

rule check_shortreads_files:
    input:
        forward_reads=lambda wildcards: get_shortreads(wildcards.sample)["forward"],
        reverse_reads=lambda wildcards: get_shortreads(wildcards.sample)["reverse"]
    output:
        r1 = temp("resources/{sample}/rawreads/{sample}_R1.fastq.gz"),
        r2 = temp("resources/{sample}/rawreads/{sample}_R2.fastq.gz")
    shell:
        """
        cp {input.forward_reads} {output.r1}
        cp {input.reverse_reads} {output.r2}
        """

rule check_fastq_longreads:
    input:
        long_reads=lambda wildcards: get_longreads(wildcards.sample)["long_reads"],
    output:
        temp("resources/{sample}/rawreads/{sample}.fastq.gz")
    shell:
        """
        if [ -f {input} ]; then
            cp {input} {output}
        else
            touch {output}
        fi
        """

rule check_fasta_longreads:
    input:
        long_reads=lambda wildcards: get_longreads(wildcards.sample)["long_reads"],
    output:
        temp("resources/{sample}/rawreads/{sample}.fasta")
    shell:
        """
        if [ -f {input} ]; then
            cp {input} {output}
        else
            touch {output}
        fi
        """

rule check_adapters:
    input:
        lambda wildcards: f"resources/{wildcards.sample}/adapters.fasta" if config["samples"][wildcards.sample].get("adapters", "") else "test_data/adapters.fasta"
    output:
        "resources/{sample}/adapters.fasta"
    shell:
        """
        if [ -f {input} ]; then
            cp {input} {output}
        else
            touch {output}
        fi
        """

rule check_nhmmer_db:
    output:
        "resources/nhmmer_db.hmm"
    shell:
        """
        NHMMER_DB=$(grep 'nhmmer_db:' config/snakemake_config.yaml | grep -v "''" | sed -E 's/.*: *//; s/[[:space:]]*$//' | uniq) && \
        if [[ ! -f "{output}" ]]; then \
            cp "$NHMMER_DB" {output}
        fi
        """

rule extract_table2asn:
    input:
        "resources/table2asn.linux64.gz"
    output:
        "resources/table2asn.linux64"
    shell:
        """
        cd resources && gunzip table2asn.linux64.gz && chmod +x table2asn.linux64
        """

rule index_hmmer_db:
    input:
        "resources/nhmmer_db.hmm"
    output:
        expand("resources/nhmmer_db.hmm.{ext}", ext=['h3f', 'h3i', 'h3m', 'h3p'])
    singularity:
        f"{config["sif_dir"]}/hmmer.sif"
    shell:
        """
        hmmpress {input}
        """

rule check_reference:
    input:
        lambda wildcards: config["samples"][wildcards.sample]["reference"] if config["samples"][wildcards.sample].get("reference", "") else "test_data/reference.fasta"
    output:
        "resources/{sample}/reference.fasta"
    shell:
        """
        if [ -f {input} ]; then
            cp {input} {output}
        else
            touch {output}
        fi
        """
