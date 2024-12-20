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

    single_reads = [file for file in files_list if sample in file]
    return {"single_reads": single_reads[0]}

rule check_shortreads_files:
    input: 
        forward_reads=lambda wildcards: get_shortreads(wildcards.sample)["forward"],
        reverse_reads=lambda wildcards: get_shortreads(wildcards.sample)["reverse"]
    output:
        r1 = temp("resources/{sample}/rawreads/{sample}_R1.fastq.gz"),
        r2 = temp("resources/{sample}/rawreads/{sample}_R2.fastq.gz")
    threads: 1
    shell:
        """
        cp {input.forward_reads} {output.r1}
        cp {input.reverse_reads} {output.r2}
        """

rule check_longreads_files:
    input: 
        single_reads=lambda wildcards: get_longreads(wildcards.sample)["single_reads"],
    output:
        temp("resources/{sample}/rawreads/{sample}.fasta")
    threads: 1
    shell:
        """
        cp {input.single_reads} {output}
        """

rule check_adapters:
    input: 
        lambda wildcards: config["samples"][wildcards.sample]["adapters"]
    output:
        "resources/{sample}/adapters.txt"
    threads: 1
    shell:
        """
        cp {input} {output}
        """

rule download_cpgavas_sif:
    output:
        "resources/cpgavas2.sif"
    threads: 1
    shell:
        """
        wget https://xxjzna.bn.files.1drv.com/y4miznP-svjCz0-AfKjAUT9u6C4Idv4mDba2TiWspJt6vaFLDHo-D9hiZ3wGotg2nh3dGdnWELOtuQxE6eTmyVOFXLSd1HqvG1ANEpzAj_kxtEbKISeF_Kx4SMvhT4HEpCSM8IU3AOe-Iiw0SpNcU6PL5OLqq65u6JdDhJpvLa5WB29_6PaXYKaDO5kQZFSfKjLdlJCHgPodkgcZZ7KwgCAWg -O {output}
        """


def build_mitohifi():
    command = "singularity build resources/mitohifi.sif docker://ghcr.io/marcelauliano/mitohifi:master"
    subprocess.run(command, shell=True, check=True)

rule build_mitohifi_sif:
    output:
        "resources/mitohifi.sif"
    threads: 1
    run:
        build_mitohifi()

rule check_nhmmer_db:
    output:
        "resources/nhmmer_db.hmm"
    threads: 1
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
    threads: 1
    shell:
        """
        cd resources && gunzip table2asn.linux64.gz && chmod +x table2asn.linux64
        """

rule index_hmmer_db:
    input:
        "resources/nhmmer_db.hmm"
    output:
        expand("resources/nhmmer_db.hmm.{ext}", ext=['h3f', 'h3i', 'h3m', 'h3p'])
    threads: 1
    singularity:
        "docker://staphb/hmmer:3.4"
    shell:
        """
        hmmpress {input}
        """