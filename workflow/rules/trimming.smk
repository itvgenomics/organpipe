rule run_fastp:
    input:
        r1 = "resources/{sample}/rawreads/{sample}_R1.fastq.gz",
        r2 = "resources/{sample}/rawreads/{sample}_R2.fastq.gz",
        adapters = "resources/{sample}/adapters.fasta"
    output:
        r1 = temp("resources/{sample}/rawreads/{sample}.R1.trimmed.gz"),
        r2 = temp("resources/{sample}/rawreads/{sample}.R2.trimmed.gz")
    log:
        "logs/{sample}/trimming/{sample}_run_fastp.log"
    threads: 4
    benchmark:
        "benchmarks/{sample}/trimming/{sample}_run_fastp.benchmark"
    params:
        minquality=lambda wildcards: config["samples"][wildcards.sample]["minquality"],
        minlength=lambda wildcards: config["samples"][wildcards.sample]["minlength"],
        adapters=lambda wildcards: f"--adapter_fasta resources/{wildcards.sample}/adapters.fasta" if config["samples"][wildcards.sample].get("adapters", "") else " --detect_adapter_for_pe "
    singularity:
        f"{config["sif_dir"]}/fastp.sif"
    shell:
        """
        fastp --in1 {input.r1} --in2 {input.r2} {params.adapters} \
        --length_required {params.minlength} --cut_mean_quality {params.minquality} --thread {threads} \
        --html resources/{wildcards.sample}/rawreads/fastp.html \
        --json resources/{wildcards.sample}/rawreads/fastp.json \
        --out1 {output.r1} --out2 {output.r2} >> {log} 2>&1
        """

rule run_cutadapt:
    input:
        "resources/{sample}/rawreads/{sample}.fastq.gz"
    output:
        temp("resources/{sample}/rawreads/{sample}.trimmed.fastq.gz")
    log:
        "logs/{sample}.run_cutadapt.log"
    benchmark:
        "benchmarks/{sample}.run_cutadapt.txt"
    threads: 8
    params:
        adapters=lambda wildcards: config["samples"][wildcards.sample]["pacbio_adapters"]
    singularity:
        f"{config["sif_dir"]}/cutadapt.sif"
    shell:
        """
        cutadapt -j {threads} --discard -O 35 --rc \
            -e 0.1 {params} -o {output} {input} >> {log} 2>&1
        """

rule run_seqtk:
    input:
        lambda wildcards:
            "resources/{sample}/rawreads/{sample}.trimmed.fastq.gz" if config["samples"][wildcards.sample].get("run_trimming", "").lower() == "yes"
            else "resources/{sample}/rawreads/{sample}.fastq.gz",
    output:
        temp("resources/{sample}/rawreads/{sample}.trimmed.fasta")
    log:
        "logs/{sample}.run_seqtk.log"
    benchmark:
        "benchmarks/{sample}.run_seqtk.txt"
    threads: 1
    params:
        adapters=lambda wildcards: config["samples"][wildcards.sample]["pacbio_adapters"]
    singularity:
        f"{config["sif_dir"]}/seqtk.sif"
    shell:
        """
        seqtk seq -a {input} > {output} 2>&1
        """
