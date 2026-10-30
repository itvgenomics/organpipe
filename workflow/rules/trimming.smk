rule run_fastp:
    input:
        r1 = "resources/{sample}/rawreads/{sample}_R1.fastq.gz",
        r2 = "resources/{sample}/rawreads/{sample}_R2.fastq.gz",
        adapters = "resources/{sample}/adapters.fasta"
    output:
        r1 = temp("resources/{sample}/rawreads/{sample}.R1.trimmed.gz"),
        r2 = temp("resources/{sample}/rawreads/{sample}.R2.trimmed.gz"),
        html = "resources/{sample}/rawreads/fastp.html",
        json = "resources/{sample}/rawreads/fastp.json"
    log:
        "logs/{sample}/trimming/{sample}_run_fastp.log"
    benchmark:
        "benchmarks/{sample}/trimming/{sample}_run_fastp.benchmark"
    params:
        minquality=sample_param("minquality"),
        minlength=sample_param("minlength"),
        adapters=lambda wildcards: f"--adapter_fasta resources/{wildcards.sample}/adapters.fasta" if config["samples"][wildcards.sample].get("adapters", "") else " --detect_adapter_for_pe "
    singularity:
        sif("fastp")
    shell:
        """
        fastp --in1 {input.r1} --in2 {input.r2} {params.adapters} \
        --length_required {params.minlength} --cut_front --cut_tail --cut_mean_quality {params.minquality} --thread {threads} \
        --html {output.html} \
        --json {output.json} \
        --out1 {output.r1} --out2 {output.r2} >> {log} 2>&1
        """

rule run_cutadapt:
    input:
        "resources/{sample}/rawreads/{sample}.fastq.gz"
    output:
        reads=temp("resources/{sample}/rawreads/{sample}.trimmed.fastq.gz"),
        check="resources/{sample}/rawreads/{sample}.trimmed.check"
    log:
        "logs/{sample}/trimming/{sample}_run_cutadapt.log"
    benchmark:
        "benchmarks/{sample}/trimming/{sample}_run_cutadapt.benchmark"
    params:
        adapters=sample_param("pacbio_adapters")
    singularity:
        sif("cutadapt")
    shell:
        """
        cutadapt -j {threads} --discard -O 35 --rc \
            -e 0.1 {params} -o {output.reads} {input} > {log} 2>&1

        touch {output.check}
        """

rule run_seqtk:
    input:
        lambda wildcards:
            "resources/{sample}/rawreads/{sample}.trimmed.fastq.gz" if is_yes(wildcards.sample, "run_trimming")
            else "resources/{sample}/rawreads/{sample}.fastq.gz",
    output:
        temp("resources/{sample}/rawreads/{sample}.converted.fasta")
    log:
        "logs/{sample}/trimming/{sample}_run_seqtk.log"
    benchmark:
        "benchmarks/{sample}/trimming/{sample}_run_seqtk.benchmark"
    singularity:
        sif("seqtk")
    shell:
        """
        seqtk seq -a {input} > {output} 2> {log}
        """
