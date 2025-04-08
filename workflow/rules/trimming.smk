rule run_adapterremoval2:
    input:
        r1 = "resources/{sample}/rawreads/{sample}_R1.fastq.gz",
        r2 = "resources/{sample}/rawreads/{sample}_R2.fastq.gz",
        adapters="resources/{sample}/adapters.txt",
    output:
        r1 = temp("resources/{sample}/rawreads/{sample}.pair1.truncated.gz"),
        r2 = temp("resources/{sample}/rawreads/{sample}.pair2.truncated.gz")
    log:
        "logs/{sample}/trimming/{sample}_run_adapterremoval2.log"
    threads: 4
    benchmark:
        "benchmarks/{sample}/trimming/{sample}_run_adapterremoval2.benchmark"
    params:
        minquality=lambda wildcards: config["samples"][wildcards.sample]["minquality"],
        minlength=lambda wildcards: config["samples"][wildcards.sample]["minlength"]
    singularity:
        f"{config["sif_dir"]}/pimba_adapterremoval.sif"
    shell:
        """
        AdapterRemoval --file1 {input.r1} --file2 {input.r2} \
        --threads {threads} --adapter-list {input.adapters} \
        --trimwindows 10 --minquality {params.minquality} --minlength {params.minlength} \
        --qualitymax 42 --mm 5 --basename resources/{wildcards.sample}/rawreads/{wildcards.sample} --gzip >> {log} 2>&1
        """
