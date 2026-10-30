def get_shortreads(sample):
    return find_short_reads(sample, config["samples"][sample]["reads_path"])

def get_longreads(sample):
    return find_long_reads(sample, config["samples"][sample]["reads_path"])

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
        long_reads=lambda wildcards: get_longreads(wildcards.sample),
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
        long_reads=lambda wildcards: get_longreads(wildcards.sample),
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
        lambda wildcards: config["samples"][wildcards.sample]["adapters"] or "test_data/adapters.fasta"
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
    run:
        db_path = config.get("nhmmer_db")

        if db_path:
            shell(f"cp {db_path} {output}")
        else:
            raise ValueError("The 'nhmmer_db' path is not specified in the configuration file.")

rule extract_table2asn:
    input:
        "resources/table2asn.linux64.gz"
    output:
        "resources/table2asn.linux64"
    shell:
        """
        gunzip -c {input} > {output} && chmod +x {output}
        """

rule index_hmmer_db:
    input:
        "resources/nhmmer_db.hmm"
    output:
        expand("resources/nhmmer_db.hmm.{ext}", ext=['h3f', 'h3i', 'h3m', 'h3p'])
    singularity:
        sif("hmmer")
    shell:
        """
        hmmpress {input}
        """

rule check_reference:
    input:
        lambda wildcards: config["samples"][wildcards.sample]["reference"] or "test_data/reference.fasta"
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

rule download_getorganelle_db:
    output:
        "resources/getorganelle_db/getorganelle_db.check"
    log:
        "logs/getorganelle/download_getorganelle_db.log"
    singularity:
        sif("getorganelle")
    shell:
        """
        export PATH=/opt/conda/bin:$PATH && \
        get_organelle_config.py -a all --config-dir resources/getorganelle_db/ >> {log} 2>&1 && \
        touch {output}
        """
