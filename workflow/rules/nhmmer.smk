def annotation_done(wildcards):
    """Output of the step that wrote the annotation used in the nhmmscan search."""
    if wildcards.assembler == "mitohifi":
        return f"results/{wildcards.sample}/mitohifi/{wildcards.run}/contigs_stats.tsv"
    if config["samples"][wildcards.sample]["organelle"] == "mito":
        return f"results/{wildcards.sample}/mitos2/{wildcards.assembler}/{wildcards.run}/{wildcards.assembly}/mitos2.fasta"
    return f"results/{wildcards.sample}/genbanks/{wildcards.assembler}/{wildcards.run}/{wildcards.assembly}.cpgavas2.genbank.check"

def annotation_source(wildcards):
    """Argument of the sequence extraction scripts with the annotation to read."""
    if wildcards.assembler == "mitohifi":
        return f"--genbank results/{wildcards.sample}/mitohifi/{wildcards.run}/{MITOHIFI_ASSEMBLY}.gb"
    if config["samples"][wildcards.sample]["organelle"] == "mito":
        return f"--mitos_dir results/{wildcards.sample}/mitos2/{wildcards.assembler}/{wildcards.run}/{wildcards.assembly}"
    return f"--genbank {genbank_path(wildcards, 'cpgavas2')}"

rule get_ncRNA_sequences:
    input:
        annotation_done
    output:
        "results/{sample}/nhmmer/{assembler}/{run}/{assembly}/rRNA-tRNA.fasta"
    log:
        "logs/{sample}/nhmmer/{assembler}/{run}/{assembly}_get_ncRNA_sequences.log"
    benchmark:
        "benchmarks/{sample}/nhmmer/{assembler}/{run}/{assembly}_get_ncRNA_sequences.benchmark"
    params:
        source=annotation_source,
        dir=output_dir
    shell:
        """
        python workflow/scripts/ncRNA_nhmmer.py {params.source} --output_dir {params.dir} >> {log} 2>&1
        """

rule get_intergenes_sequences:
    input:
        annotation_done
    output:
        "results/{sample}/nhmmer/{assembler}/{run}/{assembly}/intergenes_filter.fasta"
    log:
        "logs/{sample}/nhmmer/{assembler}/{run}/{assembly}_get_intergenes_sequences.log"
    benchmark:
        "benchmarks/{sample}/nhmmer/{assembler}/{run}/{assembly}_get_intergenes_sequences.benchmark"
    params:
        source=annotation_source,
        dir=output_dir
    shell:
        """
        python workflow/scripts/intergenes_nhmmer.py {params.source} --output_dir {params.dir} >> {log} 2>&1
        """

rule run_nhmmer:
    input:
        fasta = "results/{sample}/nhmmer/{assembler}/{run}/{assembly}/{region}.fasta",
        db = "resources/nhmmer_db.hmm",
        db_index = expand("resources/nhmmer_db.hmm.{ext}", ext=['h3f', 'h3i', 'h3m', 'h3p'])
    output:
        "results/{sample}/nhmmer/{assembler}/{run}/{assembly}/{region}.nhmmer.check"
    log:
        "logs/{sample}/nhmmer/{assembler}/{run}/{assembly}_{region}_run_nhmmer.log"
    benchmark:
        "benchmarks/{sample}/nhmmer/{assembler}/{run}/{assembly}_{region}_run_nhmmer.benchmark"
    singularity:
        sif("hmmer")
    params:
        prefix=lambda wildcards, input: input.fasta[: -len(".fasta")]
    shell:
        """
        if [ -s {input.fasta} ]; then
            nhmmscan --cpu {threads} --noali \
            -o {params.prefix}.nhmmer.out \
            --tblout {params.prefix}.tblout.out \
            --dfamtblout {params.prefix}.dfamtblout.out \
            {input.db} {input.fasta} 2>> {log}
        else
            echo "Input fasta file is empty. Might be an annotation error." >> {log} 2>&1
        fi
        touch {output}
        """
