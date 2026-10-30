rule run_bwa_index:
    input:
        "results/{sample}/assemblies/{assembler}/{run}/{assembly}.fasta"
    output:
        fasta = "results/{sample}/pilon/{assembler}/{run}/{assembly}/draft.fasta",
        index = multiext("results/{sample}/pilon/{assembler}/{run}/{assembly}/draft.fasta", ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac")
    log:
        "logs/{sample}/pilon/{assembler}/{run}/{assembly}_run_bwa_index.log"
    benchmark:
        "benchmarks/{sample}/pilon/{assembler}/{run}/{assembly}_run_bwa_index.benchmark"
    singularity:
        sif("hic_mapping")
    shell:
        """
        cp {input} {output.fasta} && \
        bwa-mem2.avx index {output.fasta} >> {log} 2>&1
        """

rule run_bwa_mem:
    input:
        unpack(assembled_reads),
        fasta = "results/{sample}/pilon/{assembler}/{run}/{assembly}/draft.fasta",
        index = multiext("results/{sample}/pilon/{assembler}/{run}/{assembly}/draft.fasta", ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac")
    output:
        "results/{sample}/pilon/{assembler}/{run}/{assembly}/{assembly}_mapping.bam"
    log:
        "logs/{sample}/pilon/{assembler}/{run}/{assembly}_run_bwa_mem.log"
    benchmark:
        "benchmarks/{sample}/pilon/{assembler}/{run}/{assembly}_run_bwa_mem.benchmark"
    singularity:
        sif("hic_mapping")
    shell:
        """
        {{ bwa-mem2.avx mem -t {threads} {input.fasta} {input.r1} {input.r2} 2>> {log}
        }} | samtools view - -Sb | samtools sort - -@ {threads} -o {output} >> {log} 2>&1
        """

rule run_samtools_index:
    input:
        "results/{sample}/pilon/{assembler}/{run}/{assembly}/{assembly}_mapping.bam"
    output:
        "results/{sample}/pilon/{assembler}/{run}/{assembly}/{assembly}_mapping.bam.bai"
    log:
        "logs/{sample}/pilon/{assembler}/{run}/{assembly}_run_samtools_index.log"
    benchmark:
        "benchmarks/{sample}/pilon/{assembler}/{run}/{assembly}_run_samtools_index.benchmark"
    singularity:
        sif("hic_mapping")
    shell:
        """
        samtools index {input} >> {log} 2>&1
        """

rule run_pilon:
    input:
        fasta = "results/{sample}/pilon/{assembler}/{run}/{assembly}/draft.fasta",
        bam = "results/{sample}/pilon/{assembler}/{run}/{assembly}/{assembly}_mapping.bam",
        bai = "results/{sample}/pilon/{assembler}/{run}/{assembly}/{assembly}_mapping.bam.bai"
    output:
        fasta = "results/{sample}/pilon/{assembler}/{run}/{assembly}/{assembly}.fasta",
        changes = "results/{sample}/pilon/{assembler}/{run}/{assembly}/{assembly}.changes",
        pilon_log = "results/{sample}/pilon/{assembler}/{run}/{assembly}/pilon.log"
    log:
        "logs/{sample}/pilon/{assembler}/{run}/{assembly}_run_pilon.log"
    benchmark:
        "benchmarks/{sample}/pilon/{assembler}/{run}/{assembly}_run_pilon.benchmark"
    singularity:
        sif("pilon")
    shell:
        """
        chmod +x resources/pilon.sh && \
        export PARALLEL_GC_THREADS={threads} && \
        export JAVA_TOOL_OPTIONS='-XX:ParallelGCThreads={threads}' && \
        ./resources/pilon.sh --genome {input.fasta} \
            --fix all --changes \
            --bam {input.bam} \
            --output results/{wildcards.sample}/pilon/{wildcards.assembler}/{wildcards.run}/{wildcards.assembly}/{wildcards.assembly} --threads {threads} \
            > {output.pilon_log} 2>> {log}

        # Keep the assembly name in the header: Pilon adds a _pilon suffix to it
        awk '/^>/{{sub(/_pilon$/,"",$0)}}1' {output.fasta} > {output.fasta}.tmp && \
        mv {output.fasta}.tmp {output.fasta}
        """
