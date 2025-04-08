rule run_bwa_index:
    input:
        "results/{sample}/assemblies/{seed}_kmer{kmer}.fasta"
    output:
        "results/{sample}/pilon/{seed}_kmer{kmer}.bwa_index.check"
    log:
        "logs/{sample}/pilon/{kmer}_{seed}_run_bwa_index.log"
    threads: 1
    benchmark:
        "benchmarks/{sample}/pilon/{kmer}_{seed}_run_bwa_index.benchmark"
    singularity:
        f"{config["sif_dir"]}/hic_mapping.sif"
    shell:
        """
        python workflow/scripts/split_fasta.py --fasta_file {input} \
        --output_dir results/{wildcards.sample}/assemblies/{wildcards.seed}_kmer{wildcards.kmer} >> {log} 2>&1 && \
        for fasta_file in results/{wildcards.sample}/assemblies/{wildcards.seed}_kmer{wildcards.kmer}/*.fasta; do
            fasta_header=$(awk '/^>/ {{print; exit}}' "$fasta_file" | sed 's/^>//') && \
            if [ "$fasta_header" != "INVALIDSEED_1" ]; then
                bwa-mem2.avx index $fasta_file >> {log} 2>&1
            fi
        done
        touch {output}
        """

rule run_bwa_mem:
    input:
        check="results/{sample}/pilon/{seed}_kmer{kmer}.bwa_index.check",
    output:
        "results/{sample}/pilon/{seed}_kmer{kmer}.bwa_mem.check"
    log:
        "logs/{sample}/pilon/{kmer}_{seed}_run_bwa_mem.log"
    threads: 4
    benchmark:
        "benchmarks/{sample}/pilon/{kmer}_{seed}_run_bwa_mem.benchmark"
    singularity:
        f"{config["sif_dir"]}/hic_mapping.sif"
    shell:
        """
        for fasta_file in results/{wildcards.sample}/assemblies/{wildcards.seed}_kmer{wildcards.kmer}/*.fasta; do
            fasta_header=$(awk '/^>/ {{print; exit}}' "$fasta_file" | sed 's/^>//') && \
            mkdir -p results/{wildcards.sample}/pilon/$fasta_header && \
            if [ "$fasta_header" != "INVALIDSEED_1" ]; then
                bwa-mem2.avx mem -t {threads} "$fasta_file" \
                results/{wildcards.sample}/novoplasty/{wildcards.seed}/kmer{wildcards.kmer}/Assembled_reads_{wildcards.sample}_R1.fasta \
                results/{wildcards.sample}/novoplasty/{wildcards.seed}/kmer{wildcards.kmer}/Assembled_reads_{wildcards.sample}_R2.fasta | samtools view - -Sb | samtools sort - -@ {threads} \
                -o results/{wildcards.sample}/pilon/$fasta_header/"$fasta_header"_mapping.bam >> {log} 2>&1
            fi
        done
        touch {output}
        """

rule run_samtools_index:
    input:
        "results/{sample}/pilon/{seed}_kmer{kmer}.bwa_mem.check"
    output:
        "results/{sample}/pilon/{seed}_kmer{kmer}.samtools_index.check"
    log:
        "logs/{sample}/pilon/{kmer}_{seed}_run_samtools_index.log"
    threads: 1
    benchmark:
        "benchmarks/{sample}/pilon/{kmer}_{seed}_run_samtools_index.benchmark"
    singularity:
        f"{config["sif_dir"]}/hic_mapping.sif"
    shell:
        """
        for fasta_file in results/{wildcards.sample}/assemblies/{wildcards.seed}_kmer{wildcards.kmer}/*.fasta; do
            fasta_header=$(awk '/^>/ {{print; exit}}' "$fasta_file" | sed 's/^>//') && \
            if [ "$fasta_header" != "INVALIDSEED_1" ]; then
                samtools index results/{wildcards.sample}/pilon/$fasta_header/"$fasta_header"_mapping.bam >> {log} 2>&1
            fi
        done
        touch {output}
        """

rule run_pilon:
    input:
        "results/{sample}/pilon/{seed}_kmer{kmer}.samtools_index.check"
    output:
        "results/{sample}/pilon/{seed}_kmer{kmer}.pilon.check"
    log:
        "logs/{sample}/pilon/{kmer}_{seed}_run_pilon.log"
    threads: 2
    benchmark:
        "benchmarks/{sample}/pilon/{kmer}_{seed}_run_pilon.benchmark"
    singularity:
        f"{config["sif_dir"]}/pilon.sif"
    shell:
        """
        export PARALLEL_GC_THREADS={threads} && \
		export JAVA_TOOL_OPTIONS='-XX:ParallelGCThreads={threads}' && \
        for fasta_file in results/{wildcards.sample}/assemblies/{wildcards.seed}_kmer{wildcards.kmer}/*.fasta; do
            fasta_header=$(awk '/^>/ {{print; exit}}' "$fasta_file" | sed 's/^>//') && \
            if [ "$fasta_header" != "INVALIDSEED_1" ]; then
                pilon --genome "$fasta_file" \
                --fix all --changes \
                --bam results/{wildcards.sample}/pilon/$fasta_header/"$fasta_header"_mapping.bam \
                --output results/{wildcards.sample}/pilon/$fasta_header/$fasta_header --threads {threads} \
                > results/{wildcards.sample}/pilon/$fasta_header/pilon.log 2>{log}
            fi
        done
        touch {output}
        """
