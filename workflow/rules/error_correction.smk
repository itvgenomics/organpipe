rule run_novoplasty_bwa_index:
    input:
        "results/{sample}/assemblies/novoplasty/{seed}_kmer{kmer}.fasta"
    output:
        "results/{sample}/pilon/novoplasty/{seed}_kmer{kmer}.bwa_index.check"
    log:
        "logs/{sample}/pilon/novoplasty/{kmer}_{seed}_run_bwa_index.log"
    benchmark:
        "benchmarks/{sample}/pilon/novoplasty/{kmer}_{seed}_run_bwa_index.benchmark"
    singularity:
        f"{config["sif_dir"]}/hic_mapping.sif"
    shell:
        """
        python workflow/scripts/split_fasta.py --fasta_file {input} --assembler 'novoplasty' \
        --output_dir results/{wildcards.sample}/assemblies/novoplasty/{wildcards.seed}_kmer{wildcards.kmer} >> {log} 2>&1 && \
        for fasta_file in results/{wildcards.sample}/assemblies/novoplasty/{wildcards.seed}_kmer{wildcards.kmer}/*.fasta; do
            fasta_header=$(awk '/^>/ {{print; exit}}' "$fasta_file" | sed 's/^>//') && \
            if [ "$fasta_header" != "INVALIDSEED_1" ]; then
                bwa-mem2.avx index $fasta_file >> {log} 2>&1
            fi
        done
        touch {output}
        """

rule run_novoplasty_bwa_mem:
    input:
        check="results/{sample}/pilon/novoplasty/{seed}_kmer{kmer}.bwa_index.check",
    output:
        "results/{sample}/pilon/novoplasty/{seed}_kmer{kmer}.bwa_mem.check"
    log:
        "logs/{sample}/pilon/novoplasty/{kmer}_{seed}_run_bwa_mem.log"
    benchmark:
        "benchmarks/{sample}/pilon/novoplasty/{kmer}_{seed}_run_bwa_mem.benchmark"
    singularity:
        f"{config["sif_dir"]}/hic_mapping.sif"
    shell:
        """
        for fasta_file in results/{wildcards.sample}/assemblies/novoplasty/{wildcards.seed}_kmer{wildcards.kmer}/*.fasta; do
            fasta_header=$(awk '/^>/ {{print; exit}}' "$fasta_file" | sed 's/^>//') && \
            mkdir -p results/{wildcards.sample}/pilon/novoplasty/$fasta_header && \
            if [ "$fasta_header" != "INVALIDSEED_1" ]; then
                {{ bwa-mem2.avx mem -t {threads} "$fasta_file" \
                    results/{wildcards.sample}/novoplasty/{wildcards.seed}/kmer{wildcards.kmer}/Assembled_reads_{wildcards.sample}_R1.fasta \
                    results/{wildcards.sample}/novoplasty/{wildcards.seed}/kmer{wildcards.kmer}/Assembled_reads_{wildcards.sample}_R2.fasta 2>> {log}
                }} | samtools view - -Sb | samtools sort - -@ {threads} \
                -o results/{wildcards.sample}/pilon/novoplasty/$fasta_header/"$fasta_header"_mapping.bam >> {log} 2>&1
            fi
        done
        touch {output}
        """

rule run_novoplasty_samtools_index:
    input:
        "results/{sample}/pilon/novoplasty/{seed}_kmer{kmer}.bwa_mem.check"
    output:
        "results/{sample}/pilon/novoplasty/{seed}_kmer{kmer}.samtools_index.check"
    log:
        "logs/{sample}/pilon/novoplasty/{kmer}_{seed}_run_samtools_index.log"
    benchmark:
        "benchmarks/{sample}/pilon/novoplasty/{kmer}_{seed}_run_samtools_index.benchmark"
    singularity:
        f"{config["sif_dir"]}/hic_mapping.sif"
    shell:
        """
        for fasta_file in results/{wildcards.sample}/assemblies/novoplasty/{wildcards.seed}_kmer{wildcards.kmer}/*.fasta; do
            fasta_header=$(awk '/^>/ {{print; exit}}' "$fasta_file" | sed 's/^>//') && \
            if [ "$fasta_header" != "INVALIDSEED_1" ]; then
                samtools index results/{wildcards.sample}/pilon/novoplasty/$fasta_header/"$fasta_header"_mapping.bam >> {log} 2>&1
            fi
        done
        touch {output}
        """

rule run_novoplasty_pilon:
    input:
        "results/{sample}/pilon/novoplasty/{seed}_kmer{kmer}.samtools_index.check"
    output:
        "results/{sample}/pilon/novoplasty/{seed}_kmer{kmer}.pilon.check"
    log:
        "logs/{sample}/pilon/novoplasty/{kmer}_{seed}_run_pilon.log"
    benchmark:
        "benchmarks/{sample}/pilon/novoplasty/{kmer}_{seed}_run_pilon.benchmark"
    singularity:
        f"{config["sif_dir"]}/pilon.sif"
    shell:
        """
        chmod +x resources/pilon.sh && \
        export PARALLEL_GC_THREADS={threads} && \
		export JAVA_TOOL_OPTIONS='-XX:ParallelGCThreads={threads}' && \
        for fasta_file in results/{wildcards.sample}/assemblies/novoplasty/{wildcards.seed}_kmer{wildcards.kmer}/*.fasta; do
            fasta_header=$(awk '/^>/ {{print; exit}}' "$fasta_file" | sed 's/^>//') && \
            if [ "$fasta_header" != "INVALIDSEED_1" ]; then
                ./resources/pilon.sh --genome "$fasta_file" \
                --fix all --changes \
                --bam results/{wildcards.sample}/pilon/novoplasty/$fasta_header/"$fasta_header"_mapping.bam \
                --output results/{wildcards.sample}/pilon/novoplasty/$fasta_header/$fasta_header --threads {threads} \
                > results/{wildcards.sample}/pilon/novoplasty/$fasta_header/pilon.log 2>{log}
            fi
        done
        touch {output}
        """


rule run_getorganelle_bwa_index:
    input:
        "results/{sample}/getorganelle/sequences.fasta"
    output:
        "results/{sample}/pilon/getorganelle/getorganelle_bwa_index.check"
    log:
        "logs/{sample}/pilon/getorganelle/run_bwa_index.log"
    benchmark:
        "benchmarks/{sample}/pilon/getorganelle/run_bwa_index.benchmark"
    singularity:
        f"{config["sif_dir"]}/hic_mapping.sif"
    shell:
        """
        python workflow/scripts/split_fasta.py --fasta_file {input} --assembler 'getorganelle' \
        --output_dir results/{wildcards.sample}/assemblies/getorganelle >> {log} 2>&1 && \
        for fasta_file in results/{wildcards.sample}/assemblies/getorganelle/*.fasta; do
            bwa-mem2.avx index $fasta_file >> {log} 2>&1
        done
        touch {output}
        """

rule run_getorganelle_bwa_mem:
    input:
        check="results/{sample}/pilon/getorganelle/getorganelle_bwa_index.check",
    output:
        "results/{sample}/pilon/getorganelle/getorganelle_bwa_mem.check"
    log:
        "logs/{sample}/pilon/getorganelle/run_bwa_mem.log"
    benchmark:
        "benchmarks/{sample}/pilon/getorganelle/run_bwa_mem.benchmark"
    singularity:
        f"{config["sif_dir"]}/hic_mapping.sif"
    shell:
        """
        for fasta_file in results/{wildcards.sample}/assemblies/getorganelle/*.fasta; do
            fasta_header=$(awk '/^>/ {{print; exit}}' "$fasta_file" | sed 's/^>//') && \
            mkdir -p results/{wildcards.sample}/pilon/getorganelle/$fasta_header && \
            {{ bwa-mem2.avx mem -t {threads} "$fasta_file" \
                results/{wildcards.sample}/getorganelle/extended_1_paired.fq \
                results/{wildcards.sample}/getorganelle/extended_2_paired.fq 2>> {log}
            }} | samtools view - -Sb | samtools sort - -@ {threads} \
            -o results/{wildcards.sample}/pilon/getorganelle/$fasta_header/"$fasta_header"_mapping.bam >> {log} 2>&1
        done
        touch {output}
        """

rule run_getorganelle_samtools_index:
    input:
        "results/{sample}/pilon/getorganelle/getorganelle_bwa_mem.check"
    output:
        "results/{sample}/pilon/getorganelle/getorganelle_samtools_index.check"
    log:
        "logs/{sample}/pilon/getorganelle/run_samtools_index.log"
    benchmark:
        "benchmarks/{sample}/pilon/getorganelle/run_samtools_index.benchmark"
    singularity:
        f"{config["sif_dir"]}/hic_mapping.sif"
    shell:
        """
        for fasta_file in results/{wildcards.sample}/assemblies/getorganelle/*.fasta; do
            fasta_header=$(awk '/^>/ {{print; exit}}' "$fasta_file" | sed 's/^>//') && \
            if [ "$fasta_header" != "INVALIDSEED_1" ]; then
                samtools index results/{wildcards.sample}/pilon/getorganelle/$fasta_header/"$fasta_header"_mapping.bam >> {log} 2>&1
            fi
        done
        touch {output}
        """

rule run_getorganelle_pilon:
    input:
        "results/{sample}/pilon/getorganelle/getorganelle_samtools_index.check"
    output:
        "results/{sample}/pilon/getorganelle/getorganelle_pilon.check"
    log:
        "logs/{sample}/pilon/getorganelle/{run_pilon.log"
    benchmark:
        "benchmarks/{sample}/pilon/getorganelle/{run_pilon.benchmark"
    singularity:
        f"{config["sif_dir"]}/pilon.sif"
    shell:
        """
        chmod +x resources/pilon.sh && \
        export PARALLEL_GC_THREADS={threads} && \
		export JAVA_TOOL_OPTIONS='-XX:ParallelGCThreads={threads}' && \
        for fasta_file in results/{wildcards.sample}/assemblies/getorganelle/*.fasta; do
        fasta_header=$(awk '/^>/ {{print; exit}}' "$fasta_file" | sed 's/^>//') && \
            ./resources/pilon.sh --genome "$fasta_file" \
            --fix all --changes \
            --bam results/{wildcards.sample}/pilon/getorganelle/$fasta_header/"$fasta_header"_mapping.bam \
            --output results/{wildcards.sample}/pilon/getorganelle/$fasta_header/$fasta_header --threads {threads} \
            > results/{wildcards.sample}/pilon/getorganelle/$fasta_header/pilon.log 2>{log}
        done
        touch {output}
        """
