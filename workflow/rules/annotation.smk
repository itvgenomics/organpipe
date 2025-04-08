rule run_mitos2:
    input:
        "results/{sample}/pilon/{seed}_kmer{kmer}.pilon.check"
    output:
        "results/{sample}/mitos2/{seed}_kmer{kmer}/{seed}_kmer{kmer}.mitos2.check"
    log:
        "logs/{sample}/mitos2/{kmer}_{seed}_run_mitos2.log"
    threads: 1
    benchmark:
        "benchmarks/{sample}/mitos2/{kmer}_{seed}_run_mitos2.benchmark"
    singularity:
        f"{config["sif_dir"]}/mitos.sif"
    params:
        genetic_code=lambda wildcards: config["samples"][wildcards.sample]["genetic_code"],
        refseq_dir="resources/refseq89m"
    shell:
        """
        for fasta_file in results/{wildcards.sample}/assemblies/{wildcards.seed}_kmer{wildcards.kmer}/*.fasta; do
            fasta_header=$(awk '/^>/ {{print; exit}}' "$fasta_file" | sed 's/^>//') && \
            if [ "$fasta_header" != "INVALIDSEED_1" ]; then
                awk '/^>/{{sub(/_pilon$/,"",$0)}}1' results/{wildcards.sample}/pilon/$fasta_header/$fasta_header.fasta > results/{wildcards.sample}/pilon/$fasta_header/temp.fasta && \
                mv results/{wildcards.sample}/pilon/$fasta_header/temp.fasta results/{wildcards.sample}/pilon/$fasta_header/$fasta_header.fasta && \
                mkdir -p results/{wildcards.sample}/mitos2/{wildcards.seed}_kmer{wildcards.kmer}/$fasta_header && \
                runmitos.py --code {params.genetic_code} \
                --input results/{wildcards.sample}/pilon/$fasta_header/$fasta_header.fasta \
                --outdir results/{wildcards.sample}/mitos2/{wildcards.seed}_kmer{wildcards.kmer}/$fasta_header \
                --refdir {params.refseq_dir} --noplots --best \
                > results/{wildcards.sample}/mitos2/{wildcards.seed}_kmer{wildcards.kmer}/$fasta_header/mitos.log 2>{log}
            fi
        done
        touch {output}
        """


rule run_cpgavas2:
    input:
        fasta="results/{sample}/pilon/{seed}_kmer{kmer}.pilon.check",
    output:
        "results/{sample}/cpgavas2/{seed}_kmer{kmer}/{seed}_kmer{kmer}.cpgavas2.check"
    log:
        "logs/{sample}/cpgavas2/{kmer}_{seed}_run_cpgavas2.log"
    threads: 1
    benchmark:
        "benchmarks/{sample}/cpgavas2/{kmer}_{seed}_run_cpgavas2.benchmark"
    singularity:
        f"{config["sif_dir"]}/pimba_adapterremoval.sif"
    shell:
        """
        sed -i '/maker/s/-quiet/--ignore_nfs_tmp -quiet/' /apps/cpgavas2C/modules/plasAnno/bin/Annotation_Chloroplast_King.py && \
        for fasta_file in results/{wildcards.sample}/assemblies/{wildcards.seed}_kmer{wildcards.kmer}/*.fasta; do
            fasta_header=$(awk '/^>/ {{print; exit}}' "$fasta_file" | sed 's/^>//') && \
            if [ "$fasta_header" != "INVALIDSEED_1" ]; then
                random_pid=$(( (RANDOM + RANDOM * 32768 + RANDOM * 32768 * 32768) % 999999999 + 1 )) && \
                awk '/^>/{{sub(/_pilon$/,"",$0)}}1' results/{wildcards.sample}/pilon/$fasta_header/$fasta_header.fasta > results/{wildcards.sample}/pilon/$fasta_header/temp.fasta && \
                mv results/{wildcards.sample}/pilon/$fasta_header/temp.fasta results/{wildcards.sample}/pilon/$fasta_header/$fasta_header.fasta && \
                mkdir -p results/{wildcards.sample}/cpgavas2/{wildcards.seed}_kmer{wildcards.kmer}/$fasta_header && \
                run-cpgavas2 -pid $random_pid -in results/{wildcards.sample}/pilon/$fasta_header/$fasta_header.fasta -db 2 >> {log} 2>1 && \
                mv /tmp/dir_$random_pid/* results/{wildcards.sample}/cpgavas2/{wildcards.seed}_kmer{wildcards.kmer}/$fasta_header/
            fi
        done
        touch {output}
        """

rule run_chloe:
    input:
        "results/{sample}/pilon/{seed}_kmer{kmer}.pilon.check"
    output:
        "results/{sample}/chloe/{seed}_kmer{kmer}.chloe.check"
    log:
        "logs/{sample}/chloe/{kmer}_{seed}_run_chloe.log"
    threads: 1
    benchmark:
        "benchmarks/{sample}/chloe/{kmer}_{seed}_run_chloe.benchmark"
    singularity:
        f"{config["sif_dir"]}/chloe.sif"
    shell:
        """
        mkdir -p results/{wildcards.sample}/chloe && \
        for fasta_file in results/{wildcards.sample}/assemblies/{wildcards.seed}_kmer{wildcards.kmer}/*.fasta; do
            fasta_header=$(awk '/^>/ {{print; exit}}' "$fasta_file" | sed 's/^>//') && \
            if [ "$fasta_header" != "INVALIDSEED_1" ]; then
                awk '/^>/{{sub(/_pilon$/,"",$0)}}1' results/{wildcards.sample}/pilon/$fasta_header/$fasta_header.fasta > results/{wildcards.sample}/pilon/$fasta_header/temp.fasta && \
                mv results/{wildcards.sample}/pilon/$fasta_header/temp.fasta results/{wildcards.sample}/pilon/$fasta_header/$fasta_header.fasta && \
                julia --project=/opt/chloe /opt/chloe/chloe.jl annotate -o results/{wildcards.sample}/chloe results/{wildcards.sample}/pilon/$fasta_header/$fasta_header.fasta >> {log} 2>&1
            fi
        done
        touch {output}
        """
