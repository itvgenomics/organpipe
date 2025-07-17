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
            original_header=$(awk '/^>/ {{print; exit}}' "$fasta_file" | sed 's/^>//') && \
            if [ "$original_header" != "INVALIDSEED_1" ]; then
                pilon_dir=results/{wildcards.sample}/pilon/$original_header

                # Remove the _pilon sufix from the header
                awk '/^>/{{sub(/_pilon$/,"",$0)}}1' $pilon_dir/$original_header.fasta > $pilon_dir/temp.fasta && \
                mv $pilon_dir/temp.fasta $pilon_dir/$original_header.fasta && \

                mitos_outdir=results/{wildcards.sample}/mitos2/{wildcards.seed}_kmer{wildcards.kmer}/$original_header

                # Rewrite header to temporary short name
                awk -v new_header=">temp_header" '/^>/{{$0=new_header}}1' "$pilon_dir/$original_header.fasta" > "$pilon_dir/temp.fasta"

                # Run MITOS2
                mkdir -p "$mitos_outdir" && \
                runmitos.py --code {params.genetic_code} \
                    --input "$pilon_dir/temp.fasta" \
                    --outdir "$mitos_outdir" \
                    --refdir {params.refseq_dir} --noplots --best \
                    > "$mitos_outdir/mitos.log" 2>{log}

                # Remove the temp.fasta file
                rm "$pilon_dir/temp.fasta"
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
        f"{config["sif_dir"]}/cpgavas2.sif"
    shell:
        """
        sed -i '/maker/s/-quiet/--ignore_nfs_tmp -quiet/' /apps/cpgavas2C/modules/plasAnno/bin/Annotation_Chloroplast_King.py && \
        for fasta_file in results/{wildcards.sample}/assemblies/{wildcards.seed}_kmer{wildcards.kmer}/*.fasta; do
            fasta_header=$(awk '/^>/ {{print; exit}}' "$fasta_file" | sed 's/^>//') && \
            if [ "$fasta_header" != "INVALIDSEED_1" ]; then
                random_pid=$(( (RANDOM + RANDOM * 32768 + RANDOM * 32768 * 32768) % 999999999 + 1 )) && \
                awk '/^>/{{sub(/_pilon$/,"",$0)}}1' results/{wildcards.sample}/pilon/$fasta_header/$fasta_header.fasta > results/{wildcards.sample}/pilon/$fasta_header/temp.fasta && \
                mv results/{wildcards.sample}/pilon/$fasta_header/temp.fasta results/{wildcards.sample}/pilon/$fasta_header/$fasta_header.fasta && \
                rm -rf results/{wildcards.sample}/cpgavas2/{wildcards.seed}_kmer{wildcards.kmer}/$fasta_header && \
                mkdir -p results/{wildcards.sample}/cpgavas2/{wildcards.seed}_kmer{wildcards.kmer}/$fasta_header && \
                run-cpgavas2 -pid $random_pid -in results/{wildcards.sample}/pilon/$fasta_header/$fasta_header.fasta -db 2 >> {log} 2>&1 && \
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
