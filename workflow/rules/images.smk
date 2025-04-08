rule get_genbank_fastas:
    input:
        lambda wildcards: [
            f"results/{wildcards.sample}/genbanks/{wildcards.seed}_kmer{wildcards.kmer}.mitos2.genbank.rotated.check"
        ] if config["samples"][wildcards.sample]["organelle"] == "mito" else [
            f"results/{wildcards.sample}/genbanks/{wildcards.seed}_kmer{wildcards.kmer}.cpgavas2.genbank.rotated.check",
            f"results/{wildcards.sample}/genbanks/{wildcards.seed}_kmer{wildcards.kmer}.chloe.genbank.rotated.check"
        ]
    output:
        "results/{sample}/images/{seed}_kmer{kmer}/{seed}_kmer{kmer}.get_genbank_fastas.check"
    threads: 1
    log:
        "logs/{sample}/images/{seed}_{kmer}_get_genbank_fastas.log"
    benchmark:
        "benchmarks/{sample}/images/{seed}_{kmer}_get_genbank_fastas.benchmark"
    params:
        organelle=lambda wildcards: config["samples"][wildcards.sample]["organelle"]
    shell:
        """
        python workflow/scripts/images.py --sample {wildcards.sample} \
            --kmer {wildcards.kmer} --seed {wildcards.seed} --parse_gb --organelle {params.organelle} >> {log} 2>&1 && \
        touch {output}
        """

rule run_blastn:
    input:
        "results/{sample}/images/{seed}_kmer{kmer}/{seed}_kmer{kmer}.get_genbank_fastas.check"
    output:
        "results/{sample}/images/{seed}_kmer{kmer}/{seed}_kmer{kmer}.blastn.check"
    threads: 1
    log:
        "logs/{sample}/images/{seed}_{kmer}_run_blastn.log"
    benchmark:
        "benchmarks/{sample}/images/{seed}_{kmer}_run_blastn.benchmark"
    singularity:
        "docker://pegi3s/blast:2.13.0"
    shell:
        """
        cat results/{wildcards.sample}/novoplasty/{wildcards.seed}/kmer{wildcards.kmer}/Assembled_reads_{wildcards.sample}_R1.fasta \
        results/{wildcards.sample}/novoplasty/{wildcards.seed}/kmer{wildcards.kmer}/Assembled_reads_{wildcards.sample}_R2.fasta >> \
        results/{wildcards.sample}/images/{wildcards.seed}_kmer{wildcards.kmer}/Assembled_reads_{wildcards.seed}_kmer{wildcards.kmer}.fa && \
        for fasta_file in results/{wildcards.sample}/images/{wildcards.seed}_kmer{wildcards.kmer}/*.fasta; do
            fasta_header=$(awk '/^>/ {{print; exit}}' "$fasta_file" | sed 's/^>//') && \
            makeblastdb -in "$fasta_file" -dbtype nucl >> {log} 2>&1 && \
            blastn -db "$fasta_file" -max_target_seqs 1 \
            -query results/{wildcards.sample}/images/{wildcards.seed}_kmer{wildcards.kmer}/Assembled_reads_{wildcards.seed}_kmer{wildcards.kmer}.fa \
            -out results/{wildcards.sample}/images/{wildcards.seed}_kmer{wildcards.kmer}/"$fasta_header".blastn.out \
            -outfmt 6 -evalue 0.00001 -task blastn -num_threads {threads} >> {log} 2>&1
        done
        rm results/{wildcards.sample}/images/{wildcards.seed}_kmer{wildcards.kmer}/Assembled_reads_{wildcards.seed}_kmer{wildcards.kmer}.fa && \
        touch {output}
        """

rule run_recruitment_plot:
    input:
        "results/{sample}/images/{seed}_kmer{kmer}/{seed}_kmer{kmer}.blastn.check"
    output:
        "results/{sample}/images/{seed}_kmer{kmer}/{seed}_kmer{kmer}.recruitment_plot.check"
    threads: 1
    log:
        "logs/{sample}/images/{seed}_{kmer}_run_recruitment_plot.log"
    benchmark:
        "benchmarks/{sample}/images/{seed}_{kmer}_run_recruitment_plot.benchmark"
    params:
        organelle=lambda wildcards: config["samples"][wildcards.sample]["organelle"]
    shell:
        """
        for fasta_file in results/{wildcards.sample}/images/{wildcards.seed}_kmer{wildcards.kmer}/*.fasta; do
            fasta_header=$(awk '/^>/ {{print; exit}}' "$fasta_file" | sed 's/^>//') && \
            python workflow/scripts/images.py --sample {wildcards.sample} \
            --kmer {wildcards.kmer} --seed {wildcards.seed} --organelle {params.organelle} \
            --recruitment_plot --blastn_fasta "$fasta_file" \
            --blastn_out results/{wildcards.sample}/images/{wildcards.seed}_kmer{wildcards.kmer}/"$fasta_header".blastn.out >> {log} 2>&1
        done
        touch {output}
        """

rule run_bwa_index_rotated:
    input:
        "results/{sample}/images/{seed}_kmer{kmer}/{seed}_kmer{kmer}.get_genbank_fastas.check"
    output:
        "results/{sample}/images/{seed}_kmer{kmer}/{seed}_kmer{kmer}.bwa_index.check"
    log:
        "logs/{sample}/images/{kmer}_{seed}_run_bwa_index_rotated.log"
    threads: 1
    benchmark:
        "benchmarks/{sample}/images/{kmer}_{seed}_run_bwa_index_rotated.benchmark"
    singularity:
        f"{config["sif_dir"]}/hic_mapping.sif"
    shell:
        """
        for fasta_file in results/{wildcards.sample}/images/{wildcards.seed}_kmer{wildcards.kmer}/*.fasta; do
            fasta_header=$(awk '/^>/ {{print; exit}}' "$fasta_file" | sed 's/^>//') && \
            if [[ "$fasta_header" == *.rotated ]]; then
                bwa-mem2.avx index $fasta_file >> {log} 2>&1
            fi
        done
        touch {output}
        """

rule run_bwa_mem_rotated:
    input:
        "results/{sample}/images/{seed}_kmer{kmer}/{seed}_kmer{kmer}.bwa_index.check"
    output:
        "results/{sample}/images/{seed}_kmer{kmer}/{seed}_kmer{kmer}.bwa_mem.check"
    log:
        "logs/{sample}/images/{kmer}_{seed}_run_bwa_mem_rotated.log"
    threads: 4
    benchmark:
        "benchmarks/{sample}/images/{kmer}_{seed}_run_bwa_mem_rotated.benchmark"
    singularity:
        f"{config["sif_dir"]}/hic_mapping.sif"
    shell:
        """
        for fasta_file in results/{wildcards.sample}/images/{wildcards.seed}_kmer{wildcards.kmer}/*.fasta; do
            fasta_header=$(awk '/^>/ {{print; exit}}' "$fasta_file" | sed 's/^>//') && \
            if [[ "$fasta_header" == *.rotated ]]; then
                bwa-mem2.avx mem -t {threads} "$fasta_file" \
                results/{wildcards.sample}/novoplasty/{wildcards.seed}/kmer{wildcards.kmer}/Assembled_reads_{wildcards.sample}_R1.fasta \
                results/{wildcards.sample}/novoplasty/{wildcards.seed}/kmer{wildcards.kmer}/Assembled_reads_{wildcards.sample}_R2.fasta | samtools view - -Sb | samtools sort - -@ {threads} \
                -o results/{wildcards.sample}/images/{wildcards.seed}_kmer{wildcards.kmer}/"$fasta_header"_mapping.bam >> {log} 2>&1
            fi
        done
        touch {output}
        """

rule run_samtools_depth:
    input:
        "results/{sample}/pilon/{seed}_kmer{kmer}.bwa_mem.check"
    output:
        "results/{sample}/images/{seed}_kmer{kmer}/{seed}_kmer{kmer}.samtools_depth.check"
    log:
        "logs/{sample}/images/{kmer}_{seed}_run_samtools_depth.log"
    threads: 1
    benchmark:
        "benchmarks/{sample}/images/{kmer}_{seed}_run_samtools_depth.benchmark"
    singularity:
        f"{config["sif_dir"]}/hic_mapping.sif"
    shell:
        """
        for fasta_file in results/{wildcards.sample}/assemblies/{wildcards.seed}_kmer{wildcards.kmer}/*.fasta; do
            fasta_header=$(awk '/^>/ {{print; exit}}' "$fasta_file" | sed 's/^>//') && \
            if [ "$fasta_header" != "INVALIDSEED_1" ]; then
                samtools depth results/{wildcards.sample}/pilon/$fasta_header/"$fasta_header"_mapping.bam \
                > results/{wildcards.sample}/images/{wildcards.seed}_kmer{wildcards.kmer}/"$fasta_header".depth 2>{log}
            fi
        done
        touch {output}
        """

rule run_samtools_depth_rotated:
    input:
        "results/{sample}/images/{seed}_kmer{kmer}/{seed}_kmer{kmer}.bwa_mem.check"
    output:
        "results/{sample}/images/{seed}_kmer{kmer}/{seed}_kmer{kmer}.samtools_depth_rotated.check"
    log:
        "logs/{sample}/images/{kmer}_{seed}_run_samtools_depth_rotated.log"
    threads: 1
    benchmark:
        "benchmarks/{sample}/images/{kmer}_{seed}_run_samtools_depth_rotated.benchmark"
    singularity:
        f"{config["sif_dir"]}/hic_mapping.sif"
    shell:
        """
        for fasta_file in results/{wildcards.sample}/images/{wildcards.seed}_kmer{wildcards.kmer}/*.fasta; do
            fasta_header=$(awk '/^>/ {{print; exit}}' "$fasta_file" | sed 's/^>//') && \
            if [[ "$fasta_header" == *.rotated ]]; then
                samtools depth results/{wildcards.sample}/images/{wildcards.seed}_kmer{wildcards.kmer}/"$fasta_header"_mapping.bam \
                > results/{wildcards.sample}/images/{wildcards.seed}_kmer{wildcards.kmer}/"$fasta_header".depth 2>{log}
            fi
        done
        touch {output}
        """

rule run_depth_plot:
    input:
        "results/{sample}/images/{seed}_kmer{kmer}/{seed}_kmer{kmer}.samtools_depth.check",
        "results/{sample}/images/{seed}_kmer{kmer}/{seed}_kmer{kmer}.samtools_depth_rotated.check"
    output:
        "results/{sample}/images/{seed}_kmer{kmer}/{seed}_kmer{kmer}.depth_plot.check"
    threads: 1
    log:
        "logs/{sample}/images/{seed}_{kmer}_run_depth_plot.log"
    benchmark:
        "benchmarks/{sample}/images/{seed}_{kmer}_run_depth_plot.benchmark"
    params:
        organelle=lambda wildcards: config["samples"][wildcards.sample]["organelle"]
    shell:
        """
        for depth_file in results/{wildcards.sample}/images/{wildcards.seed}_kmer{wildcards.kmer}/*.depth; do
            python workflow/scripts/images.py --depth --depth_bam "$depth_file" --organelle {params.organelle} >> {log} 2>&1
        done
        touch {output}
        """

rule run_ogdraw_mito:
    input:
        lambda wildcards: expand("results/{{sample}}/genbanks/{seed}_kmer{kmer}.mitos2.genbank.check",
            kmer=[kmer for kmer in config["samples"][wildcards.sample]["kmers"]],
            seed=[seed for seed in config["samples"][wildcards.sample]["seeds"]]),
        lambda wildcards: expand("results/{{sample}}/genbanks/{seed}_kmer{kmer}.mitos2.genbank.rotated.check",
            kmer=[kmer for kmer in config["samples"][wildcards.sample]["kmers"]],
            seed=[seed for seed in config["samples"][wildcards.sample]["seeds"]]),
    output:
        "results/{sample}/images/{sample}.mito.ogdraw.check"
    threads: 1
    log:
        "logs/{sample}/images/{sample}_run_ogdraw.log"
    benchmark:
        "benchmarks/{sample}/images/{sample}_run_ogdraw.benchmark"
    singularity:
        f"{config["sif_dir"]}/ogdraw.sif"
    shell:
        """
        for gb_file in results/{wildcards.sample}/genbanks/*.gb; do
            sed -i 's/DEFINITION  ./DEFINITION  Mitochondrion, complete genome./g' $gb_file && \
            base_name=$(basename "$gb_file" .gb) && \
            drawgenemap --infile "$gb_file" --format png --outfile results/{wildcards.sample}/images/"$base_name".png --tidy --useconfig resources/ogd_xml_mitochondrion.xml >> {log} 2>&1
        done
        touch {output}
        """

rule run_ogdraw_chloro:
    input:
        lambda wildcards: expand("results/{{sample}}/genbanks/{seed}_kmer{kmer}.chloe.genbank.check",
            kmer=[kmer for kmer in config["samples"][wildcards.sample]["kmers"]],
            seed=[seed for seed in config["samples"][wildcards.sample]["seeds"]]),
        lambda wildcards: expand("results/{{sample}}/genbanks/{seed}_kmer{kmer}.chloe.genbank.rotated.check",
            kmer=[kmer for kmer in config["samples"][wildcards.sample]["kmers"]],
            seed=[seed for seed in config["samples"][wildcards.sample]["seeds"]]),
        lambda wildcards: expand("results/{{sample}}/genbanks/{seed}_kmer{kmer}.cpgavas2.genbank.check",
            kmer=[kmer for kmer in config["samples"][wildcards.sample]["kmers"]],
            seed=[seed for seed in config["samples"][wildcards.sample]["seeds"]]),
        lambda wildcards: expand("results/{{sample}}/genbanks/{seed}_kmer{kmer}.cpgavas2.genbank.rotated.check",
            kmer=[kmer for kmer in config["samples"][wildcards.sample]["kmers"]],
            seed=[seed for seed in config["samples"][wildcards.sample]["seeds"]]),
    output:
        "results/{sample}/images/{sample}.chloro.ogdraw.check"
    threads: 1
    log:
        "logs/{sample}/images/{sample}_run_ogdraw.log"
    benchmark:
        "benchmarks/{sample}/images/{sample}_run_ogdraw.benchmark"
    singularity:
        f"{config["sif_dir"]}/ogdraw.sif"
    shell:
        """
        for gb_file in results/{wildcards.sample}/genbanks/*.gb; do
            sed -i 's/DEFINITION  Arabidopsis thaliana chloroplast, complete genome./DEFINITION  Chloroplast, complete genome./g' $gb_file >> {log} 2>&1 && \
            base_name=$(basename "$gb_file" .gb) && \
            drawgenemap --infile "$gb_file" --format png --outfile results/{wildcards.sample}/images/"$base_name".png --tidy --useconfig resources/ogd_xml_plastid.xml >> {log} 2>&1
        done
        touch {output}
        """
