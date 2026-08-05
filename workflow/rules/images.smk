rule get_genbank_fastas_novoplasty:
    input:
        lambda wildcards: [
            f"results/{wildcards.sample}/genbanks/novoplasty/{wildcards.seed}_kmer{wildcards.kmer}.mitos2.genbank.rotated.check"
        ] if config["samples"][wildcards.sample]["organelle"] == "mito" else [
            f"results/{wildcards.sample}/genbanks/novoplasty/{wildcards.seed}_kmer{wildcards.kmer}.cpgavas2.genbank.rotated.check",
            f"results/{wildcards.sample}/genbanks/novoplasty/{wildcards.seed}_kmer{wildcards.kmer}.chloe.genbank.rotated.check"
        ]
    output:
        "results/{sample}/images/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.get_genbank_fastas_novoplasty.check"
    log:
        "logs/{sample}/images/novoplasty/{seed}_{kmer}_get_genbank_fastas_novoplasty.log"
    benchmark:
        "benchmarks/{sample}/images/novoplasty/{seed}_{kmer}_get_genbank_fastas_novoplasty.benchmark"
    params:
        organelle=lambda wildcards: config["samples"][wildcards.sample]["organelle"]
    shell:
        """
        python workflow/scripts/images.py --assembler 'novoplasty' --sample {wildcards.sample} \
            --kmer {wildcards.kmer} --seed {wildcards.seed} --parse_gb --organelle {params.organelle} >> {log} 2>&1 && \
        touch {output}
        """

rule run_blastn_novoplasty:
    input:
        "results/{sample}/images/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.get_genbank_fastas_novoplasty.check"
    output:
        "results/{sample}/images/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.blastn.check"
    log:
        "logs/{sample}/images/novoplasty/{seed}_{kmer}_run_blastn_novoplasty.log"
    benchmark:
        "benchmarks/{sample}/images/novoplasty/{seed}_{kmer}_run_blastn_novoplasty.benchmark"
    singularity:
        f"{config["sif_dir"]}/blast.sif"
    shell:
        """
        cat results/{wildcards.sample}/novoplasty/{wildcards.seed}/kmer{wildcards.kmer}/Assembled_reads_{wildcards.sample}_R1.fasta \
        results/{wildcards.sample}/novoplasty/{wildcards.seed}/kmer{wildcards.kmer}/Assembled_reads_{wildcards.sample}_R2.fasta >> \
        results/{wildcards.sample}/images/novoplasty/{wildcards.seed}_kmer{wildcards.kmer}/Assembled_reads_{wildcards.seed}_kmer{wildcards.kmer}.fa && \
        for fasta_file in results/{wildcards.sample}/images/novoplasty/{wildcards.seed}_kmer{wildcards.kmer}/*.fasta; do
            fasta_header=$(awk '/^>/ {{print; exit}}' "$fasta_file" | sed 's/^>//') && \
            makeblastdb -in "$fasta_file" -dbtype nucl >> {log} 2>&1 && \
            blastn -db "$fasta_file" -max_target_seqs 1 \
            -query results/{wildcards.sample}/images/novoplasty/{wildcards.seed}_kmer{wildcards.kmer}/Assembled_reads_{wildcards.seed}_kmer{wildcards.kmer}.fa \
            -out results/{wildcards.sample}/images/novoplasty/{wildcards.seed}_kmer{wildcards.kmer}/"$fasta_header".blastn.out \
            -outfmt 6 -evalue 0.00001 -task blastn -num_threads {threads} >> {log} 2>&1
        done
        rm results/{wildcards.sample}/images/novoplasty/{wildcards.seed}_kmer{wildcards.kmer}/Assembled_reads_{wildcards.seed}_kmer{wildcards.kmer}.fa && \
        touch {output}
        """

rule run_recruitment_plot_novoplasty:
    input:
        "results/{sample}/images/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.blastn.check"
    output:
        "results/{sample}/images/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.recruitment_plot.check"
    log:
        "logs/{sample}/images/novoplasty/{seed}_{kmer}_run_recruitment_plot_novoplasty.log"
    benchmark:
        "benchmarks/{sample}/images/novoplasty/{seed}_{kmer}_run_recruitment_plot_novoplasty.benchmark"
    params:
        organelle=lambda wildcards: config["samples"][wildcards.sample]["organelle"]
    shell:
        """
        for fasta_file in results/{wildcards.sample}/images/novoplasty/{wildcards.seed}_kmer{wildcards.kmer}/*.fasta; do
            fasta_header=$(awk '/^>/ {{print; exit}}' "$fasta_file" | sed 's/^>//') && \
            python workflow/scripts/images.py --assembler 'novoplasty' --sample {wildcards.sample} \
            --kmer {wildcards.kmer} --seed {wildcards.seed} --organelle {params.organelle} \
            --recruitment_plot --blastn_fasta "$fasta_file" \
            --blastn_out results/{wildcards.sample}/images/novoplasty/{wildcards.seed}_kmer{wildcards.kmer}/"$fasta_header".blastn.out >> {log} 2>&1
        done
        touch {output}
        """

rule run_bwa_index_rotated_novoplasty:
    input:
        "results/{sample}/images/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.get_genbank_fastas_novoplasty.check"
    output:
        "results/{sample}/images/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.bwa_index.check"
    log:
        "logs/{sample}/images/novoplasty/{kmer}_{seed}_run_bwa_index_rotated_novoplasty.log"
    benchmark:
        "benchmarks/{sample}/images/novoplasty/{kmer}_{seed}_run_bwa_index_rotated_novoplasty.benchmark"
    singularity:
        f"{config["sif_dir"]}/hic_mapping.sif"
    shell:
        """
        for fasta_file in results/{wildcards.sample}/images/novoplasty/{wildcards.seed}_kmer{wildcards.kmer}/*.fasta; do
            fasta_header=$(awk '/^>/ {{print; exit}}' "$fasta_file" | sed 's/^>//') && \
            if [[ "$fasta_header" == *.rotated ]]; then
                bwa-mem2.avx index $fasta_file >> {log} 2>&1
            fi
        done
        touch {output}
        """

rule run_bwa_mem_rotated_novoplasty:
    input:
        "results/{sample}/images/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.bwa_index.check"
    output:
        "results/{sample}/images/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.bwa_mem.check"
    log:
        "logs/{sample}/images/novoplasty/{kmer}_{seed}_run_bwa_mem_rotated_novoplasty.log"
    benchmark:
        "benchmarks/{sample}/images/novoplasty/{kmer}_{seed}_run_bwa_mem_rotated_novoplasty.benchmark"
    singularity:
        f"{config["sif_dir"]}/hic_mapping.sif"
    shell:
        """
        for fasta_file in results/{wildcards.sample}/images/novoplasty/{wildcards.seed}_kmer{wildcards.kmer}/*.fasta; do
            fasta_header=$(awk '/^>/ {{print; exit}}' "$fasta_file" | sed 's/^>//') && \
            if [[ "$fasta_header" == *.rotated ]]; then
                {{ bwa-mem2.avx mem -t {threads} "$fasta_file" \
                results/{wildcards.sample}/novoplasty/{wildcards.seed}/kmer{wildcards.kmer}/Assembled_reads_{wildcards.sample}_R1.fasta \
                results/{wildcards.sample}/novoplasty/{wildcards.seed}/kmer{wildcards.kmer}/Assembled_reads_{wildcards.sample}_R2.fasta 2>> {log}
                }} | samtools view - -Sb | samtools sort - -@ {threads} \
                -o results/{wildcards.sample}/images/novoplasty/{wildcards.seed}_kmer{wildcards.kmer}/"$fasta_header"_mapping.bam >> {log} 2>&1
            fi
        done
        touch {output}
        """

rule run_samtools_depth_novoplasty:
    input:
        "results/{sample}/pilon/novoplasty/{seed}_kmer{kmer}.bwa_mem.check"
    output:
        "results/{sample}/images/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.samtools_depth.check"
    log:
        "logs/{sample}/images/novoplasty/{kmer}_{seed}_run_samtools_depth_novoplasty.log"
    benchmark:
        "benchmarks/{sample}/images/novoplasty/{kmer}_{seed}_run_samtools_depth_novoplasty.benchmark"
    singularity:
        f"{config["sif_dir"]}/hic_mapping.sif"
    shell:
        """
        for fasta_file in results/{wildcards.sample}/assemblies/novoplasty/{wildcards.seed}_kmer{wildcards.kmer}/*.fasta; do
            fasta_header=$(awk '/^>/ {{print; exit}}' "$fasta_file" | sed 's/^>//') && \
            if [ "$fasta_header" != "INVALIDSEED_1" ]; then
                samtools depth results/{wildcards.sample}/pilon/novoplasty/$fasta_header/"$fasta_header"_mapping.bam \
                > results/{wildcards.sample}/images/novoplasty/{wildcards.seed}_kmer{wildcards.kmer}/"$fasta_header".depth 2>{log}
            fi
        done
        touch {output}
        """

rule run_samtools_depth_novoplasty_rotated:
    input:
        "results/{sample}/images/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.bwa_mem.check"
    output:
        "results/{sample}/images/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.samtools_depth_rotated.check"
    log:
        "logs/{sample}/images/novoplasty/{kmer}_{seed}_run_samtools_depth_novoplasty_rotated.log"
    benchmark:
        "benchmarks/{sample}/images/novoplasty/{kmer}_{seed}_run_samtools_depth_novoplasty_rotated.benchmark"
    singularity:
        f"{config["sif_dir"]}/hic_mapping.sif"
    shell:
        """
        for fasta_file in results/{wildcards.sample}/images/novoplasty/{wildcards.seed}_kmer{wildcards.kmer}/*.fasta; do
            fasta_header=$(awk '/^>/ {{print; exit}}' "$fasta_file" | sed 's/^>//') && \
            if [[ "$fasta_header" == *.rotated ]]; then
                samtools depth results/{wildcards.sample}/images/novoplasty/{wildcards.seed}_kmer{wildcards.kmer}/"$fasta_header"_mapping.bam \
                > results/{wildcards.sample}/images/novoplasty/{wildcards.seed}_kmer{wildcards.kmer}/"$fasta_header".depth 2>{log}
            fi
        done
        touch {output}
        """

rule run_depth_plot_novoplasty:
    input:
        "results/{sample}/images/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.samtools_depth.check",
        "results/{sample}/images/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.samtools_depth_rotated.check"
    output:
        "results/{sample}/images/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.depth_plot.check"
    log:
        "logs/{sample}/images/novoplasty/{seed}_{kmer}_run_depth_plot_novoplasty.log"
    benchmark:
        "benchmarks/{sample}/images/novoplasty/{seed}_{kmer}_run_depth_plot_novoplasty.benchmark"
    params:
        organelle=lambda wildcards: config["samples"][wildcards.sample]["organelle"]
    shell:
        """
        for depth_file in results/{wildcards.sample}/images/novoplasty/{wildcards.seed}_kmer{wildcards.kmer}/*.depth; do
            python workflow/scripts/images.py --assembler 'novoplasty' --depth --depth_bam "$depth_file" --organelle {params.organelle} >> {log} 2>&1
        done
        touch {output}
        """

rule run_ogdraw_mito_novoplasty:
    input:
        lambda wildcards: expand("results/{{sample}}/genbanks/novoplasty/{seed}_kmer{kmer}.mitos2.genbank.check",
            kmer=[kmer for kmer in config["samples"][wildcards.sample]["kmers"]],
            seed=[seed for seed in config["samples"][wildcards.sample]["seeds"]]),
        lambda wildcards: expand("results/{{sample}}/genbanks/novoplasty/{seed}_kmer{kmer}.mitos2.genbank.rotated.check",
            kmer=[kmer for kmer in config["samples"][wildcards.sample]["kmers"]],
            seed=[seed for seed in config["samples"][wildcards.sample]["seeds"]]),
    output:
        "results/{sample}/images/novoplasty/{sample}.mito.ogdraw.check"
    log:
        "logs/{sample}/images/novoplasty/{sample}_run_ogdraw.log"
    benchmark:
        "benchmarks/{sample}/images/novoplasty/{sample}_run_ogdraw.benchmark"
    singularity:
        f"{config["sif_dir"]}/ogdraw.sif"
    shell:
        """
        for gb_file in results/{wildcards.sample}/genbanks/novoplasty/*.gb; do
            sed -i 's/DEFINITION  ./DEFINITION  Mitochondrion, complete genome./g' $gb_file && \
            base_name=$(basename "$gb_file" .gb) && \
            drawgenemap --infile "$gb_file" --format png --outfile results/{wildcards.sample}/images/novoplasty/"$base_name".png --tidy --useconfig resources/ogd_xml_mitochondrion.xml >> {log} 2>&1
        done
        touch {output}
        """

rule run_ogdraw_chloro_novoplasty:
    input:
        lambda wildcards: expand("results/{{sample}}/genbanks/novoplasty/{seed}_kmer{kmer}.chloe.genbank.check",
            kmer=[kmer for kmer in config["samples"][wildcards.sample]["kmers"]],
            seed=[seed for seed in config["samples"][wildcards.sample]["seeds"]]),
        lambda wildcards: expand("results/{{sample}}/genbanks/novoplasty/{seed}_kmer{kmer}.chloe.genbank.rotated.check",
            kmer=[kmer for kmer in config["samples"][wildcards.sample]["kmers"]],
            seed=[seed for seed in config["samples"][wildcards.sample]["seeds"]]),
        lambda wildcards: expand("results/{{sample}}/genbanks/novoplasty/{seed}_kmer{kmer}.cpgavas2.genbank.check",
            kmer=[kmer for kmer in config["samples"][wildcards.sample]["kmers"]],
            seed=[seed for seed in config["samples"][wildcards.sample]["seeds"]]),
        lambda wildcards: expand("results/{{sample}}/genbanks/novoplasty/{seed}_kmer{kmer}.cpgavas2.genbank.rotated.check",
            kmer=[kmer for kmer in config["samples"][wildcards.sample]["kmers"]],
            seed=[seed for seed in config["samples"][wildcards.sample]["seeds"]]),
    output:
        "results/{sample}/images/novoplasty/{sample}.chloro.ogdraw.check"
    log:
        "logs/{sample}/images/novoplasty/{sample}_run_ogdraw.log"
    benchmark:
        "benchmarks/{sample}/images/novoplasty/{sample}_run_ogdraw.benchmark"
    singularity:
        f"{config["sif_dir"]}/ogdraw.sif"
    shell:
        """
        for gb_file in results/{wildcards.sample}/genbanks/novoplasty/*.gb; do
            sed -i 's/DEFINITION  Arabidopsis thaliana chloroplast, complete genome./DEFINITION  Chloroplast, complete genome./g' $gb_file >> {log} 2>&1 && \
            base_name=$(basename "$gb_file" .gb) && \
            drawgenemap --infile "$gb_file" --format png --outfile results/{wildcards.sample}/images/novoplasty/"$base_name".png --tidy --useconfig resources/ogd_xml_plastid.xml >> {log} 2>&1
        done
        touch {output}
        """


rule get_genbank_fastas_getorganelle:
    input:
        lambda wildcards: [
            f"results/{wildcards.sample}/genbanks/getorganelle/mitos2.genbank.rotated.check"
        ] if config["samples"][wildcards.sample]["organelle"] == "mito" else [
            f"results/{wildcards.sample}/genbanks/getorganelle/cpgavas2.genbank.rotated.check",
            f"results/{wildcards.sample}/genbanks/getorganelle/chloe.genbank.rotated.check"
        ]
    output:
        "results/{sample}/images/getorganelle/get_genbank_fastas_getorganelle.check"
    log:
        "logs/{sample}/images/getorganelle/get_genbank_fastas_getorganelle.log"
    benchmark:
        "benchmarks/{sample}/images/getorganelle/get_genbank_fastas_getorganelle.benchmark"
    params:
        organelle=lambda wildcards: config["samples"][wildcards.sample]["organelle"]
    shell:
        """
        python workflow/scripts/images.py --sample {wildcards.sample} \
            --assembler 'getorganelle' --parse_gb --organelle {params.organelle} >> {log} 2>&1 && \
        touch {output}
        """

rule run_blastn_getorganelle:
    input:
        "results/{sample}/images/getorganelle/get_genbank_fastas_getorganelle.check"
    output:
        "results/{sample}/images/getorganelle/blastn.check"
    log:
        "logs/{sample}/images/getorganelle/run_blastn_getorganelle.log"
    benchmark:
        "benchmarks/{sample}/images/getorganelle/run_blastn_getorganelle.benchmark"
    singularity:
        f"{config["sif_dir"]}/blast.sif"
    shell:
        """
        perl workflow/scripts/get_fastq.pl -1 results/{wildcards.sample}/getorganelle/extended_1_paired.fq \
        -2 results/{wildcards.sample}/getorganelle/extended_2_paired.fq \
        -out results/{wildcards.sample}/images/getorganelle/concat_reads.fa && \
        for fasta_file in results/{wildcards.sample}/images/getorganelle/*.fasta; do
            fasta_header=$(awk '/^>/ {{print; exit}}' "$fasta_file" | sed 's/^>//') && \
            makeblastdb -in "$fasta_file" -dbtype nucl >> {log} 2>&1 && \
            blastn -db "$fasta_file" -max_target_seqs 1 \
            -query results/{wildcards.sample}/images/getorganelle/concat_reads.fa \
            -out results/{wildcards.sample}/images/getorganelle/"$fasta_header".blastn.out \
            -outfmt 6 -evalue 0.00001 -task blastn -num_threads {threads} >> {log} 2>&1
        done
        rm results/{wildcards.sample}/images/getorganelle/concat_reads.fa && \
        touch {output}
        """

rule run_recruitment_plot_getorganelle:
    input:
        "results/{sample}/images/getorganelle/blastn.check"
    output:
        "results/{sample}/images/getorganelle/recruitment_plot.check"
    log:
        "logs/{sample}/images/getorganelle/run_recruitment_plot_getorganelle.log"
    benchmark:
        "benchmarks/{sample}/images/getorganelle/run_recruitment_plot_getorganelle.benchmark"
    params:
        organelle=lambda wildcards: config["samples"][wildcards.sample]["organelle"]
    shell:
        """
        for fasta_file in results/{wildcards.sample}/images/getorganelle/*.fasta; do
            fasta_header=$(awk '/^>/ {{print; exit}}' "$fasta_file" | sed 's/^>//') && \
            python workflow/scripts/images.py --sample {wildcards.sample} \
            --assembler 'getorganelle' --organelle {params.organelle} \
            --recruitment_plot --blastn_fasta "$fasta_file" \
            --blastn_out results/{wildcards.sample}/images/getorganelle/"$fasta_header".blastn.out >> {log} 2>&1
        done
        touch {output}
        """

rule run_bwa_index_rotated_getorganelle:
    input:
        "results/{sample}/images/getorganelle/get_genbank_fastas_getorganelle.check"
    output:
        "results/{sample}/images/getorganelle/bwa_index.check"
    log:
        "logs/{sample}/images/getorganelle/run_bwa_index_rotated_getorganelle.log"
    benchmark:
        "benchmarks/{sample}/images/getorganelle/run_bwa_index_rotated_getorganelle.benchmark"
    singularity:
        f"{config["sif_dir"]}/hic_mapping.sif"
    shell:
        """
        for fasta_file in results/{wildcards.sample}/images/getorganelle/*.fasta; do
            fasta_header=$(awk '/^>/ {{print; exit}}' "$fasta_file" | sed 's/^>//') && \
            if [[ "$fasta_header" == *.rotated ]]; then
                bwa-mem2.avx index $fasta_file >> {log} 2>&1
            fi
        done
        touch {output}
        """

rule run_bwa_mem_rotated_getorganelle:
    input:
        "results/{sample}/images/getorganelle/bwa_index.check"
    output:
        "results/{sample}/images/getorganelle/bwa_mem.check"
    log:
        "logs/{sample}/images/getorganelle/run_bwa_mem_rotated_getorganelle.log"
    benchmark:
        "benchmarks/{sample}/images/getorganelle/run_bwa_mem_rotated_getorganelle.benchmark"
    singularity:
        f"{config["sif_dir"]}/hic_mapping.sif"
    shell:
        """
        for fasta_file in results/{wildcards.sample}/images/getorganelle/*.fasta; do
            fasta_header=$(awk '/^>/ {{print; exit}}' "$fasta_file" | sed 's/^>//') && \
            if [[ "$fasta_header" == *.rotated ]]; then
                {{ bwa-mem2.avx mem -t {threads} "$fasta_file" \
                results/{wildcards.sample}/getorganelle/extended_1_paired.fq \
                results/{wildcards.sample}/getorganelle/extended_2_paired.fq 2>> {log}
                }} | samtools view - -Sb | samtools sort - -@ {threads} \
                -o results/{wildcards.sample}/images/getorganelle/"$fasta_header"_mapping.bam >> {log} 2>&1
            fi
        done
        touch {output}
        """

rule run_samtools_depth_getorganelle:
    input:
        "results/{sample}/pilon/getorganelle/getorganelle_bwa_mem.check"
    output:
        "results/{sample}/images/getorganelle/samtools_depth.check"
    log:
        "logs/{sample}/images/getorganelle/run_samtools_depth_getorganelle.log"
    benchmark:
        "benchmarks/{sample}/images/getorganelle/run_samtools_depth_getorganelle.benchmark"
    singularity:
        f"{config["sif_dir"]}/hic_mapping.sif"
    shell:
        """
        for fasta_file in results/{wildcards.sample}/assemblies/getorganelle/*.fasta; do
            fasta_header=$(awk '/^>/ {{print; exit}}' "$fasta_file" | sed 's/^>//') && \
            if [ "$fasta_header" != "INVALIDSEED_1" ]; then
                samtools depth results/{wildcards.sample}/pilon/getorganelle/$fasta_header/"$fasta_header"_mapping.bam \
                > results/{wildcards.sample}/images/getorganelle/"$fasta_header".depth 2>{log}
            fi
        done
        touch {output}
        """

rule run_samtools_depth_getorganelle_rotated:
    input:
        "results/{sample}/images/getorganelle/bwa_mem.check"
    output:
        "results/{sample}/images/getorganelle/samtools_depth_rotated.check"
    log:
        "logs/{sample}/images/getorganelle/run_samtools_depth_getorganelle_rotated.log"
    benchmark:
        "benchmarks/{sample}/images/getorganelle/run_samtools_depth_getorganelle_rotated.benchmark"
    singularity:
        f"{config["sif_dir"]}/hic_mapping.sif"
    shell:
        """
        for fasta_file in results/{wildcards.sample}/images/getorganelle/*.fasta; do
            fasta_header=$(awk '/^>/ {{print; exit}}' "$fasta_file" | sed 's/^>//') && \
            if [[ "$fasta_header" == *.rotated ]]; then
                samtools depth results/{wildcards.sample}/images/getorganelle/"$fasta_header"_mapping.bam \
                > results/{wildcards.sample}/images/getorganelle/"$fasta_header".depth 2>{log}
            fi
        done
        touch {output}
        """

rule run_depth_plot_getorganelle:
    input:
        "results/{sample}/images/getorganelle/samtools_depth.check",
        "results/{sample}/images/getorganelle/samtools_depth_rotated.check"
    output:
        "results/{sample}/images/getorganelle/depth_plot.check"
    log:
        "logs/{sample}/images/getorganelle/run_depth_plot_getorganelle.log"
    benchmark:
        "benchmarks/{sample}/images/getorganelle/run_depth_plot_getorganelle.benchmark"
    params:
        organelle=lambda wildcards: config["samples"][wildcards.sample]["organelle"]
    shell:
        """
        for depth_file in results/{wildcards.sample}/images/getorganelle/*.depth; do
            python workflow/scripts/images.py --assembler 'getorganelle' --depth --depth_bam "$depth_file" --organelle {params.organelle} >> {log} 2>&1
        done
        touch {output}
        """

rule run_ogdraw_mito_getorganelle:
    input:
        lambda wildcards: expand("results/{{sample}}/genbanks/getorganelle/mitos2.genbank.check",
            kmer=[kmer for kmer in config["samples"][wildcards.sample]["kmers"]],
            seed=[seed for seed in config["samples"][wildcards.sample]["seeds"]]),
        lambda wildcards: expand("results/{{sample}}/genbanks/getorganelle/mitos2.genbank.rotated.check",
            kmer=[kmer for kmer in config["samples"][wildcards.sample]["kmers"]],
            seed=[seed for seed in config["samples"][wildcards.sample]["seeds"]]),
    output:
        "results/{sample}/images/getorganelle/mito.ogdraw.check"
    log:
        "logs/{sample}/images/getorganelle/run_ogdraw.log"
    benchmark:
        "benchmarks/{sample}/images/getorganelle/run_ogdraw.benchmark"
    singularity:
        f"{config["sif_dir"]}/ogdraw.sif"
    shell:
        """
        for gb_file in results/{wildcards.sample}/genbanks/getorganelle/*.gb; do
            sed -i 's/DEFINITION  ./DEFINITION  Mitochondrion, complete genome./g' $gb_file && \
            base_name=$(basename "$gb_file" .gb) && \
            drawgenemap --infile "$gb_file" --format png --outfile results/{wildcards.sample}/images/getorganelle/"$base_name".png --tidy --useconfig resources/ogd_xml_mitochondrion.xml >> {log} 2>&1
        done
        touch {output}
        """

rule run_ogdraw_chloro_getorganelle:
    input:
        lambda wildcards: expand("results/{{sample}}/genbanks/getorganelle/chloe.genbank.check",
            kmer=[kmer for kmer in config["samples"][wildcards.sample]["kmers"]],
            seed=[seed for seed in config["samples"][wildcards.sample]["seeds"]]),
        lambda wildcards: expand("results/{{sample}}/genbanks/getorganelle/chloe.genbank.rotated.check",
            kmer=[kmer for kmer in config["samples"][wildcards.sample]["kmers"]],
            seed=[seed for seed in config["samples"][wildcards.sample]["seeds"]]),
        lambda wildcards: expand("results/{{sample}}/genbanks/getorganelle/cpgavas2.genbank.check",
            kmer=[kmer for kmer in config["samples"][wildcards.sample]["kmers"]],
            seed=[seed for seed in config["samples"][wildcards.sample]["seeds"]]),
        lambda wildcards: expand("results/{{sample}}/genbanks/getorganelle/cpgavas2.genbank.rotated.check",
            kmer=[kmer for kmer in config["samples"][wildcards.sample]["kmers"]],
            seed=[seed for seed in config["samples"][wildcards.sample]["seeds"]]),
    output:
        "results/{sample}/images/getorganelle/chloro.ogdraw.check"
    log:
        "logs/{sample}/images/getorganelle/run_ogdraw.log"
    benchmark:
        "benchmarks/{sample}/images/getorganelle/run_ogdraw.benchmark"
    singularity:
        f"{config["sif_dir"]}/ogdraw.sif"
    shell:
        """
        for gb_file in results/{wildcards.sample}/genbanks/getorganelle/*.gb; do
            sed -i 's/DEFINITION  Arabidopsis thaliana chloroplast, complete genome./DEFINITION  Chloroplast, complete genome./g' $gb_file >> {log} 2>&1 && \
            base_name=$(basename "$gb_file" .gb) && \
            drawgenemap --infile "$gb_file" --format png --outfile results/{wildcards.sample}/images/getorganelle/"$base_name".png --tidy --useconfig resources/ogd_xml_plastid.xml >> {log} 2>&1
        done
        touch {output}
        """
