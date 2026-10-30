# The plots use the MITOS2 GenBank files of mitochondria and the CPGAVAS2 ones of chloroplasts.
def plot_annotator(wildcards):
    return annotators(wildcards.sample)[0]

def plot_base(wildcards):
    return wildcards.assembly + ANNOTATOR_SUFFIX[plot_annotator(wildcards)]

def ogdraw_genbanks(wildcards):
    genbanks = []
    for annotator in annotators(wildcards.sample):
        genbanks += [genbank_path(wildcards, annotator), genbank_path(wildcards, annotator, rotated=True)]
    return genbanks

rule get_genbank_fastas:
    input:
        lambda wildcards: f"results/{wildcards.sample}/genbanks/{wildcards.assembler}/{wildcards.run}/{wildcards.assembly}.{plot_annotator(wildcards)}.genbank.rotated.check"
    output:
        "results/{sample}/images/{assembler}/{run}/{assembly}/get_genbank_fastas.check"
    log:
        "logs/{sample}/images/{assembler}/{run}/{assembly}_get_genbank_fastas.log"
    benchmark:
        "benchmarks/{sample}/images/{assembler}/{run}/{assembly}_get_genbank_fastas.benchmark"
    params:
        genbanks=lambda wildcards: [genbank_path(wildcards, plot_annotator(wildcards)), genbank_path(wildcards, plot_annotator(wildcards), rotated=True)],
        dir=output_dir
    shell:
        """
        python workflow/scripts/images.py --parse_gb --genbanks {params.genbanks} --fasta_dir {params.dir} >> {log} 2>&1 && \
        touch {output}
        """

rule run_blastn:
    input:
        unpack(assembled_reads),
        check = "results/{sample}/images/{assembler}/{run}/{assembly}/get_genbank_fastas.check"
    output:
        "results/{sample}/images/{assembler}/{run}/{assembly}/blastn.check"
    log:
        "logs/{sample}/images/{assembler}/{run}/{assembly}_run_blastn.log"
    benchmark:
        "benchmarks/{sample}/images/{assembler}/{run}/{assembly}_run_blastn.benchmark"
    singularity:
        sif("blast")
    params:
        dir=output_dir
    shell:
        """
        shopt -s nullglob
        reads={params.dir}/assembled_reads.fa
        if [ {wildcards.assembler} = getorganelle ]; then
            perl workflow/scripts/get_fastq.pl -1 {input.r1} -2 {input.r2} -out "$reads"
        else
            cat {input.r1} {input.r2} > "$reads"
        fi
        for fasta_file in {params.dir}/*.fasta; do
            fasta_header=$(awk '/^>/ {{print; exit}}' "$fasta_file" | sed 's/^>//')
            makeblastdb -in "$fasta_file" -dbtype nucl >> {log} 2>&1
            blastn -db "$fasta_file" -max_target_seqs 1 \
            -query "$reads" \
            -out {params.dir}/"$fasta_header".blastn.out \
            -outfmt 6 -evalue 0.00001 -task blastn -num_threads {threads} >> {log} 2>&1
        done
        rm "$reads"
        touch {output}
        """

rule run_recruitment_plot:
    input:
        "results/{sample}/images/{assembler}/{run}/{assembly}/blastn.check"
    output:
        "results/{sample}/images/{assembler}/{run}/{assembly}/recruitment_plot.check"
    log:
        "logs/{sample}/images/{assembler}/{run}/{assembly}_run_recruitment_plot.log"
    benchmark:
        "benchmarks/{sample}/images/{assembler}/{run}/{assembly}_run_recruitment_plot.benchmark"
    params:
        dir=output_dir
    shell:
        """
        shopt -s nullglob
        for fasta_file in {params.dir}/*.fasta; do
            fasta_header=$(awk '/^>/ {{print; exit}}' "$fasta_file" | sed 's/^>//')
            python workflow/scripts/images.py --recruitment_plot --blastn_fasta "$fasta_file" \
            --blastn_out {params.dir}/"$fasta_header".blastn.out >> {log} 2>&1
        done
        touch {output}
        """

rule run_bwa_index_rotated:
    input:
        "results/{sample}/images/{assembler}/{run}/{assembly}/get_genbank_fastas.check"
    output:
        "results/{sample}/images/{assembler}/{run}/{assembly}/bwa_index.check"
    log:
        "logs/{sample}/images/{assembler}/{run}/{assembly}_run_bwa_index_rotated.log"
    benchmark:
        "benchmarks/{sample}/images/{assembler}/{run}/{assembly}_run_bwa_index_rotated.benchmark"
    singularity:
        sif("hic_mapping")
    params:
        rotated=lambda wildcards, output: f"{os.path.dirname(output[0])}/{plot_base(wildcards)}.rotated"
    shell:
        """
        if [ -f {params.rotated}.fasta ]; then
            bwa-mem2.avx index {params.rotated}.fasta >> {log} 2>&1
        fi
        touch {output}
        """

rule run_bwa_mem_rotated:
    input:
        unpack(assembled_reads),
        check = "results/{sample}/images/{assembler}/{run}/{assembly}/bwa_index.check"
    output:
        "results/{sample}/images/{assembler}/{run}/{assembly}/bwa_mem.check"
    log:
        "logs/{sample}/images/{assembler}/{run}/{assembly}_run_bwa_mem_rotated.log"
    benchmark:
        "benchmarks/{sample}/images/{assembler}/{run}/{assembly}_run_bwa_mem_rotated.benchmark"
    singularity:
        sif("hic_mapping")
    params:
        rotated=lambda wildcards, output: f"{os.path.dirname(output[0])}/{plot_base(wildcards)}.rotated"
    shell:
        """
        if [ -f {params.rotated}.fasta ]; then
            {{ bwa-mem2.avx mem -t {threads} {params.rotated}.fasta {input.r1} {input.r2} 2>> {log}
            }} | samtools view - -Sb | samtools sort - -@ {threads} \
            -o {params.rotated}_mapping.bam >> {log} 2>&1
        fi
        touch {output}
        """

rule run_samtools_depth:
    input:
        "results/{sample}/pilon/{assembler}/{run}/{assembly}/{assembly}_mapping.bam"
    output:
        "results/{sample}/images/{assembler}/{run}/{assembly}/{assembly}.depth"
    log:
        "logs/{sample}/images/{assembler}/{run}/{assembly}_run_samtools_depth.log"
    benchmark:
        "benchmarks/{sample}/images/{assembler}/{run}/{assembly}_run_samtools_depth.benchmark"
    singularity:
        sif("hic_mapping")
    shell:
        """
        samtools depth {input} > {output} 2>> {log}
        """

rule run_samtools_depth_rotated:
    input:
        "results/{sample}/images/{assembler}/{run}/{assembly}/bwa_mem.check"
    output:
        "results/{sample}/images/{assembler}/{run}/{assembly}/samtools_depth_rotated.check"
    log:
        "logs/{sample}/images/{assembler}/{run}/{assembly}_run_samtools_depth_rotated.log"
    benchmark:
        "benchmarks/{sample}/images/{assembler}/{run}/{assembly}_run_samtools_depth_rotated.benchmark"
    singularity:
        sif("hic_mapping")
    params:
        rotated=lambda wildcards, output: f"{os.path.dirname(output[0])}/{plot_base(wildcards)}.rotated"
    shell:
        """
        if [ -f {params.rotated}_mapping.bam ]; then
            samtools depth {params.rotated}_mapping.bam > {params.rotated}.depth 2>> {log}
        fi
        touch {output}
        """

rule run_depth_plot:
    input:
        depth = "results/{sample}/images/{assembler}/{run}/{assembly}/{assembly}.depth",
        rotated = "results/{sample}/images/{assembler}/{run}/{assembly}/samtools_depth_rotated.check"
    output:
        "results/{sample}/images/{assembler}/{run}/{assembly}/depth_plot.check"
    log:
        "logs/{sample}/images/{assembler}/{run}/{assembly}_run_depth_plot.log"
    benchmark:
        "benchmarks/{sample}/images/{assembler}/{run}/{assembly}_run_depth_plot.benchmark"
    params:
        base=lambda wildcards, output: f"{os.path.dirname(output[0])}/{plot_base(wildcards)}"
    shell:
        """
        python workflow/scripts/images.py --depth --depth_bam {input.depth} --depth_fasta {params.base}.fasta >> {log} 2>&1
        if [ -f {params.base}.rotated.depth ]; then
            python workflow/scripts/images.py --depth --depth_bam {params.base}.rotated.depth --depth_fasta {params.base}.rotated.fasta >> {log} 2>&1
        fi
        touch {output}
        """

rule run_ogdraw:
    input:
        lambda wildcards: expand("results/{{sample}}/genbanks/{{assembler}}/{{run}}/{{assembly}}.{annotator}.genbank.rotated.check", annotator=annotators(wildcards.sample))
    output:
        "results/{sample}/images/{assembler}/{run}/{assembly}/ogdraw.check"
    log:
        "logs/{sample}/images/{assembler}/{run}/{assembly}_run_ogdraw.log"
    benchmark:
        "benchmarks/{sample}/images/{assembler}/{run}/{assembly}_run_ogdraw.benchmark"
    singularity:
        sif("ogdraw")
    params:
        genbanks=ogdraw_genbanks,
        ogdraw_config=lambda wildcards: "resources/ogd_xml_mitochondrion.xml" if config["samples"][wildcards.sample]["organelle"] == "mito" else "resources/ogd_xml_plastid.xml",
        definition=lambda wildcards: "'s/DEFINITION  ./DEFINITION  Mitochondrion, complete genome./g'" if config["samples"][wildcards.sample]["organelle"] == "mito" else "'s/DEFINITION  Arabidopsis thaliana chloroplast, complete genome./DEFINITION  Chloroplast, complete genome./g'",
        dir=output_dir
    shell:
        """
        for gb_file in {params.genbanks}; do
            if [ -f "$gb_file" ]; then
                sed -i {params.definition} "$gb_file"
                base_name=$(basename "$gb_file" .gb)
                drawgenemap --infile "$gb_file" --format png --outfile {params.dir}/"$base_name".png --tidy --useconfig {params.ogdraw_config} >> {log} 2>&1
            fi
        done
        touch {output}
        """
