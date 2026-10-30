rule run_mitos2:
    input:
        "results/{sample}/pilon/{assembler}/{run}/{assembly}/{assembly}.fasta"
    output:
        fasta = "results/{sample}/mitos2/{assembler}/{run}/{assembly}/mitos2.fasta",
        mitos_log = "results/{sample}/mitos2/{assembler}/{run}/{assembly}/mitos.log"
    log:
        "logs/{sample}/mitos2/{assembler}/{run}/{assembly}_run_mitos2.log"
    benchmark:
        "benchmarks/{sample}/mitos2/{assembler}/{run}/{assembly}_run_mitos2.benchmark"
    singularity:
        sif("mitos")
    params:
        genetic_code=sample_param("genetic_code"),
        refseq_dir=config["mitos_db"],
        dir=output_dir
    shell:
        """
        # MITOS2 needs a short sequence name
        awk -v new_header=">temp_header" '/^>/{{$0=new_header}}1' {input} > {output.fasta}

        runmitos --code {params.genetic_code} \
            --input {output.fasta} \
            --outdir {params.dir} \
            -r {params.refseq_dir} --noplots --best \
            > {output.mitos_log} 2>> {log}
        """

rule run_cpgavas2:
    input:
        "results/{sample}/pilon/{assembler}/{run}/{assembly}/{assembly}.fasta"
    output:
        "results/{sample}/cpgavas2/{assembler}/{run}/{assembly}/cpgavas2.check"
    log:
        "logs/{sample}/cpgavas2/{assembler}/{run}/{assembly}_run_cpgavas2.log"
    benchmark:
        "benchmarks/{sample}/cpgavas2/{assembler}/{run}/{assembly}_run_cpgavas2.benchmark"
    singularity:
        sif("cpgavas2")
    params:
        dir=output_dir
    shell:
        """
        sed -i '/maker/s/-quiet/--ignore_nfs_tmp -quiet/' /apps/cpgavas2C/modules/plasAnno/bin/Annotation_Chloroplast_King.py && \
        random_pid=$(( (RANDOM + RANDOM * 32768 + RANDOM * 32768 * 32768) % 999999999 + 1 )) && \
        rm -rf {params.dir} && \
        mkdir -p {params.dir} && \
        run-cpgavas2 -pid $random_pid -in {input} -db 2 >> {log} 2>&1 && \
        mv /tmp/dir_$random_pid/* {params.dir}/ && \
        touch {output}
        """

rule run_chloe:
    input:
        "results/{sample}/pilon/{assembler}/{run}/{assembly}/{assembly}.fasta"
    output:
        "results/{sample}/chloe/{assembler}/{run}/{assembly}/chloe.check"
    log:
        "logs/{sample}/chloe/{assembler}/{run}/{assembly}_run_chloe.log"
    benchmark:
        "benchmarks/{sample}/chloe/{assembler}/{run}/{assembly}_run_chloe.benchmark"
    singularity:
        sif("chloe")
    params:
        dir=output_dir
    shell:
        """
        julia --project=/opt/chloe /opt/chloe/chloe.jl annotate -o {params.dir} {input} >> {log} 2>&1 && \
        touch {output}
        """
