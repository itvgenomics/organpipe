rule write_gb_mitos2:
    input:
        "results/{sample}/mitos2/{assembler}/{run}/{assembly}/mitos2.fasta"
    output:
        "results/{sample}/genbanks/{assembler}/{run}/{assembly}.mitos2.genbank.check"
    log:
        "logs/{sample}/genbanks/{assembler}/{run}/{assembly}_write_gb_mitos2.log"
    params:
        genetic_code=sample_param("genetic_code"),
        genbank=lambda wildcards: genbank_path(wildcards, "mitos2"),
        annotation_dir=lambda wildcards, input: os.path.dirname(input[0])
    shell:
        """
        python workflow/scripts/mito_gff2genbank.py --code {params.genetic_code} \
            --annotation_dir {params.annotation_dir} --output {params.genbank} >> {log} 2>&1 && \
        touch {output}
        """

rule write_gb_chloe:
    input:
        check = "results/{sample}/chloe/{assembler}/{run}/{assembly}/chloe.check",
        table2asn = "resources/table2asn.linux64"
    output:
        "results/{sample}/genbanks/{assembler}/{run}/{assembly}.chloe.genbank.check"
    log:
        "logs/{sample}/genbanks/{assembler}/{run}/{assembly}_write_gb_chloe.log"
    params:
        genbank=lambda wildcards: genbank_path(wildcards, "chloe"),
        annotation_dir=lambda wildcards, input: os.path.dirname(input.check)
    shell:
        """
        python workflow/scripts/chloro_gff2genbank.py --software chloe --gene2product resources/gene2product.txt \
            --annotation_dir {params.annotation_dir} --output {params.genbank} >> {log} 2>&1 && \
        touch {output}
        """

rule write_gb_cpgavas2:
    input:
        "results/{sample}/cpgavas2/{assembler}/{run}/{assembly}/cpgavas2.check"
    output:
        "results/{sample}/genbanks/{assembler}/{run}/{assembly}.cpgavas2.genbank.check"
    log:
        "logs/{sample}/genbanks/{assembler}/{run}/{assembly}_write_gb_cpgavas2.log"
    params:
        genbank=lambda wildcards: genbank_path(wildcards, "cpgavas2"),
        annotation_dir=lambda wildcards, input: os.path.dirname(input[0])
    shell:
        """
        python workflow/scripts/chloro_gff2genbank.py --software cpgavas2 \
            --annotation_dir {params.annotation_dir} --output {params.genbank} >> {log} 2>&1 && \
        touch {output}
        """

rule rotate_genbank:
    input:
        "results/{sample}/genbanks/{assembler}/{run}/{assembly}.{annotator}.genbank.check"
    output:
        "results/{sample}/genbanks/{assembler}/{run}/{assembly}.{annotator}.genbank.rotated.check"
    log:
        "logs/{sample}/genbanks/{assembler}/{run}/{assembly}_{annotator}_rotate_genbank.log"
    params:
        organelle=sample_param("organelle"),
        start_gene=lambda wildcards: "rrnS" if wildcards.annotator == "mitos2" else "psbA",
        genbank=lambda wildcards: genbank_path(wildcards, wildcards.annotator),
        rotated=lambda wildcards: genbank_path(wildcards, wildcards.annotator, rotated=True)
    shell:
        """
        if [ -f {params.genbank} ]; then
            python workflow/scripts/rotate_genbank.py --organelle {params.organelle} --start_gene {params.start_gene} \
                --input {params.genbank} --output {params.rotated} >> {log} 2>&1
        fi
        touch {output}
        """
