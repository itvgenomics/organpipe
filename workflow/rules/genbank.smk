rule write_genbank_mitos2:
    input:
        "results/{sample}/mitos2/{seed}_kmer{kmer}/{seed}_kmer{kmer}.mitos2.check"
    output:
        "results/{sample}/genbanks/{seed}_kmer{kmer}.mitos2.genbank.check"
    threads: 1
    log:
        "logs/{sample}/genbanks/{sample}_{kmer}_{seed}_write_genbank_mitos2.log"
    params:
        genetic_code=lambda wildcards: config["samples"][wildcards.sample]["genetic_code"]
    shell:
        """
        python workflow/scripts/mito_gff2genbank.py --code {params.genetic_code} \
            --sample {wildcards.sample} --seed {wildcards.seed} --kmer {wildcards.kmer} >> {log} 2>&1 && \
        touch {output}
        """

rule write_genbank_chloe:
    input:
        "results/{sample}/chloe/{seed}_kmer{kmer}.chloe.check"
    output:
        "results/{sample}/genbanks/{seed}_kmer{kmer}.chloe.genbank.check"
    threads: 1
    log:
        "logs/{sample}/genbanks/{sample}_{kmer}_{seed}_write_genbank_chloe.log"
    params:
        genetic_code=lambda wildcards: config["samples"][wildcards.sample]["genetic_code"]
    shell:
        """
        python workflow/scripts/chloro_gff2genbank.py --code {params.genetic_code} \
            --sample {wildcards.sample} --seed {wildcards.seed} \
            --kmer {wildcards.kmer} --software chloe --gene2product resources/gene2product.txt >> {log} 2>&1 && \
        touch {output}
        """

rule write_genbank_cpgavas2:
    input:
        "results/{sample}/cpgavas2/{seed}_kmer{kmer}/{seed}_kmer{kmer}.cpgavas2.check"
    output:
        "results/{sample}/genbanks/{seed}_kmer{kmer}.cpgavas2.genbank.check"
    threads: 1
    log:
        "logs/{sample}/genbanks/{sample}_{kmer}_{seed}_write_genbank_cpgavas2.log"
    params:
        genetic_code=lambda wildcards: config["samples"][wildcards.sample]["genetic_code"]
    shell:
        """
        python workflow/scripts/chloro_gff2genbank.py --code {params.genetic_code} \
            --sample {wildcards.sample} --seed {wildcards.seed} \
            --kmer {wildcards.kmer} --software cpgavas2 >> {log} 2>&1 && \
        touch {output}
        """

rule rotate_mito_genbanks:
    input:
        "results/{sample}/genbanks/{seed}_kmer{kmer}.mitos2.genbank.check"
    output:
        "results/{sample}/genbanks/{seed}_kmer{kmer}.mitos2.genbank.rotated.check"
    threads: 1
    log:
        "logs/{sample}/genbanks/{seed}_kmer{kmer}.genbank_rotate.log"
    params:
        organelle=lambda wildcards: config["samples"][wildcards.sample]["organelle"]
    shell:
        """
        python workflow/scripts/rotate_genbank.py --organelle {params.organelle} \
            --start_gene rrnS --seed {wildcards.seed} --kmer {wildcards.kmer} \
            --sample {wildcards.sample} >> {log} 2>&1 && \
        touch {output}
        """

rule rotate_cpgavas_genbanks:
    input:
        "results/{sample}/genbanks/{seed}_kmer{kmer}.cpgavas2.genbank.check"
    output:
        "results/{sample}/genbanks/{seed}_kmer{kmer}.cpgavas2.genbank.rotated.check"
    threads: 1
    log:
        "logs/{sample}/genbanks/{seed}_kmer{kmer}.genbank.rotated.log"
    params:
        organelle=lambda wildcards: config["samples"][wildcards.sample]["organelle"]
    shell:
        """
        python workflow/scripts/rotate_genbank.py --organelle {params.organelle} \
            --start_gene psbA --seed {wildcards.seed} --kmer {wildcards.kmer} \
            --sample {wildcards.sample} --software cpgavas2 >> {log} 2>&1 && \
        touch {output}
        """

rule rotate_chloe_genbanks:
    input:
        "results/{sample}/genbanks/{seed}_kmer{kmer}.chloe.genbank.check"
    output:
        "results/{sample}/genbanks/{seed}_kmer{kmer}.chloe.genbank.rotated.check"
    threads: 1
    log:
        "logs/{sample}/genbanks/{seed}_kmer{kmer}.genbank.rotated.log"
    params:
        organelle=lambda wildcards: config["samples"][wildcards.sample]["organelle"]
    shell:
        """
        python workflow/scripts/rotate_genbank.py --organelle {params.organelle} \
            --start_gene psbA --seed {wildcards.seed} --kmer {wildcards.kmer} \
            --sample {wildcards.sample} --software chloe >> {log} 2>&1 && \
        touch {output}
        """