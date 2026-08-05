rule write_gb_novoplasty_mitos2:
    input:
        "results/{sample}/mitos2/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.mitos2.check"
    output:
        "results/{sample}/genbanks/novoplasty/{seed}_kmer{kmer}.mitos2.genbank.check"
    log:
        "logs/{sample}/genbanks/novoplasty/{sample}_{kmer}_{seed}_write_gb_novoplasty_mitos2.log"
    params:
        genetic_code=lambda wildcards: config["samples"][wildcards.sample]["genetic_code"]
    shell:
        """
        python workflow/scripts/mito_gff2genbank.py --code {params.genetic_code} --assembler 'novoplasty' \
            --sample {wildcards.sample} --seed {wildcards.seed} --kmer {wildcards.kmer} >> {log} 2>&1 && \
        touch {output}
        """

rule write_gb_novoplasty_chloe:
    input:
        "results/{sample}/chloe/novoplasty/{seed}_kmer{kmer}.chloe.check"
    output:
        "results/{sample}/genbanks/novoplasty/{seed}_kmer{kmer}.chloe.genbank.check"
    log:
        "logs/{sample}/genbanks/novoplasty/{sample}_{kmer}_{seed}_write_gb_novoplasty_chloe.log"
    params:
        genetic_code=lambda wildcards: config["samples"][wildcards.sample]["genetic_code"]
    shell:
        """
        python workflow/scripts/chloro_gff2genbank.py --code {params.genetic_code} --assembler 'novoplasty' \
            --sample {wildcards.sample} --seed {wildcards.seed} \
            --kmer {wildcards.kmer} --software chloe --gene2product resources/gene2product.txt >> {log} 2>&1 && \
        touch {output}
        """

rule write_gb_novoplasty_cpgavas2:
    input:
        "results/{sample}/cpgavas2/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.cpgavas2.check"
    output:
        "results/{sample}/genbanks/novoplasty/{seed}_kmer{kmer}.cpgavas2.genbank.check"
    log:
        "logs/{sample}/genbanks/novoplasty/{sample}_{kmer}_{seed}_write_gb_novoplasty_cpgavas2.log"
    params:
        genetic_code=lambda wildcards: config["samples"][wildcards.sample]["genetic_code"]
    shell:
        """
        python workflow/scripts/chloro_gff2genbank.py --code {params.genetic_code} --assembler 'novoplasty' \
            --sample {wildcards.sample} --seed {wildcards.seed} \
            --kmer {wildcards.kmer} --software cpgavas2 >> {log} 2>&1 && \
        touch {output}
        """

rule rotate_mito_gb_novoplasty:
    input:
        "results/{sample}/genbanks/novoplasty/{seed}_kmer{kmer}.mitos2.genbank.check"
    output:
        "results/{sample}/genbanks/novoplasty/{seed}_kmer{kmer}.mitos2.genbank.rotated.check"
    log:
        "logs/{sample}/genbanks/novoplasty/{seed}_kmer{kmer}.genbank_rotate.log"
    params:
        organelle=lambda wildcards: config["samples"][wildcards.sample]["organelle"]
    shell:
        """
        python workflow/scripts/rotate_genbank.py --assembler 'novoplasty' --organelle {params.organelle} \
            --start_gene rrnS --seed {wildcards.seed} --kmer {wildcards.kmer} \
            --sample {wildcards.sample} >> {log} 2>&1 && \
        touch {output}
        """

rule rotate_cpgavas_gb_novoplasty:
    input:
        "results/{sample}/genbanks/novoplasty/{seed}_kmer{kmer}.cpgavas2.genbank.check"
    output:
        "results/{sample}/genbanks/novoplasty/{seed}_kmer{kmer}.cpgavas2.genbank.rotated.check"
    log:
        "logs/{sample}/genbanks/novoplasty/{seed}_kmer{kmer}.genbank.rotated.log"
    params:
        organelle=lambda wildcards: config["samples"][wildcards.sample]["organelle"]
    shell:
        """
        python workflow/scripts/rotate_genbank.py --assembler 'novoplasty' --organelle {params.organelle} \
            --start_gene psbA --seed {wildcards.seed} --kmer {wildcards.kmer} \
            --sample {wildcards.sample} --software cpgavas2 >> {log} 2>&1 && \
        touch {output}
        """

rule rotate_chloe_gb_novoplasty:
    input:
        "results/{sample}/genbanks/novoplasty/{seed}_kmer{kmer}.chloe.genbank.check"
    output:
        "results/{sample}/genbanks/novoplasty/{seed}_kmer{kmer}.chloe.genbank.rotated.check"
    log:
        "logs/{sample}/genbanks/novoplasty/{seed}_kmer{kmer}.genbank.rotated.log"
    params:
        organelle=lambda wildcards: config["samples"][wildcards.sample]["organelle"]
    shell:
        """
        python workflow/scripts/rotate_genbank.py --assembler 'novoplasty' --organelle {params.organelle} \
            --start_gene psbA --seed {wildcards.seed} --kmer {wildcards.kmer} \
            --sample {wildcards.sample} --software chloe >> {log} 2>&1 && \
        touch {output}
        """

rule write_gb_getorganelle_mitos2:
    input:
        "results/{sample}/mitos2/getorganelle/getorganelle_mitos2.check"
    output:
        "results/{sample}/genbanks/getorganelle/mitos2.genbank.check"
    log:
        "logs/{sample}/genbanks/getorganelle/write_gb_getorganelle_mitos2.log"
    params:
        genetic_code=lambda wildcards: config["samples"][wildcards.sample]["genetic_code"]
    shell:
        """
        python workflow/scripts/mito_gff2genbank.py --code {params.genetic_code} --assembler 'getorganelle' \
            --sample {wildcards.sample} >> {log} 2>&1 && \
        touch {output}
        """

rule write_gb_getorganelle_chloe:
    input:
        "results/{sample}/chloe/getorganelle/getorganelle_chloe.check"
    output:
        "results/{sample}/genbanks/getorganelle/chloe.genbank.check"
    log:
        "logs/{sample}/genbanks/getorganelle/write_gb_getorganelle_chloe.log"
    params:
        genetic_code=lambda wildcards: config["samples"][wildcards.sample]["genetic_code"]
    shell:
        """
        python workflow/scripts/chloro_gff2genbank.py --code {params.genetic_code} --assembler 'getorganelle' \
            --sample {wildcards.sample} --software chloe --gene2product resources/gene2product.txt >> {log} 2>&1 && \
        touch {output}
        """

rule write_gb_getorganelle_cpgavas2:
    input:
        "results/{sample}/cpgavas2/getorganelle/getorganelle_cpgavas2.check"
    output:
        "results/{sample}/genbanks/getorganelle/cpgavas2.genbank.check"
    log:
        "logs/{sample}/genbanks/getorganelle/write_gb_getorganelle_cpgavas2.log"
    params:
        genetic_code=lambda wildcards: config["samples"][wildcards.sample]["genetic_code"]
    shell:
        """
        python workflow/scripts/chloro_gff2genbank.py --code {params.genetic_code} --assembler 'getorganelle' \
            --sample {wildcards.sample} --software cpgavas2 >> {log} 2>&1 && \
        touch {output}
        """

rule rotate_mito_gb_getorganelle:
    input:
        "results/{sample}/genbanks/getorganelle/mitos2.genbank.check"
    output:
        "results/{sample}/genbanks/getorganelle/mitos2.genbank.rotated.check"
    log:
        "logs/{sample}/genbanks/getorganelle/genbank_rotate.log"
    params:
        organelle=lambda wildcards: config["samples"][wildcards.sample]["organelle"]
    shell:
        """
        python workflow/scripts/rotate_genbank.py --assembler 'getorganelle' --organelle {params.organelle} \
            --start_gene rrnS --sample {wildcards.sample} >> {log} 2>&1 && \
        touch {output}
        """

rule rotate_cpgavas_gb_getorganelle:
    input:
        "results/{sample}/genbanks/getorganelle/cpgavas2.genbank.check"
    output:
        "results/{sample}/genbanks/getorganelle/cpgavas2.genbank.rotated.check"
    log:
        "logs/{sample}/genbanks/getorganelle/genbank_rotate.log"
    params:
        organelle=lambda wildcards: config["samples"][wildcards.sample]["organelle"]
    shell:
        """
        python workflow/scripts/rotate_genbank.py --assembler 'getorganelle' --organelle {params.organelle} \
            --start_gene psbA --sample {wildcards.sample} --software cpgavas2 >> {log} 2>&1 && \
        touch {output}
        """

rule rotate_chloe_gb_getorganelle:
    input:
        "results/{sample}/genbanks/getorganelle/chloe.genbank.check"
    output:
        "results/{sample}/genbanks/getorganelle/chloe.genbank.rotated.check"
    log:
        "logs/{sample}/genbanks/getorganelle/genbank_rotate.log"
    params:
        organelle=lambda wildcards: config["samples"][wildcards.sample]["organelle"]
    shell:
        """
        python workflow/scripts/rotate_genbank.py --assembler 'getorganelle' --organelle {params.organelle} \
            --start_gene psbA --sample {wildcards.sample} --software chloe >> {log} 2>&1 && \
        touch {output}
        """
