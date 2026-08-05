rule get_mito_ncRNA_sequences_novoplasty:
    input:
        "results/{sample}/mitos2/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.mitos2.check"
    output:
        "results/{sample}/nhmmer/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.ncRNA_sequences.check"
    log:
        "logs/{sample}/nhmmer/novoplasty/{kmer}_{seed}_get_mito_ncRNA_sequences.log"
    benchmark:
        "benchmarks/{sample}/nhmmer/novoplasty/{kmer}_{seed}_get_mito_ncRNA_sequences.txt"
    params:
        organelle=lambda wildcards: config["samples"][wildcards.sample]["organelle"],
        sequencing_type=lambda wildcards: config["samples"][wildcards.sample]["sequencing_type"]
    shell:
        """
        python workflow/scripts/ncRNA_nhmmer.py --organelle {params.organelle} \
            --sequencing_type {params.sequencing_type} --assembler 'novoplasty' \
            --sample {wildcards.sample} --seed {wildcards.seed} --kmer {wildcards.kmer} >> {log} 2>&1 && \
        touch {output}
        """

rule run_ncRNA_nhmmer_novoplasty:
    input:
        "results/{sample}/nhmmer/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.ncRNA_sequences.check",
        "resources/nhmmer_db.hmm",
        "results/{sample}/assemblies/novoplasty/{seed}_kmer{kmer}.fasta",
        expand("resources/nhmmer_db.hmm.{ext}", ext=['h3f', 'h3i', 'h3m', 'h3p'])
    output:
        "results/{sample}/nhmmer/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.ncRNA_nhmmer.check"
    log:
        "logs/{sample}/nhmmer/novoplasty/{kmer}_{seed}_run_ncRNA_nhmmer.log"
    benchmark:
        "benchmarks/{sample}/nhmmer/novoplasty/{kmer}_{seed}_run_ncRNA_nhmmer.txt"
    singularity:
        f"{config["sif_dir"]}/hmmer.sif"
    shell:
        """
        for fasta_file in results/{wildcards.sample}/assemblies/novoplasty/{wildcards.seed}_kmer{wildcards.kmer}/*.fasta; do
            fasta_header=$(awk '/^>/ {{print; exit}}' "$fasta_file" | sed 's/^>//') && \
            if [ -e "results/{wildcards.sample}/nhmmer/novoplasty/{wildcards.seed}_kmer{wildcards.kmer}/$fasta_header/rRNA-tRNA.fasta" ]; then
                hmmscan --cpu {threads} --noali \
                -o results/{wildcards.sample}/nhmmer/novoplasty/{wildcards.seed}_kmer{wildcards.kmer}/$fasta_header/rRNA-tRNA.nhmmer.out \
                --tblout results/{wildcards.sample}/nhmmer/novoplasty/{wildcards.seed}_kmer{wildcards.kmer}/$fasta_header/rRNA-tRNA.tblout.out \
                --pfamtblout results/{wildcards.sample}/nhmmer/novoplasty/{wildcards.seed}_kmer{wildcards.kmer}/$fasta_header/rRNA-tRNA.dfamtblout.out \
                resources/nhmmer_db.hmm \
                results/{wildcards.sample}/nhmmer/novoplasty/{wildcards.seed}_kmer{wildcards.kmer}/$fasta_header/rRNA-tRNA.fasta 2> {log}
            else
                echo "Input fasta file does not exist." >> {log} 2>&1
            fi
        done
        touch {output}
        """

rule get_mito_intergenes_sequences_novoplasty:
    input:
        "results/{sample}/mitos2/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.mitos2.check"
    output:
        "results/{sample}/nhmmer/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.intergenes_sequences.check"
    log:
        "logs/{sample}/nhmmer/novoplasty/{kmer}_{seed}_get_mito_intergenes_sequences.log"
    benchmark:
        "benchmarks/{sample}/nhmmer/novoplasty/{kmer}_{seed}_get_mito_intergenes_sequences.txt"
    params:
        organelle=lambda wildcards: config["samples"][wildcards.sample]["organelle"],
        sequencing_type=lambda wildcards: config["samples"][wildcards.sample]["sequencing_type"]
    shell:
        """
        python workflow/scripts/intergenes_nhmmer.py --organelle {params.organelle} \
            --sequencing_type {params.sequencing_type} --assembler 'novoplasty' \
            --sample {wildcards.sample} --seed {wildcards.seed} --kmer {wildcards.kmer} >> {log} 2>&1 && \
        touch {output}
        """

rule run_intergenes_nhmmer_novoplasty:
    input:
        "results/{sample}/nhmmer/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.intergenes_sequences.check",
        "resources/nhmmer_db.hmm",
        "results/{sample}/assemblies/novoplasty/{seed}_kmer{kmer}.fasta",
        expand("resources/nhmmer_db.hmm.{ext}", ext=['h3f', 'h3i', 'h3m', 'h3p'])
    output:
        "results/{sample}/nhmmer/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.intergenes_nhmmer.check"
    log:
        "logs/{sample}/nhmmer/novoplasty/{kmer}_{seed}_run_intergenes_nhmmer.log"
    benchmark:
        "benchmarks/{sample}/nhmmer/novoplasty/{kmer}_{seed}_run_intergenes_nhmmer.txt"
    singularity:
        f"{config["sif_dir"]}/hmmer.sif"
    shell:
        """
        for fasta_file in results/{wildcards.sample}/assemblies/novoplasty/{wildcards.seed}_kmer{wildcards.kmer}/*.fasta; do
            fasta_header=$(awk '/^>/ {{print; exit}}' "$fasta_file" | sed 's/^>//') && \
            if [ -e "results/{wildcards.sample}/nhmmer/novoplasty/{wildcards.seed}_kmer{wildcards.kmer}/$fasta_header/intergenes_filter.fasta" ]; then
                hmmscan --cpu {threads} --noali \
                -o results/{wildcards.sample}/nhmmer/novoplasty/{wildcards.seed}_kmer{wildcards.kmer}/$fasta_header/intergenes_filter.nhmmer.out \
                --tblout results/{wildcards.sample}/nhmmer/novoplasty/{wildcards.seed}_kmer{wildcards.kmer}/$fasta_header/intergenes_filter.tblout.out \
                --pfamtblout results/{wildcards.sample}/nhmmer/novoplasty/{wildcards.seed}_kmer{wildcards.kmer}/$fasta_header/intergenes_filter.dfamtblout.out \
                resources/nhmmer_db.hmm \
                results/{wildcards.sample}/nhmmer/novoplasty/{wildcards.seed}_kmer{wildcards.kmer}/$fasta_header/intergenes_filter.fasta 2> {log}
            else
                echo "Input fasta file does not exist." >> {log} 2>&1
            fi
        done
        touch {output}
        """

rule get_mito_ncRNA_sequences_long:
    input:
        "results/{sample}/mitohifi/{seed}/contigs_stats.tsv"
    output:
        "results/{sample}/nhmmer/{seed}/mitohifi/ncRNA_sequences.check"
    log:
        "logs/{sample}/nhmmer/{seed}_get_mito_ncRNA_sequences.log"
    benchmark:
        "benchmarks/{sample}/nhmmer/{seed}_get_mito_ncRNA_sequences.txt"
    params:
        organelle=lambda wildcards: config["samples"][wildcards.sample]["organelle"],
        sequencing_type=lambda wildcards: config["samples"][wildcards.sample]["sequencing_type"]
    shell:
        """
        python workflow/scripts/ncRNA_nhmmer.py --organelle {params.organelle} \
            --sequencing_type {params.sequencing_type} \
            --sample {wildcards.sample} --seed {wildcards.seed} >> {log} 2>&1 && \
        touch {output}
        """

rule run_ncRNA_nhmmer_long:
    input:
        "results/{sample}/nhmmer/{seed}/mitohifi/ncRNA_sequences.check",
        "resources/nhmmer_db.hmm",
        expand("resources/nhmmer_db.hmm.{ext}", ext=['h3f', 'h3i', 'h3m', 'h3p'])
    output:
        "results/{sample}/nhmmer/{seed}/mitohifi/ncRNA_nhmmer.check"
    log:
        "logs/{sample}/nhmmer/{seed}_run_ncRNA_nhmmer.log"
    benchmark:
        "benchmarks/{sample}/nhmmer/{seed}_run_ncRNA_nhmmer.txt"
    singularity:
        f"{config["sif_dir"]}/hmmer.sif"
    shell:
        """
        if [ -e "results/{wildcards.sample}/nhmmer/novoplasty/{wildcards.seed}/rRNA-tRNA.fasta" ]; then
            hmmscan --cpu {threads} --noali \
            -o results/{wildcards.sample}/nhmmer/novoplasty/{wildcards.seed}/rRNA-tRNA.nhmmer.out \
            --tblout results/{wildcards.sample}/nhmmer/novoplasty/{wildcards.seed}/rRNA-tRNA.tblout.out \
            --pfamtblout results/{wildcards.sample}/nhmmer/novoplasty/{wildcards.seed}/rRNA-tRNA.dfamtblout.out \
            resources/nhmmer_db.hmm \
            results/{wildcards.sample}/nhmmer/novoplasty/{wildcards.seed}/rRNA-tRNA.fasta 2>{log}
        else
            echo "Input fasta file does not exist." >> {log} 2>&1
        fi
        touch {output}
        """

rule get_mito_intergenes_sequences_long:
    input:
        "results/{sample}/mitohifi/{seed}/contigs_stats.tsv"
    output:
        "results/{sample}/nhmmer/{seed}/mitohifi/intergenes_sequences.check"
    log:
        "logs/{sample}/nhmmer/{seed}_get_mito_intergenes_sequences.log"
    benchmark:
        "benchmarks/{sample}/nhmmer/{seed}_get_mito_intergenes_sequences.txt"
    params:
        organelle=lambda wildcards: config["samples"][wildcards.sample]["organelle"],
        sequencing_type=lambda wildcards: config["samples"][wildcards.sample]["sequencing_type"]
    shell:
        """
        python workflow/scripts/intergenes_nhmmer.py --organelle {params.organelle} \
            --sequencing_type {params.sequencing_type} --assembler 'novoplasty' \
            --sample {wildcards.sample} --seed {wildcards.seed} >> {log} 2>&1 && \
        touch {output}
        """

rule run_intergenes_nhmmer_long:
    input:
        "results/{sample}/nhmmer/{seed}/mitohifi/intergenes_sequences.check",
        "resources/nhmmer_db.hmm",
        expand("resources/nhmmer_db.hmm.{ext}", ext=['h3f', 'h3i', 'h3m', 'h3p'])
    output:
        "results/{sample}/nhmmer/{seed}/mitohifi/intergenes_nhmmer.check"
    log:
        "logs/{sample}/nhmmer/{seed}_run_intergenes_nhmmer.log"
    benchmark:
        "benchmarks/{sample}/nhmmer/{seed}_run_intergenes_nhmmer.txt"
    singularity:
        f"{config["sif_dir"]}/hmmer.sif"
    shell:
        """
        if [ -e "results/{wildcards.sample}/nhmmer/novoplasty/{wildcards.seed}/intergenes_filter.fasta" ]; then
            hmmscan --cpu {threads} --noali \
            -o results/{wildcards.sample}/nhmmer/novoplasty/{wildcards.seed}/intergenes_filter.nhmmer.out \
            --tblout results/{wildcards.sample}/nhmmer/novoplasty/{wildcards.seed}/intergenes_filter.tblout.out \
            --pfamtblout results/{wildcards.sample}/nhmmer/novoplasty/{wildcards.seed}/intergenes_filter.dfamtblout.out \
            resources/nhmmer_db.hmm \
            results/{wildcards.sample}/nhmmer/novoplasty/{wildcards.seed}/intergenes_filter.fasta 2> {log}
        else
            echo "Input fasta file does not exist." >> {log} 2>&1
        fi
        touch {output}
        """

rule get_chloro_intergenes_sequences_novoplasty:
    input:
        "results/{sample}/genbanks/novoplasty/{seed}_kmer{kmer}.cpgavas2.genbank.check"
    output:
        "results/{sample}/nhmmer/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.chloro.intergenes_sequences.check"
    log:
        "logs/{sample}/nhmmer/{seed}_kmer{kmer}_chloro_intergenes_sequences.log"
    benchmark:
        "benchmarks/{sample}/nhmmer/{seed}_kmer{kmer}_chloro_intergenes_sequences.txt"
    params:
        organelle=lambda wildcards: config["samples"][wildcards.sample]["organelle"],
        sequencing_type=lambda wildcards: config["samples"][wildcards.sample]["sequencing_type"]
    shell:
        """
        python workflow/scripts/intergenes_nhmmer.py --organelle {params.organelle} \
            --sequencing_type {params.sequencing_type} --sample {wildcards.sample} --seed {wildcards.seed} \
            --kmer {wildcards.kmer} >> {log} 2>&1 && \
        touch {output}
        """

rule run_intergenes_nhmmer_chloro_novoplasty:
    input:
        "results/{sample}/nhmmer/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.chloro.intergenes_sequences.check",
        "resources/nhmmer_db.hmm",
        expand("resources/nhmmer_db.hmm.{ext}", ext=['h3f', 'h3i', 'h3m', 'h3p'])
    output:
        "results/{sample}/nhmmer/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.chloro.intergenes_nhmmer.check",
    log:
        "logs/{sample}/nhmmer/{seed}_kmer{kmer}_run_intergenes_nhmmer.log"
    benchmark:
        "benchmarks/{sample}/nhmmer/{seed}_kmer{kmer}_run_intergenes_nhmmer.txt"
    singularity:
        f"{config["sif_dir"]}/hmmer.sif"
    shell:
        """
        for fasta_file in results/{wildcards.sample}/assemblies/novoplasty/{wildcards.seed}_kmer{wildcards.kmer}/*.fasta; do
            fasta_header=$(awk '/^>/ {{print; exit}}' "$fasta_file" | sed 's/^>//') && \
            if [ -e "results/{wildcards.sample}/nhmmer/novoplasty/{wildcards.seed}_kmer{wildcards.kmer}/$fasta_header/intergenes_filter.fasta" ]; then
                hmmscan --cpu {threads} --noali \
                -o results/{wildcards.sample}/nhmmer/novoplasty/{wildcards.seed}_kmer{wildcards.kmer}/$fasta_header/intergenes_filter.nhmmer.out \
                --tblout results/{wildcards.sample}/nhmmer/novoplasty/{wildcards.seed}_kmer{wildcards.kmer}/$fasta_header/intergenes_filter.tblout.out \
                --pfamtblout results/{wildcards.sample}/nhmmer/novoplasty/{wildcards.seed}_kmer{wildcards.kmer}/$fasta_header/intergenes_filter.dfamtblout.out \
                resources/nhmmer_db.hmm \
                results/{wildcards.sample}/nhmmer/novoplasty/{wildcards.seed}_kmer{wildcards.kmer}/$fasta_header/intergenes_filter.fasta 2> {log}
            else
                echo "Input fasta file does not exist." >> {log} 2>&1
            fi
        done
        touch {output}
        """

rule get_chloro_ncRNA_sequences_novoplasty:
    input:
        "results/{sample}/genbanks/novoplasty/{seed}_kmer{kmer}.cpgavas2.genbank.check"
    output:
        "results/{sample}/nhmmer/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.chloro.ncRNA_sequences.check"
    log:
        "logs/{sample}/nhmmer/{seed}_kmer{kmer}_chloro_ncRNA_sequences.log"
    benchmark:
        "benchmarks/{sample}/nhmmer/{seed}_kmer{kmer}_chloro_ncRNA_sequences.txt"
    params:
        organelle=lambda wildcards: config["samples"][wildcards.sample]["organelle"],
        sequencing_type=lambda wildcards: config["samples"][wildcards.sample]["sequencing_type"]
    shell:
        """
        python workflow/scripts/ncRNA_nhmmer.py --organelle {params.organelle} \
            --sequencing_type {params.sequencing_type} --sample {wildcards.sample} --seed {wildcards.seed} \
            --kmer {wildcards.kmer} >> {log} 2>&1 && \
        touch {output}
        """

rule run_ncRNA_nhmmer_chloro_novoplasty:
    input:
        "results/{sample}/nhmmer/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.chloro.ncRNA_sequences.check",
        "resources/nhmmer_db.hmm",
        expand("resources/nhmmer_db.hmm.{ext}", ext=['h3f', 'h3i', 'h3m', 'h3p'])
    output:
        "results/{sample}/nhmmer/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.chloro.ncRNA_nhmmer.check",
    log:
        "logs/{sample}/nhmmer/{seed}_kmer{kmer}_run_ncRNA_nhmmer.log"
    benchmark:
        "benchmarks/{sample}/nhmmer/{seed}_kmer{kmer}_run_ncRNA_nhmmer.txt"
    singularity:
        f"{config["sif_dir"]}/hmmer.sif"
    shell:
        """
        for fasta_file in results/{wildcards.sample}/assemblies/novoplasty/{wildcards.seed}_kmer{wildcards.kmer}/*.fasta; do
            fasta_header=$(awk '/^>/ {{print; exit}}' "$fasta_file" | sed 's/^>//') && \
            if [ -e "results/{wildcards.sample}/nhmmer/novoplasty/{wildcards.seed}_kmer{wildcards.kmer}/$fasta_header/rRNA-tRNA.fasta" ]; then
                hmmscan --cpu {threads} --noali \
                -o results/{wildcards.sample}/nhmmer/novoplasty/{wildcards.seed}_kmer{wildcards.kmer}/$fasta_header/rRNA-tRNA.nhmmer.out \
                --tblout results/{wildcards.sample}/nhmmer/novoplasty/{wildcards.seed}_kmer{wildcards.kmer}/$fasta_header/rRNA-tRNA.tblout.out \
                --pfamtblout results/{wildcards.sample}/nhmmer/novoplasty/{wildcards.seed}_kmer{wildcards.kmer}/$fasta_header/rRNA-tRNA.dfamtblout.out \
                resources/nhmmer_db.hmm \
                results/{wildcards.sample}/nhmmer/novoplasty/{wildcards.seed}_kmer{wildcards.kmer}/$fasta_header/rRNA-tRNA.fasta 2> {log}
            else
                echo "Input fasta file does not exist." >> {log} 2>&1
            fi
        done
        touch {output}
        """


rule get_mito_ncRNA_sequences_getorganelle:
    input:
        "results/{sample}/mitos2/getorganelle/getorganelle_mitos2.check"
    output:
        "results/{sample}/nhmmer/getorganelle/ncRNA_sequences.check"
    log:
        "logs/{sample}/nhmmer/getorganelle/get_mito_ncRNA_sequences.log"
    benchmark:
        "benchmarks/{sample}/nhmmer/getorganelle/get_mito_ncRNA_sequences.txt"
    params:
        organelle=lambda wildcards: config["samples"][wildcards.sample]["organelle"],
        sequencing_type=lambda wildcards: config["samples"][wildcards.sample]["sequencing_type"]
    shell:
        """
        python workflow/scripts/ncRNA_nhmmer.py --organelle {params.organelle} --sample {wildcards.sample} \
            --sequencing_type {params.sequencing_type} --assembler 'getorganelle' >> {log} 2>&1 && \
        touch {output}
        """

rule get_chloro_ncRNA_sequences_novoplasty_getorganelle:
    input:
        "results/{sample}/cpgavas2/getorganelle/getorganelle_cpgavas2.check"
    output:
        "results/{sample}/nhmmer/getorganelle/chloro.ncRNA_sequences.check"
    log:
        "logs/{sample}/nhmmer/getorganelle/chloro_ncRNA_sequences.log"
    benchmark:
        "benchmarks/{sample}/nhmmer/getorganelle/chloro_ncRNA_sequences.txt"
    params:
        organelle=lambda wildcards: config["samples"][wildcards.sample]["organelle"],
        sequencing_type=lambda wildcards: config["samples"][wildcards.sample]["sequencing_type"]
    shell:
        """
        python workflow/scripts/ncRNA_nhmmer.py --organelle {params.organelle} --sample {wildcards.sample} \
            --sequencing_type {params.sequencing_type} --assembler 'getorganelle' >> {log} 2>&1 && \
        touch {output}
        """

rule run_ncRNA_nhmmer_getorganelle:
    input:
        "results/{sample}/nhmmer/getorganelle/ncRNA_sequences.check",
        "resources/nhmmer_db.hmm",
        "results/{sample}/getorganelle/sequences.fasta",
        expand("resources/nhmmer_db.hmm.{ext}", ext=['h3f', 'h3i', 'h3m', 'h3p'])
    output:
        "results/{sample}/nhmmer/getorganelle/ncRNA_nhmmer.check"
    log:
        "logs/{sample}/nhmmer/getorganelle/run_ncRNA_nhmmer.log"
    benchmark:
        "benchmarks/{sample}/nhmmer/getorganelle/run_ncRNA_nhmmer.txt"
    singularity:
        f"{config["sif_dir"]}/hmmer.sif"
    shell:
        """
        for fasta_file in results/{wildcards.sample}/assemblies/getorganelle/*.fasta; do
            fasta_header=$(awk '/^>/ {{print; exit}}' "$fasta_file" | sed 's/^>//') && \
            if [ -e "results/{wildcards.sample}/nhmmer/getorganelle/$fasta_header/rRNA-tRNA.fasta" ]; then
                hmmscan --cpu {threads} --noali \
                -o results/{wildcards.sample}/nhmmer/getorganelle/$fasta_header/rRNA-tRNA.nhmmer.out \
                --tblout results/{wildcards.sample}/nhmmer/getorganelle/$fasta_header/rRNA-tRNA.tblout.out \
                --pfamtblout results/{wildcards.sample}/nhmmer/getorganelle/$fasta_header/rRNA-tRNA.dfamtblout.out \
                resources/nhmmer_db.hmm \
                results/{wildcards.sample}/nhmmer/getorganelle/$fasta_header/rRNA-tRNA.fasta 2> {log}
            else
                echo "Input fasta file does not exist." >> {log} 2>&1
            fi
        done
        touch {output}
        """

rule run_ncRNA_nhmmer_chloro_getorganelle:
    input:
        "results/{sample}/nhmmer/getorganelle/chloro.ncRNA_sequences.check",
        "resources/nhmmer_db.hmm",
        "results/{sample}/getorganelle/sequences.fasta",
        expand("resources/nhmmer_db.hmm.{ext}", ext=['h3f', 'h3i', 'h3m', 'h3p'])
    output:
        "results/{sample}/nhmmer/getorganelle/chloro.ncRNA_nhmmer.check",
    log:
        "logs/{sample}/nhmmer/getorganelle/run_ncRNA_nhmmer.log"
    benchmark:
        "benchmarks/{sample}/nhmmer/getorganelle/run_ncRNA_nhmmer.txt"
    singularity:
        f"{config["sif_dir"]}/hmmer.sif"
    shell:
        """
        for fasta_file in results/{wildcards.sample}/assemblies/getorganelle/*.fasta; do
            fasta_header=$(awk '/^>/ {{print; exit}}' "$fasta_file" | sed 's/^>//') && \
            if [ -e "results/{wildcards.sample}/nhmmer/getorganelle/$fasta_header/rRNA-tRNA.fasta" ]; then
                hmmscan --cpu {threads} --noali \
                -o results/{wildcards.sample}/nhmmer/getorganelle/$fasta_header/rRNA-tRNA.nhmmer.out \
                --tblout results/{wildcards.sample}/nhmmer/getorganelle/$fasta_header/rRNA-tRNA.tblout.out \
                --pfamtblout results/{wildcards.sample}/nhmmer/getorganelle/$fasta_header/rRNA-tRNA.dfamtblout.out \
                resources/nhmmer_db.hmm \
                results/{wildcards.sample}/nhmmer/getorganelle/$fasta_header/rRNA-tRNA.fasta 2> {log}
            else
                echo "Input fasta file does not exist." >> {log} 2>&1
            fi
        done
        touch {output}
        """

rule get_mito_intergenes_sequences_getorganelle:
    input:
        "results/{sample}/mitos2/getorganelle/getorganelle_mitos2.check"
    output:
        "results/{sample}/nhmmer/getorganelle/intergenes_sequences.check"
    log:
        "logs/{sample}/nhmmer/getorganelle/get_mito_intergenes_sequences.log"
    benchmark:
        "benchmarks/{sample}/nhmmer/getorganelle/get_mito_intergenes_sequences.txt"
    params:
        organelle=lambda wildcards: config["samples"][wildcards.sample]["organelle"],
        sequencing_type=lambda wildcards: config["samples"][wildcards.sample]["sequencing_type"]
    shell:
        """
        python workflow/scripts/intergenes_nhmmer.py --organelle {params.organelle} \
            --sequencing_type {params.sequencing_type} --sample {wildcards.sample}\
            --assembler 'getorganelle' >> {log} 2>&1 && \
        touch {output}
        """

rule get_chloro_intergenes_sequences_getorganelle:
    input:
        "results/{sample}/cpgavas2/getorganelle/getorganelle_cpgavas2.check"
    output:
        "results/{sample}/nhmmer/getorganelle/chloro.intergenes_sequences.check"
    log:
        "logs/{sample}/nhmmer/getorganelle/chloro_intergenes_sequences.log"
    benchmark:
        "benchmarks/{sample}/nhmmer/getorganelle/chloro_intergenes_sequences.txt"
    params:
        organelle=lambda wildcards: config["samples"][wildcards.sample]["organelle"],
        sequencing_type=lambda wildcards: config["samples"][wildcards.sample]["sequencing_type"]
    shell:
        """
        python workflow/scripts/intergenes_nhmmer.py --organelle {params.organelle} \
            --sequencing_type {params.sequencing_type} --sample {wildcards.sample} \
            --assembler 'getorganelle' >> {log} 2>&1 && \
        touch {output}
        """

rule run_intergenes_nhmmer_getorganelle:
    input:
        "results/{sample}/nhmmer/getorganelle/intergenes_sequences.check",
        "resources/nhmmer_db.hmm",
        "results/{sample}/getorganelle/sequences.fasta",
        expand("resources/nhmmer_db.hmm.{ext}", ext=['h3f', 'h3i', 'h3m', 'h3p'])
    output:
        "results/{sample}/nhmmer/getorganelle/intergenes_nhmmer.check"
    log:
        "logs/{sample}/nhmmer/getorganelle/run_intergenes_nhmmer.log"
    benchmark:
        "benchmarks/{sample}/nhmmer/getorganelle/run_intergenes_nhmmer.txt"
    singularity:
        f"{config["sif_dir"]}/hmmer.sif"
    shell:
        """
        for fasta_file in results/{wildcards.sample}/assemblies/getorganelle/*.fasta; do
            fasta_header=$(awk '/^>/ {{print; exit}}' "$fasta_file" | sed 's/^>//') && \
            if [ -e "results/{wildcards.sample}/nhmmer/getorganelle/$fasta_header/intergenes_filter.fasta" ]; then
                hmmscan --cpu {threads} --noali \
                -o results/{wildcards.sample}/nhmmer/getorganelle/$fasta_header/intergenes_filter.nhmmer.out \
                --tblout results/{wildcards.sample}/nhmmer/getorganelle/$fasta_header/intergenes_filter.tblout.out \
                --pfamtblout results/{wildcards.sample}/nhmmer/getorganelle/$fasta_header/intergenes_filter.dfamtblout.out \
                resources/nhmmer_db.hmm \
                results/{wildcards.sample}/nhmmer/getorganelle/$fasta_header/intergenes_filter.fasta 2> {log}
            else
                echo "Input fasta file does not exist." >> {log} 2>&1
            fi
        done
        touch {output}
        """

rule run_intergenes_nhmmer_chloro_getorganelle:
    input:
        "results/{sample}/nhmmer/getorganelle/chloro.intergenes_sequences.check",
        "resources/nhmmer_db.hmm",
        expand("resources/nhmmer_db.hmm.{ext}", ext=['h3f', 'h3i', 'h3m', 'h3p'])
    output:
        "results/{sample}/nhmmer/getorganelle/chloro.intergenes_nhmmer.check",
    log:
        "logs/{sample}/nhmmer/getorganelle/run_intergenes_nhmmer.log"
    benchmark:
        "benchmarks/{sample}/nhmmer/getorganelle/run_intergenes_nhmmer.txt"
    singularity:
        f"{config["sif_dir"]}/hmmer.sif"
    shell:
        """
        for fasta_file in results/{wildcards.sample}/assemblies/getorganelle/*.fasta; do
            fasta_header=$(awk '/^>/ {{print; exit}}' "$fasta_file" | sed 's/^>//') && \
            if [ -e "results/{wildcards.sample}/nhmmer/getorganelle/$fasta_header/intergenes_filter.fasta" ]; then
                hmmscan --cpu {threads} --noali \
                -o results/{wildcards.sample}/nhmmer/getorganelle/$fasta_header/intergenes_filter.nhmmer.out \
                --tblout results/{wildcards.sample}/nhmmer/getorganelle/$fasta_header/intergenes_filter.tblout.out \
                --pfamtblout results/{wildcards.sample}/nhmmer/getorganelle/$fasta_header/intergenes_filter.dfamtblout.out \
                resources/nhmmer_db.hmm \
                results/{wildcards.sample}/nhmmer/getorganelle/$fasta_header/intergenes_filter.fasta 2> {log}
            else
                echo "Input fasta file does not exist." >> {log} 2>&1
            fi
        done
        touch {output}
        """
