rule get_mito_ncRNA_sequences_short:
    input:
        "results/{sample}/mitos2/{seed}_kmer{kmer}/{seed}_kmer{kmer}.mitos2.check"
    output:
        "results/{sample}/nhmmer/{seed}_kmer{kmer}/{seed}_kmer{kmer}.ncRNA_sequences.check"
    log:
        "logs/{sample}/nhmmer/{kmer}_{seed}_get_mito_ncRNA_sequences.log"
    threads: 1
    benchmark:
        "benchmarks/{sample}/nhmmer/{kmer}_{seed}_get_mito_ncRNA_sequences.txt"
    params:
        organelle=lambda wildcards: config["samples"][wildcards.sample]["organelle"],
        sequencing_type=lambda wildcards: config["samples"][wildcards.sample]["sequencing_type"]
    shell:
        """
        python workflow/scripts/ncRNA_nhmmer.py --organelle {params.organelle} \
            --sequencing_type {params.sequencing_type} \
            --sample {wildcards.sample} --seed {wildcards.seed} --kmer {wildcards.kmer} >> {log} 2>&1 && \
        touch {output}
        """

rule run_ncRNA_nhmmer_short:
    input:
        "results/{sample}/nhmmer/{seed}_kmer{kmer}/{seed}_kmer{kmer}.ncRNA_sequences.check",
        "resources/nhmmer_db.hmm",
        "results/{sample}/assemblies/{seed}_kmer{kmer}.fasta",
        expand("resources/nhmmer_db.hmm.{ext}", ext=['h3f', 'h3i', 'h3m', 'h3p'])
    output:
        "results/{sample}/nhmmer/{seed}_kmer{kmer}/{seed}_kmer{kmer}.ncRNA_nhmmer.check"
    threads: 4
    log:
        "logs/{sample}/nhmmer/{kmer}_{seed}_run_ncRNA_nhmmer.log"
    benchmark:
        "benchmarks/{sample}/nhmmer/{kmer}_{seed}_run_ncRNA_nhmmer.txt"
    singularity:
        "docker://staphb/hmmer:3.3"
    shell:
        """
        for fasta_file in results/{wildcards.sample}/assemblies/{wildcards.seed}_kmer{wildcards.kmer}/*.fasta; do
            fasta_header=$(awk '/^>/ {{print; exit}}' "$fasta_file" | sed 's/^>//') && \
            if [ -e "results/{wildcards.sample}/nhmmer/{wildcards.seed}_kmer{wildcards.kmer}/$fasta_header/rRNA-tRNA.fasta" ]; then
                hmmscan --cpu {threads} --noali \
                -o results/{wildcards.sample}/nhmmer/{wildcards.seed}_kmer{wildcards.kmer}/$fasta_header/rRNA-tRNA.nhmmer.out \
                --tblout results/{wildcards.sample}/nhmmer/{wildcards.seed}_kmer{wildcards.kmer}/$fasta_header/rRNA-tRNA.tblout.out \
                --pfamtblout results/{wildcards.sample}/nhmmer/{wildcards.seed}_kmer{wildcards.kmer}/$fasta_header/rRNA-tRNA.dfamtblout.out \
                resources/nhmmer_db.hmm \
                results/{wildcards.sample}/nhmmer/{wildcards.seed}_kmer{wildcards.kmer}/$fasta_header/rRNA-tRNA.fasta 2> {log}
            else
                echo "Input fasta file does not exist." >> {log} 2>&1
            fi
        done
        touch {output}
        """

rule get_mito_intergenes_sequences_short:
    input:
        "results/{sample}/mitos2/{seed}_kmer{kmer}/{seed}_kmer{kmer}.mitos2.check"
    output:
        "results/{sample}/nhmmer/{seed}_kmer{kmer}/{seed}_kmer{kmer}.intergenes_sequences.check"
    log:
        "logs/{sample}/nhmmer/{kmer}_{seed}_get_mito_intergenes_sequences.log"
    threads: 1
    benchmark:
        "benchmarks/{sample}/nhmmer/{kmer}_{seed}_get_mito_intergenes_sequences.txt"
    params:
        organelle=lambda wildcards: config["samples"][wildcards.sample]["organelle"],
        sequencing_type=lambda wildcards: config["samples"][wildcards.sample]["sequencing_type"]
    shell:
        """
        python workflow/scripts/intergenes_nhmmer.py --organelle {params.organelle} \
            --sequencing_type {params.sequencing_type} \
            --sample {wildcards.sample} --seed {wildcards.seed} --kmer {wildcards.kmer} >> {log} 2>&1 && \
        touch {output}
        """

rule run_intergenes_nhmmer_short:
    input:
        "results/{sample}/nhmmer/{seed}_kmer{kmer}/{seed}_kmer{kmer}.intergenes_sequences.check",
        "resources/nhmmer_db.hmm",
        "results/{sample}/assemblies/{seed}_kmer{kmer}.fasta",
        expand("resources/nhmmer_db.hmm.{ext}", ext=['h3f', 'h3i', 'h3m', 'h3p'])
    output:
        "results/{sample}/nhmmer/{seed}_kmer{kmer}/{seed}_kmer{kmer}.intergenes_nhmmer.check"
    threads: 4
    log:
        "logs/{sample}/nhmmer/{kmer}_{seed}_run_intergenes_nhmmer.log"
    benchmark:
        "benchmarks/{sample}/nhmmer/{kmer}_{seed}_run_intergenes_nhmmer.txt"
    singularity:
        "docker://staphb/hmmer:3.3"
    shell:
        """
        for fasta_file in results/{wildcards.sample}/assemblies/{wildcards.seed}_kmer{wildcards.kmer}/*.fasta; do
            fasta_header=$(awk '/^>/ {{print; exit}}' "$fasta_file" | sed 's/^>//') && \
            if [ -e "results/{wildcards.sample}/nhmmer/{wildcards.seed}_kmer{wildcards.kmer}/$fasta_header/intergenes_filter.fasta" ]; then
                hmmscan --cpu {threads} --noali \
                -o results/{wildcards.sample}/nhmmer/{wildcards.seed}_kmer{wildcards.kmer}/$fasta_header/intergenes_filter.nhmmer.out \
                --tblout results/{wildcards.sample}/nhmmer/{wildcards.seed}_kmer{wildcards.kmer}/$fasta_header/intergenes_filter.tblout.out \
                --pfamtblout results/{wildcards.sample}/nhmmer/{wildcards.seed}_kmer{wildcards.kmer}/$fasta_header/intergenes_filter.dfamtblout.out \
                resources/nhmmer_db.hmm \
                results/{wildcards.sample}/nhmmer/{wildcards.seed}_kmer{wildcards.kmer}/$fasta_header/intergenes_filter.fasta 2> {log}
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
        "results/{sample}/nhmmer/{seed}/ncRNA_sequences.check"
    log:
        "logs/{sample}/nhmmer/{seed}_get_mito_ncRNA_sequences.log"
    threads: 1
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
        "results/{sample}/nhmmer/{seed}/ncRNA_sequences.check",
        "resources/nhmmer_db.hmm",
        expand("resources/nhmmer_db.hmm.{ext}", ext=['h3f', 'h3i', 'h3m', 'h3p'])
    output:
        "results/{sample}/nhmmer/{seed}/ncRNA_nhmmer.check"
    threads: 4
    log:
        "logs/{sample}/nhmmer/{seed}_run_ncRNA_nhmmer.log"
    benchmark:
        "benchmarks/{sample}/nhmmer/{seed}_run_ncRNA_nhmmer.txt"
    singularity:
        "docker://staphb/hmmer:3.3"
    shell:
        """
        if [ -e "results/{wildcards.sample}/nhmmer/{wildcards.seed}/rRNA-tRNA.fasta" ]; then
            hmmscan --cpu {threads} --noali \
            -o results/{wildcards.sample}/nhmmer/{wildcards.seed}/rRNA-tRNA.nhmmer.out \
            --tblout results/{wildcards.sample}/nhmmer/{wildcards.seed}/rRNA-tRNA.tblout.out \
            --pfamtblout results/{wildcards.sample}/nhmmer/{wildcards.seed}/rRNA-tRNA.dfamtblout.out \
            resources/nhmmer_db.hmm \
            results/{wildcards.sample}/nhmmer/{wildcards.seed}/rRNA-tRNA.fasta 2>{log}
        else
            echo "Input fasta file does not exist." >> {log} 2>&1
        fi
        touch {output}
        """

rule get_mito_intergenes_sequences_long:
    input:
        "results/{sample}/mitohifi/{seed}/contigs_stats.tsv"
    output:
        "results/{sample}/nhmmer/{seed}/intergenes_sequences.check"
    log:
        "logs/{sample}/nhmmer/{seed}_get_mito_intergenes_sequences.log"
    threads: 1
    benchmark:
        "benchmarks/{sample}/nhmmer/{seed}_get_mito_intergenes_sequences.txt"
    params:
        organelle=lambda wildcards: config["samples"][wildcards.sample]["organelle"],
        sequencing_type=lambda wildcards: config["samples"][wildcards.sample]["sequencing_type"]
    shell:
        """
        python workflow/scripts/intergenes_nhmmer.py --organelle {params.organelle} \
            --sequencing_type {params.sequencing_type} \
            --sample {wildcards.sample} --seed {wildcards.seed} >> {log} 2>&1 && \
        touch {output}
        """

rule run_intergenes_nhmmer_long:
    input:
        "results/{sample}/nhmmer/{seed}/intergenes_sequences.check",
        "resources/nhmmer_db.hmm",
        expand("resources/nhmmer_db.hmm.{ext}", ext=['h3f', 'h3i', 'h3m', 'h3p'])
    output:
        "results/{sample}/nhmmer/{seed}/intergenes_nhmmer.check"
    threads: 4
    log:
        "logs/{sample}/nhmmer/{seed}_run_intergenes_nhmmer.log"
    benchmark:
        "benchmarks/{sample}/nhmmer/{seed}_run_intergenes_nhmmer.txt"
    singularity:
        "docker://staphb/hmmer:3.3"
    shell:
        """
        if [ -e "results/{wildcards.sample}/nhmmer/{wildcards.seed}/intergenes_filter.fasta" ]; then
            hmmscan --cpu {threads} --noali \
            -o results/{wildcards.sample}/nhmmer/{wildcards.seed}/intergenes_filter.nhmmer.out \
            --tblout results/{wildcards.sample}/nhmmer/{wildcards.seed}/intergenes_filter.tblout.out \
            --pfamtblout results/{wildcards.sample}/nhmmer/{wildcards.seed}/intergenes_filter.dfamtblout.out \
            resources/nhmmer_db.hmm \
            results/{wildcards.sample}/nhmmer/{wildcards.seed}/intergenes_filter.fasta 2> {log}
        else
            echo "Input fasta file does not exist." >> {log} 2>&1
        fi
        touch {output}
        """

rule get_chloro_intergenes_sequences:
    input:
        "results/{sample}/genbanks/{seed}_kmer{kmer}.cpgavas2.genbank.check"
    output:
        "results/{sample}/nhmmer/{seed}_kmer{kmer}/{seed}_kmer{kmer}.chloro.intergenes_sequences.check"
    log:
        "logs/{sample}/nhmmer/{seed}_kmer{kmer}_chloro_intergenes_sequences.log"
    threads: 1
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

rule run_intergenes_nhmmer_chloro:
    input:
        "results/{sample}/nhmmer/{seed}_kmer{kmer}/{seed}_kmer{kmer}.chloro.intergenes_sequences.check",
        "resources/nhmmer_db.hmm",
        expand("resources/nhmmer_db.hmm.{ext}", ext=['h3f', 'h3i', 'h3m', 'h3p'])
    output:
        "results/{sample}/nhmmer/{seed}_kmer{kmer}/{seed}_kmer{kmer}.chloro.intergenes_nhmmer.check",
    threads: 4
    log:
        "logs/{sample}/nhmmer/{seed}_kmer{kmer}_run_intergenes_nhmmer.log"
    benchmark:
        "benchmarks/{sample}/nhmmer/{seed}_kmer{kmer}_run_intergenes_nhmmer.txt"
    singularity:
        "docker://staphb/hmmer:3.4"
    shell:
        """
        for fasta_file in results/{wildcards.sample}/assemblies/{wildcards.seed}_kmer{wildcards.kmer}/*.fasta; do
            fasta_header=$(awk '/^>/ {{print; exit}}' "$fasta_file" | sed 's/^>//') && \
            if [ -e "results/{wildcards.sample}/nhmmer/{wildcards.seed}_kmer{wildcards.kmer}/$fasta_header/intergenes_filter.fasta" ]; then
                hmmscan --cpu {threads} --noali \
                -o results/{wildcards.sample}/nhmmer/{wildcards.seed}_kmer{wildcards.kmer}/$fasta_header/intergenes_filter.nhmmer.out \
                --tblout results/{wildcards.sample}/nhmmer/{wildcards.seed}_kmer{wildcards.kmer}/$fasta_header/intergenes_filter.tblout.out \
                --pfamtblout results/{wildcards.sample}/nhmmer/{wildcards.seed}_kmer{wildcards.kmer}/$fasta_header/intergenes_filter.dfamtblout.out \
                resources/nhmmer_db.hmm \
                results/{wildcards.sample}/nhmmer/{wildcards.seed}_kmer{wildcards.kmer}/$fasta_header/intergenes_filter.fasta 2> {log}
            else
                echo "Input fasta file does not exist." >> {log} 2>&1
            fi
        done
        touch {output}
        """

rule get_chloro_ncRNA_sequences:
    input:
        "results/{sample}/genbanks/{seed}_kmer{kmer}.cpgavas2.genbank.check"
    output:
        "results/{sample}/nhmmer/{seed}_kmer{kmer}/{seed}_kmer{kmer}.chloro.ncRNA_sequences.check"
    log:
        "logs/{sample}/nhmmer/{seed}_kmer{kmer}_chloro_ncRNA_sequences.log"
    threads: 1
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

rule run_ncRNA_nhmmer_chloro:
    input:
        "results/{sample}/nhmmer/{seed}_kmer{kmer}/{seed}_kmer{kmer}.chloro.ncRNA_sequences.check",
        "resources/nhmmer_db.hmm",
        expand("resources/nhmmer_db.hmm.{ext}", ext=['h3f', 'h3i', 'h3m', 'h3p'])
    output:
        "results/{sample}/nhmmer/{seed}_kmer{kmer}/{seed}_kmer{kmer}.chloro.ncRNA_nhmmer.check",
    threads: 4
    log:
        "logs/{sample}/nhmmer/{seed}_kmer{kmer}_run_ncRNA_nhmmer.log"
    benchmark:
        "benchmarks/{sample}/nhmmer/{seed}_kmer{kmer}_run_ncRNA_nhmmer.txt"
    singularity:
        "docker://staphb/hmmer:3.4"
    shell:
        """
        for fasta_file in results/{wildcards.sample}/assemblies/{wildcards.seed}_kmer{wildcards.kmer}/*.fasta; do
            fasta_header=$(awk '/^>/ {{print; exit}}' "$fasta_file" | sed 's/^>//') && \
            if [ -e "results/{wildcards.sample}/nhmmer/{wildcards.seed}_kmer{wildcards.kmer}/$fasta_header/rRNA-tRNA.fasta" ]; then
                hmmscan --cpu {threads} --noali \
                -o results/{wildcards.sample}/nhmmer/{wildcards.seed}_kmer{wildcards.kmer}/$fasta_header/rRNA-tRNA.nhmmer.out \
                --tblout results/{wildcards.sample}/nhmmer/{wildcards.seed}_kmer{wildcards.kmer}/$fasta_header/rRNA-tRNA.tblout.out \
                --pfamtblout results/{wildcards.sample}/nhmmer/{wildcards.seed}_kmer{wildcards.kmer}/$fasta_header/rRNA-tRNA.dfamtblout.out \
                resources/nhmmer_db.hmm \
                results/{wildcards.sample}/nhmmer/{wildcards.seed}_kmer{wildcards.kmer}/$fasta_header/rRNA-tRNA.fasta 2> {log}
            else
                echo "Input fasta file does not exist." >> {log} 2>&1
            fi
        done
        touch {output}
        """
