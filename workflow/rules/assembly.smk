rule create_hash:
    input:
        r1=lambda wildcards:
            "resources/{sample}/rawreads/{sample}.R1.trimmed.gz" if config["samples"][wildcards.sample].get("run_trimming", "").lower() == "yes"
            else "resources/{sample}/rawreads/{sample}_R1.fastq.gz",
        r2=lambda wildcards:
            "resources/{sample}/rawreads/{sample}.R2.trimmed.gz" if config["samples"][wildcards.sample].get("run_trimming", "").lower() == "yes"
            else "resources/{sample}/rawreads/{sample}_R2.fastq.gz"
    output:
        "results/{sample}/hashtable/kmer{kmer}/hash_config.txt",
        temp("results/{sample}/hashtable/kmer{kmer}/HASH2B_{sample}.txt"),
        temp("results/{sample}/hashtable/kmer{kmer}/HASH2C_{sample}.txt"),
        temp("results/{sample}/hashtable/kmer{kmer}/HASH_{sample}.txt"),
    log:
        "logs/{sample}/novoplasty/{sample}_{kmer}_create_hash.log"
    benchmark:
        "benchmarks/{sample}/novoplasty/{sample}_{kmer}_create_hash.benchmark"
    singularity:
        f"{config["sif_dir"]}/novoplasty.sif"
    params:
        organelle=lambda wildcards: config["samples"][wildcards.sample]["organelle"],
        genome_range=lambda wildcards: config["samples"][wildcards.sample]["genome_range"],
        reads_length=lambda wildcards: config["samples"][wildcards.sample]["reads_length"],
        insert_size=lambda wildcards: config["samples"][wildcards.sample]["insert_size"],
        max_memory=lambda wildcards: config["samples"][wildcards.sample]["max_memory"],
    shell:
        """
        mkdir -p results/{wildcards.sample}/hashtable/kmer{wildcards.kmer}/ && \
        python workflow/scripts/initial_config.py --organelle {params.organelle} \
        --sample {wildcards.sample} --genome_range {params.genome_range} --reads_length {params.reads_length} \
        --insert_size {params.insert_size} --forward {input.r1} \
        --reverse {input.r2} --max_memory {params.max_memory} --kmer {wildcards.kmer} && \
        NOVOPlasty.pl -c results/{wildcards.sample}/hashtable/kmer{wildcards.kmer}/hash_config.txt >> {log} 2>&1
        """

rule run_novoplasty:
    input:
        hash2b = "results/{sample}/hashtable/kmer{kmer}/HASH2B_{sample}.txt",
        hash2c = "results/{sample}/hashtable/kmer{kmer}/HASH2C_{sample}.txt",
        hashtable = "results/{sample}/hashtable/kmer{kmer}/HASH_{sample}.txt",
        seed = "resources/{sample}/seeds/{seed}.fasta",
        reference = "resources/{sample}/reference.fasta"
    output:
        "results/{sample}/novoplasty/{seed}/kmer{kmer}/log_{sample}.txt"
    singularity:
        f"{config["sif_dir"]}/novoplasty.sif"
    log:
        "logs/{sample}/novoplasty/{sample}_{kmer}_{seed}_novoplasty.log"
    benchmark:
        "benchmarks/{sample}/novoplasty/{sample}_{kmer}_{seed}_novoplasty.benchmark"
    params:
        organelle=lambda wildcards: config["samples"][wildcards.sample]["organelle"],
        genome_range=lambda wildcards: config["samples"][wildcards.sample]["genome_range"],
        reads_length=lambda wildcards: config["samples"][wildcards.sample]["reads_length"],
        insert_size=lambda wildcards: config["samples"][wildcards.sample]["insert_size"],
        max_memory=lambda wildcards: config["samples"][wildcards.sample]["max_memory"],
        reference=lambda wildcards: f"resources/{wildcards.sample}/reference.fasta" if config["samples"][wildcards.sample].get("reference", "") else ""
    shell:
        """
        python workflow/scripts/create_novoplasty_config.py \
            --organelle {params.organelle} \
            --sample {wildcards.sample} \
            --genome_range {params.genome_range} \
            --reads_length {params.reads_length} \
            --insert_size {params.insert_size} \
            --max_memory {params.max_memory} \
            --kmer {wildcards.kmer} \
            --seed {wildcards.seed} \
            --hashtable {input.hashtable} \
            --hash2b {input.hash2b} \
            --hash2c {input.hash2c} \
            --reference {params.reference} && \
        NOVOPlasty.pl -c results/{wildcards.sample}/novoplasty/{wildcards.seed}/kmer{wildcards.kmer}/config.txt >> {log} 2>&1
        """

rule run_mitohifi:
    input:
        reads = lambda wildcards:
            "resources/{sample}/rawreads/{sample}.trimmed.fasta" if config["samples"][wildcards.sample].get("run_trimming", "").lower() == "yes"
            else "resources/{sample}/rawreads/{sample}.fasta",
        reference_fasta="resources/{sample}/seeds/{seed}.fasta",
        reference_gb="resources/{sample}/seeds/{seed}.gb",
    output:
        "results/{sample}/mitohifi/{seed}/contigs_stats.tsv"
    log:
        "logs/{sample}/mitohifi/{seed}_run_mitohifi.log"
    benchmark:
       "benchmarks/{sample}/mitohifi/{seed}_run_mitohifi.txt"
    params:
        genetic_code=lambda wildcards: config["samples"][wildcards.sample]["genetic_code"]
    singularity:
        f"{config["sif_dir"]}/mitohifi.sif"
    shell:
        """
        mkdir -p results/{wildcards.sample}/mitohifi/{wildcards.seed} && \
        cd results/{wildcards.sample}/mitohifi/{wildcards.seed} && \
        mitohifi.py -t {threads} -r ../../../../{input.reads} \
            -f ../../../../{input.reference_fasta} \
            -g ../../../../{input.reference_gb} \
            -o {params.genetic_code} >> ../../../../{log} 2>&1
        """

rule run_getorganelle:
    input:
        r1 = lambda wildcards:
            "resources/{sample}/rawreads/{sample}.R1.trimmed.gz" if config["samples"][wildcards.sample].get("run_trimming", "").lower() == "yes"
            else "resources/{sample}/rawreads/{sample}_R1.fastq.gz",
        r2 = lambda wildcards:
            "resources/{sample}/rawreads/{sample}.R2.trimmed.gz" if config["samples"][wildcards.sample].get("run_trimming", "").lower() == "yes"
            else "resources/{sample}/rawreads/{sample}_R2.fastq.gz",
        database = "resources/getorganelle_db/getorganelle_db.check"
    output:
        "results/{sample}/getorganelle/get_org.log.txt"
    log:
        "logs/{sample}/getorganelle/run_getorganelle.log"
    benchmark:
       "benchmarks/{sample}/getorganelle/run_getorganelle.txt"
    params:
        database=lambda wildcards: config["samples"][wildcards.sample]["database"],
        n_rounds=lambda wildcards: config["samples"][wildcards.sample]["n_rounds"],
        target_size=lambda wildcards: config["samples"][wildcards.sample]["target_size"],
        extra_flags=lambda wildcards: config["samples"][wildcards.sample]["extra_flags"],
        spades_kmers=lambda wildcards: config["samples"][wildcards.sample]["spades_kmers"]
    singularity:
        f"{config["sif_dir"]}/getorganelle.sif"
    shell:
        """
        mkdir -p results/{wildcards.sample}/getorganelle && \
        export PATH=/opt/conda/bin:$PATH && \
        get_organelle_from_reads.py -1 {input.r1} -2 {input.r2} \
            -o results/{wildcards.sample}/getorganelle -t {threads} \
            --config-dir resources/getorganelle_db/ --overwrite \
            -F {params.database} -R {params.n_rounds} \
            -k {params.spades_kmers} --target-genome-size {params.target_size} \
            {params.extra_flags} >> {log} 2>&1
        """

rule get_novoplasty_assemblies:
    input:
        "results/{sample}/novoplasty/{seed}/kmer{kmer}/log_{sample}.txt"
    output:
        "results/{sample}/assemblies/novoplasty/{seed}_kmer{kmer}.fasta"
    log:
        "logs/{sample}/novoplasty/{sample}_{kmer}_{seed}_get_assemblies.log"
    params:
        organelle=lambda wildcards: config["samples"][wildcards.sample]["organelle"],
        genome_range=lambda wildcards: config["samples"][wildcards.sample]["genome_range"],
    shell:
        """
        python workflow/scripts/get_novoplasty_assemblies.py \
            --organelle {params.organelle} --genome_range {params.genome_range} \
            --sample {wildcards.sample} --kmer {wildcards.kmer} --seed {wildcards.seed} >> {log} 2>&1
        """

rule get_getorganelle_assemblies:
    input:
        "results/{sample}/getorganelle/get_org.log.txt"
    output:
        "results/{sample}/getorganelle/sequences.fasta"
    log:
        "logs/{sample}/getorganelle/get_assemblies.log"
    shell:
        """
        python workflow/scripts/get_getorganelle_assemblies.py \
            -i results/{wildcards.sample}/getorganelle \
            -o results/{wildcards.sample}/getorganelle/sequences.fasta >> {log} 2>&1
        """
