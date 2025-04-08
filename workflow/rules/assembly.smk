rule create_hash:
    input:
        r1=lambda wildcards:
            "resources/{sample}/rawreads/{sample}.pair1.truncated.gz" if config["samples"][wildcards.sample].get("adapterremoval", "").lower() == "yes"
            else "resources/{sample}/rawreads/{sample}_R1.fastq.gz",
        r2=lambda wildcards:
            "resources/{sample}/rawreads/{sample}.pair2.truncated.gz" if config["samples"][wildcards.sample].get("adapterremoval", "").lower() == "yes"
            else "resources/{sample}/rawreads/{sample}_R2.fastq.gz"
    output:
        "results/{sample}/hashtable/kmer{kmer}/hash_config.txt",
        "results/{sample}/hashtable/kmer{kmer}/HASH2B_{sample}.txt",
        "results/{sample}/hashtable/kmer{kmer}/HASH2C_{sample}.txt",
        "results/{sample}/hashtable/kmer{kmer}/HASH_{sample}.txt",
    threads: 1
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
        max_memory=lambda wildcards: config["samples"][wildcards.sample]["max_memory"]
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
        seed = "resources/{sample}/seeds/{seed}.fasta"
    output:
        "results/{sample}/novoplasty/{seed}/kmer{kmer}/log_{sample}.txt"
    threads: 1
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
        reference=lambda wildcards: config["samples"][wildcards.sample]["reference"]
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
        reads= "resources/{sample}/rawreads/{sample}.fasta",
        reference_fasta="resources/{sample}/seeds/{seed}.fasta",
        reference_gb="resources/{sample}/seeds/{seed}.gb",
    output:
        "results/{sample}/mitohifi/{seed}/contigs_stats.tsv"
    threads: 8
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

rule get_assemblies:
    input:
        "results/{sample}/novoplasty/{seed}/kmer{kmer}/log_{sample}.txt"
    output:
        "results/{sample}/assemblies/{seed}_kmer{kmer}.fasta"
    threads: 1
    log:
        "logs/{sample}/novoplasty/{sample}_{kmer}_{seed}_get_assemblies.log"
    params:
        organelle=lambda wildcards: config["samples"][wildcards.sample]["organelle"],
        genome_range=lambda wildcards: config["samples"][wildcards.sample]["genome_range"],
    shell:
        """
        python workflow/scripts/get_assemblies.py \
            --organelle {params.organelle} --genome_range {params.genome_range} \
            --sample {wildcards.sample} --kmer {wildcards.kmer} --seed {wildcards.seed} >> {log} 2>&1
        """
