def novoplasty_mem_mb(wildcards):
    """NOVOPlasty uses up to max_memory GB; the extra 10% keeps the job within its memory limit."""
    return int(config["samples"][wildcards.sample]["max_memory"]) * 1100

def mitohifi_reads(wildcards):
    """MitoHiFi reads FASTA: the file as provided, or the reads converted by run_seqtk."""
    sample = wildcards.sample
    if get_longreads(sample).endswith(".fasta") and not is_yes(sample, "run_trimming"):
        return f"resources/{sample}/rawreads/{sample}.fasta"
    return f"resources/{sample}/rawreads/{sample}.converted.fasta"

def getorganelle_options(wildcards):
    """Optional GetOrganelle flags: a blank field keeps the default of the organelle type."""
    sample_config = config["samples"][wildcards.sample]
    flags = {"-R": "n_rounds", "-k": "spades_kmers", "--target-genome-size": "target_size"}
    return " ".join(f"{flag} {sample_config[field]}" for flag, field in flags.items() if sample_config[field] != "")

rule create_hash:
    input:
        unpack(short_reads)
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
        sif("novoplasty")
    resources:
        mem_mb=novoplasty_mem_mb
    params:
        organelle=sample_param("organelle"),
        genome_range=sample_param("genome_range"),
        reads_length=sample_param("reads_length"),
        insert_size=sample_param("insert_size"),
        max_memory=sample_param("max_memory"),
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
        sif("novoplasty")
    log:
        "logs/{sample}/novoplasty/{sample}_{kmer}_{seed}_novoplasty.log"
    benchmark:
        "benchmarks/{sample}/novoplasty/{sample}_{kmer}_{seed}_novoplasty.benchmark"
    resources:
        mem_mb=novoplasty_mem_mb
    params:
        organelle=sample_param("organelle"),
        genome_range=sample_param("genome_range"),
        reads_length=sample_param("reads_length"),
        insert_size=sample_param("insert_size"),
        max_memory=sample_param("max_memory"),
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
        reads = mitohifi_reads,
        reference_fasta="resources/{sample}/seeds/{seed}.fasta",
        reference_gb="resources/{sample}/seeds/{seed}.gb",
    output:
        "results/{sample}/mitohifi/{seed}/contigs_stats.tsv"
    log:
        "logs/{sample}/mitohifi/{seed}_run_mitohifi.log"
    benchmark:
       "benchmarks/{sample}/mitohifi/{seed}_run_mitohifi.txt"
    params:
        genetic_code=sample_param("genetic_code")
    singularity:
        sif("mitohifi")
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
        unpack(short_reads),
        database = "resources/getorganelle_db/getorganelle_db.check"
    output:
        "results/{sample}/getorganelle/get_org.log.txt"
    log:
        "logs/{sample}/getorganelle/run_getorganelle.log"
    benchmark:
       "benchmarks/{sample}/getorganelle/run_getorganelle.txt"
    params:
        database=sample_param("database"),
        options=getorganelle_options,
        extra_flags=sample_param("extra_flags"),
    singularity:
        sif("getorganelle")
    shell:
        """
        mkdir -p results/{wildcards.sample}/getorganelle && \
        export PATH=/opt/conda/bin:$PATH && \
        get_organelle_from_reads.py -1 {input.r1} -2 {input.r2} \
            -o results/{wildcards.sample}/getorganelle -t {threads} \
            --config-dir resources/getorganelle_db/ --overwrite \
            -F {params.database} {params.options} \
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
        organelle=sample_param("organelle"),
        genome_range=sample_param("genome_range"),
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
        "results/{sample}/assemblies/getorganelle/all.fasta"
    log:
        "logs/{sample}/getorganelle/get_assemblies.log"
    shell:
        """
        python workflow/scripts/get_getorganelle_assemblies.py \
            -i results/{wildcards.sample}/getorganelle \
            -o {output} >> {log} 2>&1
        """

# One FASTA per assembly. Each assembly then has its own jobs, so a failed job only
# repeats that assembly. The number of assemblies is only known after this step, and
# that is why a dry run does not list the jobs after it until the assemblies exist.
checkpoint split_assemblies:
    input:
        "results/{sample}/assemblies/{assembler}/{run}.fasta"
    output:
        directory("results/{sample}/assemblies/{assembler}/{run}")
    log:
        "logs/{sample}/assemblies/{assembler}/{run}_split_assemblies.log"
    shell:
        """
        python workflow/scripts/split_fasta.py --fasta_file {input} --assembler {wildcards.assembler} \
            --output_dir {output} >> {log} 2>&1
        """
