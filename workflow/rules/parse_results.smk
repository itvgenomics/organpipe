def assembly_targets(sample, assembler, run, assembly):
    """Last outputs of the steps enabled for one assembly. The earlier steps follow from them."""
    sample_config = config["samples"][sample]
    prefix = f"{assembler}/{run}/{assembly}"
    targets = [f"results/{sample}/pilon/{prefix}/{assembly}.fasta"]

    if not is_yes(sample, "annotation"):
        return targets

    targets += [
        f"results/{sample}/genbanks/{prefix}.{annotator}.genbank.rotated.check"
        for annotator in annotators(sample)
    ]

    if is_yes(sample, "run_images"):
        targets += [
            f"results/{sample}/images/{prefix}/recruitment_plot.check",
            f"results/{sample}/images/{prefix}/depth_plot.check",
            f"results/{sample}/images/{prefix}/ogdraw.check",
        ]

    # Chloroplast assemblies of NOVOPlasty are isomers of the same genome: only the first is searched
    first_isomer = re.match(r"(CA|OPT)_1_", assembly)
    if is_yes(sample, "run_nhmmer") and (
        sample_config["organelle"] == "mito" or assembler == "getorganelle" or first_isomer
    ):
        targets += [f"results/{sample}/nhmmer/{prefix}/{region}.nhmmer.check" for region in NHMMER_REGIONS]

    return targets

def run_targets(wildcards):
    """Every output of a run. It waits for split_assemblies to know the assemblies."""
    split_dir = checkpoints.split_assemblies.get(**wildcards).output[0]
    assemblies = glob_wildcards(os.path.join(split_dir, "{assembly}.fasta")).assembly
    return [
        target
        for assembly in assemblies
        for target in assembly_targets(wildcards.sample, wildcards.assembler, wildcards.run, assembly)
    ]

def get_output_files(wildcards):
    sample = wildcards.sample_id
    sample_config = config["samples"][sample]
    sample_output = []

    if sample_config["sequencing_type"] == "long":
        if is_yes(sample, "run_trimming"):
            sample_output.append(f"resources/{sample}/rawreads/{sample}.trimmed.check")

        sample_output += expand("results/{sample}/mitohifi/{seed}/contigs_stats.tsv", sample=sample, seed=sample_config["seeds"])

        if is_yes(sample, "run_nhmmer"):
            sample_output += expand(
                "results/{sample}/nhmmer/mitohifi/{seed}/{assembly}/{region}.nhmmer.check",
                sample=sample, seed=sample_config["seeds"], assembly=MITOHIFI_ASSEMBLY, region=NHMMER_REGIONS,
            )

        return sample_output

    if is_yes(sample, "run_trimming"):
        sample_output.append(f"resources/{sample}/rawreads/fastp.html")

    # One file per run, so the runs of a sample do not wait for each other's assemblies
    if is_yes(sample, "run_novoplasty"):
        sample_output += expand(
            "results/{sample}/finished_runs/novoplasty/{seed}_kmer{kmer}.done",
            sample=sample, seed=sample_config["seeds"], kmer=sample_config["kmers"],
        )

    if is_yes(sample, "run_getorganelle"):
        sample_output.append(f"results/{sample}/finished_runs/getorganelle/{GETORGANELLE_RUN}.done")

    return sample_output

rule finish_run:
    input:
        run_targets
    output:
        "results/{sample}/finished_runs/{assembler}/{run}.done"
    shell:
        """
        touch {output}
        """

rule parse_results:
    input:
        get_output_files
    output:
        "workflow/reports/{sample_id}/parse_results.done"
    log:
        "logs/{sample_id}.parse_results.log"
    shell:
        """
        python workflow/scripts/parse_results.py --sample {wildcards.sample_id} >> {log} 2>&1 && \
        touch {output}
        """
