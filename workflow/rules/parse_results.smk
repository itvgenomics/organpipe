def get_output_files(wildcards):
    sample_id = wildcards.sample_id

    sample_output = []

    sample_config = config["samples"][sample_id]

    if config["samples"][sample_id].get("sequencing_type", "").lower() == "short" and config["samples"][sample_id].get("run_novoplasty", "").lower() == "yes":
        kmers = [kmer for kmer in config["samples"][sample_id]["kmers"]]
        seeds = [seed for seed in config["samples"][sample_id]["seeds"]]

    if config["samples"][sample_id].get("sequencing_type", "").lower() == "short":

        if config["samples"][sample_id].get("run_trimming", "").lower() == "yes":
            sample_output.extend(expand("resources/{sample_id}/rawreads/fastp.html", sample_id=sample_id))
            sample_output.extend(expand("resources/{sample_id}/rawreads/fastp.json", sample_id=sample_id))

        if config["samples"][sample_id].get("run_novoplasty", "").lower() == "yes":
            sample_output.extend(expand("results/{sample_id}/novoplasty/{seed}/kmer{kmer}/log_{sample_id}.txt", sample_id=sample_id, seed=seeds, kmer=kmers))
            sample_output.extend(expand("results/{sample_id}/assemblies/novoplasty/{seed}_kmer{kmer}.fasta", sample_id=sample_id, seed=seeds, kmer=kmers))

            sample_output.extend(expand("results/{sample_id}/pilon/novoplasty/{seed}_kmer{kmer}.bwa_index.check", sample_id=sample_id, seed=seeds, kmer=kmers))
            sample_output.extend(expand("results/{sample_id}/pilon/novoplasty/{seed}_kmer{kmer}.bwa_mem.check", sample_id=sample_id, seed=seeds, kmer=kmers))
            sample_output.extend(expand("results/{sample_id}/pilon/novoplasty/{seed}_kmer{kmer}.samtools_index.check", sample_id=sample_id, seed=seeds, kmer=kmers))
            sample_output.extend(expand("results/{sample_id}/pilon/novoplasty/{seed}_kmer{kmer}.pilon.check", sample_id=sample_id, seed=seeds, kmer=kmers))

        if config["samples"][sample_id].get("run_getorganelle", "").lower() == "yes":
            sample_output.extend(expand("results/{sample_id}/getorganelle/get_org.log.txt", sample_id=sample_id))
            sample_output.extend(expand("results/{sample_id}/getorganelle/sequences.fasta", sample_id=sample_id))

            sample_output.extend(expand("results/{sample_id}/pilon/getorganelle/getorganelle_bwa_index.check", sample_id=sample_id))
            sample_output.extend(expand("results/{sample_id}/pilon/getorganelle/getorganelle_bwa_mem.check", sample_id=sample_id))
            sample_output.extend(expand("results/{sample_id}/pilon/getorganelle/getorganelle_samtools_index.check", sample_id=sample_id))
            sample_output.extend(expand("results/{sample_id}/pilon/getorganelle/getorganelle_pilon.check", sample_id=sample_id))

        if config["samples"][sample_id].get("annotation", "").lower() == "yes":
            if config["samples"][sample_id].get("run_novoplasty", "").lower() == "yes":
                if config["samples"][sample_id].get("organelle", "").lower() == "mito":
                    sample_output.extend(expand("results/{sample_id}/mitos2/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.mitos2.check", sample_id=sample_id, seed=seeds, kmer=kmers))
                    sample_output.extend(expand("results/{sample_id}/genbanks/novoplasty/{seed}_kmer{kmer}.mitos2.genbank.check", sample_id=sample_id, seed=seeds, kmer=kmers))
                    sample_output.extend(expand("results/{sample_id}/genbanks/novoplasty/{seed}_kmer{kmer}.mitos2.genbank.rotated.check", sample_id=sample_id, seed=seeds, kmer=kmers))

                    if config["samples"][sample_id].get("run_images", "").lower() == "yes":
                        sample_output.extend(expand("results/{sample_id}/images/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.get_genbank_fastas_novoplasty.check", sample_id=sample_id, seed=seeds, kmer=kmers))
                        sample_output.extend(expand("results/{sample_id}/images/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.blastn.check", sample_id=sample_id, seed=seeds, kmer=kmers))
                        sample_output.extend(expand("results/{sample_id}/images/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.recruitment_plot.check", sample_id=sample_id, seed=seeds, kmer=kmers))
                        sample_output.extend(expand("results/{sample_id}/images/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.bwa_index.check", sample_id=sample_id, seed=seeds, kmer=kmers))
                        sample_output.extend(expand("results/{sample_id}/images/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.bwa_mem.check", sample_id=sample_id, seed=seeds, kmer=kmers))
                        sample_output.extend(expand("results/{sample_id}/images/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.samtools_depth.check", sample_id=sample_id, seed=seeds, kmer=kmers))
                        sample_output.extend(expand("results/{sample_id}/images/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.samtools_depth_rotated.check", sample_id=sample_id, seed=seeds, kmer=kmers))
                        sample_output.extend(expand("results/{sample_id}/images/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.depth_plot.check", sample_id=sample_id, seed=seeds, kmer=kmers))
                        sample_output.extend(expand("results/{sample_id}/images/novoplasty/{sample_id}.mito.ogdraw.check", sample_id=sample_id))

                elif config["samples"][sample_id].get("organelle", "").lower() == "chloro":
                    sample_output.extend(expand("results/{sample_id}/cpgavas2/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.cpgavas2.check", sample_id=sample_id, seed=seeds, kmer=kmers))
                    sample_output.extend(expand("results/{sample_id}/genbanks/novoplasty/{seed}_kmer{kmer}.cpgavas2.genbank.check", sample_id=sample_id, seed=seeds, kmer=kmers))
                    sample_output.extend(expand("results/{sample_id}/genbanks/novoplasty/{seed}_kmer{kmer}.cpgavas2.genbank.rotated.check", sample_id=sample_id, seed=seeds, kmer=kmers))

                    sample_output.extend(expand("results/{sample_id}/chloe/novoplasty/{seed}_kmer{kmer}.chloe.check", sample_id=sample_id, seed=seeds, kmer=kmers))
                    sample_output.extend(expand("results/{sample_id}/genbanks/novoplasty/{seed}_kmer{kmer}.chloe.genbank.check", sample_id=sample_id, seed=seeds, kmer=kmers))
                    sample_output.extend(expand("results/{sample_id}/genbanks/novoplasty/{seed}_kmer{kmer}.chloe.genbank.rotated.check", sample_id=sample_id, seed=seeds, kmer=kmers))

                    if config["samples"][sample_id].get("run_images", "").lower() == "yes":
                        sample_output.extend(expand("results/{sample_id}/images/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.get_genbank_fastas_novoplasty.check", sample_id=sample_id, seed=seeds, kmer=kmers))
                        sample_output.extend(expand("results/{sample_id}/images/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.blastn.check", sample_id=sample_id, seed=seeds, kmer=kmers))
                        sample_output.extend(expand("results/{sample_id}/images/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.recruitment_plot.check", sample_id=sample_id, seed=seeds, kmer=kmers))
                        sample_output.extend(expand("results/{sample_id}/images/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.bwa_index.check", sample_id=sample_id, seed=seeds, kmer=kmers))
                        sample_output.extend(expand("results/{sample_id}/images/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.bwa_mem.check", sample_id=sample_id, seed=seeds, kmer=kmers))
                        sample_output.extend(expand("results/{sample_id}/images/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.samtools_depth.check", sample_id=sample_id, seed=seeds, kmer=kmers))
                        sample_output.extend(expand("results/{sample_id}/images/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.samtools_depth_rotated.check", sample_id=sample_id, seed=seeds, kmer=kmers))
                        sample_output.extend(expand("results/{sample_id}/images/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.depth_plot.check", sample_id=sample_id, seed=seeds, kmer=kmers))
                        sample_output.extend(expand("results/{sample_id}/images/novoplasty/{sample_id}.chloro.ogdraw.check", sample_id=sample_id))

                if config["samples"][sample_id].get("run_nhmmer", "").lower() == "yes":
                    sample_output.extend(expand("resources/nhmmer_db.hmm.{ext}", ext=['h3f', 'h3i', 'h3m', 'h3p']))

                    if config["samples"][sample_id].get("organelle", "").lower() == "mito":
                        sample_output.extend(expand("results/{sample_id}/nhmmer/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.ncRNA_sequences.check", sample_id=sample_id, seed=seeds, kmer=kmers))
                        sample_output.extend(expand("results/{sample_id}/nhmmer/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.ncRNA_nhmmer.check", sample_id=sample_id, seed=seeds, kmer=kmers))
                        sample_output.extend(expand("results/{sample_id}/nhmmer/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.intergenes_sequences.check", sample_id=sample_id, seed=seeds, kmer=kmers))
                        sample_output.extend(expand("results/{sample_id}/nhmmer/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.intergenes_nhmmer.check", sample_id=sample_id, seed=seeds, kmer=kmers))

                    elif config["samples"][sample_id].get("organelle", "").lower() == "chloro":
                        sample_output.extend(expand("results/{sample_id}/nhmmer/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.chloro.intergenes_sequences.check", sample_id=sample_id, seed=seeds, kmer=kmers))
                        sample_output.extend(expand("results/{sample_id}/nhmmer/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.chloro.intergenes_nhmmer.check", sample_id=sample_id, seed=seeds, kmer=kmers))
                        sample_output.extend(expand("results/{sample_id}/nhmmer/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.chloro.ncRNA_sequences.check", sample_id=sample_id, seed=seeds, kmer=kmers))
                        sample_output.extend(expand("results/{sample_id}/nhmmer/novoplasty/{seed}_kmer{kmer}/{seed}_kmer{kmer}.chloro.ncRNA_nhmmer.check", sample_id=sample_id, seed=seeds, kmer=kmers))

            if config["samples"][sample_id].get("run_getorganelle", "").lower() == "yes":
                if config["samples"][sample_id].get("organelle", "").lower() == "mito":
                    sample_output.extend(expand("results/{sample_id}/mitos2/getorganelle/getorganelle_mitos2.check", sample_id=sample_id))
                    sample_output.extend(expand("results/{sample_id}/genbanks/getorganelle/mitos2.genbank.check", sample_id=sample_id))
                    sample_output.extend(expand("results/{sample_id}/genbanks/getorganelle/mitos2.genbank.rotated.check", sample_id=sample_id))

                    if config["samples"][sample_id].get("run_images", "").lower() == "yes":
                        sample_output.extend(expand("results/{sample_id}/images/getorganelle/get_genbank_fastas_getorganelle.check", sample_id=sample_id))
                        sample_output.extend(expand("results/{sample_id}/images/getorganelle/blastn.check", sample_id=sample_id))
                        sample_output.extend(expand("results/{sample_id}/images/getorganelle/recruitment_plot.check", sample_id=sample_id))
                        sample_output.extend(expand("results/{sample_id}/images/getorganelle/bwa_index.check", sample_id=sample_id))
                        sample_output.extend(expand("results/{sample_id}/images/getorganelle/bwa_mem.check", sample_id=sample_id))
                        sample_output.extend(expand("results/{sample_id}/images/getorganelle/samtools_depth.check", sample_id=sample_id))
                        sample_output.extend(expand("results/{sample_id}/images/getorganelle/samtools_depth_rotated.check", sample_id=sample_id))
                        sample_output.extend(expand("results/{sample_id}/images/getorganelle/depth_plot.check", sample_id=sample_id))
                        sample_output.extend(expand("results/{sample_id}/images/getorganelle/mito.ogdraw.check", sample_id=sample_id))

                elif config["samples"][sample_id].get("organelle", "").lower() == "chloro":
                    sample_output.extend(expand("results/{sample_id}/cpgavas2/getorganelle/getorganelle_cpgavas2.check", sample_id=sample_id))
                    sample_output.extend(expand("results/{sample_id}/chloe/getorganelle/getorganelle_chloe.check", sample_id=sample_id))
                    sample_output.extend(expand("results/{sample_id}/genbanks/getorganelle/chloe.genbank.check", sample_id=sample_id))
                    sample_output.extend(expand("results/{sample_id}/genbanks/getorganelle/cpgavas2.genbank.check", sample_id=sample_id))
                    sample_output.extend(expand("results/{sample_id}/genbanks/getorganelle/cpgavas2.genbank.rotated.check", sample_id=sample_id))
                    sample_output.extend(expand("results/{sample_id}/genbanks/getorganelle/chloe.genbank.rotated.check", sample_id=sample_id))

                    if config["samples"][sample_id].get("run_images", "").lower() == "yes":
                        sample_output.extend(expand("results/{sample_id}/images/getorganelle/get_genbank_fastas_getorganelle.check", sample_id=sample_id))
                        sample_output.extend(expand("results/{sample_id}/images/getorganelle/blastn.check", sample_id=sample_id))
                        sample_output.extend(expand("results/{sample_id}/images/getorganelle/recruitment_plot.check", sample_id=sample_id))
                        sample_output.extend(expand("results/{sample_id}/images/getorganelle/bwa_index.check", sample_id=sample_id))
                        sample_output.extend(expand("results/{sample_id}/images/getorganelle/bwa_mem.check", sample_id=sample_id))
                        sample_output.extend(expand("results/{sample_id}/images/getorganelle/samtools_depth.check", sample_id=sample_id))
                        sample_output.extend(expand("results/{sample_id}/images/getorganelle/samtools_depth_rotated.check", sample_id=sample_id))
                        sample_output.extend(expand("results/{sample_id}/images/getorganelle/depth_plot.check", sample_id=sample_id))
                        sample_output.extend(expand("results/{sample_id}/images/getorganelle/chloro.ogdraw.check", sample_id=sample_id))

                if config["samples"][sample_id].get("run_nhmmer", "").lower() == "yes":
                    sample_output.extend(expand("resources/nhmmer_db.hmm.{ext}", ext=['h3f', 'h3i', 'h3m', 'h3p']))

                    if config["samples"][sample_id].get("organelle", "").lower() == "mito":
                        sample_output.extend(expand("results/{sample_id}/nhmmer/getorganelle/ncRNA_sequences.check", sample_id=sample_id))
                        sample_output.extend(expand("results/{sample_id}/nhmmer/getorganelle/ncRNA_nhmmer.check", sample_id=sample_id))
                        sample_output.extend(expand("results/{sample_id}/nhmmer/getorganelle/intergenes_sequences.check", sample_id=sample_id))
                        sample_output.extend(expand("results/{sample_id}/nhmmer/getorganelle/intergenes_nhmmer.check", sample_id=sample_id))


                    if config["samples"][sample_id].get("organelle", "").lower() == "chloro":
                        sample_output.extend(expand("results/{sample_id}/nhmmer/getorganelle/chloro.ncRNA_sequences.check", sample_id=sample_id))
                        sample_output.extend(expand("results/{sample_id}/nhmmer/getorganelle/chloro.ncRNA_nhmmer.check", sample_id=sample_id))
                        sample_output.extend(expand("results/{sample_id}/nhmmer/getorganelle/chloro.intergenes_sequences.check", sample_id=sample_id))
                        sample_output.extend(expand("results/{sample_id}/nhmmer/getorganelle/chloro.intergenes_nhmmer.check", sample_id=sample_id))

    elif config["samples"][sample_id].get("sequencing_type", "").lower() == "long":
        seeds = [seed for seed in config["samples"][sample_id]["seeds"]]

        if config["samples"][sample_id].get("run_trimming", "").lower() == "yes":
            sample_output.extend(expand("resources/{sample_id}/rawreads/{sample_id}.trimmed.check", sample_id=sample_id))

        sample_output.extend(expand("results/{sample_id}/mitohifi/{seed}/contigs_stats.tsv", sample_id=sample_id, seed=seeds))

        if config["samples"][sample_id].get("run_nhmmer", "").lower() == "yes":
            sample_output.extend(expand("results/{sample_id}/nhmmer/mitohifi/{seed}/intergenes_sequences.check", sample_id=sample_id, seed=seeds))
            sample_output.extend(expand("results/{sample_id}/nhmmer/mitohifi/{seed}/intergenes_nhmmer.check", sample_id=sample_id, seed=seeds))
            sample_output.extend(expand("results/{sample_id}/nhmmer/mitohifi/{seed}/ncRNA_sequences.check", sample_id=sample_id, seed=seeds))
            sample_output.extend(expand("results/{sample_id}/nhmmer/mitohifi/{seed}/ncRNA_nhmmer.check", sample_id=sample_id, seed=seeds))


    return sample_output

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
