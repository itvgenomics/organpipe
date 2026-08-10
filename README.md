# OrganPipe: A pipeline for the assembly, annotation, and curation of mitochondrial and chloroplast genomes

![OrganPipe Logo](OrganPipe.png)

---

> 🛠️ **Generate your config file online**
> Use our interactive web tool to easily create the configuration file for your run:
> **[https://itvgenomics.github.io/config_generator/](https://itvgenomics.github.io/config_generator/)**

---

## How to cite OrganPipe (preprint)

 ```
Moreira-Oliveira, R. R., Silva, B. M., Molina, M., Oliveira-Lima, M., Leão, T. F., Vasconcelos, S., & Nunes, G. L. (2025). OrganPipe: An automated tool to facilitate the assembly, annotation, and curation of mitochondrial and chloroplast genomes. https://doi.org/10.21203/RS.3.RS-5686696/V1
 ```
[https://doi.org/10.21203/RS.3.RS-5686696/V1](https://doi.org/10.21203/RS.3.RS-5686696/V1)

We encourage the users to also cite the tools in this pipeline:

```
Nicolas Dierckxsens, Patrick Mardulyn, Guillaume Smits, NOVOPlasty: de novo assembly of organelle genomes from whole genome data, Nucleic Acids Research, Volume 45, Issue 4, 28 February 2017, Page e18, https://doi.org/10.1093/nar/gkw955
```

```
Uliano-Silva, M., Ferreira, J.G.R.N., Krasheninnikova, K. et al. MitoHiFi: a python pipeline for mitochondrial genome assembly from PacBio high fidelity reads. BMC Bioinformatics 24, 288 (2023). https://doi.org/10.1186/s12859-023-05385-y
```

```
Alexander Donath, Frank Jühling, Marwa Al-Arab, Stephan H Bernhart, Franziska Reinhardt, Peter F Stadler, Martin Middendorf, Matthias Bernt, Improved annotation of protein-coding genes boundaries in metazoan mitochondrial genomes, Nucleic Acids Research, Volume 47, Issue 20, 18 November 2019, Pages 10543–10552, https://doi.org/10.1093/nar/gkz833
```

```
Walker BJ, Abeel T, Shea T, Priest M, Abouelliel A, Sakthikumar S, Cuomo CA, Zeng Q, Wortman J, Young SK, Earl AM. Pilon: an integrated tool for comprehensive microbial variant detection and genome assembly improvement. PLoS One. 2014 Nov 19;9(11):e112963. doi: 10.1371/journal.pone.0112963. PMID: 25409509; PMCID: PMC4237348.
```

```
Linchun Shi, Haimei Chen, Mei Jiang, Liqiang Wang, Xi Wu, Linfang Huang, Chang Liu, CPGAVAS2, an integrated plastome sequence annotator and analyzer, Nucleic Acids Research, Volume 47, Issue W1, 02 July 2019, Pages W65–W73, https://doi.org/10.1093/nar/gkz345
```

```
https://github.com/ian-small/Chloe.jl
```

```
Stephan Greiner, Pascal Lehwark, Ralph Bock, OrganellarGenomeDRAW (OGDRAW) version 1.3.1: expanded toolkit for the graphical visualization of organellar genomes, Nucleic Acids Research, Volume 47, Issue W1, 02 July 2019, Pages W59–W64, https://doi.org/10.1093/nar/gkz238
```

## How to Install OrganPipe Environment

### Prerequisites
Before installing the required software, make sure you have the following:
- A Linux-based operating system (e.g., Ubuntu, CentOS, Fedora)
>⚠️ **Windows Users**: If you are using Windows, OrganPipe must be run in a WSL2 local folder (not on a Windows folder, e.g /mnt/c) to ensure proper file system and performance compatibility.
- Python (version 3.5 or later) installed on your system
- [Git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) installed on your system


### Snakemake/Singularity/Docker Installation
1. **Clone the OrganPipe Repository**:
   Begin by cloning the OrganPipe project repository to your local machine. This will provide you with the necessary files, including an `environment.yaml` file to simplify the installation process.

   ```bash
   git clone https://github.com/itvgenomics/organpipe.git
   cd organpipe

2. **Install Conda or Mamba**:
    The OrganPipe pipeline requires a working Conda or Mamba installation to manage dependencies. You can find the installation instructions for these tools here:

    - [Conda Installation Guide](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)
    - [Mamba Installation Guide](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html)

    If you don't already have Conda or Mamba installed, follow the guides above to set them up on your system.

**Make sure you are installing Singularity, Snakemake and Docker compatible versions.**

3. **Install Docker**: Docker installation can be found [here](https://docs.docker.com/engine/install/). If you install our environment recommended versions, you can install the Docker 24.0.7 version.

4. **Create the OrganPipe Environment**:
    After setting up Conda or Mamba, activate your base environment and create a dedicated environment for OrganPipe using the provided conda_env.yaml file:

    ```bash
    conda env create -n organpipe -f conda_env.yaml

5. **Activate the OrganPipe Environment**:
    Once the environment is created, activate it to begin using OrganPipe:

    ```bash
    conda activate organpipe

### Additional Resources
- [Snakemake Documentation](https://snakemake.readthedocs.io/en/stable/)
- [Snakemake GitHub Repository](https://github.com/snakemake/snakemake)

## How to run OrganPipe

OrganPipe uses different assembly tools depending on your data and configuration. Here is a brief overview of how each tool works within the pipeline:

- **[NOVOPlasty](https://github.com/ndierckx/NOVOPlasty)** (Short Reads): A *de novo* assembler that requires a **seed sequence** (a fasta or genbank file) to initiate the assembly. It extends this seed iteratively using the provided short reads to assemble the circular organelle genome.
- **[GetOrganelle](https://github.com/Kinggerm/GetOrganelle)** (Short Reads): Uses a pre-compiled **database** (e.g., `animal_mt`, `embplant_pt`) rather than a single seed. It recruits reads mapping to the target database and performs a graph-based assembly to resolve the complete organelle genome.
- **[MitoHiFi](https://github.com/marcelauliano/MitoHiFi)** (Long Reads): Designed for PacBio HiFi or ONT reads. It queries NCBI to **download reference FASTA and GenBank files** based on a specified species name, and then uses these references to identify, filter, and circularize the mitochondrial contigs from your long-read assembly.

>⚠️ **Warning about parallelization**:
Snakemake uses the specified number of threads (CPUs) in local mode to check job availability. If the product of `-t` × `max_memory` exceeds your system's available RAM, it can cause the system to crash. Make sure your -t value and max_memory settings are compatible with your machine’s memory.

1. **Configure the Pipeline**:
    - Before running the pipeline, you must edit the configuration file to provide information specific to your samples and variables. We provide this file in two formats: YAML and CSV at the `config` directory. **If you don't need to use a specific variable, leave the value blank rather than removing the key/column. (e. g. reference: '')**

   - **YAML Format**: Each key in the file represents a variable, and each value corresponds to the parameters or settings for that variable.

   - **CSV Format**: Each column corresponds to a variable, and each row represents a sample with specific values for the variables.

The fields to be edited are the following:

| **Field**           | **Example**          | **Description/Comment**                                                                                                                                                              | **Field Requirement**                         |
|----------------------|----------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------------------------------|
| `sample`            | `"sample1"`           | Specify the sample name to process. Use "all" to process all samples in `reads_path`. You can also set a comma separated list of samples (e.g. "sample1,sample2,sample3"). All samples will be running with the same parameters, if you do want different values for each sample, use the .csv format | Required.                                     |
| `reads_path`        | `"/path/to/reads"`    | Directory containing sequencing read files (_R1/_R2 or _pair1/_pair2) in fastq.gz format. For short reads, ensure that only the sequencing read files (in fastq.gz format) are present. For long reads, in fasta or in fastq.gz format.                                                                                         | Required.                                     |
| `organelle`         | `"mito"`              | Type of organelle genome to assemble: "mito" for mitochondrial or "chloro" for chloroplast.                                                                                         | Required.                                     |
| `genetic_code`      | `2`                 | Genetic code for mitochondrial genome annotation (e.g., 2 for The Vertebrate Mitochondrial Code).                                                                                                       | Required.                                     |
| `reference`         | `"/path/to/ref.fasta"`| Path to a reference genome in fasta format. Leave blank if unavailable. This is used by NOVOPlasty as a guide to resolve duplicated regions in the plant mitochondria or the inverted repeat in the chloroplast                                                                                                            | Optional.                                     |
| `sequencing_type`   | `"Short"`             | Type of sequencing data: "Short" for Illumina or "Long" for PacBio/ONT reads.                                                                                                       | Required.                                     |
| `genome_range`      | `"12000-22000"`       | Expected genome size range in kilobases. Every sequence between this range will be annotated if `annotation` is set to "yes"                                                                                                                                             | Required.                                     |
| `run_trimming`    | `"Yes"`               | Run Fastp for adapter anda quality trimming: "Yes" to perform trimming, "No" to skip. For long reads, only .fastq.gz files are supported.                                                                                                  | Optional, default is "No".                   |
| `adapters`          | `"/path/to/adapters.fasta"` | Path to a text file listing adapter sequences to remove. Use this field if Fastp do not support your adapters.                                                                                                                             | Optional        |
| `pacbio_adapters`          | `"-b ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT"` | Set the adapters used at Cutadapt to perform the trimming. You can specify multiple adapters like: `"-b ATTGCGCGCGTA -b TGCATGGCTAGT"`                                                                                                                           | Optional        |
| `minlength`         | `50`                | Minimum read length to retain after trimming.                                                                                                                                       | Required if `run_trimming` is "Yes".        |
| `minquality`        | `20`                | Minimum base quality score threshold for trimming.                                                                                                                                 | Required if `run_trimming` is "Yes".        |
| `seed_format`       | `"fasta"`             | Format of the seed file: "fasta" or "genbank".                                                                                                                                       | Optional.                                     |
| `seed_file`         | `"/path/to/seed.fasta"`| Path to the seed file to initialize assembly.                                                                                                                                         | Optional.                                     |
| `feature`           | `"CDS"`               | Feature type for assembly if GenBank seed file is used: "CDS", "rRNA", or "tRNA".                                                                                                   | Required if `seed_format` is "genbank".       |
| `search_ncbi`       | `"Yes"`               | Search NCBI for reference sequences: "Yes" to search, "No" to skip.                                                                                                                 | Optional.                                     |
| `search_genes`      | "`COI,16S,ATP6"`      | Genes to search for on NCBI if `search_ncbi` is enabled.                                                                                                                            | Required if `search_ncbi` is "Yes".           |
| `search_term`       | `"Amphisbaena"`       | If no matching record is found, an error will be raised.                                                                                                                                             | Required if `search_ncbi` is "Yes".           |
| `max_references`    | `5`                | Maximum reference sequences to download per gene.                                                                                                                                  | Required if `search_ncbi` is "Yes".                                     |
| `kmers`             | `"19,23,33,39"`   | List of k-mer sizes for genome assembly.                                                                                                                                             | Required for Short Reads assembly.               |
| `max_memory`        | `4`                | Maximum memory (in GB) for genome assembly.                                                                                                                                         | Required for Short Reads assembly.                  |
| `reads_length`      | `150`               | Read length (in base pairs) of sequencing data.                                                                                                                                      | Required for Short Reads.                  |
| `insert_size`       | `300`               | Average insert size (in base pairs) of sequencing data.                                                                                                                             | Required for Short Reads.                  |
| `annotation`        | `"Yes"`               | Run annotation pipeline: "Yes" to annotate, "No" to skip.                                                                                                                           | Optional. Default is "No".                    |
| `run_nhmmer`        | `"No"`                | Run nhmmer to identify ncRNA and intergenic regions. **Note:** Enabling this can  slow down the pipeline.                                                                                                                               | Optional. Default is "No".                    |
| `run_images`        | `"Yes"`               | Generate visualizations like OGDraw diagrams and depth plots. **Note:** Enabling this can slow down the pipeline.                                                                                                                      | Optional. Default is "Yes".                   |
| `run_novoplasty`    | `"Yes"`               | Run NOVOPlasty for genome assembly: "Yes" to run, "No" to skip. You can run both assemblers simultaneously by setting both `run_novoplasty` and `run_getorganelle` to "Yes".                                                          | Optional. Default is "Yes".                   |
| `run_getorganelle`  | `"No"`                | Run GetOrganelle for genome assembly: "Yes" to run, "No" to skip. You can run both assemblers simultaneously by setting both `run_novoplasty` and `run_getorganelle` to "Yes".                                                         | Optional. Default is "No".                    |
| `database`          | `"animal_mt"`         | GetOrganelle organelle type (e.g., `embplant_pt`, `other_pt`, `embplant_mt`, `embplant_nr`, `animal_mt`, `fungus_mt`, `fungus_nr`). Multiple types can be combined with commas.                                                       | Required if `run_getorganelle` is "Yes".      |
| `n_rounds`          | `10`                  | Maximum number of extending rounds for GetOrganelle (suggested: ≥ 2). Defaults vary by organelle type (e.g., 15 for `embplant_pt`, 10 for `animal_mt`).                                                                              | Optional (GetOrganelle).                      |
| `target_size`       | `13000`               | Hypothetical target genome size used by GetOrganelle to estimate word size. Defaults vary by organelle type. Should be a comma-separated list of integers in multi-organelle mode.                                                     | Optional (GetOrganelle).                      |
| `spades_kmers`      | `"21,55,85,115"`      | SPAdes k-mer settings passed to GetOrganelle. Use the same format as SPAdes (e.g., `21,55,85,115`).                                                                                                                                   | Optional (GetOrganelle).                      |
| `extra_flags`       | `"--overwrite"`       | Any additional flags to pass directly to the GetOrganelle command line.                                                                                                                                                               | Optional (GetOrganelle).                      |
| `search_species`    | `"Amphisbaena"`       | Taxon name used when searching NCBI for complete mitogenome references (MitoHiFi / long reads only).                                                                                                                                  | Required for Long Reads.                      |
| `n_references`      | `5`                   | Maximum number of reference sequences to download during the NCBI search (MitoHiFi / long reads only).                                                                                                                               | Required for Long Reads.                      |

>If search_ncbi is set to "Yes" and no matching sequences are found for the provided `search_term`, OrganPipe will raise an error and stop the pipeline. Make sure the term you are searching for has available sequences in NCBI before running the pipeline.

2. **Run OrganPipe**:
    - Run the pipeline with default parameters:

    ```
    bash OrganPipe.sh -d </path/to/work/dir> -t <n_threads> -c </path/to/configfile>
    ```

    - **Flags**:
        - **-d** </path/to/work/dir> (Required) = Path to your working directory where all the workflow file are
        - **-c** </path/to/config.yaml> (Required) = Overwrite the default configuration file with all needed parameters (e.g. config/config.yaml/csv)
        - **-t** {int} (Required) = Number of threads to use. If running in **SLURM** mode, this value determines how many jobs will be submitted to the queue.
        - **-np** (Optional) = Perform a dry run to see what jobs will be executed without actually running them.
        - **-unlock** (Optional) = Unlock the working directory if Snakemake has somehow locked it.
        - **-batch** (Optional) = If you are running a large number of samples, or number of rules executed > 3000, consider using this flag. This slightly improves the DAG resolution time from Snakemake. You can set the number with `-nbatch` (Default = 15)
        - **-sifdir** (Optional) = Choose a directory to build all singularity image files used in the pipeline. If the path already contains the images, they will not be pulled. Default: resources/sif_dir
        - **-rerun** (Optional) = Delete previous results and temporary files for the specified sample(s) to ensure a clean re-run with updated configurations. Use this when reprocessing samples with different parameters.
        - **-nhmmer_db** {path} (Optional) = Path to the HMM database used by nhmmer. If not specified, defaults to `resources/rfam.hmm`. You can use any HMMER database compatible with HMMER version 3.4 for improved accuracy.
        - **-slurm** (Optional) = Use the `config/slurm_params.yaml` file to run the workflow with SLURM job submission using Snakemake's profile system. This enables use of SLURM-specific resource configuration, submission rules, and cluster-specific options. If you want to change any default SLURM settings, such as the partition: Edit `config/slurm_params.yaml` and set the appropriate value for the `slurm_partition` variable.
        - **-j** {int} (Optional, **-slurm** mode only) = Controls how many jobs are submitted to the SLURM queue at once.
        - **-partition** {string} (Required when **-slurm** is used) — Specifies the SLURM partition (queue) to which the jobs will be submitted.

    - If you want to change any default settings, such as the threads number and memory usage: Edit config/local_params.yaml. **DO NOT CHANGE THE `{WORKDIR}`, `{THREADS}` and `{PARTITION}` VARIABLES**. Change `{PARTITION}` only if you want to redirect specific rules to diffent partitions.

    - **OOM (Out Of Memory) Errors**: If a specific rule fails due to an Out-Of-Memory error during pipeline execution, you can increase the `mem_mb` directive for that specific rule in the `config/local_params.yaml` file.

    - We recommend initially running the pipeline with the -np (dry run) flag. This will allow you to verify that all paths and configurations are correct and that the pipeline will execute as expected. It's a good way to ensure everything is set up properly before running the actual workflow.

    ```
    bash OrganPipe.sh -d </path/to/work/dir> -t 1 -c </path/to/configfile> -np
    ```

3. **Testing with Sample Data**:
    - To test if everything is set up correctly, you can run the pipeline using a test dataset provided in the test_data directory. This ensures that the pipeline is functioning as expected before working with your own data.

    - Make sure you are in the directory where you cloned the Git repository, then execute the following command:

    ```
    bash OrganPipe.sh -d . -t 4 -c test_data/config.csv
    ```

4. **Checking the Results**:

All results will be compiled in the `workflow/reports` directory. If you want to check the raw output files from the software used in the pipeline, you can find them in the `workflow/results` directory.

Each sample directory contains the following folders:

| Folder | Description |
|:--------|:-------------|
| `fastas/` | Contains all assembly FASTA files. |
| `files/` | Includes important files to assist with manual curation. |
| `genbanks/` | Contains GenBank format output files. |
| `genes/` | Holds all assembled genes for each seed/k-mer. |
| `mitos2/` | Contains mitochondrial-specific results (mito only). |
| `novoplasty/` | Results from NOVOPlasty (short reads only). |
| `getorganelle/` | Results from GetOrganelle (short reads only, when enabled). |
| `pilon/` | Polishing results from Pilon (short reads only). |
| `mitohifi/` | Results from MitoHiFi (long reads only). |
| `nhmmer/` | HMMER-based search results. |
| `cpgavas2/` | Chloroplast-specific results (chloro only). |

Each sample folder also includes the following key files:

| File | Description |
|:------|:-------------|
| `summary.csv` | Summary of all runs and results. |
| `mitos2.csv` | Mitochondrial assembly metrics (mito only). |
| `novoplasty.csv` | NOVOPlasty results (short reads only). |
| `pilon.csv` | Pilon results (short reads only). |
| `mitohifi.csv` | MitoHiFi results (long reads only). |
| `nhmmer_intergenes.csv` | HMMER intergenic region annotations. |
| `nhmmer_ncRNA.csv` | HMMER ncRNA annotations. |
| `cpgavas2_codon_usage.csv` | Codon usage statistics (chloro only). |
| `cpgavas2_gene_composition.csv` | Gene composition summary (chloro only). |
| `cpgavas2_intron_exon.csv` | Intron–exon structure data (chloro only). |
| `cpgavas2_problems.csv` | Quality or annotation issues (chloro only). |

Entries in the `summary.csv` file follow this format:

**NOVOPlasty assemblies:**

```
{CA/Option/Contig}_{sample}_{seed}_{TaxID}_{kmer}_1
```

- **CA** → Circularized Assembly (NOVOPlasty successfully assembled one complete genome)
- **Option** → Multiple circularizations detected for the same sample
- **Contig** → Non-circularized assembly (the genome could not be closed)
- **sample** → Sample identifier (e.g., ITV00872)
- **seed** → Seed number used in the assembly
- **TaxID** → NCBI taxonomic identifier
- **kmer** → K-mer size used

For the entry: `CA_1_ITV00872_25_1-COI_941666_19_1`

- **CA** → Circularized Assembly
- **ITV00872_25** → Sample name
- **1-COI** → Seed used in assembly
- **941666** → TaxID from the seed`s organism
- **19** → K-mer size used
- **_1** → Assembly Number

**GetOrganelle assemblies:**

```
{database}_{status}_{number}
```

- **database** → GetOrganelle database used, as set in the config file (e.g., `animal_mt`, `embplant_pt`)
- **status** → Assembly status: `complete` for a closed circular genome, `scaffold` for a non-closed sequence
- **number** → Assembly number (incremental, starting at 1)

For the entry: `animal_mt_complete_1`

- **animal_mt** → Database used (animal mitogenome)
- **complete** → Genome was successfully circularized
- **1** → First assembly produced

# Changelog

## OrganPipe 1.2.0 – Changelog

### New Features
- Introduced a new `-j` flag in SLURM mode to control how many jobs are submitted to the queue.
- Converted result parsing into a Snakemake rule:
  - Parsing now runs automatically when each sample completes all its jobs.
  - Reports are generated incrementally during pipeline execution instead of waiting until the end.
- Added **GetOrganelle** to the pipeline as an alternative (or complementary) assembler:
  - New config directives: `database`, `n_rounds`, `target_size`, `spades_kmers`, `extra_flags`.
  - New config directive `run_getorganelle` to enable/disable GetOrganelle.
- Users can now choose which assembler to use (`run_novoplasty`, `run_getorganelle`), or run both simultaneously.
- Added `-nhmmer_db` flag: if not specified, defaults to `resources/rfam.hmm`;
    - Removed the `nhmmer_db` directive from the config file.
- New long-read config directives: `search_species` and `n_references` (replaces earlier naming).

### Improvements
- Updated MitoS to version 2.1.10.
- Improved robustness when handling MitoS BLAST outputs with fewer than 12 columns.
- Updated the all test data and its config files.

### Bug Fixes
- Fixed parsing issues when `result.geneorder` is missing from MitoS2 output files.
- Attempted to resolve failures during Singularity `.sif` image builds.
- Fixed when `novoplasty.csv` returned empty columns.

## OrganPipe v1.1 — Change Log

### Bug Fixes
- Fixed several incorrect variable examples in configuration files.
- Fixed `stderr` output handling for **bwa-mem2**.
- Fixed how the reference FASTA is handled in **NOVOPlasty**.
- Fixed `parse_results.py` when **CpGAVAS2** fails to annotate.
- Fixed rule `all` inputs so Snakemake does not re-run finished samples when re-running with an updated config file.
- Fixed issues when using `all` as the sample name.
- Added proper help output when an invalid flag is passed to `OrganPipe.sh`.

### Read Trimming Improvements
- Added **fastp** as an alternative to AdapterRemoval:
  - The config variable `adapterremoval` was renamed to `run_trimming`.
  - fastp now **auto-detects adapters**, removing the need for an `adapters.txt` file.
  - A custom adapters file can still be used if fastp fails to detect adapters automatically.

### Re-run Support
- Added `--rerun` flag:
  - Deletes previous results and temporary files for selected sample(s).
  - Ensures a clean re-run when changing parameters or configurations.
  - Recommended when reprocessing samples with updated settings.

### Long Reads Support
- Improved handling of long-read sequencing data.
- Added test data for long reads.
- Added an example configuration file for long-read workflows.

### Documentation & Config Improvements
- Updated the **README** to document all new features.
- Improved comments in `config/config.yaml` for better readability.
- Added a **config check script** to verify that all required parameters are correctly set.

### SLURM Support
- OrganPipe now supports **SLURM**:
  - Default SLURM parameters can be modified in `profiles/slurm/config.yaml`.
  - Users must set their partition using the `slurm_partition` variable.

### Output & Reporting Changes
- Replaced **CIRCOS** with **OGDRAW** for generating circular mitogenome/plastome images.
- Renamed `abstract.csv` to `summary.csv`.
- Chloroplast short-read assemblies now include a summary.
- Reports generation was reworked — see the *Results* section in the README for details.

# Output Structure

All outputs are located in: `workflow/reports/<sample_name>/`


## Directories
- `fastas` — All assembly FASTA files
- `files` — Important files for manual curation
- `genbanks`
- `genes` — Gene assemblies for each seed/k-mer
- `mitos2` — Mitochondrial assemblies only
- `novoplasty` — Short reads only
- `pilon` — Short reads only
- `mitohifi` — Long reads only
- `nhmmer`
- `cpgavas2` — Chloroplast only

## Output Files
- `summary.csv`
- `mitos2.csv` (mitochondria only)
- `novoplasty.csv` (short reads only)
- `getorganelle.csv` (short reads only, when enabled)
- `pilon.csv` (short reads only)
- `mitohifi.csv` (long reads only)
- `nhmmer_intergenes.csv`
- `nhmmer_ncRNA.csv`
- `cpgavas2_codon_usage.csv` (chloroplast only)
- `cpgavas2_gene_composition.csv` (chloroplast only)
- `cpgavas2_intron_exon.csv` (chloroplast only)
- `cpgavas2_problems.csv` (chloroplast only)
