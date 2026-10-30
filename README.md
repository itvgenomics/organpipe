# OrganPipe: A pipeline for the assembly, annotation, and curation of mitochondrial and chloroplast genomes

![OrganPipe Logo](OrganPipe.png)

---

> 🛠️ **Generate your config file online**
> Use our interactive web tool to easily create the configuration file for your run:
> **[https://itvgenomics.github.io/config_generator/](https://itvgenomics.github.io/config_generator/)**

---

OrganPipe is a Snakemake workflow that assembles, polishes, annotates and reports mitochondrial and chloroplast genomes from short (Illumina) or long (PacBio HiFi / ONT) reads. Every tool runs inside a Singularity container, so the only software you install is the conda environment.

## Contents

- [Quick start](#quick-start)
- [Installation](#installation)
- [How OrganPipe works](#how-organpipe-works)
- [Input reads](#input-reads)
- [Configuration](#configuration)
- [Running OrganPipe](#running-organpipe)
- [Threads, memory and SLURM resources](#threads-memory-and-slurm-resources)
- [Container images](#container-images)
- [Outputs](#outputs)
- [Logs and benchmarks](#logs-and-benchmarks)
- [Troubleshooting](#troubleshooting)
- [How to cite](#how-to-cite)
- [Changelog](CHANGELOG.md)

## Quick start

```bash
git clone https://github.com/itvgenomics/organpipe.git && cd organpipe
conda env create -n organpipe -f conda_env.yaml && conda activate organpipe
bash OrganPipe.sh -d . -t 4 -c test_data/config.csv
```

The last command runs the bundled test dataset (see [Testing with the sample data](#testing-with-the-sample-data)). To run your own samples, fill in a [configuration file](#configuration) and point `-c` to it.

## Installation

### Requirements

- **Linux** (e.g. Ubuntu, CentOS, Fedora).
  > ⚠️ **Windows users**: run OrganPipe inside a WSL2 local folder, not on a mounted Windows drive (e.g. `/mnt/c`), to avoid file system and performance problems.
- **[Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) or [Mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html)**. The conda environment provides Python 3.12, Snakemake 8.5.1, the Snakemake SLURM plugin and Singularity 3.8.5.
- **[Git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)**.
- **Internet access** during the run, to:
  - pull the Singularity images on the first run;
  - search NCBI for seeds (`search_ncbi: "Yes"`) and for MitoHiFi references (long reads);
  - download the GetOrganelle databases (`run_getorganelle: "Yes"`).
- **Disk space** for the Singularity images (stored once in `resources/sif_dir`, or in the directory given with `-sifdir`), plus the intermediate files of your samples.
- **Enough RAM** for the assembly step. See [Threads, memory and SLURM resources](#threads-memory-and-slurm-resources).

Docker is **not** required to run OrganPipe. Singularity pulls the images directly from Docker Hub. Docker is only needed if you want to rebuild the images from the `dockerfiles/` directory.

### Install

1. **Clone the repository**

   ```bash
   git clone https://github.com/itvgenomics/organpipe.git
   cd organpipe
   ```

2. **Create and activate the conda environment**

   ```bash
   conda env create -n organpipe -f conda_env.yaml
   conda activate organpipe
   ```

### Additional resources

- [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/)
- [Snakemake GitHub repository](https://github.com/snakemake/snakemake)

## How OrganPipe works

OrganPipe chooses the assembler from your data and configuration:

- **[NOVOPlasty](https://github.com/ndierckx/NOVOPlasty)** (short reads): a *de novo* assembler that starts from a **seed sequence** (FASTA or GenBank) and extends it with the reads until the circular genome is closed. Seeds come from your own file or from an NCBI search. Every seed is assembled with every k-mer in `kmers`.
- **[GetOrganelle](https://github.com/Kinggerm/GetOrganelle)** (short reads): recruits reads with a pre-compiled **database** (e.g. `animal_mt`, `embplant_pt`) instead of a single seed, and resolves the genome from the assembly graph.
- **[MitoHiFi](https://github.com/marcelauliano/MitoHiFi)** (long reads, mitochondria only): downloads **reference FASTA and GenBank files** from NCBI for `search_species`, then uses them to identify, filter and circularize the mitochondrial contigs of the long-read assembly.

For short reads you can run NOVOPlasty, GetOrganelle or both. The main steps after assembly are:

| Step | Tools | When it runs |
|:-----|:------|:-------------|
| Trimming | fastp (adapters and low-quality bases, short reads), cutadapt (long reads) | `run_trimming: "Yes"` |
| Error correction | bwa-mem2, samtools, Pilon | Short reads, always |
| Annotation | MITOS2 (mito), CPGAVAS2 and Chloë (chloro) | `annotation: "Yes"` |
| Images | OGDRAW maps, read recruitment and depth plots | `run_images: "Yes"` (short reads only) |
| ncRNA / intergenic search | HMMER `nhmmscan` (`-nhmmer_db`, default Rfam) | `run_nhmmer: "Yes"` |
| Reports | `parse_results.py` | Always, as soon as each sample finishes |

## Input reads

All read files go directly in the directory set in `reads_path`, and their names must contain the sample name as a whole name: `sample1` matches `sample1_R1.fastq.gz` and `sample1_S1_R1_001.fastq.gz`, but not `sample10_R1.fastq.gz`. A file named exactly after the sample is preferred. If a sample still matches more than one file (e.g. one per lane), the run stops: merge the files first.

- **Short reads**: paired-end files in `.fastq.gz` (or `.fasta.gz`), with `_R1`/`_R2` or `_pair1`/`_pair2` marking the read pair, e.g. `sample1_R1.fastq.gz` and `sample1_R2.fastq.gz`. Keep only read files in this directory. Files whose name contains `Undetermined` are ignored when `sample: "all"` is used.
- **Long reads**: one file per sample, in `.fastq.gz` or uncompressed `.fasta`, e.g. `sample1.fastq.gz`. Trimming is only supported for `.fastq.gz`, and the run stops if `run_trimming` is `Yes` for a `.fasta` file.

## Configuration

Before running the pipeline, describe your samples in a configuration file. Templates are provided in `config/` in two formats:

- **YAML** (`config/config.yaml`): one set of parameters applied to every sample listed in `sample`.
- **CSV** (`config/config.csv`): one row per sample, so each sample can have its own parameters. Columns may be separated by `,` or `;`.

> **Keep every key/column.** If you don't need a field, leave its value blank instead of removing it (e.g. `reference: ''`).

`Yes`/`No` values are case-insensitive. The configuration is validated by `workflow/scripts/check_config.py` before every run, and the run stops with an explanatory error if a required field is missing or invalid. The online [config generator](https://itvgenomics.github.io/config_generator/) creates a valid file for you.

### General

| Field | Example | Description | Requirement |
|:------|:--------|:------------|:------------|
| `sample` | `"sample1"` | Sample name, with only letters, numbers and underscores. Use `"all"` to process every sample in `reads_path`, or a comma-separated list (`"sample1,sample2"`). In the YAML format all listed samples share the same parameters; use the CSV format to set different values per sample. | Required |
| `reads_path` | `"/path/to/reads"` | Directory with the read files (see [Input reads](#input-reads)). | Required |
| `organelle` | `"mito"` | `"mito"` for mitochondrial or `"chloro"` for chloroplast genomes. Chloroplast assembly is only supported for short reads. | Required |
| `genetic_code` | `2` | NCBI genetic code used for annotation (e.g. `2` for the vertebrate mitochondrial code, `11` for plastids). | Required |
| `sequencing_type` | `"Short"` | `"Short"` for Illumina or `"Long"` for PacBio/ONT reads. | Required |
| `genome_range` | `"12000-22000"` | Expected genome size range, **in base pairs** (`min-max`). Assemblies within this range are kept and annotated. | Required |

### Pipeline steps

For short reads, each of these fields must be set to `Yes` or `No`. At least one of `run_novoplasty` and `run_getorganelle` must be `Yes`.

| Field | Example | Description | Requirement |
|:------|:--------|:------------|:------------|
| `run_trimming` | `"Yes"` | Trim adapters and low-quality bases: fastp for short reads, cutadapt (with `pacbio_adapters`) for long reads. | Required (short reads) |
| `annotation` | `"Yes"` | Annotate the assemblies. Must be `"Yes"` for long reads. | Required |
| `run_nhmmer` | `"No"` | Search ncRNA and intergenic regions with HMMER. **Note:** this can slow down the pipeline. | Required (short reads) |
| `run_images` | `"No"` | Generate OGDRAW maps, read recruitment and depth plots. Not supported for long reads. **Note:** this can slow down the pipeline. | Required (short reads) |
| `run_novoplasty` | `"Yes"` | Assemble with NOVOPlasty. Can be combined with `run_getorganelle`. | Required (short reads) |
| `run_getorganelle` | `"No"` | Assemble with GetOrganelle. Can be combined with `run_novoplasty`. | Required (short reads) |

### Trimming

| Field | Example | Description | Requirement |
|:------|:--------|:------------|:------------|
| `adapters` | `"/path/to/adapters.fasta"` | FASTA file with adapter sequences for fastp. Leave blank to let fastp detect the adapters automatically. | Optional (short reads) |
| `minlength` | `50` | Minimum read length kept after trimming. | Required if `run_trimming` is `Yes` (short reads) |
| `minquality` | `20` | Quality threshold of the fastp sliding window: low-quality bases are cut from both read ends (`--cut_front --cut_tail --cut_mean_quality`). | Required if `run_trimming` is `Yes` (short reads) |
| `pacbio_adapters` | `"-b ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT"` | Adapters removed by cutadapt. Repeat `-b` for several adapters: `"-b ATTGCGCGCGTA -b TGCATGGCTAGT"`. | Required if `run_trimming` is `Yes` (long reads) |

> The NOVOPlasty manual recommends removing only the adapters, without quality trimming. OrganPipe keeps quality trimming because it helps with degraded data, such as museum samples. With good-quality reads and NOVOPlasty, consider `run_trimming: "No"` or a low `minquality`.

### Seeds and reference (NOVOPlasty)

Provide a seed file **or** let OrganPipe search NCBI for seeds.

| Field | Example | Description | Requirement |
|:------|:--------|:------------|:------------|
| `seed_format` | `"fasta"` | Format of `seed_file`: `"fasta"` (`.fasta`, `.fa`, `.fna`) or `"genbank"` (`.gb`, `.gbk`, `.genbank`). Leave blank to search NCBI instead. Must be blank for long reads. | Optional |
| `seed_file` | `"/path/to/seed.fasta"` | Seed file used to start the assembly. | Required if `seed_format` is set |
| `feature` | `"CDS"` | Features extracted from a GenBank seed: `"CDS"`, `"rRNA"` or `"tRNA"` (case-sensitive). | Required if `seed_format` is `"genbank"` |
| `search_ncbi` | `"Yes"` | Download seeds from NCBI. Required when no seed file is given and NOVOPlasty is enabled. | Required (short reads) |
| `search_genes` | `"COI,16S,ATP6"` | Genes to search on NCBI. | Required if `search_ncbi` is `Yes` |
| `search_term` | `"Amphisbaena"` | Taxon name used in the NCBI search. If no record is found, the run stops with an error. | Required if `search_ncbi` is `Yes` |
| `max_references` | `5` | Maximum number of seeds downloaded per gene. | Required if `search_ncbi` is `Yes` |
| `reference` | `"/path/to/ref.fasta"` | Reference genome (FASTA) used by NOVOPlasty to resolve duplicated regions in plant mitochondria or the inverted repeat in chloroplasts. Leave blank if unavailable. | Optional |

> If `search_ncbi` is `Yes` and no sequence matches `search_term`, OrganPipe stops with an error. Check that NCBI has sequences for your taxon before running, or use a broader taxon (e.g. the genus or family).

### NOVOPlasty

| Field | Example | Description | Requirement |
|:------|:--------|:------------|:------------|
| `kmers` | `"19,23,33,39"` | K-mer sizes. Each seed is assembled once per k-mer. | Required if `run_novoplasty` is `Yes` |
| `max_memory` | `40` | Maximum memory (GB) used by NOVOPlasty per job. In SLURM mode, these jobs request this value plus 10%. | Required if `run_novoplasty` is `Yes` |
| `reads_length` | `150` | Read length (bp). | Required if `run_novoplasty` is `Yes` |
| `insert_size` | `300` | Average insert size (bp). | Required if `run_novoplasty` is `Yes` |

### GetOrganelle

| Field | Example | Description | Requirement |
|:------|:--------|:------------|:------------|
| `database` | `"animal_mt"` | GetOrganelle organelle type (`embplant_pt`, `other_pt`, `embplant_mt`, `embplant_nr`, `animal_mt`, `fungus_mt`, `fungus_nr`). Several types can be combined with commas. | Required if `run_getorganelle` is `Yes` |
| `n_rounds` | `10` | Maximum number of extension rounds (suggested ≥ 2). Defaults depend on the organelle type (e.g. 15 for `embplant_pt`, 10 for `animal_mt`). | Optional |
| `target_size` | `13000` | Expected genome size used to estimate the word size. Defaults depend on the organelle type. Comma-separated list in multi-organelle mode. | Optional |
| `spades_kmers` | `"21,55,85,115"` | SPAdes k-mers passed to GetOrganelle. | Optional |
| `extra_flags` | `"--overwrite"` | Additional flags passed directly to GetOrganelle. | Optional |

### Long reads (MitoHiFi)

| Field | Example | Description | Requirement |
|:------|:--------|:------------|:------------|
| `search_species` | `"Amphisbaena"` | Taxon used to search NCBI for complete mitogenome references. | Required for long reads |
| `n_references` | `5` | Maximum number of references downloaded from NCBI. | Required for long reads |

## Running OrganPipe

```bash
bash OrganPipe.sh -d </path/to/work/dir> -c </path/to/configfile> -t <cores>
```

We recommend a dry run (`-np`) first. It checks the paths and the configuration and lists the jobs that will run, without running them:

```bash
bash OrganPipe.sh -d </path/to/work/dir> -c </path/to/configfile> -t 1 -np
```

> The number of assemblies is only known after each assembler finishes, so a dry run only lists the jobs up to the assembly step. Once the assemblies exist, it also lists the jobs that follow them.

### Flags

| Flag | Mode | Description |
|:-----|:-----|:------------|
| `-d <path>` | All (required) | Working directory, i.e. the OrganPipe directory with the `workflow/` and `config/` folders. `OrganPipe.sh` can be called from any directory; relative paths in the configuration file are relative to this directory. |
| `-c <file>` | All (required) | Configuration file (`.yaml` or `.csv`). |
| `-t <int>` | Local (required) | Number of CPU cores Snakemake may use at once. Not allowed with `-slurm`. |
| `-slurm` | SLURM | Submit each job to SLURM, using the settings in `config/slurm_params.yaml`, `config/threads.yaml` and `config/resources.yaml`. |
| `-j <int>` | SLURM (required) | Maximum number of jobs queued or running at once. |
| `-partition <name>` | SLURM (required) | Default SLURM partition (queue) for all jobs. |
| `-np` | All | Dry run: show the jobs that would run without running them. |
| `-unlock` | All | Unlock the working directory after an interrupted run (see [Troubleshooting](#troubleshooting)). |
| `-rerun` | All | Delete `results/`, `resources/` and `workflow/reports/` of the samples from the **previous** run before starting. See [Re-running samples](#re-running-samples-with-new-parameters). |
| `-batch` | All | Split the DAG into batches. Useful for many samples (over ~3000 jobs), as it reduces the DAG resolution time. |
| `-nbatch <int>` | All | Number of batches used with `-batch` (default: 3). |
| `-sifdir <path>` | All | Directory for the Singularity images. Images already in it are not pulled again (default: `resources/sif_dir`). |
| `-nhmmer_db <path>` | All | Nucleotide HMM database used by HMMER (default: `resources/rfam.hmm`). Any database compatible with HMMER 3.4 can be used. |
| `-mitos_db <path>` | All | MITOS2 reference directory, inside the working directory (default: `resources/refseq89m`, for Metazoa). Download another MITOS2 reference to annotate other groups, e.g. `resources/refseq89f` for fungi. |
| `-mem <int>` | Local | Memory (GB) that the NOVOPlasty jobs may use at once. Snakemake then runs only as many of them as fit in this value. |
| `-h`, `--help` | | Show the usage and exit. |
| `-v`, `--version` | | Show the OrganPipe version and exit. |

Example in SLURM mode:

```bash
bash OrganPipe.sh -d . -c config/config.csv -slurm -j 50 -partition <your_partition>
```

### Testing with the sample data

From the repository directory, run:

```bash
bash OrganPipe.sh -d . -t 4 -c test_data/config.csv
```

`test_data/config.csv` tests each feature once, so the run stays short:

| Sample | What it tests |
|:-------|:--------------|
| `hs1` | The whole mitochondrial chain with NOVOPlasty: FASTA seed, `reference`, two k-mers (which give identical assemblies), annotation, images and HMMER. |
| `hs2` | fastp with an adapters file and quality trimming, the NCBI seed search and `annotation: "No"`. |
| `hs3` | A GenBank seed (`feature: "rRNA"`, two seeds) and fastp detecting the adapters. |
| `hs4` | GetOrganelle on a mitochondrion, with the optional fields blank. |
| `SRR1328237` | The whole chloroplast chain with NOVOPlasty and GetOrganelle, and `extra_flags`. |
| `hifi` | Long reads: cutadapt, MitoHiFi and HMMER. |

It needs internet access for the NCBI searches and the image downloads. When it finishes, each test sample has its own folder in `workflow/reports/`. Other test configurations are available in `test_data/` (single sample, two samples, `sample: all` and long reads).

### Automated tests

The tests need only the conda environment, not the containers or the test reads:

```bash
python -m unittest discover -s tests
```

`tests/test_scripts.py` checks the workflow scripts. `tests/test_workflow.py` uses stand-ins for the tools to check that every combination of assembler, organelle and pipeline steps builds the expected jobs. GitHub Actions runs them on every push and pull request.

### Re-running samples with new parameters

Snakemake skips jobs whose outputs already exist. To reprocess a sample with different parameters, use `-rerun`.

> ⚠️ `-rerun` deletes `results/<sample>`, `resources/<sample>` and `workflow/reports/<sample>` for **every sample in `config/snakemake_config.yaml`**, i.e. the samples of the **previous** run, before the new configuration is read. Make sure that is what you want before using it.

## Threads, memory and SLURM resources

> ⚠️ **Memory in local mode**: NOVOPlasty can use up to `max_memory` GB per job, and Snakemake may run several assembly jobs at the same time. Use `-mem <GB>` to limit the memory of these jobs to what your machine has. Without it, keep `-t` × `max_memory` within your RAM, or the system may run out of memory.

As a reference, in a run of 24 amphibian mitochondrial samples with `max_memory: 40`, the NOVOPlasty steps (`create_hash`, `run_novoplasty`) peaked at ~38 GB, while all other steps stayed below 2.5 GB. You can measure your own runs with the [benchmark report](#benchmark-report).

Threads and resources are set in Snakemake profile **templates** in `config/`, with the rules grouped by pipeline step:

| File | Used in | What you can set |
|:-----|:--------|:-----------------|
| `config/threads.yaml` | Local and SLURM | `set-threads`: threads (CPUs) per rule |
| `config/local_params.yaml` | Local | General Snakemake settings |
| `config/slurm_params.yaml` | SLURM | General Snakemake settings and `default-resources`: partition, memory (2 GB) and runtime (72 h) for all rules |
| `config/resources.yaml` | SLURM | `set-resources`: memory, runtime or partition of the rules that need more than the defaults. `create_hash` and `run_novoplasty` get their memory from `max_memory` |

- Rules not listed in `config/threads.yaml` run with 1 thread. In local mode a job never uses more than `-t` cores.
- In SLURM mode, each job is submitted with `--cpus-per-task` = threads, `--mem` = `mem_mb` (MB) and `--time` = `runtime` (minutes).
- `OrganPipe.sh` joins these files into `profiles/local/config.yaml` or `profiles/slurm/config.yaml` on every run, so **always edit the files in `config/`**. Changes made in `profiles/` are overwritten.
- **Do not change the `{WORKDIR}`, `{THREADS}` and `{PARTITION}` placeholders**. `OrganPipe.sh` fills them in from `-d`, `-j` and `-partition`.

To send a specific rule to another partition, or to give it more memory or time, add or edit its line in `config/resources.yaml`:

```yaml
set-resources:
  run_getorganelle: {mem_mb: 100000, runtime: 1440, slurm_partition: "'highmem'"}
```

String values such as partition names must be wrapped in single quotes inside double quotes (`"'highmem'"`), so Snakemake reads them as text.

## Container images

`workflow/containers.yaml` lists the image and version of every tool, in one line per tool:

```yaml
hmmer: staphb/hmmer:3.4
```

Before each run, `OrganPipe.sh` pulls the missing images into `resources/sif_dir` (or the directory given with `-sifdir`) as `<name>_<tag>.sif`, e.g. `hmmer_3.4.sif`.

To update a tool, change its tag. On the next run the new image is pulled, and Snakemake re-runs the finished jobs that used the old image and every job that depends on them. The old `.sif` file is kept; delete it when no run needs it.

Without internet access, `OrganPipe.sh` stops and prints the `singularity pull` commands to run on another machine. Copy the resulting files to the paths shown. A failed pull of `cpgavas2` only stops the jobs that use it (chloroplast annotation), not the whole run.

### NCBI API key

Seed searches on NCBI work without an API key, at up to 3 requests per second. If you have a key, set it to allow up to 10 requests per second:

```bash
export NCBI_API_KEY=<your_key>
```

## Outputs

OrganPipe writes three output trees in the working directory:

| Directory | Content |
|:----------|:--------|
| `workflow/reports/<sample>/` | Curated results and reports: **start here**. |
| `results/<sample>/` | Raw outputs of every tool in the pipeline, per assembly: `results/<sample>/<tool>/<assembler>/<run>/<assembly>/`. A run is `<seed>_kmer<kmer>` for NOVOPlasty and `all` for GetOrganelle. |
| `logs/` and `benchmarks/` | Per-job logs and resource usage (see [Logs and benchmarks](#logs-and-benchmarks)). |

### Report folders

Each `workflow/reports/<sample>/` folder contains the subfolders below. Only the ones that apply to your data and configuration are created.

| Folder | Description |
|:-------|:------------|
| `fastas/` | All assembly FASTA files. |
| `genbanks/` | GenBank files of the annotated assemblies. |
| `genes/` | Genes extracted from each assembly (per seed/k-mer). |
| `files/` | Files to help with manual curation. |
| `novoplasty/` | NOVOPlasty results (short reads). |
| `pilon/` | Pilon error-correction results (short reads). |
| `mitos2/` | MITOS2 annotation results (mitochondria). |
| `cpgavas2/` | CPGAVAS2 annotation results (chloroplast). |
| `mitohifi/` | MitoHiFi results (long reads). |
| `nhmmer/` | HMMER results (`run_nhmmer: "Yes"`). |
| `images/` | OGDRAW maps and read recruitment/depth plots (`run_images: "Yes"`). |

### Report files

| File | Description |
|:-----|:------------|
| `summary.csv` | Summary of all assemblies and their results. |
| `warnings.txt` | Problems that did not stop the run, e.g. a seed without an assembly in `genome_range`, an assembly that could not be rotated or a plot that could not be drawn. Only present if there is a warning. |
| `novoplasty.csv` | NOVOPlasty assembly results (short reads). |
| `pilon.csv` | Pilon error-correction results (short reads). |
| `mitos2.csv` | MITOS2 annotation metrics (mitochondria). |
| `mitohifi.csv` | MitoHiFi results (long reads). |
| `nhmmer_intergenes.csv` | HMMER intergenic region annotations. |
| `nhmmer_ncRNA.csv` | HMMER ncRNA annotations. |
| `cpgavas2_codon_usage.csv` | Codon usage statistics (chloroplast). |
| `cpgavas2_gene_composition.csv` | Gene composition summary (chloroplast). |
| `cpgavas2_intron_exon.csv` | Intron–exon structure (chloroplast). |
| `cpgavas2_problems.csv` | Annotation issues reported by CPGAVAS2 (chloroplast). |

### Identical assemblies

Different seeds and k-mers often produce the same genome. The `identical_group` column of `summary.csv`, `pilon.csv` and `mitohifi.csv` marks these assemblies: rows with the same label (`group_1`, `group_2`...) have the same sequence, regardless of the start position and the strand. The column is blank for assemblies that are different from all the others.

### Rotated GenBank files

`<assembly>.rotated.gb` starts at `rrnS`, on the forward strand (mitochondria), and `<assembly>.cpgavas2.rotated.gb` and `<assembly>.chloe.rotated.gb` start at `psbA`, on the reverse strand (chloroplast). All the assemblies of a sample then have the same start and orientation, which makes them easier to compare. An assembly without the gene has no rotated file, and it is listed in `warnings.txt`.

### Assembly names

**NOVOPlasty assemblies** are named:

```
{status}_{n}_{sample}_{seed}_{kmer}_{index}
```

- **status**: `CA` (circularized assembly, the genome was closed), `OPT` (option: several possible circularizations were found) or `CTG` (contig: the genome could not be closed).
- **n**: assembly number reported by NOVOPlasty.
- **sample**: sample name.
- **seed**: seed identifier. Seeds downloaded from NCBI are named `{number}-{gene}_{TaxID}`.
- **kmer**: k-mer size.
- **index**: sequence number within the assembly.

For example, `CA_1_ITV00872_25_1-COI_941666_19_1` is the circularized assembly 1 of sample `ITV00872_25`, built from seed `1-COI_941666` (a COI sequence from TaxID 941666) with k-mer 19. It is the first sequence of that assembly.

**GetOrganelle assemblies** are named:

```
{database}_{status}_{number}
```

- **database**: GetOrganelle database set in the configuration (e.g. `animal_mt`, `embplant_pt`).
- **status**: `complete` for a closed circular genome, otherwise the GetOrganelle status of the sequence (e.g. scaffolds).
- **number**: sequence number, starting at 1.

For example, `animal_mt_complete_1` is the first complete genome assembled with the `animal_mt` database.

## Logs and benchmarks

| Path | Content |
|:-----|:--------|
| `logs/<sample>/<step>/` | Log of every job, e.g. `logs/<sample>/pilon/<assembler>/<run>/<assembly>_run_pilon.log` |
| `logs/<sample>.parse_results.log` | Log of the report generation for each sample |
| `logs/benchmark_report.log` | Log of the benchmark report |
| `benchmarks/<sample>/<step>/` | Runtime, CPU and memory usage of every job |
| `.snakemake/log/` | Full Snakemake log of each run |

When a job fails, Snakemake prints its log (`show-failed-logs`) and keeps running the jobs that do not depend on it (`keep-going`), so other samples still finish.

### Benchmark report

At the end of every run (but not on dry runs), OrganPipe writes an HTML report of runtime, CPU and memory usage to `workflow/reports/benchmark.html`. The report is also written when some jobs failed, so it always covers the benchmarks written so far. It breaks the benchmarks down by stage, rule, sample and k-mer, lists the longest jobs, and flags samples with missing jobs. If the report cannot be generated, the run is not affected and the reason is written to `logs/benchmark_report.log`.

To build the report manually, e.g. for another configuration, run from the OrganPipe working directory:

```bash
python workflow/scripts/benchmark_report.py -c config/config.csv
```

If `-c` is omitted, `config/snakemake_config.yaml` from the last run is used when present. Samples in `benchmarks/` that are not in the config are left out unless `--all-samples` is set. See `--help` for the other options.

## Troubleshooting

**The run stopped or was interrupted. How do I resume it?**
Run the same command again. Snakemake skips the jobs that already finished and re-runs the incomplete ones. Each assembly has its own jobs, so a failed job only repeats that assembly.

**"Directory cannot be locked" / "LockException".**
A previous run was interrupted. Run the same command with `-unlock`, then run it again without `-unlock`.

**NCBI search error: no seed found.**
No NCBI record matched `search_term` and `search_genes`. Use a broader taxon (genus or family), other genes, or provide your own `seed_file`.

**A job failed with an out-of-memory error.**
- Local mode: use `-mem`, or lower `-t` or `max_memory` so that `-t` × `max_memory` fits in your RAM.
- SLURM mode: increase `mem_mb` for that rule in `config/resources.yaml`. For NOVOPlasty, increase `max_memory`.

**A SLURM job was killed for exceeding its time limit.**
Increase `runtime` (minutes) for that rule in `config/resources.yaml`.

**A Singularity image failed to download or build.**
Run the pipeline again; images already in the `-sifdir` directory are not pulled again. Check your internet connection and the free space in that directory. Without internet access, `OrganPipe.sh` prints the `singularity pull` commands to run on another machine (see [Container images](#container-images)).

**A pull failed with `unauthorized: incorrect username or password`.**
Singularity is using an expired Docker Hub login stored in `~/.docker/config.json`. Run `docker login` again, or remove that login with `docker logout`. The images are public and need no login.

**I changed the parameters, but the sample was not re-processed.**
Snakemake does not re-run jobs whose outputs already exist. Use `-rerun` (read the warning in [Re-running samples](#re-running-samples-with-new-parameters)), or delete `results/<sample>`, `resources/<sample>` and `workflow/reports/<sample>` manually.

**A sample finished without some results.**
Check `workflow/reports/<sample>/warnings.txt`, which lists the problems that did not stop the run.

**The DAG takes too long to build.**
With many samples, use `-batch` (and `-nbatch` to set the number of batches).

## How to cite

If you use OrganPipe, please cite the preprint:

```
Moreira-Oliveira, R. R., Silva, B. M., Molina, M., Oliveira-Lima, M., Leão, T. F., Vasconcelos, S., & Nunes, G. L. (2025). OrganPipe: An automated tool to facilitate the assembly, annotation, and curation of mitochondrial and chloroplast genomes. https://doi.org/10.21203/RS.3.RS-5686696/V1
```

```bibtex
@misc{moreira-oliveira2025organpipe,
  title     = {OrganPipe: An automated tool to facilitate the assembly, annotation, and curation of mitochondrial and chloroplast genomes},
  author    = {Moreira-Oliveira, R. R. and Silva, B. M. and Molina, M. and Oliveira-Lima, M. and Le{\~a}o, T. F. and Vasconcelos, S. and Nunes, G. L.},
  year      = {2025},
  publisher = {Research Square},
  doi       = {10.21203/RS.3.RS-5686696/V1}
}
```

We encourage users to also cite the tools used in the pipeline:

- **NOVOPlasty**: Dierckxsens, N., Mardulyn, P., & Smits, G. (2017). NOVOPlasty: de novo assembly of organelle genomes from whole genome data. *Nucleic Acids Research*, 45(4), e18. https://doi.org/10.1093/nar/gkw955
- **GetOrganelle**: Jin, J.-J., Yu, W.-B., Yang, J.-B., Song, Y., dePamphilis, C. W., Yi, T.-S., & Li, D.-Z. (2020). GetOrganelle: a fast and versatile toolkit for accurate de novo assembly of organelle genomes. *Genome Biology*, 21, 241. https://doi.org/10.1186/s13059-020-02154-5
- **MitoHiFi**: Uliano-Silva, M., Ferreira, J. G. R. N., Krasheninnikova, K., et al. (2023). MitoHiFi: a python pipeline for mitochondrial genome assembly from PacBio high fidelity reads. *BMC Bioinformatics*, 24, 288. https://doi.org/10.1186/s12859-023-05385-y
- **MITOS2**: Donath, A., Jühling, F., Al-Arab, M., Bernhart, S. H., Reinhardt, F., Stadler, P. F., Middendorf, M., & Bernt, M. (2019). Improved annotation of protein-coding genes boundaries in metazoan mitochondrial genomes. *Nucleic Acids Research*, 47(20), 10543–10552. https://doi.org/10.1093/nar/gkz833
- **Pilon**: Walker, B. J., Abeel, T., Shea, T., et al. (2014). Pilon: an integrated tool for comprehensive microbial variant detection and genome assembly improvement. *PLoS One*, 9(11), e112963. https://doi.org/10.1371/journal.pone.0112963
- **CPGAVAS2**: Shi, L., Chen, H., Jiang, M., Wang, L., Wu, X., Huang, L., & Liu, C. (2019). CPGAVAS2, an integrated plastome sequence annotator and analyzer. *Nucleic Acids Research*, 47(W1), W65–W73. https://doi.org/10.1093/nar/gkz345
- **Chloë**: https://github.com/ian-small/Chloe.jl
- **OGDRAW**: Greiner, S., Lehwark, P., & Bock, R. (2019). OrganellarGenomeDRAW (OGDRAW) version 1.3.1: expanded toolkit for the graphical visualization of organellar genomes. *Nucleic Acids Research*, 47(W1), W59–W64. https://doi.org/10.1093/nar/gkz238

## License

OrganPipe is distributed under the GNU General Public License v3.0. See [LICENSE](LICENSE).
