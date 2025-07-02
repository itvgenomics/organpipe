# OrganPipe: A pipeline for the assembly, annotation, and curation of mitochondrial and chloroplast genomes

## How to cite OrganPipe (preprint)

 ```
Moreira-Oliveira, R. R., Silva, B. M., Molina, M., Oliveira-Lima, M., Leão, T. F., Vasconcelos, S., & Nunes, G. L. (2025). OrganPipe: An automated tool to facilitate the assembly, annotation, and curation of mitochondrial and chloroplast genomes. https://doi.org/10.21203/RS.3.RS-5686696/V1
 ```
[https://doi.org/10.21203/RS.3.RS-5686696/V1](https://doi.org/10.21203/RS.3.RS-5686696/V1)


## How to Install OrganPipe Environment

### Prerequisites
Before installing the required software, make sure you have the following:
- A Linux-based operating system (e.g., Ubuntu, CentOS, Fedora)
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

1. **Configure the Pipeline**:
    - Before running the pipeline, you must edit the configuration file to provide information specific to your samples and variables. We provide this file in two formats: YAML and CSV at the `config` directory. **If you don't need to use a specific variable, leave the value blank rather than removing the key/column. (e. g. reference: '')**

   - **YAML Format**: Each key in the file represents a variable, and each value corresponds to the parameters or settings for that variable.

   - **CSV Format**: Each column corresponds to a variable, and each row represents a sample with specific values for the variables.

The fields to be edited are the following:

| **Field**           | **Example**          | **Description/Comment**                                                                                                                                                              | **Field Requirement**                         |
|----------------------|----------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------------------------------|
| `sample`            | `sample1`           | Specify the sample name to process. Use "all" to process all samples in `reads_path`. You can also set a comma separated list of samples (e.g. "sample1,sample2,sample3"). All samples will be running with the same parameters, if you do want different values for each sample, use the .csv format | Required.                                     |
| `reads_path`        | `/path/to/reads`    | Directory containing sequencing read files (_R1/_R2 or _pair1/_pair2) in fastq.gz format.                                                                                          | Required.                                     |
| `organelle`         | `mito`              | Type of organelle genome to assemble: "mito" for mitochondrial or "chloro" for chloroplast.                                                                                         | Required.                                     |
| `genetic_code`      | `2`                 | Genetic code for mitochondrial genome annotation (e.g., 2 for The Vertebrate Mitochondrial Code).                                                                                                       | Required.                                     |
| `reference`         | `/path/to/ref.fasta`| Path to a reference genome in fasta format. Leave blank if unavailable. This is used by NOVOPlasty as a guide to resolve duplicated regions in the plant mitochondria or the inverted repeat in the chloroplast                                                                                                            | Optional.                                     |
| `sequencing_type`   | `Short`             | Type of sequencing data: "Short" for Illumina or "Long" for PacBio/ONT reads.                                                                                                       | Required.                                     |
| `genome_range`      | `12000-22000`       | Expected genome size range in kilobases. Every sequence between this range will be annotated if `annotation` is set to "yes"                                                                                                                                             | Required.                                     |
| `adapterremoval`    | `Yes`               | Run AdapterRemoval for adapter trimming: "Yes" to perform trimming, "No" to skip.                                                                                                   | Optional, default is "No".                   |
| `adapters`          | `/path/to/adapters.txt` | Path to a text file listing adapter sequences to remove.                                                                                                                             | Required if `adapterremoval` is "Yes".        |
| `minlength`         | `50`                | Minimum read length to retain after trimming.                                                                                                                                       | Required if `adapterremoval` is "Yes".        |
| `minquality`        | `20`                | Minimum base quality score threshold for trimming.                                                                                                                                 | Required if `adapterremoval` is "Yes".        |
| `seed_format`       | `fasta`             | Format of the seed file: "fasta" or "genbank".                                                                                                                                       | Optional.                                     |
| `seed_file`         | `/path/to/seed.fasta`| Path to the seed file to initialize assembly.                                                                                                                                         | Optional.                                     |
| `feature`           | `CDS`               | Feature type for assembly if GenBank seed file is used: "CDS", "rRNA", or "tRNA".                                                                                                   | Required if `seed_format` is "genbank".       |
| `search_ncbi`       | `Yes`               | Search NCBI for reference sequences: "Yes" to search, "No" to skip.                                                                                                                 | Optional.                                     |
| `search_genes`      | `COI,16S,ATP6`      | Genes to search for on NCBI if `search_ncbi` is enabled.                                                                                                                            | Required if `search_ncbi` is "Yes".           |
| `search_term`       | `Amphisbaena`       | If no matching record is found, an error will be raised.                                                                                                                                             | Required if `search_ncbi` is "Yes".           |
| `max_references`    | `5`                | Maximum reference sequences to download per gene.                                                                                                                                  | Required if `search_ncbi` is "Yes".                                     |
| `kmers`             | `19,23,33,39`   | List of k-mer sizes for genome assembly.                                                                                                                                             | Required for Short Reads assembly.               |
| `max_memory`        | `4`                | Maximum memory (in GB) for genome assembly.                                                                                                                                         | Required for Short Reads assembly.                  |
| `reads_length`      | `150`               | Read length (in base pairs) of sequencing data.                                                                                                                                      | Required for Short Reads.                  |
| `insert_size`       | `300`               | Average insert size (in base pairs) of sequencing data.                                                                                                                             | Required for Short Reads.                  |
| `annotation`        | `Yes`               | Run annotation pipeline: "Yes" to annotate, "No" to skip.                                                                                                                           | Optional. Default is "No".                    |
| `run_nhmmer`        | `No`                | Run nhmmer to identify ncRNA and intergenic regions. **Note:** Enabling this can  slow down the pipeline.                                                                                                                               | Optional. Default is "No".                    |
| `nhmmer_db`         | `resources/rfam.hmm`| Path to the HMM database for nhmmer. We are providing a simple database in the `resources` directory, but you can use any HMMER database compatible with HMMER version 3.4 for improved accuracy.                                                                                                                                                | Required if `run_nhmmer` is "Yes".            |
| `run_images`        | `Yes`               | Generate visualizations like OGDraw diagrams and depth plots. **Note:** Enabling this can slow down the pipeline.                                                                                                                      | Optional. Default is "Yes".                   |

2. **Run OrganPipe**:
    - Run the pipeline with default parameters:

    ```
    bash OrganPipe.sh -d </path/to/work/dir> -t <n_threads> -c </path/to/configfile>
    ```

    - **Flags**:
        - **-d** </path/to/work/dir> (Required) = Path to your working directory where all the workflow file are
        - **-c** </path/to/config.yaml> (Required) = Overwrite the default configuration file with all needed parameters (e.g. config/config.yaml/csv)
        - **-t** {int} (Required) = Number of threads to use
        - **-np** (Optional) = Perform a dry run to see what jobs will be executed without actually running them.
        - **-unlock** (Optional) = Unlock the working directory if Snakemake has somehow locked it.
        - **-batch** (Optional) = If you are running a large number of samples, or number of rules executed > 3000, consider using this flag. This slightly improves the DAG resolution time from Snakemake. You can set the number with -nbatch (Default = 15)
        - **-sifdir** (Optional) = Choose a directory to build all singularity image files used in the pipeline. If the path already contains the images, they will not be pulled. Default: resources/sif_dir


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
