# Changelog

## OrganPipe 1.3.0

### New Features
- Added a benchmark report: an HTML report of runtime, CPU and memory usage (`workflow/reports/benchmark.html`) built from the Snakemake benchmark files at the end of every run, including runs where some jobs failed. It reads the rule definitions, so it works with any combination of samples, assemblers and steps. It can also be built manually with `workflow/scripts/benchmark_report.py`.
- Added `-h`/`--help` and `-v`/`--version` to `OrganPipe.sh`.
- Added `workflow/containers.yaml`, which sets the image and version of every tool in one place:
  - Images are now stored as `<name>_<tag>.sif`, so changing a tag pulls the new version and keeps the old file. Images from earlier versions are pulled once more under the new names; the old `<name>.sif` files can be deleted.
  - Snakemake re-runs the finished jobs of a tool when its image changes. This includes the first run after this update in a working directory with finished jobs.
  - Interrupted pulls no longer leave an incomplete `.sif` file. When a pull fails, `OrganPipe.sh` prints the `singularity pull` commands to run on another machine.
- Added the `identical_group` column to `summary.csv`, `pilon.csv` and `mitohifi.csv`, which marks assemblies with the same sequence regardless of their start position and strand.
- Added `workflow/reports/<sample>/warnings.txt`, which lists the problems that did not stop the run.
- Added `-mitos_db` to choose the MITOS2 reference, since `resources/refseq89m` is for Metazoa.
- Added `-mem` to limit the memory of the NOVOPlasty jobs in local mode.
- Added automated tests (`tests/`) and a GitHub Actions workflow.
- Rewrote `test_data/config.csv` to test each feature once, with fewer seeds and k-mers. It also covers what was missing: an adapters file, `reference`, blank GetOrganelle fields, `extra_flags` and HMMER for long reads.

### Improvements
- Reorganized the workflow around one set of rules for all assemblers, with one set of jobs per assembly:
  - A failed job now repeats only that assembly, assemblies run in parallel, and the rules reduced from 86 to 45.
  - Raw results moved to `results/<sample>/<tool>/<assembler>/<run>/<assembly>/`. The structure of `workflow/reports/` did not change. Working directories from earlier versions need `-rerun`.
  - A dry run now lists the jobs only up to the assembly step until the assemblies exist.
  - Rule names changed, e.g. `run_novoplasty_pilon` and `run_getorganelle_pilon` are now `run_pilon`. Update your custom settings in `config/threads.yaml` and `config/resources.yaml`.
- The memory of `create_hash` and `run_novoplasty` now follows `max_memory` (plus 10%), instead of a fixed 50 GB.
- `run_blastn` now uses 4 threads.
- Rotated GenBank files now have the same orientation in all the assemblies of a sample: `rrnS` on the forward strand (mitochondria) and `psbA` on the reverse strand (chloroplast).
- GenBank files of mitochondria now include every coding gene annotated by MITOS2, not only the 13 usual ones.
- The ncRNA and intergenic search now uses `nhmmscan`, the HMMER tool for nucleotide databases, instead of `hmmscan`. The columns of `nhmmer_ncRNA.csv` and `nhmmer_intergenes.csv` changed accordingly.
- The NCBI API key is now read from the `NCBI_API_KEY` environment variable, and is optional.
- `OrganPipe.sh` can now be called from any directory.
- `sample: "all"` now also works for long reads and ignores `Undetermined` files in `.fastq.gz`, as documented. The run stops if it finds no samples.
- Completed `.gitignore`, so a run leaves the repository clean. Removed the compiled Python files from the repository and `workflow/scripts/create_hmmdb.py`, which no step used.
- Reorganized `config/local_params.yaml` and `config/slurm_params.yaml`, with comments explaining each setting:
  - Threads moved to `config/threads.yaml`, shared by the local and SLURM modes. It lists only the rules that use more than 1 thread.
  - Memory per rule moved to `config/resources.yaml`. It lists only the rules that need more than the defaults.
  - SLURM partition, memory and runtime are now set once in `default-resources`. The default memory is now 2 GB (it was 1 GB).
  - `run_cutadapt` now uses 4 threads, like `run_fastp`.
  - `max-jobs-per-second` is fixed at 10 instead of the value of `-j`, to avoid flooding the scheduler.
- Rewrote the README: quick start, table of contents, configuration fields grouped by topic, resource and output sections, troubleshooting, and the changelog moved to this file.

### Bug Fixes
- `minquality` had no effect: fastp only uses `--cut_mean_quality` together with a cut option. Low-quality bases are now cut from both read ends (`--cut_front --cut_tail`), so results of samples with `run_trimming` enabled can change.
- A sample could get the reads of another one, when its name was part of the other's (e.g. `hs1` got `hs10_R2.fastq.gz`). A sample that still matches several files now stops the run with an error.
- `sample: "small_x"` processed every sample in `reads_path`, because the name contains `all`.
- Bases other than A, C, G and T were deleted from NOVOPlasty assemblies, which shifted the coordinates. They are now replaced with `N` and reported in `warnings.txt`.
- The ncRNA and intergenic search found nothing for chloroplast assemblies of NOVOPlasty, because it looked for the old assembly names.
- `run_nhmmer: "Yes"` failed for long reads.
- `sequencing_type: "Short"` and other values with capital letters failed in the ncRNA and intergenic search. `organelle` and `sequencing_type` are now case-insensitive in every step.
- CPGAVAS2 could no longer be downloaded, because the Sylabs library (`library.sylabs.io`) is offline. The image now comes from Docker Hub (`brunomsilva/cpgavas2:latest`).
- Chloroplast samples failed on a new installation, because `resources/table2asn.linux64` was never extracted. Extracting it no longer deletes the `.gz` file.
- CPGAVAS2 and Chloë could overwrite each other's temporary files when running at the same time.
- `adapters` failed with a cyclic dependency error whenever it was set.
- Long reads in `.fastq.gz` without trimming were given to MitoHiFi as a compressed file named `.fasta`.
- Genes with the same product (e.g. two `tRNA-Leu`) were lost in the ncRNA search.
- A MITOS2 run without annotation stopped the report of the sample. It is now listed in `warnings.txt`.
- A NOVOPlasty run could mix circularized assemblies with options, depending on the order of the files in its folder.
- Messages of `seqtk` could end up in the FASTA file given to MitoHiFi, and repeating a `run_blastn` job duplicated its reads.
- `kmers` with a single value, a list in `target_size` and blank GetOrganelle fields no longer fail.
- `cpgavas2_codon_usage.csv`, `cpgavas2_gene_composition.csv`, `cpgavas2_intron_exon.csv` and `cpgavas2_problems.csv` were never written.
- A long-read sample without references on NCBI now stops with an error, instead of finishing without results.
- Seeds that were counted without being written no longer reduce the number of seeds downloaded from NCBI.
- `check_config.py` now stops the run for an invalid YAML file, invalid `pacbio_adapters` or a long-read `.fasta` file with trimming, and checks `max_memory` and the sample names.
- Logs are no longer overwritten when a step processes several files, and rules no longer share log files.
- `OrganPipe.sh` now prints a clear error for unknown options, options missing their value and missing required options (`-c`, `-d`, `-t`), instead of exiting silently. The SLURM option checks run before the banner is printed.
- Fixed the benchmark path of `run_getorganelle_pilon`, which was written to a file literally named `{run_pilon.benchmark`.
- Fixed rule names in the profiles that did not match any rule (`check_longreads_files`, `get_chloro_ncRNA_sequences_getorganelle`), so their settings were ignored.
- `run_cutadapt`, `run_seqtk`, `check_fastq_longreads` and `get_chloro_ncRNA_sequences_novoplasty_getorganelle` were submitted to SLURM without a partition or runtime.
- The SLURM partition name is now quoted, so names that are also Python built-ins (e.g. `all`) are passed correctly.
- Removed the unused `node` resource from the SLURM profile.

## OrganPipe 1.2.0

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

## OrganPipe 1.1

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
- Added `-rerun` flag:
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
  - Default SLURM parameters can be modified in `config/slurm_params.yaml`.
  - Users must set their partition using the `slurm_partition` variable.

### Output & Reporting Changes
- Replaced **CIRCOS** with **OGDRAW** for generating circular mitogenome/plastome images.
- Renamed `abstract.csv` to `summary.csv`.
- Chloroplast short-read assemblies now include a summary.
- Reports generation was reworked — see the *Outputs* section in the README for details.
