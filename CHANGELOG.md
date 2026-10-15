# Changelog

## OrganPipe 1.3.0

### New Features
- Added a benchmark report: an HTML report of runtime, CPU and memory usage (`workflow/reports/benchmark.html`) built from the Snakemake benchmark files at the end of every run, including runs where some jobs failed. It reads the rule definitions, so it works with any combination of samples, assemblers and steps. It can also be built manually with `workflow/scripts/benchmark_report.py`.
- Added `-h`/`--help` and `-v`/`--version` to `OrganPipe.sh`.

### Improvements
- Reorganized `config/local_params.yaml` and `config/slurm_params.yaml`, with comments explaining each setting:
  - SLURM partition, memory and runtime are now set once in `default-resources`; `set-resources` lists only the rules that need more memory.
  - `set-threads` lists only the rules that use more than 1 thread.
  - `run_cutadapt` now uses 4 threads, like `run_fastp`.
  - `max-jobs-per-second` is fixed at 10 instead of the value of `-j`, to avoid flooding the scheduler.
- Rewrote the README: quick start, table of contents, configuration fields grouped by topic, resource and output sections, troubleshooting, and the changelog moved to this file.

### Bug Fixes
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
