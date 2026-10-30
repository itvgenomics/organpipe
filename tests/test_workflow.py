"""Checks that the workflow builds a complete chain of jobs for each kind of run.

No container or real tool is used. The tools that produce the assemblies are simple
stand-ins, which lets Snakemake get past the assembly step and list every later job.
Run with: python -m unittest discover tests
"""

import csv
import os
import re
import shutil
import stat
import subprocess
import tempfile
import unittest
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent

GENOME = "ACGT" * 4000

TOOLS = {
    # Runs the job command directly: singularity [options] exec [options] image bash -c 'command'
    "singularity": """
case "$*" in
    *--version*|version) echo "singularity version 3.8.5"; exit 0 ;;
    pull*) touch "${@: -2:1}"; exit 0 ;;
esac
# The image must be the one the pull step wrote
image="${@: -4:1}"
[ -f "$image" ] || { echo "image not found: $image" >&2; exit 1; }
"${@: -3}"
""",
    "NOVOPlasty.pl": f"""
if [ -n "${{MOCK_EVENTS:-}}" ]; then
    echo "start" >> "$MOCK_EVENTS"
    sleep 4
    echo "end" >> "$MOCK_EVENTS"
fi
config="$2"
value() {{ sed -n "s/^ *$1 *= *//p" "$config"; }}
project=$(value "Project name"); out=$(value "Output path")
if [ "$(value "Store Hash")" = "yes" ]; then
    touch "$out/HASH2B_$project.txt" "$out/HASH2C_$project.txt" "$out/HASH_$project.txt"
else
    echo "Total contigs : 1" > "$out/log_$project.txt"
    if [ -n "${{MOCK_CONTIGS:-}}" ]; then
        # A sample that does not circularize: two contigs in the genome range and a short one
        printf '>Contig1\\n{GENOME}\\n>Contig2\\n{GENOME}\\n>Contig3\\nACGTACGT\\n' > "$out/Contigs_1_$project.fasta"
    elif [ -z "${{MOCK_NO_ASSEMBLY:-}}" ]; then
        printf '>Contig1\\n{GENOME}\\n' > "$out/Circularized_assembly_1_$project.fasta"
    fi
    printf '>r\\nACGT\\n' | tee "$out/Assembled_reads_${{project}}_R1.fasta" > "$out/Assembled_reads_${{project}}_R2.fasta"
fi
""",
    "get_organelle_config.py": "exit 0",
    "get_organelle_from_reads.py": f"""
while [ $# -gt 0 ]; do
    if [ "$1" = "-o" ]; then out="$2"; fi
    case "$1" in
        -R|-k|--target-genome-size) [[ "${{2:-}}" =~ ^[0-9] ]] || {{ echo "$1 without a value" >&2; exit 2; }} ;;
    esac
    shift
done
echo done > "$out/get_org.log.txt"
printf '>c\\n{GENOME}\\n' > "$out/animal_mt.K115.complete.graph1.1.path_sequence.fasta"
touch "$out/extended_1_paired.fq" "$out/extended_2_paired.fq"
""",
    "fastp": """
while [ $# -gt 0 ]; do
    case "$1" in --out1|--out2|--html|--json) touch "$2" ;; esac
    shift
done
""",
    "cutadapt": """
while [ $# -gt 0 ]; do
    if [ "$1" = "-o" ]; then touch "$2"; fi
    shift
done
""",
    "seqtk": "echo '>r'; echo ACGT",
    "bwa-mem2.avx": """
if [ "$1" = "index" ]; then
    for extension in 0123 amb ann bwt.2bit.64 pac; do touch "$2.$extension"; done
fi
""",
    "samtools": """
case "$1" in
    view) cat ;;
    sort) while [ $# -gt 0 ]; do if [ "$1" = "-o" ]; then cat > "$2"; fi; shift; done ;;
    index) touch "$2.bai" ;;
esac
""",
    "mitohifi.py": "touch contigs_stats.tsv",
    "hmmpress": 'touch "$1.h3f" "$1.h3i" "$1.h3m" "$1.h3p"',
}

# Replaces resources/pilon.sh, which needs the Pilon container.
PILON = """#!/bin/bash
while [ $# -gt 0 ]; do
    case "$1" in
        --genome) genome="$2" ;;
        --output) output="$2" ;;
    esac
    shift
done
sed '/^>/s/$/_pilon/' "$genome" > "$output.fasta"
touch "$output.changes"
"""

SHORT = "test_data/config_1sample.yaml"
LONG = "test_data/config_longreads.yaml"

NOVOPLASTY_ONLY = {"run_novoplasty": "Yes", "run_getorganelle": "No"}
GETORGANELLE_ONLY = {
    "run_novoplasty": "No",
    "run_getorganelle": "Yes",
    "database": "animal_mt",
    "seed_format": None,
    "seed_file": None,
    "reference": None,
}
BOTH_ASSEMBLERS = {"run_novoplasty": "Yes", "run_getorganelle": "Yes", "database": "embplant_pt"}
CHLORO = {"organelle": "chloro", "genetic_code": 11, "genome_range": "10000-20000"}
TRIMMING = {"run_trimming": "Yes", "minlength": 50, "minquality": 20}
ALL_STEPS = {"annotation": "Yes", "run_images": "Yes", "run_nhmmer": "Yes"}
NO_STEPS = {"annotation": "No", "run_images": "No", "run_nhmmer": "No"}

MITO_RULES = {"run_bwa_index", "run_bwa_mem", "run_samtools_index", "run_pilon", "run_mitos2", "write_gb_mitos2", "rotate_genbank"}
CHLORO_RULES = {"run_pilon", "run_cpgavas2", "run_chloe", "write_gb_cpgavas2", "write_gb_chloe", "rotate_genbank", "extract_table2asn"}
IMAGE_RULES = {
    "get_genbank_fastas",
    "run_blastn",
    "run_recruitment_plot",
    "run_bwa_index_rotated",
    "run_bwa_mem_rotated",
    "run_samtools_depth",
    "run_samtools_depth_rotated",
    "run_depth_plot",
    "run_ogdraw",
}
NHMMER_RULES = {"get_ncRNA_sequences", "get_intergenes_sequences", "run_nhmmer", "index_hmmer_db"}
ANNOTATION_RULES = {"run_mitos2", "run_cpgavas2", "run_chloe", "rotate_genbank"}

# name: (base config, fields to change, rules that must run, rules that must not run)
SCENARIOS = {
    "mito_novoplasty_all_steps": (
        SHORT,
        {**NOVOPLASTY_ONLY, **ALL_STEPS},
        MITO_RULES | IMAGE_RULES | NHMMER_RULES,
        {"run_cpgavas2", "run_chloe", "run_fastp", "run_getorganelle"},
    ),
    "mito_novoplasty_without_annotation": (
        SHORT,
        {**NOVOPLASTY_ONLY, **NO_STEPS},
        {"run_pilon"},
        ANNOTATION_RULES | IMAGE_RULES | NHMMER_RULES,
    ),
    "mito_getorganelle_blank_options": (
        SHORT,
        {**GETORGANELLE_ONLY, **ALL_STEPS, "n_rounds": None, "target_size": None, "spades_kmers": None},
        MITO_RULES | IMAGE_RULES | NHMMER_RULES,
        {"run_novoplasty", "create_hash"},
    ),
    "mito_trimming_with_adapters": (
        SHORT,
        {**NOVOPLASTY_ONLY, **TRIMMING, "adapters": "test_data/adapters.fasta"},
        {"run_fastp", "check_adapters", "run_pilon"},
        set(),
    ),
    "chloro_both_assemblers_all_steps": (
        SHORT,
        {**BOTH_ASSEMBLERS, **CHLORO, **TRIMMING, **ALL_STEPS, "n_rounds": 15, "target_size": 160000, "spades_kmers": "21,55"},
        CHLORO_RULES | IMAGE_RULES | NHMMER_RULES | {"run_fastp", "run_getorganelle", "run_novoplasty"},
        {"run_mitos2", "write_gb_mitos2"},
    ),
    "long_reads_trimmed_with_nhmmer": (
        LONG,
        {"run_trimming": "Yes", "run_nhmmer": "Yes"},
        {"run_cutadapt", "run_seqtk", "run_mitohifi", "run_nhmmer", "get_ncRNA_sequences", "get_intergenes_sequences"},
        {"run_pilon", "check_fasta_longreads"},
    ),
    "long_reads_fastq_without_trimming": (
        LONG,
        {"run_trimming": "No", "run_nhmmer": "No"},
        {"check_fastq_longreads", "run_seqtk", "run_mitohifi"},
        {"run_cutadapt", "check_fasta_longreads", "run_nhmmer"},
    ),
    "long_reads_fasta_without_trimming": (
        LONG,
        {"run_trimming": "No", "run_nhmmer": "No", "reads_path": "fasta_reads"},
        {"check_fasta_longreads", "run_mitohifi"},
        {"run_cutadapt", "check_fastq_longreads", "run_seqtk"},
    ),
}


def job_counts(log):
    """Return {rule: count} from the last job table of a Snakemake log."""
    tables = re.findall(r"^job +count\n-+ +-+\n(.*?)^total", log, flags=re.MULTILINE | re.DOTALL)
    if not tables:
        return {}
    return {line.split()[0]: int(line.split()[1]) for line in tables[-1].strip().splitlines()}


@unittest.skipUnless(shutil.which("snakemake"), "snakemake is not installed")
class WorkflowScenarios(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bin = tempfile.TemporaryDirectory()
        for name, body in TOOLS.items():
            path = Path(cls.bin.name) / name
            path.write_text("#!/bin/bash\n" + body + "\n")
            path.chmod(path.stat().st_mode | stat.S_IEXEC)
        cls.env = {**os.environ, "PATH": f"{cls.bin.name}:{os.environ['PATH']}"}

    @classmethod
    def tearDownClass(cls):
        cls.bin.cleanup()

    def prepare(self, workdir, base, fields):
        # Only the files of the project: what a real run left in the repository can be very large
        listing = subprocess.run(
            ["git", "ls-files", "--cached", "--others", "--exclude-standard"],
            cwd=ROOT, capture_output=True, text=True, check=True,
        ).stdout
        # Large inputs that the jobs only read are linked, so each scenario takes little disk space
        linked = ["resources/refseq89m", "resources/table2asn.linux64.gz", "test_data/rawreads"]
        for path in listing.splitlines():
            if (ROOT / path).exists() and not any(path == link or path.startswith(f"{link}/") for link in linked):
                (workdir / path).parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(ROOT / path, workdir / path)
        for path in linked:
            (workdir / path).parent.mkdir(parents=True, exist_ok=True)
            (workdir / path).symlink_to(ROOT / path)

        data = yaml.safe_load((workdir / base).read_text())
        data.update(fields)
        config = workdir / "scenario.yaml"
        config.write_text(yaml.safe_dump(data))

        if data["sequencing_type"].lower() == "long":
            seeds = workdir / "resources" / "hifi" / "seeds"
            seeds.mkdir(parents=True)
            (seeds / "REF.fasta").write_text(">REF\nACGT\n")
            (seeds / "REF.gb").touch()
            fasta_reads = workdir / "fasta_reads"
            fasta_reads.mkdir()
            (fasta_reads / "hifi.fasta").write_text(">r\nACGT\n")

        return config

    def run_tool(self, workdir, *command, **env):
        return subprocess.run(command, cwd=workdir, env={**self.env, **env}, capture_output=True, text=True)

    def check_scenario(self, name):
        base, fields, expected, unexpected = SCENARIOS[name]
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp) / "organpipe"
            config = self.prepare(workdir, base, fields)

            result = self.run_tool(workdir, "bash", "OrganPipe.sh", "-d", ".", "-c", str(config), "-t", "2", "-np")
            self.assertEqual(result.returncode, 0, result.stdout[-3000:] + result.stderr[-3000:])
            jobs = set(job_counts(result.stdout))

            # Run up to the assemblies, so the later jobs can be listed
            snakemake = ["snakemake", "-d", ".", "-s", "workflow/Snakefile", "--profile", "profiles/local", "--cores", "2"]
            until = "split_assemblies" if base == SHORT else "run_mitohifi"
            result = self.run_tool(workdir, *snakemake, "--until", until)
            self.assertEqual(result.returncode, 0, result.stdout[-3000:] + result.stderr[-3000:])

            result = self.run_tool(workdir, *snakemake, "-n")
            self.assertEqual(result.returncode, 0, result.stdout[-3000:] + result.stderr[-3000:])
            jobs |= set(job_counts(result.stdout))

            self.assertEqual(expected - jobs, set(), f"missing jobs in {name}")
            self.assertEqual(unexpected & jobs, set(), f"unexpected jobs in {name}")

    def run_to_the_report(self, workdir, **env):
        base, fields, _, _ = SCENARIOS["mito_novoplasty_without_annotation"]
        config = self.prepare(workdir, base, fields)
        (workdir / "resources" / "pilon.sh").write_text(PILON)
        result = self.run_tool(workdir, "bash", "OrganPipe.sh", "-d", ".", "-c", str(config), "-t", "2", **env)
        self.assertEqual(result.returncode, 0, result.stdout[-3000:] + result.stderr[-3000:])
        return workdir / "workflow" / "reports" / "hs1"

    def test_run_to_the_report(self):
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp) / "organpipe"
            report = self.run_to_the_report(workdir)
            self.assertTrue((report / "parse_results.done").exists())
            self.assertFalse((report / "warnings.txt").exists())

            # All the k-mers of the stand-in assembler give the same genome
            with open(report / "pilon.csv", newline="") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 3)
            self.assertEqual({row["identical_group"] for row in rows}, {"group_1"})

            # The _pilon suffix must not reach the next steps
            for row in rows:
                run = f"{row['seed']}_kmer{row['kmer']}"
                fasta = workdir / "results" / "hs1" / "pilon" / "novoplasty" / run / row["assembly"] / f"{row['assembly']}.fasta"
                self.assertEqual(fasta.read_text().splitlines()[0], f">{row['assembly']}")

    def test_called_from_another_directory(self):
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp) / "organpipe"
            other = Path(tmp) / "other"
            other.mkdir()
            base, fields, _, _ = SCENARIOS["mito_novoplasty_without_annotation"]
            config = self.prepare(workdir, base, fields)
            result = self.run_tool(
                other, "bash", str(workdir / "OrganPipe.sh"), "-d", str(workdir), "-c", str(config), "-t", "1", "-np"
            )
            self.assertEqual(result.returncode, 0, result.stdout[-3000:] + result.stderr[-3000:])
            self.assertEqual(list(other.iterdir()), [])

    def test_fasta_reads_with_trimming_stop_the_run(self):
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp) / "organpipe"
            config = self.prepare(workdir, LONG, {"run_trimming": "Yes", "reads_path": "fasta_reads"})
            result = self.run_tool(workdir, "bash", "OrganPipe.sh", "-d", ".", "-c", str(config), "-t", "1", "-np")
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Trimming is only supported for .fastq.gz", result.stdout + result.stderr)

    def test_mitos_db_outside_the_working_directory_is_an_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp) / "organpipe"
            base, fields, _, _ = SCENARIOS["mito_novoplasty_without_annotation"]
            config = self.prepare(workdir, base, fields)
            result = self.run_tool(
                workdir, "bash", "OrganPipe.sh", "-d", ".", "-c", str(config), "-t", "1", "-np", "-mitos_db", tmp
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("-mitos_db must be inside", result.stderr)

    def peak_novoplasty_jobs(self, *options):
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp) / "organpipe"
            events = Path(tmp) / "events"
            events.touch()
            base, fields, _, _ = SCENARIOS["mito_novoplasty_without_annotation"]
            config = self.prepare(workdir, base, fields)
            (workdir / "resources" / "pilon.sh").write_text(PILON)
            result = self.run_tool(
                workdir, "bash", "OrganPipe.sh", "-d", ".", "-c", str(config), "-t", "8", *options, MOCK_EVENTS=str(events)
            )
            self.assertEqual(result.returncode, 0, result.stdout[-3000:] + result.stderr[-3000:])

            # The events are in the order they happened: the clock of a virtual machine can go back
            running = peak = 0
            for event in events.read_text().splitlines():
                running += 1 if event == "start" else -1
                peak = max(peak, running)
            return peak

    def test_memory_limit(self):
        # config_1sample.yaml has 3 k-mers and max_memory: 2, so each job takes 2.2 GB
        self.assertEqual(self.peak_novoplasty_jobs("-mem", "7"), 3)
        self.assertEqual(self.peak_novoplasty_jobs("-mem", "3"), 1)

    def jobs_for_contigs(self, fields):
        """Return the jobs after the assembly step, and the assemblies, of a sample that does not circularize."""
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp) / "organpipe"
            config = self.prepare(workdir, SHORT, {**NOVOPLASTY_ONLY, **ALL_STEPS, "kmers": "33", **fields})
            result = self.run_tool(workdir, "bash", "OrganPipe.sh", "-d", ".", "-c", str(config), "-t", "2", "-np")
            self.assertEqual(result.returncode, 0, result.stdout[-3000:] + result.stderr[-3000:])

            snakemake = ["snakemake", "-d", ".", "-s", "workflow/Snakefile", "--profile", "profiles/local", "--cores", "2"]
            result = self.run_tool(workdir, *snakemake, "--until", "split_assemblies", MOCK_CONTIGS="1")
            self.assertEqual(result.returncode, 0, result.stdout[-3000:] + result.stderr[-3000:])
            result = self.run_tool(workdir, *snakemake, "-n")
            self.assertEqual(result.returncode, 0, result.stdout[-3000:] + result.stderr[-3000:])

            fastas = (workdir / "results" / "hs1" / "assemblies" / "novoplasty").glob("*/*.fasta")
            return job_counts(result.stdout), sorted(fasta.name for fasta in fastas)

    def test_contigs_in_the_range_are_annotated_on_mitochondria(self):
        jobs, assemblies = self.jobs_for_contigs({})
        self.assertEqual(len(assemblies), 2, assemblies)
        self.assertTrue(all(name.startswith("CTG_1_hs1_") for name in assemblies), assemblies)
        for rule in ["run_pilon", "run_mitos2", "write_gb_mitos2", "rotate_genbank", "run_ogdraw"]:
            self.assertEqual(jobs.get(rule), 2, rule)
        self.assertEqual(jobs.get("run_nhmmer"), 4)

    def test_contigs_in_the_range_are_annotated_on_chloroplasts(self):
        jobs, assemblies = self.jobs_for_contigs(CHLORO)
        self.assertEqual(len(assemblies), 2, assemblies)
        for rule in ["run_pilon", "run_cpgavas2", "run_chloe", "write_gb_cpgavas2", "write_gb_chloe"]:
            self.assertEqual(jobs.get(rule), 2, rule)
        self.assertEqual(jobs.get("rotate_genbank"), 4)

    def test_seed_without_assembly_is_reported(self):
        with tempfile.TemporaryDirectory() as tmp:
            report = self.run_to_the_report(Path(tmp) / "organpipe", MOCK_NO_ASSEMBLY="1")
            self.assertTrue((report / "parse_results.done").exists())
            self.assertIn("produced no assembly", (report / "warnings.txt").read_text())
            self.assertFalse((report / "pilon.csv").exists())

    def test_mito_novoplasty_all_steps(self):
        self.check_scenario("mito_novoplasty_all_steps")

    def test_mito_novoplasty_without_annotation(self):
        self.check_scenario("mito_novoplasty_without_annotation")

    def test_mito_getorganelle_blank_options(self):
        self.check_scenario("mito_getorganelle_blank_options")

    def test_mito_trimming_with_adapters(self):
        self.check_scenario("mito_trimming_with_adapters")

    def test_chloro_both_assemblers_all_steps(self):
        self.check_scenario("chloro_both_assemblers_all_steps")

    def test_long_reads_trimmed_with_nhmmer(self):
        self.check_scenario("long_reads_trimmed_with_nhmmer")

    def test_long_reads_fastq_without_trimming(self):
        self.check_scenario("long_reads_fastq_without_trimming")

    def test_long_reads_fasta_without_trimming(self):
        self.check_scenario("long_reads_fasta_without_trimming")


if __name__ == "__main__":
    unittest.main()
