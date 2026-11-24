"""Unit tests of the workflow scripts. Run with: python -m unittest discover tests"""

import csv
import os
import random
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = ROOT / "workflow" / "scripts"
# Config files as written by the config generator (https://itvgenomics.github.io/config_generator/)
GENERATOR_CONFIGS = ROOT / "tests" / "data"
sys.path.insert(0, str(SCRIPTS))

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

import get_novoplasty_assemblies
import naming
import parse_results
import rotate_genbank


def touch_files(directory, names):
    for name in names:
        (Path(directory) / name).touch()


def run_script(name, *args, cwd=None):
    return subprocess.run(
        [sys.executable, str(SCRIPTS / name), *args],
        cwd=cwd,
        capture_output=True,
        text=True,
    )


class FindShortReads(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.path = self.tmp.name

    def names(self, sample):
        reads = naming.find_short_reads(sample, self.path)
        return os.path.basename(reads["forward"]), os.path.basename(reads["reverse"])

    def test_prefix_of_another_sample(self):
        touch_files(self.path, ["hs10_R1.fastq.gz", "hs10_R2.fastq.gz", "hs1_R1.fastq.gz", "hs1_R2.fastq.gz"])
        self.assertEqual(self.names("hs1"), ("hs1_R1.fastq.gz", "hs1_R2.fastq.gz"))

    def test_only_another_sample_has_files(self):
        touch_files(self.path, ["hs10_R1.fastq.gz", "hs10_R2.fastq.gz"])
        with self.assertRaisesRegex(ValueError, "No forward read file"):
            naming.find_short_reads("hs1", self.path)

    def test_longer_file_names(self):
        touch_files(self.path, ["s1_S3_L001_R1_001.fastq.gz", "s1_S3_L001_R2_001.fastq.gz"])
        self.assertEqual(self.names("s1"), ("s1_S3_L001_R1_001.fastq.gz", "s1_S3_L001_R2_001.fastq.gz"))

    def test_pair_markers(self):
        touch_files(self.path, ["s1_pair1.fastq.gz", "s1_pair2.fastq.gz"])
        self.assertEqual(self.names("s1"), ("s1_pair1.fastq.gz", "s1_pair2.fastq.gz"))

    def test_exact_name_wins(self):
        touch_files(self.path, ["s1_R1.fastq.gz", "s1_R2.fastq.gz", "s1_old_R1.fastq.gz", "s1_old_R2.fastq.gz"])
        self.assertEqual(self.names("s1"), ("s1_R1.fastq.gz", "s1_R2.fastq.gz"))

    def test_several_lanes_are_an_error(self):
        touch_files(
            self.path,
            ["s1_L001_R1.fastq.gz", "s1_L001_R2.fastq.gz", "s1_L002_R1.fastq.gz", "s1_L002_R2.fastq.gz"],
        )
        with self.assertRaisesRegex(ValueError, "More than one"):
            naming.find_short_reads("s1", self.path)

    def test_other_files_are_ignored(self):
        touch_files(self.path, ["s1_R1.fastq.gz", "s1_R2.fastq.gz", "s1_R1.fastq.gz.md5"])
        self.assertEqual(self.names("s1"), ("s1_R1.fastq.gz", "s1_R2.fastq.gz"))


class FindLongReads(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.path = self.tmp.name

    def test_prefix_of_another_sample(self):
        touch_files(self.path, ["hifi1.fastq.gz", "hifi10.fastq.gz"])
        self.assertEqual(os.path.basename(naming.find_long_reads("hifi1", self.path)), "hifi1.fastq.gz")

    def test_two_formats_are_an_error(self):
        touch_files(self.path, ["hifi.fastq.gz", "hifi.fasta"])
        with self.assertRaisesRegex(ValueError, "More than one"):
            naming.find_long_reads("hifi", self.path)


class FindSingleReads(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.path = self.tmp.name

    def name(self, sample):
        return os.path.basename(naming.find_single_reads(sample, self.path))

    def test_merged_reads_next_to_paired_files(self):
        touch_files(self.path, ["s1_R1.fastq.gz", "s1_R2.fastq.gz", "s1_M.fastq.gz", "s10_M.fastq.gz"])
        self.assertEqual(self.name("s1"), "s1_M.fastq.gz")
        self.assertEqual(self.name("s1_M"), "s1_M.fastq.gz")

    def test_long_read_formats_are_ignored(self):
        touch_files(self.path, ["s1.fasta", "s1.fastq.gz"])
        self.assertEqual(self.name("s1"), "s1.fastq.gz")

    def test_only_paired_files_are_an_error(self):
        touch_files(self.path, ["s1_R1.fastq.gz", "s1_R2.fastq.gz"])
        with self.assertRaisesRegex(ValueError, "No read file"):
            naming.find_single_reads("s1", self.path)


class CleanSequence(unittest.TestCase):
    def test_keeps_the_length(self):
        cleaned = get_novoplasty_assemblies.clean_sequence(Seq("ACGTNRYacgt*"), "contig")
        self.assertEqual(str(cleaned), "ACGTNNNACGTN")


class IdenticalAssemblies(unittest.TestCase):
    SEQUENCE = "ATGCGTACGTTAGCCTAGGATCCATTGCA"

    def test_rotation_and_strand(self):
        rotated = self.SEQUENCE[7:] + self.SEQUENCE[:7]
        reverse = str(Seq(self.SEQUENCE).reverse_complement())
        groups = parse_results.find_identical_groups(
            {"a": self.SEQUENCE, "b": rotated, "c": reverse, "d": reverse[3:] + reverse[:3]}
        )
        self.assertEqual(groups, {"a": "group_1", "b": "group_1", "c": "group_1", "d": "group_1"})

    def test_different_sequences(self):
        different = "T" + self.SEQUENCE[1:]
        groups = parse_results.find_identical_groups(
            {"a": self.SEQUENCE, "b": different, "c": self.SEQUENCE + "A"}
        )
        self.assertEqual(groups, {})

    def test_two_groups_are_numbered_in_order(self):
        other = "GGGATTACACCCGATTTAGCAT"
        groups = parse_results.find_identical_groups(
            {"a": self.SEQUENCE, "b": self.SEQUENCE, "c": other, "d": other[5:] + other[:5], "e": "ACGT"}
        )
        self.assertEqual(groups, {"a": "group_1", "b": "group_1", "c": "group_2", "d": "group_2"})

    def test_least_rotation_matches_brute_force(self):
        rng = random.Random(1)
        sequences = ["A", "AAAA", "ATAT", "ACAC", "TTTA", "GATTACAGATTACA"]
        sequences += ["".join(rng.choice("ACGTN") for _ in range(rng.randint(1, 40))) for _ in range(500)]
        sequences += ["".join(rng.choice("AT") for _ in range(rng.randint(1, 40))) for _ in range(500)]
        for sequence in sequences:
            start = parse_results.least_rotation(sequence)
            rotations = [sequence[i:] + sequence[:i] for i in range(len(sequence))]
            self.assertEqual(sequence[start:] + sequence[:start], min(rotations), sequence)

    def test_one_difference_separates_the_groups(self):
        changed = self.SEQUENCE[:10] + ("A" if self.SEQUENCE[10] != "A" else "C") + self.SEQUENCE[11:]
        self.assertEqual(parse_results.find_identical_groups({"a": self.SEQUENCE, "b": changed}), {})

    def test_repeated_sequence(self):
        self.assertEqual(parse_results.canonical_sequence("ATAT"), parse_results.canonical_sequence("TATA"))


class RotateGenbank(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.input = os.path.join(self.tmp.name, "in.gb")
        self.output = os.path.join(self.tmp.name, "out.gb")

    def write_record(self, strand, feature_type="rRNA", qualifiers=None):
        sequence = "AAAAACCCCCGGGGGTTTTTACGTACGTAC"
        feature = SeqFeature(
            FeatureLocation(10, 20, strand=strand),
            type=feature_type,
            qualifiers=qualifiers or {"product": "rrnS"},
        )
        cds = SeqFeature(FeatureLocation(21, 27, strand=1), type="CDS", qualifiers={"gene": "nad1"})
        record = SeqRecord(
            Seq(sequence),
            id="contig",
            name="contig",
            description="",
            annotations={"molecule_type": "DNA", "topology": "circular"},
            features=[feature, cds],
        )
        SeqIO.write(record, self.input, "genbank")
        return feature.extract(record.seq), cds.extract(record.seq)

    def find(self, record, feature_type):
        return next(f for f in record.features if f.type == feature_type)

    def test_mito_forward_strand(self):
        gene, cds = self.write_record(strand=1)
        rotate_genbank.rotate_genbank(self.input, self.output, "rrnS", "mito")
        record = SeqIO.read(self.output, "genbank")
        feature = self.find(record, "rRNA")
        self.assertEqual(int(feature.location.start), 0)
        self.assertEqual(feature.location.strand, 1)
        self.assertEqual(feature.extract(record.seq), gene)
        self.assertEqual(self.find(record, "CDS").extract(record.seq), cds)

    def test_mito_reverse_strand_is_flipped(self):
        gene, cds = self.write_record(strand=-1)
        rotate_genbank.rotate_genbank(self.input, self.output, "rrnS", "mito")
        record = SeqIO.read(self.output, "genbank")
        feature = self.find(record, "rRNA")
        self.assertEqual(int(feature.location.start), 0)
        self.assertEqual(feature.location.strand, 1)
        self.assertEqual(feature.extract(record.seq), gene)
        self.assertEqual(self.find(record, "CDS").extract(record.seq), cds)
        self.assertEqual(record.annotations["topology"], "circular")

    def test_chloro_forward_strand_is_flipped(self):
        gene, _ = self.write_record(strand=1, feature_type="gene", qualifiers={"gene": "psbA"})
        rotate_genbank.rotate_genbank(self.input, self.output, "psbA", "chloro")
        record = SeqIO.read(self.output, "genbank")
        feature = self.find(record, "gene")
        self.assertEqual(int(feature.location.start), 0)
        self.assertEqual(feature.location.strand, -1)
        self.assertEqual(feature.extract(record.seq), gene)

    def test_missing_gene_is_a_warning(self):
        self.write_record(strand=1, qualifiers={"product": "rrnL"})
        result = run_script(
            "rotate_genbank.py", "--organelle", "mito", "--start_gene", "rrnS", "--input", self.input, "--output", self.output
        )
        self.assertEqual(result.returncode, 0)
        self.assertIn(naming.WARNING_PREFIX, result.stdout)
        self.assertFalse(os.path.exists(self.output))

    def test_name_and_definition(self):
        self.write_record(strand=1)
        rotate_genbank.rotate_genbank(self.input, self.output, "rrnS", "mito", name="s1_x_rot", definition="rotated")
        record = SeqIO.read(self.output, "genbank")
        self.assertEqual((record.name, record.description), ("s1_x_rot", "rotated"))

    def test_name_is_kept_without_one(self):
        self.write_record(strand=1)
        rotate_genbank.rotate_genbank(self.input, self.output, "rrnS", "mito")
        self.assertEqual(SeqIO.read(self.output, "genbank").name, "contig")


class NcRNA(unittest.TestCase):
    def test_genes_with_the_same_product_are_kept(self):
        with tempfile.TemporaryDirectory() as tmp:
            genbank = os.path.join(tmp, "in.gb")
            features = [
                SeqFeature(FeatureLocation(0, 4, strand=1), type="tRNA", qualifiers={"product": "tRNA-Leu"}),
                SeqFeature(FeatureLocation(8, 12, strand=1), type="tRNA", qualifiers={"product": "tRNA-Leu"}),
                SeqFeature(FeatureLocation(4, 8, strand=1), type="rRNA", qualifiers={"product": "16S"}),
            ]
            record = SeqRecord(
                Seq("AAAACCCCGGGG"), id="x", name="x", description="",
                annotations={"molecule_type": "DNA"}, features=features,
            )
            SeqIO.write(record, genbank, "genbank")

            result = run_script("ncRNA_nhmmer.py", "--genbank", genbank, "--output_dir", tmp)
            self.assertEqual(result.returncode, 0, result.stderr)

            sequences = {r.id: str(r.seq) for r in SeqIO.parse(os.path.join(tmp, "rRNA-tRNA.fasta"), "fasta")}
            self.assertEqual(sequences, {"tRNA-Leu": "AAAA", "tRNA-Leu_2": "GGGG", "16S": "CCCC"})

    def test_no_annotation_writes_an_empty_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = run_script("ncRNA_nhmmer.py", "--genbank", os.path.join(tmp, "absent.gb"), "--output_dir", tmp)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn(naming.WARNING_PREFIX, result.stdout)
            self.assertEqual(os.path.getsize(os.path.join(tmp, "rRNA-tRNA.fasta")), 0)


class ReadNhmmerTable(unittest.TestCase):
    def test_description_with_spaces(self):
        with tempfile.TemporaryDirectory() as tmp:
            table = os.path.join(tmp, "tblout.out")
            with open(table, "w") as handle:
                handle.write("# target name\n")
                handle.write("RF00177 - q1 - 1 100 5 104 1 110 1500 + 1e-10 50.0 0.1 small subunit rRNA\n")
            df = parse_results.read_nhmmer_table(table)
            self.assertEqual(list(df.columns), parse_results.HMMER_COLUMNS)
            self.assertEqual(df.loc[0, "description of target"], "small subunit rRNA")
            self.assertEqual(df.loc[0, "strand"], "+")


class CreateSnakemakeConfig(unittest.TestCase):
    """Runs the script on configs that do not need NCBI."""

    BASE = {
        "reads_path": "reads", "organelle": "Mito", "genetic_code": "2", "sequencing_type": "Short",
        "run_trimming": "No", "annotation": "Yes", "run_nhmmer": "No", "run_images": "No", "adapters": "",
        "minlength": "", "minquality": "", "pacbio_adapters": "", "run_novoplasty": "Yes",
        "genome_range": "14000-18000", "reference": "", "seed_format": "fasta", "seed_file": "seed.fasta",
        "feature": "", "search_ncbi": "No", "search_genes": "", "search_term": "", "max_references": "",
        "kmers": "19,23", "max_memory": "2", "reads_length": "150", "insert_size": "250",
        "run_getorganelle": "No", "database": "", "n_rounds": "", "target_size": "", "spades_kmers": "",
        "extra_flags": "", "search_species": "", "n_references": "",
    }

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.dir = Path(self.tmp.name)
        (self.dir / "config").mkdir()
        (self.dir / "reads").mkdir()
        (self.dir / "seed.fasta").write_text(">s1\nACGT\n")

    def create(self, filename):
        result = run_script("create_snakemake_config.py", "--configfile", filename, cwd=self.dir)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        import yaml

        return yaml.safe_load((self.dir / "config" / "snakemake_config.yaml").read_text())["samples"]

    def write_yaml(self, **fields):
        import yaml

        data = {**self.BASE, **fields}
        data = {key: (None if value == "" else value) for key, value in data.items()}
        (self.dir / "c.yaml").write_text(yaml.safe_dump(data))

    def write_csv(self, *rows):
        with open(self.dir / "c.csv", "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["sample", *self.BASE], delimiter=";")
            writer.writeheader()
            for row in rows:
                writer.writerow({**self.BASE, **row})

    def test_sample_name_containing_all(self):
        touch_files(self.dir / "reads", ["small_x_R1.fastq.gz", "small_x_R2.fastq.gz", "other_R1.fastq.gz", "other_R2.fastq.gz"])
        self.write_yaml(sample="small_x")
        self.assertEqual(list(self.create("c.yaml")), ["small_x"])

    def test_sample_all(self):
        touch_files(
            self.dir / "reads",
            ["a_R1.fastq.gz", "a_R2.fastq.gz", "b_S1_R1_001.fastq.gz", "b_S1_R2_001.fastq.gz",
             "Undetermined_R1.fastq.gz", "Undetermined_R2.fastq.gz"],
        )
        self.write_yaml(sample="all")
        self.assertEqual(sorted(self.create("c.yaml")), ["a", "b_S1"])

    def test_sample_all_with_single_reads(self):
        touch_files(self.dir / "reads", ["a.fastq.gz", "b_M.fastq.gz", "c_R1.fastq.gz", "c_R2.fastq.gz", "Undetermined.fastq.gz"])
        self.write_yaml(sample="all", read_layout="Single", run_novoplasty="No", run_getorganelle="Yes", database="animal_mt")
        samples = self.create("c.yaml")
        self.assertEqual(sorted(samples), ["a", "b_M"])
        self.assertEqual(samples["a"]["read_layout"], "single")

    def test_read_layout_is_optional(self):
        # BASE has no read_layout, like the configs written before it existed
        self.write_yaml(sample="s1")
        self.assertEqual(self.create("c.yaml")["s1"]["read_layout"], "paired")
        self.write_csv({"sample": "s1"})
        self.assertEqual(self.create("c.csv")["s1"]["read_layout"], "paired")

    def test_single_kmer_as_a_number(self):
        self.write_yaml(sample="s1", kmers=33)
        self.assertEqual(self.create("c.yaml")["s1"]["kmers"], ["33"])

    def test_single_kmer_in_csv(self):
        self.write_csv({"sample": "s1", "kmers": "33"})
        self.assertEqual(self.create("c.csv")["s1"]["kmers"], ["33"])

    def test_list_of_target_sizes(self):
        self.write_yaml(sample="s1", target_size="39000,390000")
        self.assertEqual(self.create("c.yaml")["s1"]["target_size"], "39000,390000")

    def test_yaml_and_csv_produce_the_same_config(self):
        self.write_yaml(sample="s1")
        from_yaml = self.create("c.yaml")
        self.write_csv({"sample": "s1"})
        self.assertEqual(self.create("c.csv"), from_yaml)

    def test_text_fields_are_lowercase(self):
        self.write_yaml(sample="s1", run_trimming="Yes", minlength=50, minquality=20)
        entry = self.create("c.yaml")["s1"]
        self.assertEqual(
            (entry["organelle"], entry["sequencing_type"], entry["run_trimming"]), ("mito", "short", "yes")
        )
        self.assertEqual((entry["minlength"], entry["minquality"]), (50, 20))

    def test_seeds_are_only_prepared_for_novoplasty(self):
        # The seed file does not exist: preparing the seeds would stop the run
        self.write_yaml(
            sample="s1", run_novoplasty="No", run_getorganelle="Yes", database="animal_mt",
            seed_file="missing.fasta", kmers="",
        )
        entry = self.create("c.yaml")["s1"]
        self.assertEqual((entry["seeds"], entry["kmers"]), ([], []))

    def test_blank_and_unquoted_values(self):
        # YAML reads unquoted Yes/No as booleans
        self.write_yaml(sample="s1", reference="  ", run_images=False, annotation=True)
        entry = self.create("c.yaml")["s1"]
        self.assertEqual((entry["reference"], entry["run_images"], entry["annotation"]), ("", "no", "yes"))


class CheckConfig(unittest.TestCase):
    def run_check(self, text):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "c.yaml")
            with open(path, "w") as handle:
                handle.write(text)
            return run_script("check_config.py", "-c", path)

    def test_invalid_yaml_stops_the_run(self):
        self.assertNotEqual(self.run_check("sample: [unclosed").returncode, 0)

    def test_yaml_that_is_not_a_mapping_stops_the_run(self):
        self.assertNotEqual(self.run_check("- a\n- b\n").returncode, 0)

    def test_sample_names(self):
        import check_config

        self.assertEqual(check_config.check_sample_names({"sample": "a_1,b2"}), [])
        self.assertEqual(check_config.check_sample_names({"sample": "all"}), [])
        self.assertEqual(len(check_config.check_sample_names({"sample": "a-1"})), 1)

    def test_max_memory(self):
        import check_config

        config = {"sequencing_type": "short", "run_novoplasty": "yes", "reads_length": 150, "insert_size": 300, "kmers": "33"}
        for value, errors in [(40, 0), (40.0, 0), ("40", 0), ("", 1), (None, 1), ("8G", 1), (0, 1)]:
            self.assertEqual(
                len(check_config.check_assembler_requirements({**config, "max_memory": value})), errors, value
            )

    def test_pacbio_adapters_only_matter_when_trimming(self):
        import check_config

        long_reads = {"sequencing_type": "long", "pacbio_adapters": None}
        self.assertEqual(check_config.check_pacbio_adapters({**long_reads, "run_trimming": "no"}), [])
        self.assertEqual(len(check_config.check_pacbio_adapters({**long_reads, "run_trimming": "yes"})), 1)

    def test_short_reads_in_fasta_cannot_be_trimmed(self):
        # fastp reads no reads from a .fasta.gz file and still ends without error
        import check_config

        with tempfile.TemporaryDirectory() as tmp:
            touch_files(tmp, ["p_R1.fasta.gz", "p_R2.fasta.gz", "q_R1.fastq.gz", "q_R2.fastq.gz", "s.fasta.gz", "t.fastq.gz"])
            short = {"sequencing_type": "short", "reads_path": tmp}

            for sample, layout, errors in [
                ("p", "paired", 2), ("q", "paired", 0), ("all", "paired", 2),
                ("s", "single", 1), ("t", "single", 0), ("all", "single", 1),
            ]:
                config = {**short, "sample": sample, "read_layout": layout}
                found = check_config.check_reads_path({**config, "run_trimming": "Yes"})
                self.assertEqual(len(found), errors, (sample, layout, found))
                self.assertEqual(check_config.check_reads_path({**config, "run_trimming": "No"}), [], (sample, layout))

    def test_long_reads_in_fasta_cannot_be_trimmed(self):
        import check_config

        with tempfile.TemporaryDirectory() as tmp:
            # Undetermined files are not samples of `sample: all`
            touch_files(tmp, ["hifi.fasta", "ont.fastq.gz", "Undetermined.fasta"])
            long_reads = {"sequencing_type": "long", "reads_path": tmp}

            for sample, errors in [("hifi", 1), ("ont", 0), ("all", 1)]:
                config = {**long_reads, "sample": sample}
                found = check_config.check_reads_path({**config, "run_trimming": "Yes"})
                self.assertEqual(len(found), errors, (sample, found))
                self.assertEqual(check_config.check_reads_path({**config, "run_trimming": "No"}), [], sample)

    def test_configs_from_the_generator(self):
        # Blank fields are written as "", and the NOVOPlasty fields stay blank when it is disabled
        for name in ["generator_getorganelle.yaml", "generator_getorganelle.csv", "generator_longreads.yaml"]:
            result = run_script("check_config.py", "-c", str(GENERATOR_CONFIGS / name), cwd=ROOT)
            self.assertEqual(result.returncode, 0, name + result.stdout + result.stderr)

    def test_blank_values(self):
        import config_values

        for value in [None, "", "  ", float("nan")]:
            self.assertEqual(config_values.as_text(value), "", value)
        self.assertEqual([config_values.as_text(v) for v in [True, False, 2, " a"]], ["Yes", "No", "2", " a"])

    def test_unquoted_yes_and_no(self):
        import check_config

        fields = ["annotation", "run_trimming", "run_nhmmer", "run_images", "run_novoplasty"]
        config = {"sequencing_type": "short", "run_getorganelle": True, **dict.fromkeys(fields, False)}
        self.assertEqual(check_config.check_yes_no_fields(config), [])

    def test_long_reads_only_need_annotation(self):
        import check_config

        fields = ["run_trimming", "run_nhmmer", "run_images", "run_novoplasty", "run_getorganelle"]
        config = {"sequencing_type": "long", "annotation": "Yes", **dict.fromkeys(fields, None)}
        self.assertEqual(check_config.check_yes_no_fields(config), [])
        self.assertEqual(len(check_config.check_yes_no_fields({**config, "annotation": None})), 1)
        self.assertEqual(len(check_config.check_yes_no_fields({**config, "run_nhmmer": "maybe"})), 1)
        # Short reads need every step
        self.assertEqual(len(check_config.check_yes_no_fields({**config, "sequencing_type": "short"})), 5)

    def test_novoplasty_fields_are_ignored_without_novoplasty(self):
        import check_config

        config = {
            "sequencing_type": "short", "run_novoplasty": "No", "run_getorganelle": "Yes", "database": "animal_mt",
            "genome_range": "", "reference": "", "seed_format": "", "search_ncbi": "", "kmers": "", "max_memory": "",
        }
        checks = [
            check_config.check_reference_file,
            check_config.check_genome_range,
            check_config.check_seed_format_and_file,
            check_config.check_search_ncbi_requirements,
            check_config.check_assembler_requirements,
        ]
        for check in checks:
            self.assertEqual(check(config), [], check.__name__)

    def test_novoplasty_needs_kmers_and_seeds(self):
        import check_config

        config = {
            "sequencing_type": "short", "run_novoplasty": "Yes", "reads_length": "150", "insert_size": "300",
            "max_memory": "40",
        }
        for kmers, errors in [("19,23,33", 0), ("33", 0), (33, 0), ("", 1), (None, 1), ("19;23", 1), ("0", 1)]:
            self.assertEqual(
                len(check_config.check_assembler_requirements({**config, "kmers": kmers})), errors, kmers
            )

        # Without a seed file the seeds come from NCBI, and a blank search_ncbi means "No"
        for seed_format, search_ncbi, errors in [("", "", 1), ("", "No", 1), ("", "Yes", 0), ("fa", "Yes", 1)]:
            fields = {"seed_format": seed_format, "search_ncbi": search_ncbi}
            self.assertEqual(
                len(check_config.check_seed_format_and_file({**config, **fields})), errors, fields
            )

    def test_read_layout(self):
        import check_config

        getorganelle = {"sequencing_type": "short", "run_novoplasty": "No", "run_getorganelle": "Yes"}
        novoplasty = {**getorganelle, "run_novoplasty": "Yes"}
        cases = [
            (getorganelle, None, 0), (getorganelle, "", 0), (getorganelle, "Single", 0), (getorganelle, "paired", 0),
            (getorganelle, "merged", 1), (novoplasty, "single", 1), (novoplasty, "paired", 0), (novoplasty, None, 0),
        ]
        for config, layout, errors in cases:
            self.assertEqual(
                len(check_config.check_read_layout({**config, "read_layout": layout})), errors, (config, layout)
            )

    def test_reads_path_with_single_reads(self):
        import check_config

        with tempfile.TemporaryDirectory() as tmp:
            touch_files(tmp, ["s1_M.fastq.gz", "s2_R1.fastq.gz", "s2_R2.fastq.gz"])
            config = {"reads_path": tmp, "sequencing_type": "short", "read_layout": "single"}
            self.assertEqual(check_config.check_reads_path({**config, "sample": "s1"}), [])
            self.assertEqual(len(check_config.check_reads_path({**config, "sample": "s2"})), 1)
            # A paired sample needs its _R1/_R2 files
            self.assertEqual(len(check_config.check_reads_path({**config, "sample": "s1", "read_layout": ""})), 1)

    def test_genetic_code_is_only_needed_when_used(self):
        import check_config

        cases = [
            ({"sequencing_type": "short", "organelle": "mito", "annotation": "Yes"}, 1),
            ({"sequencing_type": "short", "organelle": "mito", "annotation": "No"}, 0),
            ({"sequencing_type": "short", "organelle": "chloro", "annotation": "Yes"}, 0),
            ({"sequencing_type": "long", "organelle": "mito", "annotation": "No"}, 1),
        ]
        for config, errors in cases:
            for blank in [None, ""]:
                self.assertEqual(
                    len(check_config.check_genetic_code({**config, "genetic_code": blank})), errors, config
                )

        mito = cases[0][0]
        for code, errors in [(2, 0), ("11", 0), (2.0, 0), (7, 1), ("two", 1)]:
            self.assertEqual(len(check_config.check_genetic_code({**mito, "genetic_code": code})), errors, code)


class SplitFasta(unittest.TestCase):
    def test_file_name_matches_the_header(self):
        with tempfile.TemporaryDirectory() as tmp:
            fasta = os.path.join(tmp, "in.fasta")
            with open(fasta, "w") as handle:
                handle.write(">CA_1_s_seed_19\nACGT\nACGT\n>CA_1_s_seed_19\nTTTT\n")
            output = os.path.join(tmp, "out")
            result = run_script("split_fasta.py", "--fasta_file", fasta, "--assembler", "novoplasty", "--output_dir", output)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(sorted(os.listdir(output)), ["CA_1_s_seed_19_1.fasta", "CA_1_s_seed_19_2.fasta"])
            for name in os.listdir(output):
                record = SeqIO.read(os.path.join(output, name), "fasta")
                self.assertEqual(f"{record.id}.fasta", name)

    def test_empty_assembly_creates_an_empty_directory(self):
        with tempfile.TemporaryDirectory() as tmp:
            fasta = os.path.join(tmp, "in.fasta")
            open(fasta, "w").close()
            for assembler in ["novoplasty", "getorganelle"]:
                output = os.path.join(tmp, assembler)
                result = run_script("split_fasta.py", "--fasta_file", fasta, "--assembler", assembler, "--output_dir", output)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(os.listdir(output), [])


class GenbankIdentity(unittest.TestCase):
    def identity(self, *args, **kwargs):
        return naming.genbank_identity(*args, **kwargs)

    def test_novoplasty_uses_the_seed_name(self):
        name, definition = self.identity(
            "DNA1", "chloro", "novoplasty", "1-rbcL_1858743_kmer33", "OPT_2_DNA1_1-rbcL_1858743_33_2", "chloe"
        )
        self.assertEqual(name, "DNA1_1-rbcL_1858743_k33_OPT2_2_chloe")
        self.assertEqual(
            definition,
            "DNA1 chloroplast, NOVOPlasty assembly OPT_2_DNA1_1-rbcL_1858743_33_2 "
            "(seed 1-rbcL_1858743, k-mer 33), annotated with Chloe",
        )

    def test_rotated(self):
        name, definition = self.identity(
            "s1", "mito", "novoplasty", "1-COI_9606_kmer19", "CA_1_s1_1-COI_9606_19_1", "mitos2", rotated=True
        )
        self.assertEqual(name, "s1_1-COI_9606_k19_CA1_1_rot")
        self.assertTrue(definition.endswith("annotated with MITOS2, rotated to start at rrnS"))

    def test_getorganelle(self):
        name, definition = self.identity(
            "s1", "chloro", "getorganelle", "all", "embplant_pt_complete_1", "cpgavas2", rotated=True
        )
        self.assertEqual(name, "s1_embplant_pt_complete_1_cpgavas2_rot")
        self.assertEqual(
            definition,
            "s1 chloroplast, GetOrganelle assembly embplant_pt_complete_1, annotated with CPGAVAS2, "
            "rotated to start at psbA",
        )

    def test_mitohifi(self):
        self.assertEqual(
            self.identity("s1", "mito", "mitohifi", "PQ540662.1"),
            ("s1_PQ540662.1_MitoHiFi", "s1 mitochondrion, MitoHiFi assembly (reference PQ540662.1)"),
        )

    def test_long_seed_names_are_shortened_and_stay_unique(self):
        names = set()
        for number in [1, 2]:
            seed = f"{number}-NADH_dehydrogenase_subunit_4L"
            name, definition = self.identity(
                "ITV00872_25", "mito", "novoplasty", f"{seed}_kmer19", f"CTG_1_ITV00872_25_{seed}_19_3", "mitos2",
                rotated=True,
            )
            self.assertLessEqual(len(name), naming.LOCUS_MAX_LENGTH)
            self.assertTrue(name.startswith(f"ITV00872_25_{number}-NADH"), name)
            self.assertTrue(name.endswith("_k19_CTG1_3_rot"), name)
            self.assertIn(f"seed {seed},", definition)
            names.add(name)
        self.assertEqual(len(names), 2)

    def test_only_locus_characters(self):
        name, _ = self.identity("s1", "mito", "novoplasty", "1-a:b c_kmer21", "CA_1_s1_1-a:b c_21_1", "mitos2")
        self.assertEqual(name, "s1_1-a_b_c_k21_CA1_1")

    def test_unexpected_assembly_name(self):
        name, _ = self.identity("s1", "mito", "novoplasty", "seed_kmer21", "other_name", "mitos2")
        self.assertEqual(name, "s1_other_name")

    def test_long_name_is_written_and_read_back(self):
        with tempfile.TemporaryDirectory() as tmp:
            genbank = os.path.join(tmp, "out.gb")
            record = SeqRecord(Seq("ACGT" * 10), annotations={"molecule_type": "DNA", "topology": "circular"})
            name = "s" * naming.LOCUS_MAX_LENGTH
            naming.write_genbank(naming.set_genbank_identity(record, name, "a description"), genbank)
            record = SeqIO.read(genbank, "genbank")
            self.assertEqual((record.name, record.description), (name, "a description"))
            self.assertEqual(record.annotations["topology"], "circular")


class MitosGenbankName(unittest.TestCase):
    def test_name_and_definition(self):
        with tempfile.TemporaryDirectory() as tmp:
            with open(os.path.join(tmp, "mitos2.fasta"), "w") as handle:
                handle.write(">temp_header\n" + "ACGT" * 10 + "\n")
            with open(os.path.join(tmp, "result.gff"), "w") as handle:
                handle.write("temp_header\tmitos\tgene\t1\t12\t.\t+\t.\tgene_id=nad1\n")
            output = os.path.join(tmp, "out.gb")
            result = run_script(
                "mito_gff2genbank.py", "--code", "2", "--annotation_dir", tmp, "--output", output,
                "--name", "s1_seed_k21_CA1_1", "--definition", "s1 mitochondrion",
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            record = SeqIO.read(output, "genbank")
            self.assertEqual((record.name, record.description), ("s1_seed_k21_CA1_1", "s1 mitochondrion"))
            self.assertEqual(record.features[0].qualifiers["gene"], ["nad1"])


class Cpgavas2Genbank(unittest.TestCase):
    """CPGAVAS2 writes {pid}.gbf, with the organism of its example file."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.output = os.path.join(self.tmp.name, "out.gb")

    def write_gbf(self, with_genes):
        source = SeqFeature(
            FeatureLocation(0, 40, strand=1),
            type="source",
            qualifiers={"organism": ["Arabidopsis thaliana"], "db_xref": ["taxon:3702"]},
        )
        gene = SeqFeature(FeatureLocation(0, 12, strand=1), type="gene", qualifiers={"gene": ["psbA"]})
        record = SeqRecord(
            Seq("ACGT" * 10), id="123", name="123", description="Chloroplast, complete genome",
            annotations={"molecule_type": "DNA", "topology": "circular", "organism": "Arabidopsis thaliana"},
            features=[source, gene] if with_genes else [source],
        )
        SeqIO.write(record, os.path.join(self.tmp.name, "123.gbf"), "genbank")

    def run_cpgavas2(self):
        return run_script(
            "chloro_gff2genbank.py", "--software", "cpgavas2", "--annotation_dir", self.tmp.name,
            "--output", self.output, "--name", "s1_seed_k21_CA1_1_cpgavas2", "--definition", "s1 chloroplast",
        )

    def test_name_and_organism(self):
        self.write_gbf(with_genes=True)
        result = self.run_cpgavas2()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn(naming.WARNING_PREFIX, result.stdout)
        record = SeqIO.read(self.output, "genbank")
        self.assertEqual((record.name, record.description), ("s1_seed_k21_CA1_1_cpgavas2", "s1 chloroplast"))
        self.assertEqual(record.annotations["organism"], "Unknown.")
        self.assertEqual(record.features[0].qualifiers["organism"], ["unknown"])
        self.assertNotIn("db_xref", record.features[0].qualifiers)
        self.assertEqual(record.annotations["topology"], "circular")

    def test_source_feature_only_is_a_warning(self):
        self.write_gbf(with_genes=False)
        result = self.run_cpgavas2()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn(naming.WARNING_PREFIX, result.stdout)
        self.assertIn("annotated no genes", result.stdout)
        self.assertFalse(os.path.exists(self.output))


class Cpgavas2ReportFiles(unittest.TestCase):
    def test_report_and_its_local_links(self):
        with tempfile.TemporaryDirectory() as tmp:
            html = """<html><head><link href="/static_cpgavas2018/css/plantidit.css"></head>
            <a href="123.gbf">GenBank</a> <a href="123_image4.png"><img src="123_image5.png"></a>
            <a href="123_reportTable.txt">table</a> <a href="http://tandem.bu.edu/trf/trf.html">TRF</a>
            <a href="123_missing.txt">missing</a>
            <!-- <a href="123.tar.gz">results</a> -->
            </html>"""
            Path(tmp, "123.html").write_text(html)
            touch_files(tmp, [
                "123.gbf", "123_image4.png", "123_image5.png", "123_reportTable.txt",
                "123.tar.gz", "123.fas.suf", "Cp01.maker.gff", "123.fas.2.7.7.80.10.50.500.1.html",
            ])
            os.mkdir(os.path.join(tmp, "ref"))
            self.assertEqual(
                parse_results.cpgavas2_report_files(tmp),
                ["123.html", "123.gbf", "123_image4.png", "123_image5.png", "123_reportTable.txt"],
            )


if __name__ == "__main__":
    unittest.main()
