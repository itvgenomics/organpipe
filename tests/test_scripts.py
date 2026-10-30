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

        config = {"sequencing_type": "short", "run_novoplasty": "yes", "reads_length": 150, "insert_size": 300}
        for value, errors in [(40, 0), (40.0, 0), ("40", 0), ("", 1), (None, 1), ("8G", 1), (0, 1)]:
            self.assertEqual(
                len(check_config.check_assembler_requirements({**config, "max_memory": value})), errors, value
            )

    def test_pacbio_adapters_only_matter_when_trimming(self):
        import check_config

        long_reads = {"sequencing_type": "long", "pacbio_adapters": None}
        self.assertEqual(check_config.check_pacbio_adapters({**long_reads, "run_trimming": "no"}), [])
        self.assertEqual(len(check_config.check_pacbio_adapters({**long_reads, "run_trimming": "yes"})), 1)


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


if __name__ == "__main__":
    unittest.main()
