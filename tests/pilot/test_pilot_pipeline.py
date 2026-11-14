"""End-to-end pilot preparation and verification on a tiny synthetic design.

No network, biological tool, or real genome is used. A fake ``art_illumina``
reproduces only ART's command-line and paired FASTQ output contract, and the
download step copies local gzip files.
"""

import contextlib
import csv
import gzip
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[2]


def module(name):
    spec = importlib.util.spec_from_file_location(f"{name}_e2e", ROOT / "scripts" / f"{name}.py")
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


prep = module("prepare_pilot")

KNOWN, HOLDOUT, DISTRACTOR = "GCF_000000001.1", "GCF_000000002.1", "GCF_000000003.1"

FAKE_ART = r'''
import os
import sys

args, values, i = sys.argv[1:], {}, 0
while i < len(args):
    if args[i] in ("-p", "-na"):
        values[args[i]] = True
        i += 1
    else:
        values[args[i]] = args[i + 1]
        i += 2
if os.environ.get("FAKE_ART_EXIT"):
    print("fake ART failure requested")
    sys.exit(int(os.environ["FAKE_ART_EXIT"]))
if not (values.get("-p") and values.get("-na")):
    sys.exit("paired mode without alignment output is required")
length = int(values["-l"])
pairs = max(0, int(values["-c"]) - int(os.environ.get("FAKE_ART_SHORTFALL", "0")))
seed = int(values["-rs"])
with open(values["-i"]) as stream:
    lines = stream.read().split("\n")
name, sequence = lines[0][1:].split()[0], "".join(lines[1:])
complement = str.maketrans("ACGT", "TGCA")
with open(values["-o"] + "1.fq", "w") as one, open(values["-o"] + "2.fq", "w") as two:
    for n in range(pairs):
        start = (seed + 7 * n) % (len(sequence) - length + 1)
        read = sequence[start:start + length]
        one.write(f"@{name}-{n}/1\n{read}\n+\n{'I' * length}\n")
        two.write(f"@{name}-{n}/2\n{read.translate(complement)[::-1]}\n+\n{'I' * length}\n")
print(f"fake ART wrote {pairs} pairs")
'''


def read_tsv(path):
    with Path(path).open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def build_design(root):
    """Write a three-genome, two-sample design and its offline source archives."""
    rng = random.Random(7)
    design, sources = root / "design", root / "sources"
    design.mkdir()
    sources.mkdir()
    genomes = [
        (KNOWN, "known", "562", [400, 150]),
        (HOLDOUT, "species_holdout_A", "1423", [500]),
        (DISTRACTOR, "distractor", "1299", [300]),
    ]
    rows = []
    for genome_id, role, species, lengths in genomes:
        archive = sources / f"{genome_id}.fna.gz"
        with gzip.open(archive, "wt") as out:
            for index, length in enumerate(lengths):
                sequence = "".join(rng.choice("ACGT") for _ in range(length))
                out.write(f">rec{index} synthetic record\n")
                out.writelines(sequence[i:i + 60] + "\n" for i in range(0, length, 60))
        rows.append({
            "genome_id": genome_id, "assembly_name": "ASM1v1", "organism": "synthetic",
            "taxid": species + "0", "species_taxid": species, "domain": "Bacteria", "role": role,
            "bases": sum(lengths), "records": len(lengths),
            "source_md5": hashlib.md5(archive.read_bytes()).hexdigest(),
        })
    with (design / "genomes.tsv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    (design / "composition.tsv").write_text(
        "sample_id\tcommunity\tgenome_id\tpairs\n"
        f"A1\tcommunity_A\t{KNOWN}\t30\n"
        f"A1\tcommunity_A\t{HOLDOUT}\t10\n"
        f"B1\tcommunity_B\t{KNOWN}\t40\n"
    )
    config = {
        "schema_version": 1, "seed": 11, "genomes": "genomes.tsv",
        "composition": "composition.tsv", "read_pairs_per_sample": 40,
        "simulation": {
            "tool": "art_illumina", "version": "fake", "profile": "HS25", "read_length": 30,
            "fragment_mean": 60, "fragment_sd": 5, "gzip_level": 3,
            "art_surplus_fraction": 0.01, "art_min_surplus_pairs": 2,
        },
        "resources": {
            "threads": 4, "memory_gib": 1, "scratch_gib": 3, "genome_download_limit_mib": 1,
            "simulation_jobs": 1, "sample_jobs": 2,
        },
        "catalogues": {"all_known": [], "species_holdout": [HOLDOUT]},
    }
    (design / "pilot.json").write_text(json.dumps(config, indent=2))
    return design / "pilot.json", sources


@unittest.skipIf(os.name == "nt", "the fake ART executable relies on a POSIX shebang")
class PreparationEndToEndTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.config, self.sources = build_design(self.root)
        self.art = self.root / "bin" / "art_illumina"
        self.art.parent.mkdir()
        self.art.write_text(f"#!{sys.executable}\n{FAKE_ART}")
        self.art.chmod(0o755)
        self.urls = []

    def fake_download(self, url, destination, max_bytes):
        self.assertTrue(url.startswith("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/000/"))
        self.urls.append(url)
        accession = "_".join(url.rsplit("/", 1)[1].split("_")[:2])
        data = (self.sources / f"{accession}.fna.gz").read_bytes()
        if len(data) > max_bytes:
            raise ValueError("genome download byte limit exceeded")
        destination.write_bytes(data)
        return len(data)

    @contextlib.contextmanager
    def offline(self, art=True, **environment):
        real_which = shutil.which

        def which(name, *args, **kwargs):
            if name == "art_illumina":
                return str(self.art) if art else None
            return real_which(name, *args, **kwargs)

        with contextlib.ExitStack() as stack:
            stack.enter_context(mock.patch.object(prep, "download", self.fake_download))
            stack.enter_context(mock.patch.object(prep.time, "sleep"))
            stack.enter_context(mock.patch.object(prep.shutil, "which", which))
            stack.enter_context(mock.patch.dict(os.environ, environment))
            stack.enter_context(contextlib.redirect_stdout(io.StringIO()))
            yield

    def prepared(self, name="dataset", **environment):
        destination = self.root / name
        with self.offline(**environment):
            prep.prepare(self.config, destination)
        return destination

    def verify(self, destination):
        with contextlib.redirect_stdout(io.StringIO()):
            prep.verify(self.config, destination)

    def test_prepared_dataset_preserves_truth_catalogues_and_provenance(self):
        destination = self.prepared()
        self.verify(destination)
        self.assertFalse((destination / "PREPARING").exists())
        manifest = json.loads((destination / "dataset.json").read_text())
        self.assertEqual(manifest["status"], "complete")
        self.assertEqual(len(self.urls), 3)

        # Records are qualified by assembly, so identical source names stay distinct.
        record_ids = [row["record_id"] for row in read_tsv(destination / "record-map.tsv")]
        self.assertEqual(len(record_ids), 4)
        self.assertEqual(len(set(record_ids)), 4)
        self.assertIn(f"{KNOWN}__rec0", record_ids)
        self.assertIn(f"{HOLDOUT}__rec0", record_ids)

        for case, expected in [("all_known", {KNOWN, HOLDOUT, DISTRACTOR}),
                               ("species_holdout", {KNOWN, DISTRACTOR})]:
            directory = destination / "catalogues" / case
            rows = read_tsv(directory / "catalogue.tsv")
            self.assertEqual({row["reference_id"] for row in rows}, expected, case)
            listed = (directory / "references.txt").read_text().splitlines()
            self.assertEqual(listed, [row["fasta"] for row in rows])
            for row in rows:
                fasta = (directory / row["fasta"]).resolve()
                self.assertEqual(fasta.parent, (destination / "references").resolve())
                self.assertEqual(prep.digest(fasta), row["sha256"])

        self.assertEqual(read_tsv(destination / "samples.tsv"), [
            {"sample_id": "A1", "read1": "reads/A1_R1.fastq.gz", "read2": "reads/A1_R2.fastq.gz"},
            {"sample_id": "B1", "read1": "reads/B1_R1.fastq.gz", "read2": "reads/B1_R2.fastq.gz"},
        ])
        self.assertEqual(read_tsv(destination / "metadata.tsv"), [
            {"sample_id": "A1", "community": "community_A"},
            {"sample_id": "B1", "community": "community_B"},
        ])

        # Read headers carry sample-aware, mate-aware source truth.
        expected_counts = {"A1": {KNOWN: 30, HOLDOUT: 10}, "B1": {KNOWN: 40}}
        for sample, expected in expected_counts.items():
            reads = destination / "reads"
            pairs = list(prep.paired_records(reads / f"{sample}_R1.fastq.gz", reads / f"{sample}_R2.fastq.gz"))
            counts = {}
            for ordinal, (one, two) in enumerate(pairs, 1):
                for mate, record in [(1, one), (2, two)]:
                    owner, genome, number = record[0].split("|")
                    self.assertEqual((owner, number), (sample, f"{ordinal:09d}/{mate}"))
                    self.assertEqual(len(record[1]), 30)
                counts[genome] = counts.get(genome, 0) + 1
            self.assertEqual(counts, expected)
            self.assertNotIn(DISTRACTOR, counts)

        # ART is invoked with argument lists and a distinct recorded seed per part.
        commands = [command for report in manifest["samples"] for command in report["commands"]]
        self.assertEqual(len(commands), 5)  # two KNOWN records per sample + one HOLDOUT record
        self.assertTrue(all(isinstance(command, list) for command in commands))
        seeds = [command[command.index("-rs") + 1] for command in commands]
        self.assertEqual(len(set(seeds)), len(seeds))
        self.assertEqual(len(list((destination / "logs").iterdir())), 5)

    def test_preparation_is_byte_reproducible_across_destinations(self):
        first, second = self.prepared("first"), self.prepared("second")
        for sample in ["A1", "B1"]:
            for mate in [1, 2]:
                name = f"reads/{sample}_R{mate}.fastq.gz"
                self.assertEqual(prep.digest(first / name), prep.digest(second / name), name)
        for name in ["record-map.tsv", "catalogues/species_holdout/catalogue.tsv"]:
            self.assertEqual((first / name).read_bytes(), (second / name).read_bytes(), name)

    def test_existing_or_partial_destinations_are_never_reused(self):
        destination = self.prepared()
        before = prep.digest(destination / "reads/A1_R1.fastq.gz")
        with self.assertRaisesRegex(ValueError, "output exists"), self.offline():
            prep.prepare(self.config, destination)
        self.assertEqual(prep.digest(destination / "reads/A1_R1.fastq.gz"), before)

    def test_verify_detects_changed_reads_references_and_configuration(self):
        destination = self.prepared()
        read = destination / "reads/A1_R1.fastq.gz"
        original = read.read_bytes()
        read.write_bytes(gzip.compress(gzip.decompress(original).replace(b"\nA", b"\nC", 1), mtime=0))
        with self.assertRaisesRegex(ValueError, "checksum mismatch"):
            self.verify(destination)
        read.write_bytes(original)
        self.verify(destination)

        reference = destination / "references" / f"{KNOWN}.fna"
        reference.write_text(reference.read_text().lower())
        with self.assertRaisesRegex(ValueError, "reference checksum mismatch"):
            self.verify(destination)

        config = json.loads(self.config.read_text())
        config["seed"] += 1
        self.config.write_text(json.dumps(config))
        with self.assertRaisesRegex(ValueError, "another configuration"):
            self.verify(destination)

    def test_art_shortfall_fails_and_leaves_an_unverifiable_marker(self):
        destination = self.root / "short"
        with self.assertRaisesRegex(ValueError, "below target"), self.offline(FAKE_ART_SHORTFALL="1000"):
            prep.prepare(self.config, destination)
        self.assertTrue((destination / "PREPARING").exists())
        self.assertFalse((destination / "dataset.json").exists())
        with self.assertRaises(OSError):
            self.verify(destination)

    def test_failed_art_process_is_not_reported_as_completed(self):
        destination = self.root / "failed"
        with self.assertRaises(subprocess.CalledProcessError), self.offline(FAKE_ART_EXIT="3"):
            prep.prepare(self.config, destination)
        self.assertTrue((destination / "PREPARING").exists())
        self.assertFalse((destination / "dataset.json").exists())
        self.assertEqual(list((destination / "reads").iterdir()), [])

    def test_source_checksum_mismatch_stops_before_simulation(self):
        archive = self.sources / f"{KNOWN}.fna.gz"
        archive.write_bytes(gzip.compress(b">rec0\nACGT\n", mtime=0))
        destination = self.root / "tampered-source"
        with self.assertRaisesRegex(ValueError, "MD5 mismatch"), self.offline():
            prep.prepare(self.config, destination)
        self.assertTrue((destination / "PREPARING").exists())
        self.assertEqual(list((destination / "reads").iterdir()), [])

    def test_preflight_refusals_create_no_destination(self):
        destination = self.root / "never-created"
        with self.assertRaisesRegex(ValueError, "art_illumina is missing"), self.offline(art=False):
            prep.prepare(self.config, destination)
        self.assertFalse(destination.exists())

        tiny_disk = mock.patch.object(prep.shutil, "disk_usage", return_value=SimpleNamespace(free=0))
        with self.assertRaisesRegex(ValueError, "insufficient free space"), self.offline(), tiny_disk:
            prep.prepare(self.config, destination)
        self.assertFalse(destination.exists())

        config = json.loads(self.config.read_text())
        config["resources"]["scratch_gib"] = 2
        self.config.write_text(json.dumps(config))
        with self.assertRaisesRegex(ValueError, "scratch budget"), self.offline():
            prep.prepare(self.config, destination)
        self.assertFalse(destination.exists())

    def test_verify_rejects_a_held_out_genome_in_its_catalogue(self):
        # Holdouts must be physically absent from each index input set
        # (docs/pilot-interfaces.md).
        destination = self.prepared()
        directory = destination / "catalogues" / "species_holdout"
        holdout = next(r for r in read_tsv(destination / "catalogues/all_known/catalogue.tsv")
                       if r["reference_id"] == HOLDOUT)
        with (directory / "catalogue.tsv").open("a") as stream:
            stream.write("\t".join(holdout.values()) + "\n")
        with (directory / "references.txt").open("a") as stream:
            stream.write(holdout["fasta"] + "\n")
        with self.assertRaisesRegex(ValueError, f"held-out references \\['{HOLDOUT}'\\]"):
            self.verify(destination)

    def test_verify_rejects_catalogues_that_drift_from_the_design(self):
        def tamper(edit):
            destination = self.prepared(f"dataset-{len(list(self.root.iterdir()))}")
            edit(destination / "catalogues")
            return destination

        def drop_known(root):
            path = root / "species_holdout" / "catalogue.tsv"
            lines = path.read_text().splitlines(keepends=True)
            path.write_text("".join(line for line in lines if not line.startswith(KNOWN)))

        def change_taxid(root):
            path = root / "all_known" / "catalogue.tsv"
            path.write_text(path.read_text().replace("\t5620\t562\t", "\t5621\t562\t"))

        def reorder_list(root):
            path = root / "all_known" / "references.txt"
            path.write_text("".join(reversed(path.read_text().splitlines(keepends=True))))

        for label, edit, message in [
            ("known reference dropped", drop_known, "catalogue.tsv differs"),
            ("taxid changed", change_taxid, "catalogue.tsv differs"),
            ("index list reordered", reorder_list, "references.txt differs"),
            ("unrecorded case", lambda root: (root / "unrecorded_case").mkdir(), "catalogue cases"),
            ("case removed", lambda root: shutil.rmtree(root / "species_holdout"), "catalogue cases"),
        ]:
            with self.subTest(label):
                with self.assertRaisesRegex(ValueError, message):
                    self.verify(tamper(edit))

    def test_verify_checks_samples_and_references_against_the_design(self):
        destination = self.prepared()
        path = destination / "dataset.json"
        original = json.loads(path.read_text())
        for field, keep, message in [
            ("samples", lambda s: s["sample_id"] != "B1", "manifest samples"),
            ("references", lambda g: g["genome_id"] != DISTRACTOR, "manifest references"),
        ]:
            with self.subTest(field=field):
                manifest = json.loads(json.dumps(original))
                manifest[field] = [item for item in manifest[field] if keep(item)]
                path.write_text(json.dumps(manifest))
                with self.assertRaisesRegex(ValueError, message):
                    self.verify(destination)
        duplicated = json.loads(json.dumps(original))
        duplicated["samples"].append(duplicated["samples"][0])
        path.write_text(json.dumps(duplicated))
        with self.assertRaisesRegex(ValueError, "manifest samples"):
            self.verify(destination)
        path.write_text(json.dumps(original))
        self.verify(destination)


class InputRobustnessTests(unittest.TestCase):
    def test_nameless_or_headerless_records_are_clean_errors(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "nameless.fq").write_text("@\nACGT\n+\nIIII\n")
            with self.assertRaisesRegex(ValueError, "without a name"):
                list(prep.fastq_records(root / "nameless.fq"))
            (root / "nameless.fna").write_text(">\nACGT\n")
            with self.assertRaisesRegex(ValueError, "without a name"):
                list(prep.fasta_records(root / "nameless.fna"))
            (root / "headerless.fna").write_text("ACGT\n>r1\nACGT\n")
            with self.assertRaisesRegex(ValueError, "before the first header"):
                list(prep.fasta_records(root / "headerless.fna"))
            (root / "leading-blank.fna").write_text("\n>r1 description\nac\ngt\n")
            self.assertEqual(list(prep.fasta_records(root / "leading-blank.fna")), [("r1", "ACGT")])

    def test_every_resource_limit_must_be_positive_before_any_download(self):
        pilot = ROOT / "experiments" / "pilot"
        config = json.loads((pilot / "pilot.json").read_text())
        config["genomes"] = str(pilot / config["genomes"])
        config["composition"] = str(pilot / config["composition"])
        for key in config["resources"]:
            with self.subTest(key), tempfile.TemporaryDirectory() as temp:
                changed = json.loads(json.dumps(config))
                changed["resources"][key] = 0
                path = Path(temp) / "pilot.json"
                path.write_text(json.dumps(changed))
                with self.assertRaises(ValueError):
                    prep.load_design(path)


class DownloadTests(unittest.TestCase):
    def test_download_streams_with_a_byte_limit_without_network(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            source = root / "source.fna.gz"
            source.write_bytes(os.urandom(4096))
            target = root / "copy.fna.gz"
            self.assertEqual(prep.download(source.as_uri(), target, 4096), 4096)
            self.assertEqual(target.read_bytes(), source.read_bytes())
            with self.assertRaisesRegex(ValueError, "byte limit"):
                prep.download(source.as_uri(), root / "limited.fna.gz", 1024)

    def test_source_urls_follow_the_ncbi_assembly_layout(self):
        genome = {"genome_id": "GCF_000005845.2", "assembly_name": "ASM584v2"}
        self.assertEqual(
            prep.source_url(genome),
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/"
            "GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz",
        )


class CommandLineTests(unittest.TestCase):
    def run_script(self, *args):
        return subprocess.run([sys.executable, str(ROOT / "scripts/prepare_pilot.py"), *args],
                              capture_output=True, text=True, timeout=60)

    def test_validate_action_checks_the_frozen_design_offline(self):
        result = self.run_script("validate")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("9 references, 8 samples", result.stdout)

    def test_verify_of_missing_dataset_is_a_clean_error(self):
        with tempfile.TemporaryDirectory() as temp:
            result = self.run_script("verify", "--output-dir", str(Path(temp) / "absent"))
        self.assertEqual(result.returncode, 1)
        self.assertTrue(result.stderr.startswith("error:"), result.stderr)
        self.assertNotIn("Traceback", result.stderr)


if __name__ == "__main__":
    unittest.main()
