"""Consistency of the recorded M1 decisions, frozen design, and analysis helpers.

These checks read only repository files. They fail when an owner-confirmed
choice changes in one place but not in the documents or records that cite it.
"""

import csv
import importlib.util
import json
import math
from pathlib import Path
import re
import struct
import subprocess
import sys
import tempfile
import unittest
import zipfile


ROOT = Path(__file__).resolve().parents[2]
PILOT = ROOT / "experiments" / "pilot"


def module(name):
    spec = importlib.util.spec_from_file_location(f"{name}_contracts", ROOT / "scripts" / f"{name}.py")
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


regions = module("inspect_locove")


def text(relative):
    return (ROOT / relative).read_text(encoding="utf-8")


def load_json(relative):
    return json.loads(text(relative))


def table(relative):
    with (ROOT / relative).open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def environment_pins():
    pins = {}
    for line in text("environment.yml").splitlines():
        match = re.fullmatch(r"\s*-\s*([A-Za-z0-9_.-]+)=([^\s#=]+)\s*", line)
        if match:
            pins[match.group(1)] = match.group(2)
    return pins


class RecordedToolchainTests(unittest.TestCase):
    def test_environment_pins_match_recorded_tool_versions(self):
        pins = environment_pins()
        recorded = load_json("experiments/pilot/tool-versions.json")["tools"]
        released = {name: info["version"] for name, info in recorded.items() if "source_revision" not in info}
        self.assertEqual(pins, released)

    def test_decision_and_readme_cite_the_pinned_versions(self):
        pins = environment_pins()
        decision = text("docs/decisions/0001-pilot.md")
        display = {
            "sylph": "Sylph", "kraken2": "Kraken2", "strobealign": "Strobealign",
            "minimap2": "Minimap2", "samtools": "Samtools", "megahit": "MEGAHIT",
            "spades": "SPAdes", "art": "ART", "quast": "QUAST",
        }
        for package, name in display.items():
            self.assertIn(f"{name} {pins[package]}", decision, package)
        self.assertIn(f"Rust {pins['rust']}", text("README.md"))

    def test_locove_source_revision_is_recorded_consistently(self):
        source = load_json("tools/source-tools.json")["locove"]
        revision = source["revision"]
        self.assertRegex(revision, r"\A[0-9a-f]{40}\Z", "pin a full commit, not a branch or tag")
        self.assertEqual(source["features"], [])
        self.assertIn(f"--rev {revision}", text("experiments/pilot/README.md"))
        self.assertIn(revision, text("docs/research/locove.md"))
        self.assertIn(f"LoCovE {source['package_version']}", text("docs/decisions/0001-pilot.md"))
        tool = load_json("experiments/pilot/tool-versions.json")["tools"]["locove"]
        self.assertEqual((tool["source_revision"], tool["version"]), (revision, source["package_version"]))
        self.assertEqual(load_json("experiments/pilot/verification.json")["locove"]["source_revision"], revision)
        self.assertEqual(load_json("experiments/pilot/locove-demo-benchmark.json")["source_revision"], revision)


class FrozenDesignTests(unittest.TestCase):
    def setUp(self):
        self.config = load_json("experiments/pilot/pilot.json")
        self.roles = {row["genome_id"]: row["role"] for row in table("experiments/pilot/genomes.tsv")}
        self.composition = table("experiments/pilot/composition.tsv")

    def with_role(self, prefix):
        return {genome for genome, role in self.roles.items() if role.startswith(prefix)}

    def test_catalogue_cases_exclude_exactly_the_matching_holdouts(self):
        species, strain = self.with_role("species_holdout"), self.with_role("strain_holdout")
        self.assertEqual(len(species), 2)
        self.assertEqual(len(strain), 1)
        expected = {
            "all_known": set(),
            "species_holdout": species,
            "strain_holdout": strain,
            "combined_holdout": species | strain,
        }
        self.assertEqual({k: set(v) for k, v in self.config["catalogues"].items()}, expected)
        catalogue_only = self.with_role("distractor")
        self.assertTrue(catalogue_only)
        for case, excluded in self.config["catalogues"].items():
            self.assertFalse(catalogue_only & set(excluded), f"{case} must keep distractors")

    def test_species_holdouts_are_confined_to_their_own_community(self):
        communities = {}
        for row in self.composition:
            communities.setdefault(row["community"], set()).add(row["sample_id"])
        self.assertEqual(sorted(map(len, communities.values())), [4, 4])
        for genome in self.with_role("species_holdout"):
            suffix = self.roles[genome].rsplit("_", 1)[1]
            present = {row["sample_id"] for row in self.composition if row["genome_id"] == genome}
            self.assertEqual(present, communities[f"community_{suffix}"], genome)

    def test_documented_composition_fractions_hold_for_every_sample(self):
        total = self.config["read_pairs_per_sample"]
        by_sample = {}
        for row in self.composition:
            by_sample.setdefault(row["sample_id"], {})[self.roles[row["genome_id"]]] = int(row["pairs"]) / total
        for sample, fractions in by_sample.items():
            species = [v for role, v in fractions.items() if role.startswith("species_holdout")]
            self.assertEqual(len(species), 1, sample)
            self.assertTrue(0.15 <= species[0] <= 0.18, sample)
            self.assertTrue(0.01 <= fractions["strain_holdout"] <= 0.04, sample)
            self.assertIn("known_rare", fractions, sample)

    def test_owner_confirmed_thresholds_match_the_decision_record(self):
        policies, decision = self.config["policies"], text("docs/decisions/0001-pilot.md")
        flat = " ".join(decision.split())
        assignment, recruitment = policies["reference_assignment"], policies["recruitment"]
        evaluation = policies["evaluation"]
        self.assertIn(f"ANI at least {policies['sylph_min_containment_ani_percent']}%", flat)
        self.assertIn(
            f"MAPQ {assignment['min_mapq']}, identity {round(assignment['min_identity'] * 100)}%, "
            f"aligned fraction {round(assignment['min_aligned_fraction'] * 100)}%", flat)
        self.assertIn(
            f"identity {round(recruitment['min_identity'] * 100)}% and aligned fraction "
            f"{round(recruitment['min_aligned_fraction'] * 100)}%", flat)
        self.assertEqual(recruitment["max_rounds"], 3)
        self.assertIn("three-round maximum", flat)
        self.assertIn(f"at least {round(evaluation['min_assembly_input_pair_reduction'] * 100)}% fewer", flat)
        self.assertEqual(evaluation["max_heldout_species_recovery_loss_percentage_points"], 5)
        self.assertIn("five percentage points", flat)
        self.assertIn(f"{evaluation['min_contig_length']:,} bp contig cutoff", " ".join(text("experiments/pilot/README.md").split()))

    def test_conservative_retention_and_locove_limits_cannot_be_silently_relaxed(self):
        policies = self.config["policies"]
        assignment = policies["reference_assignment"]
        self.assertTrue(assignment["require_both_mates"])
        self.assertTrue(assignment["require_same_genome"])
        self.assertTrue(assignment["retain_ambiguous_pairs"])
        self.assertIsNone(policies["additional_abundance_floor"])
        locove = policies["locove"]
        self.assertFalse(locove["read_assignments_available"])
        self.assertTrue(locove["coverage_is_pooled"])
        self.assertTrue(locove["undefined_segments"].startswith("unknown"))
        self.assertIn("not validated", policies["status"])


def write_npz(path, arrays):
    """Write NPY v1 members the way LoCovE does: 1-D, little-endian, C order."""
    with zipfile.ZipFile(path, "w") as archive:
        for name, (dtype, values, *shape) in arrays.items():
            shape = shape[0] if shape else (len(values),)
            header = repr({"descr": dtype, "fortran_order": False, "shape": shape}).encode()
            header += b" " * ((64 - (10 + len(header) + 1) % 64) % 64) + b"\n"
            code = {"<i8": "<q", "<f8": "<d", "<i4": "<i"}[dtype]
            data = b"\x93NUMPY\x01\x00" + len(header).to_bytes(2, "little") + header
            archive.writestr(name + ".npy", data + b"".join(struct.pack(code, v) for v in values))


def coverage(directory, records, segments):
    """segments: (start, end, record, anchors, breadth) rows."""
    directory = Path(directory)
    (directory / "record-map.json").write_text(json.dumps({
        "genomes": ["genome"], "records": records, "rec_genome": [0] * len(records)}))
    columns = list(zip(*segments)) if segments else [[]] * 5
    write_npz(directory / "segments.npz", {
        "seg_start": ("<i8", columns[0]), "seg_end": ("<i8", columns[1]),
        "seg_record": ("<i8", columns[2]), "n_anchors": ("<i8", columns[3]),
        "breadth_anchor": ("<f8", columns[4]),
    })


def experiment(result, threshold, flank):
    return next(e for e in result["experiments"]
                if e["breadth_strictly_greater_than"] == threshold and e["flank_bases"] == flank)


class LocoveWindowTests(unittest.TestCase):
    def test_flanks_are_clipped_per_record_and_records_are_never_merged(self):
        with tempfile.TemporaryDirectory() as temp:
            coverage(temp, ["r0", "r1"], [
                (0, 100, 0, 5, 0.0), (100, 200, 0, 5, 1.0),
                (0, 100, 1, 5, 1.0), (100, 300, 1, 5, 0.0),
            ])
            result = regions.inspect(temp)
        self.assertEqual(result["reference_bases"], 500)
        self.assertEqual(experiment(result, 0.5, 0)["candidate_bases"], 200)
        clipped = experiment(result, 0.5, 500)
        self.assertEqual(clipped["candidate_bases"], 500)
        self.assertEqual(clipped["fraction_of_reference"], 1.0)
        self.assertTrue(all(e["retain_unknown"] for e in result["experiments"]))

    def test_nan_breadth_with_anchors_is_unknown_not_absent(self):
        with tempfile.TemporaryDirectory() as temp:
            coverage(temp, ["r0"], [(0, 100, 0, 3, math.nan), (100, 200, 0, 3, 0.0)])
            result = regions.inspect(temp)
        self.assertEqual(experiment(result, 0.95, 0)["candidate_bases"], 100)

    def test_malformed_arrays_and_coordinates_are_rejected(self):
        for name, arrays in [
            ("dtype", {"seg_start": ("<i4", [0])}),
            ("shape", {"seg_start": ("<i8", [0], (1, 1))}),
            ("length", {"seg_start": ("<i8", [0], (2,))}),
        ]:
            with self.subTest(name), tempfile.TemporaryDirectory() as temp:
                path = Path(temp) / "segments.npz"
                write_npz(path, arrays)
                with zipfile.ZipFile(path) as archive, self.assertRaises(ValueError):
                    regions.read_array(archive, "seg_start")
        for segment in [(10, 5, 0, 1, 1.0), (0, 10, 1, 1, 1.0), (-1, 10, 0, 1, 1.0)]:
            with self.subTest(segment=segment), tempfile.TemporaryDirectory() as temp:
                coverage(temp, ["r0"], [segment])
                with self.assertRaises(ValueError):
                    regions.inspect(temp)

    def test_existing_analysis_output_is_not_overwritten(self):
        with tempfile.TemporaryDirectory() as temp:
            coverage(temp, ["r0"], [(0, 100, 0, 1, 1.0)])
            output = Path(temp) / "windows.json"
            output.write_text("{}")
            result = subprocess.run(
                [sys.executable, str(ROOT / "scripts/inspect_locove.py"), temp, "--output", str(output)],
                capture_output=True, text=True, timeout=60)
            self.assertNotEqual(result.returncode, 0)
            self.assertEqual(output.read_text(), "{}")

    def test_records_without_segments_or_bases_are_unknown_not_zero(self):
        # Record lengths are inferred from segment ends, so a record without
        # segments has unknown length and coverage; it must not count as zero.
        for records, segments, message in [
            (["r0", "r1"], [(0, 100, 0, 5, 1.0)], "no segments, e.g. \\['r1'\\]"),
            (["r0"], [], "no segments"),
            ([], [], "no reference bases"),
            (["r0"], [(0, 0, 0, 0, math.nan)], "no reference bases"),
        ]:
            with self.subTest(records=records, segments=segments), tempfile.TemporaryDirectory() as temp:
                coverage(temp, records, segments)
                with self.assertRaisesRegex(ValueError, message):
                    regions.inspect(temp)


if __name__ == "__main__":
    unittest.main()
