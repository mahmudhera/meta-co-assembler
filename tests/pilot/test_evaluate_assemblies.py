"""Parsing and gate logic of the assembly evaluator (MetaQUAST is not run)."""

import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("evaluate_assemblies", ROOT / "scripts" / "evaluate_assemblies.py")
evaluator = importlib.util.module_from_spec(spec)
spec.loader.exec_module(evaluator)

REPORT = """Assembly\tdiscovery\tall-reads
# contigs (>= 1000 bp)\t12\t40
Total length (>= 1000 bp)\t4100000\t30000000
N50\t250000\t180000
# misassemblies\t1\t3
Unaligned length\t0\t1500
Genome fraction (%)\t11.5\t97.25
"""

FRACTION = """Assemblies\tdiscovery\tall-reads
GCF_heldout\t96.000\t98.500
GCF_known\t0.500\t99.000
GCF_absent\t-\t-
"""


class EvaluateAssembliesTests(unittest.TestCase):
    def test_reports_and_genome_fractions_parse_per_assembly(self):
        report = evaluator.parse_report(REPORT)
        self.assertEqual(report["discovery"]["contigs_at_least_1000"], 12)
        self.assertEqual(report["all-reads"]["combined_genome_fraction_percent"], 97.25)
        self.assertIsNone(report["discovery"]["indels_per_100kbp"], "missing metrics are null")
        fraction = evaluator.parse_genome_fraction(FRACTION)
        self.assertEqual(fraction["discovery"]["GCF_heldout"], 96)
        self.assertIsNone(fraction["all-reads"]["GCF_absent"], "QUAST '-' means not aligned")

    def test_gate_requires_reduction_and_bounded_loss_for_every_heldout_genome(self):
        fraction = evaluator.parse_genome_fraction(FRACTION)
        result = evaluator.gate(400, 1000, fraction["discovery"], fraction["all-reads"], ["GCF_heldout"])
        self.assertEqual(result["input_pair_reduction"], 0.6)
        self.assertEqual(result["heldout_recovery_loss_pp"], {"GCF_heldout": 2.5})
        self.assertTrue(result["passes_reduction"] and result["passes_recovery"])
        worse = evaluator.gate(600, 1000, {"G": 90.0}, {"G": 96.0}, ["G"])
        self.assertFalse(worse["passes_reduction"])
        self.assertFalse(worse["passes_recovery"], "a 6-point loss exceeds 5 points")
        unmeasured = evaluator.gate(100, 1000, {}, {}, ["G"])
        self.assertFalse(unmeasured["passes_recovery"], "a held-out genome without a measurement cannot pass")


def stage(stage_id, seconds, peak=1000):
    return {"id": stage_id, "status": "complete", "admitted_ms": 0, "finished_ms": seconds * 1000,
            "peak_rss_kib": peak, "scratch_peak_bytes": 10}


def write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def iterate_fixture(root):
    """A coassembly of S1+S2 (round 1) and a two-round iterate run regrouped per sample."""
    coassembly = root / "coassembly"
    write(coassembly / "run.json", json.dumps({
        "status": "succeeded",
        "configuration": {"command": "coassemble", "source": {"kind": "discovery-pools"}, "grouping": {"groups": [
            {"group_id": "group-000001", "samples": ["S1", "S2"], "pairs": 100, "bases": 30000}]}},
        "stages": [stage("assemble-group-000001", 20)],
    }))
    write(coassembly / "stages/assemble-group-000001/assembly.json", json.dumps({"status": "contigs"}))
    write(coassembly / "stages/assemble-group-000001/group.json", json.dumps({"group": {"members": [
        {"sample_id": "S1", "pairs": 60}, {"sample_id": "S2", "pairs": 40}]}}))
    run = root / "iterate"
    groups = [{"group_id": "group-000001", "samples": ["S1"]}, {"group_id": "group-000002", "samples": ["S2"]}]
    write(run / "run.json", json.dumps({
        "status": "succeeded",
        "configuration": {
            "command": "iterate", "max_rounds": 2, "grouping": {"groups": groups},
            "inputs": [{"index": 1, "path": str(coassembly), "groups": [
                {"group_id": "group-000001", "samples": ["S1", "S2"]}]}],
        },
        "stages": [stage("accept-r01", 1), stage("check-r02-group-000001", 2), stage("check-r02-group-000002", 2),
                   stage("assemble-r02-group-000001", 5, peak=4000), stage("assemble-r02-group-000002", 0),
                   stage("accept-r02", 1)],
    }))
    write(run / "stages/assemble-r02-group-000001/assembly.json",
          json.dumps({"status": "contigs", "input": {"pairs": 10, "bases": 3000}}))
    write(run / "stages/check-r02-group-000001/round.json",
          json.dumps({"members": [{"sample_id": "S1", "counts": {"pairs_retained": 10}}]}))
    write(run / "stages/assemble-r02-group-000002/assembly.json",
          json.dumps({"status": "skipped", "input": {"pairs": 0, "bases": 0}}))
    write(run / "stages/accept-r01/cumulative.fa", ">r01.run01.group-000001__k1 x\nACGT\n")
    write(run / "stages/accept-r02/cumulative.fa",
          ">r01.run01.group-000001__k1 x\nACGT\n>r02.group-000001__k2\nGG\nCC\n")
    return run


class IterateArmTests(unittest.TestCase):
    def test_iterate_arms_count_every_round_and_trace_contigs_to_groups(self):
        with tempfile.TemporaryDirectory() as tmp:
            run = iterate_fixture(Path(tmp))
            arm = evaluator.load_arm(str(run))
            self.assertEqual(arm["command"], "iterate")
            self.assertEqual(arm["rounds_evaluated"], 2)
            self.assertEqual([(a["round"], a["pairs"]) for a in arm["assemblies"]], [(1, 100), (2, 10)],
                             "skipped assemblies add nothing")
            self.assertEqual([a["pairs_by_sample"] for a in arm["assemblies"]], [{"S1": 60, "S2": 40}, {"S1": 10}])
            self.assertEqual(arm["check_and_accept_seconds"], 6.0)
            origins = {tuple(s["samples"]): s["ids"] for s in arm["contig_sets"]}
            self.assertEqual(origins, {("S1", "S2"): {"r01.run01.group-000001__k1"},
                                       ("S1",): {"r02.group-000001__k2"}})
            first = evaluator.load_arm(f"{run}@1")
            self.assertEqual(len(first["assemblies"]), 1)
            self.assertEqual(first["contig_sets"][0]["path"].parent.name, "accept-r01")
            with self.assertRaises(ValueError):
                evaluator.load_arm(f"{run}@3")
            with self.assertRaises(ValueError):
                evaluator.load_arm(f"{Path(tmp) / 'coassembly'}@1")
            out = io.StringIO()
            evaluator.write_contigs(out, run / "stages/accept-r02/cumulative.fa", {"r02.group-000001__k2"})
            self.assertEqual(out.getvalue(), ">r02.group-000001__k2\nGG\nCC\n")

    def test_heldout_genomes_exclude_species_with_a_catalogued_relative(self):
        catalogue = [{"reference_id": "ecoli_mg1655", "species_taxid": "562"}]
        species = {"ecoli_mg1655": "562", "ecoli_w3110": "562", "bsub": "1423"}
        references = ["bsub", "ecoli_mg1655", "ecoli_w3110"]
        self.assertEqual(evaluator.heldout_genomes(references, catalogue, species), (["bsub"], ["ecoli_w3110"]))
        self.assertEqual(evaluator.heldout_genomes(references, catalogue, None), (["bsub", "ecoli_w3110"], []))

    def test_contigs_need_a_known_origin_and_pooled_inputs_split_by_sample(self):
        with self.assertRaises(ValueError):
            evaluator.origin_samples("k141_1", {}, {})
        communities = {"S1": "A", "S2": "A", "S3": "B"}
        self.assertEqual(evaluator.pairs_by_community({"S1": 5, "S2": 3, "S3": 7}, communities), {"A": 8, "B": 7})


if __name__ == "__main__":
    unittest.main()
