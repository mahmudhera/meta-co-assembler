"""Parsing and gate logic of the assembly evaluator (MetaQUAST is not run)."""

import importlib.util
from pathlib import Path
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


if __name__ == "__main__":
    unittest.main()
