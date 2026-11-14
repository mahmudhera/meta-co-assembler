"""Contig assignment and interval logic of the bin evaluator (minimap2 is not run)."""

import importlib.util
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("evaluate_bins", ROOT / "scripts" / "evaluate_bins.py")
evaluator = importlib.util.module_from_spec(spec)
spec.loader.exec_module(evaluator)


def paf(query, qlen, target, start, end, matches, kind="P"):
    return "\t".join(map(str, [query, qlen, 0, qlen, "+", target, 10_000, start, end, matches, qlen, 60,
                                f"tp:A:{kind}"]))


class EvaluateBinsTests(unittest.TestCase):
    def test_contigs_go_to_the_genome_with_most_matches_over_half_their_length(self):
        text = "\n".join([
            paf("c1", 1000, "G1|chr", 0, 600, 590),
            paf("c1", 1000, "G2|chr", 0, 1000, 940),   # congener: more matches overall
            paf("c1", 1000, "G1|chr", 700, 1100, 395),  # G1 total 985 wins
            paf("c2", 1000, "G2|chr", 0, 400, 400),     # under half of the contig
            paf("c3", 1000, "G3|chr", 0, 1000, 990),
            paf("c3", 1000, "G4|chr", 0, 1000, 995, "S"),  # secondary: interval only
        ])
        assigned = evaluator.assign_contigs(text)
        self.assertEqual(assigned["c1"][0], "G1")
        self.assertEqual(assigned["c1"][2]["G1"], [(0, 600), (700, 1100)])
        self.assertIsNone(assigned["c2"][0], "400 of 1,000 bases is not enough")
        self.assertEqual(assigned["c3"][0], "G3", "secondary matches do not assign")
        self.assertEqual(assigned["c3"][2]["G4"], [(0, 1000)])

    def test_intervals_merge_and_tiers_match_the_pipeline(self):
        self.assertEqual(evaluator.merged_length([(0, 10), (5, 20), (30, 40), (35, 38)]), 30)
        self.assertEqual(evaluator.merged_length([]), 0)
        self.assertEqual(evaluator.tier(90, 4.9), "high")
        self.assertEqual(evaluator.tier(50, 9.9), "medium")
        self.assertEqual(evaluator.tier(95, 12), "low")


if __name__ == "__main__":
    unittest.main()
