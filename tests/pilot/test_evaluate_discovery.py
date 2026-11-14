"""The discovery evaluator on a tiny hand-built `discover` run directory."""

import gzip
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("evaluate_discovery", ROOT / "scripts" / "evaluate_discovery.py")
evaluator = importlib.util.module_from_spec(spec)
spec.loader.exec_module(evaluator)

HEADER = "pair_index\tread_name\toutcome\treason\treference_id\tspecies_taxid\tmapq1\tmapq2\n"


def make_run(root, rows, status="succeeded"):
    """rows: (read_name, outcome). Catalogue holds KNOWN only; HELD is held out."""
    root = Path(root)
    (root / "run" / "stages" / "map-S1").mkdir(parents=True)
    catalogue = root / "case" / "catalogue.tsv"
    catalogue.parent.mkdir()
    catalogue.write_text("reference_id\tfasta\ttaxid\tspecies_taxid\tsha256\nKNOWN\tk.fna\t1\t1\t" + "0" * 64 + "\n")
    run = {
        "status": status,
        "configuration": {
            "profiler": "sylph",
            "tool": {"version": "sylph 1.0.0"},
            "mapper": {"name": "strobealign", "tool": {"version": "strobealign 0.20.0"}},
            "catalogue": {"path": str(catalogue), "sha256": "a" * 64},
            "assignment_policy": {"min_identity": 0.98},
        },
        "stages": [
            {"id": "profile-S1", "admitted_ms": 0, "finished_ms": 1000},
            {"id": "map-S1", "admitted_ms": 1000, "finished_ms": 4000},
        ],
    }
    (root / "run" / "run.json").write_text(json.dumps(run))
    stage = root / "run" / "stages" / "map-S1"
    with gzip.open(stage / "decisions.tsv.gz", "wt") as stream:
        stream.write(HEADER)
        for index, (name, outcome) in enumerate(rows, 1):
            stream.write(f"{index}\t{name}\t{outcome}\tr\t\t\t\t\n")
    (stage / "extraction.json").write_text(json.dumps({"reasons": {"unique": 3, "mate-unmapped": 2}}))
    return root / "run"


ROWS = [
    ("S1|KNOWN|000000001", "reference-assigned"),
    ("S1|KNOWN|000000002", "reference-assigned"),
    ("S1|KNOWN|000000003", "discovery-unassigned"),
    ("S1|HELD|000000004", "discovery-ambiguous"),
    ("S1|HELD|000000005", "reference-assigned"),
]


class EvaluateDiscoveryTests(unittest.TestCase):
    def test_false_removal_contamination_and_reduction(self):
        with tempfile.TemporaryDirectory() as temp:
            result = evaluator.evaluate(make_run(temp, ROWS), {"KNOWN": "known", "HELD": "species_holdout_A"})
        self.assertEqual(result["false_removal"], {"HELD": 0.5})
        self.assertEqual(result["residual_contamination"], {"KNOWN": round(1 / 3, 6)})
        self.assertEqual(result["totals"]["pairs"], 5)
        self.assertEqual(result["totals"]["retained"], 2)
        self.assertEqual(result["totals"]["input_pair_reduction"], 0.6)
        self.assertEqual(result["totals"]["pool_heldout_fraction"], 0.5)
        self.assertEqual(result["by_genome"]["HELD"]["role"], "species_holdout_A")
        self.assertEqual(result["catalogue"], "case")
        self.assertEqual(result["stage_seconds"], {"map": 3.0, "profile": 1.0})
        self.assertEqual(result["tool_versions"]["mapper"], "strobealign 0.20.0")

    def test_incomplete_runs_and_foreign_read_names_are_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            with self.assertRaisesRegex(ValueError, "not succeeded"):
                evaluator.evaluate(make_run(temp, ROWS, status="failed"), {})
        with tempfile.TemporaryDirectory() as temp:
            with self.assertRaisesRegex(ValueError, "pilot sample"):
                evaluator.evaluate(make_run(temp, [("plain-read-name", "reference-assigned")]), {})

    def test_command_line_refuses_to_overwrite_results(self):
        with tempfile.TemporaryDirectory() as temp:
            run = make_run(temp, ROWS)
            output = Path(temp) / "results.json"
            output.write_text("{}")
            result = subprocess.run(
                [sys.executable, str(ROOT / "scripts/evaluate_discovery.py"), str(run), "--output", str(output)],
                capture_output=True, text=True, timeout=60)
            self.assertNotEqual(result.returncode, 0)
            self.assertEqual(output.read_text(), "{}")


if __name__ == "__main__":
    unittest.main()
