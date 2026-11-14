import importlib.util
import json
from pathlib import Path
import struct
import tempfile
import unittest
import zipfile


ROOT = Path(__file__).resolve().parents[2]


def module(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


prep = module("prepare_pilot")
regions = module("inspect_locove")


class PreparationTests(unittest.TestCase):
    def test_frozen_design_has_exact_budgets_and_no_distractor_reads(self):
        config, genomes, samples = prep.load_design(prep.DEFAULT_CONFIG)
        self.assertEqual(len(genomes), 9)
        self.assertEqual(len(samples), 8)
        self.assertEqual({s["group"] for s in samples.values()}, {"community_A", "community_B"})
        self.assertEqual(sum(sum(s["counts"].values()) for s in samples.values()), 8_000_000)
        self.assertEqual(config["resources"]["threads"], 16)
        self.assertEqual(config["resources"]["memory_gib"], 64)

    def test_record_allocation_conserves_pairs_with_deterministic_ties(self):
        self.assertEqual(prep.allocate(11, [1, 1, 1]), [4, 4, 3])
        self.assertEqual(prep.allocate(0, [1, 2]), [0, 0])
        self.assertEqual(sum(prep.allocate(12345, [6260361, 33038])), 12345)
        with self.assertRaises(ValueError):
            prep.allocate(5, [0])

    def test_seeds_are_reproducible_and_sample_specific(self):
        one = prep.derived_seed(42, "A01", "genome", 0)
        self.assertEqual(one, prep.derived_seed(42, "A01", "genome", 0))
        self.assertNotEqual(one, prep.derived_seed(42, "B01", "genome", 0))
        self.assertGreater(one, 0)
        self.assertLess(one, 2147483647)

    def test_mate_mismatch_and_truncated_fastq_are_errors(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "one.fq").write_text("@one/1\nACGT\n+\nIIII\n")
            (root / "two.fq").write_text("@two/2\nACGT\n+\nIIII\n")
            with self.assertRaisesRegex(ValueError, "mismatched"):
                list(prep.paired_records(root / "one.fq", root / "two.fq"))
            (root / "one.fq").write_text("@one/1\nACGT\n+\n")
            with self.assertRaisesRegex(ValueError, "malformed"):
                list(prep.fastq_records(root / "one.fq"))

    def test_mixture_preserves_pairs_truth_and_reproducible_compression(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            parts = []
            for genome in ["genome_a", "genome_b"]:
                prefix = root / (genome + "-")
                for mate in [1, 2]:
                    Path(str(prefix) + f"{mate}.fq").write_text(
                        "".join(f"@raw{i}/{mate}\nACGT\n+\nIIII\n" for i in range(3))
                    )
                parts.append({"prefix": prefix, "pairs": 3, "genome": genome})
            outputs = []
            for run in [1, 2]:
                one, two = root / f"{run}-1.fastq.gz", root / f"{run}-2.fastq.gz"
                counts = prep.mix_parts(parts, one, two, "A01", 4, 42)
                self.assertEqual(counts, {"genome_a": 3, "genome_b": 3})
                pairs = list(prep.paired_records(one, two))
                self.assertEqual(len(pairs), 6)
                self.assertEqual(len({p[0][0] for p in pairs}), 6)
                self.assertTrue(all(p[0][0].startswith("A01|") for p in pairs))
                outputs.append((one.read_bytes(), two.read_bytes()))
            self.assertEqual(outputs[0], outputs[1])
            self.assertEqual(outputs[0][0][4:8], bytes(4))

    def test_existing_output_is_never_overwritten(self):
        with tempfile.TemporaryDirectory() as temp:
            with self.assertRaisesRegex(ValueError, "output exists"):
                prep.prepare(prep.DEFAULT_CONFIG, Path(temp))

    def test_verified_art_surplus_is_trimmed_to_the_exact_target(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            prefix = root / "part-"
            for mate in [1, 2]:
                Path(str(prefix) + f"{mate}.fq").write_text(
                    "".join(f"@raw{i}/{mate}\nACGT\n+\nIIII\n" for i in range(3))
                )
            part = {"prefix": prefix, "pairs": 2, "generated_pairs": 3, "genome": "genome_a"}
            one, two = root / "one.fastq.gz", root / "two.fastq.gz"
            counts = prep.mix_parts([part], one, two, "A01", 4, 42)
            self.assertEqual(counts, {"genome_a": 2})
            self.assertEqual(len(list(prep.paired_records(one, two))), 2)


class RegionTests(unittest.TestCase):
    def test_half_open_union_does_not_double_count_overlap_or_flanks(self):
        self.assertEqual(regions.merged_bases([(10, 20), (0, 15), (20, 25), (40, 50)]), 35)
        self.assertEqual(regions.merged_bases([]), 0)
        with self.assertRaises(ValueError):
            regions.merged_bases([(5, 2)])

    def test_unknown_segments_are_retained_in_analysis(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "record-map.json").write_text(json.dumps({
                "genomes": ["genome"], "records": ["record"], "rec_genome": [0]
            }))
            arrays = {
                "seg_start": ("<i8", [0, 100, 200]),
                "seg_end": ("<i8", [100, 200, 300]),
                "seg_record": ("<i8", [0, 0, 0]),
                "n_anchors": ("<i8", [10, 0, 10]),
                "breadth_anchor": ("<f8", [0.0, float("nan"), 1.0]),
            }
            with zipfile.ZipFile(root / "segments.npz", "w") as archive:
                for name, (dtype, values) in arrays.items():
                    header = repr({"descr": dtype, "fortran_order": False, "shape": (len(values),)}).encode()
                    header += b" " * ((64 - (10 + len(header) + 1) % 64) % 64) + b"\n"
                    data = b"\x93NUMPY\x01\x00" + len(header).to_bytes(2, "little") + header
                    data += b"".join(struct.pack("<q" if dtype == "<i8" else "<d", value) for value in values)
                    archive.writestr(name + ".npy", data)
            result = regions.inspect(root)
            self.assertEqual(result["experiments"][0]["candidate_bases"], 200)
            self.assertEqual(result["experiments"][1]["candidate_bases"], 300)
            self.assertFalse(result["per_read_assignments"])


if __name__ == "__main__":
    unittest.main()
