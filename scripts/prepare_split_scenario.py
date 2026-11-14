#!/usr/bin/env python3
"""Supplementary pilot scenario for pooled iteration: one novel genome spread thinly over all samples.

Copies of the eight pilot samples receive extra read pairs simulated from one
genome (by default GCF_000008645.1, a distractor in the pilot catalogues)
with the pilot's ART settings and seed scheme, and a catalogue variant drops
that genome, so its species is outside the catalogue. At the default 6,000
pairs per sample (about 1x of a 1.75 Mbp genome), each four-sample group sees
about 4x and all samples together about 8x: thin per group, assemblable when
pooled. The extra reads are appended as a second gzip member, so every
original read stays byte-identical and in place; the pilot inputs are never
modified. Writes reads, samples.tsv, metadata.tsv, composition.tsv, the
catalogue variant, and scenario.json (commands and digests) to a new
directory.
"""

import argparse
import gzip
import importlib.util
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("prepare_pilot", ROOT / "scripts" / "prepare_pilot.py")
pilot = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pilot)


def simulate(config, sample, genome, records, pairs, temp, read1, read2):
    """ART parts for one sample, written as one gzip member per mate (exactly `pairs` pairs)."""
    simulation = config["simulation"]
    lengths = [length for _, length in records]
    parts = []
    for index, ((path, _), count) in enumerate(zip(records, pilot.allocate(pairs, lengths))):
        if count == 0:
            continue
        prefix = Path(temp) / f"{sample}-{genome}-{index}-"
        surplus = max(simulation["art_min_surplus_pairs"], int(count * simulation["art_surplus_fraction"]) + 1)
        command = [shutil.which("art_illumina"), "-ss", simulation["profile"], "-p",
                   "-l", str(simulation["read_length"]), "-m", str(simulation["fragment_mean"]),
                   "-s", str(simulation["fragment_sd"]),
                   "-rs", str(pilot.derived_seed(config["seed"], "split", sample, genome, index)),
                   "-c", str(count + surplus), "-na", "-i", str(path), "-o", str(prefix)]
        parts.append(pilot.simulate_part((command, prefix, count, genome)))
    ordinal = 0
    with open(read1, "xb") as raw1, open(read2, "xb") as raw2:
        writers = [gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0,
                                 compresslevel=simulation["gzip_level"]) for raw in (raw1, raw2)]
        for part in parts:
            stream = pilot.paired_records(str(part["prefix"]) + "1.fq", str(part["prefix"]) + "2.fq")
            for _ in range(part["pairs"]):
                ordinal += 1
                for mate, ((_, sequence, quality), writer) in enumerate(zip(next(stream), writers), 1):
                    writer.write(f"@{sample}|{genome}|split{ordinal:09d}/{mate}\n{sequence}\n+\n{quality}\n".encode())
        for writer in writers:
            writer.close()
    return [part["command"] for part in parts]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pilot-dir", type=Path, default=Path("data/pilot"))
    parser.add_argument("--config", type=Path, default=pilot.DEFAULT_CONFIG)
    parser.add_argument("--output-dir", type=Path, default=Path("data/pilot-split"))
    parser.add_argument("--genome", default="GCF_000008645.1")
    parser.add_argument("--pairs-per-sample", type=int, default=6000)
    parser.add_argument("--base-catalogue", default="species_holdout")
    args = parser.parse_args()
    pilot_dir, out = args.pilot_dir.resolve(), args.output_dir.resolve()
    if out.exists():
        print(f"error: {out} exists; choose a fresh directory", file=sys.stderr)
        return 1
    if not shutil.which("art_illumina"):
        print("error: art_illumina is missing; use meta-co-assembler-env", file=sys.stderr)
        return 1
    config = json.loads(args.config.read_text(encoding="utf-8"))
    genome_rows = {row["genome_id"]: row for row in pilot.read_table(args.config.parent / "genomes.tsv")}
    genome = genome_rows[args.genome]
    record_dir = pilot_dir / "records" / args.genome
    records = []
    for path in sorted(record_dir.glob("*.fna")):
        length = sum(len(line.strip()) for line in path.read_text().splitlines() if not line.startswith(">"))
        records.append((path, length))
    base = pilot_dir / "catalogues" / args.base_catalogue
    catalogue_rows = pilot.read_table(base / "catalogue.tsv")
    if args.genome not in {row["reference_id"] for row in catalogue_rows}:
        print(f"error: {args.genome} is not in the {args.base_catalogue} catalogue", file=sys.stderr)
        return 1
    if any(row["species_taxid"] == genome["species_taxid"] and row["reference_id"] != args.genome
           for row in catalogue_rows):
        print(f"error: another catalogue genome shares the species of {args.genome}", file=sys.stderr)
        return 1
    (out / "reads").mkdir(parents=True)
    samples = pilot.read_table(pilot_dir / "samples.tsv")
    composition = pilot.read_table(args.config.parent / "composition.tsv")
    community = {row["sample_id"]: row["community"] for row in composition}
    reports, new_samples = [], []
    with tempfile.TemporaryDirectory(prefix="split-", dir=out) as temp:
        for row in samples:
            sample = row["sample_id"]
            extra = [Path(temp) / f"{sample}_extra_R{mate}.fastq.gz" for mate in (1, 2)]
            commands = simulate(config, sample, args.genome, records, args.pairs_per_sample, temp, *extra)
            outputs = []
            for mate, key in ((1, "read1"), (2, "read2")):
                target = out / "reads" / f"{sample}_R{mate}.fastq.gz"
                with target.open("xb") as stream:
                    for source in (pilot_dir / row[key], extra[mate - 1]):
                        with source.open("rb") as part:
                            shutil.copyfileobj(part, stream)
                outputs.append(target)
            new_samples.append({"sample_id": sample, "read1": f"reads/{sample}_R1.fastq.gz",
                                "read2": f"reads/{sample}_R2.fastq.gz"})
            reports.append({"sample_id": sample, "extra_pairs": args.pairs_per_sample, "commands": commands,
                            "read1_sha256": pilot.digest(outputs[0]), "read2_sha256": pilot.digest(outputs[1])})
    pilot.write_table(out / "samples.tsv", ["sample_id", "read1", "read2"], new_samples)
    shutil.copyfile(pilot_dir / "metadata.tsv", out / "metadata.tsv")
    composition += [{"sample_id": row["sample_id"], "community": community[row["sample_id"]],
                     "genome_id": args.genome, "pairs": str(args.pairs_per_sample)} for row in samples]
    pilot.write_table(out / "composition.tsv", ["sample_id", "community", "genome_id", "pairs"], composition)
    catalogue_dir = out / "catalogues" / "split_holdout"
    catalogue_dir.mkdir(parents=True)
    kept = []
    for row in catalogue_rows:
        if row["reference_id"] == args.genome:
            continue
        fasta = (base / row["fasta"]).resolve()
        kept.append({**row, "fasta": os.path.relpath(fasta, catalogue_dir)})
    pilot.write_table(catalogue_dir / "catalogue.tsv", list(catalogue_rows[0]), kept)
    genome_bases = sum(length for _, length in records)
    read_length = config["simulation"]["read_length"]
    pilot.write_json(out / "scenario.json", {
        "schema_version": 1,
        "purpose": "a novel genome too thin per group but assemblable when all samples are pooled",
        "genome": args.genome,
        "organism": genome["organism"],
        "species_taxid": genome["species_taxid"],
        "pairs_per_sample": args.pairs_per_sample,
        "coverage_per_sample": round(args.pairs_per_sample * 2 * read_length / genome_bases, 3),
        "base_catalogue": args.base_catalogue,
        "catalogue": "catalogues/split_holdout/catalogue.tsv (base minus the genome)",
        "seed": config["seed"],
        "samples": reports,
    })
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
