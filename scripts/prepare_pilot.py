#!/usr/bin/env python3
"""Explicit M1 research-data preparation, separate from the Rust runtime."""

import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import random
import re
import shutil
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from contextlib import ExitStack


DEFAULT_CONFIG = Path(__file__).resolve().parents[1] / "experiments/pilot/pilot.json"
ID = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]*\Z")


def digest(path, algorithm="sha256"):
    value = hashlib.new(algorithm)
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def read_table(path):
    with Path(path).open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError(f"{path}: missing or duplicate headers")
        rows = list(reader)
    if not rows or any(None in row or None in row.values() for row in rows):
        raise ValueError(f"{path}: empty or malformed table")
    return rows


def write_table(path, fields, rows):
    with Path(path).open("x", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path, value):
    with Path(path).open("x", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")


def allocate(total, lengths):
    """Largest remainders; exact total, deterministic ties, no floating point."""
    if total < 0 or not lengths or any(n <= 0 for n in lengths):
        raise ValueError("invalid allocation")
    denominator = sum(lengths)
    result = [total * n // denominator for n in lengths]
    order = sorted(range(len(lengths)), key=lambda i: (-(total * lengths[i] % denominator), i))
    for i in order[: total - sum(result)]:
        result[i] += 1
    return result


def derived_seed(seed, *parts):
    text = ":".join(map(str, (seed, *parts))).encode()
    return int.from_bytes(hashlib.sha256(text).digest()[:4], "big") % 2147483646 + 1


def load_design(config_path):
    config_path = Path(config_path).resolve()
    config = json.loads(config_path.read_text())
    if config.get("schema_version") != 1:
        raise ValueError("unsupported pilot schema_version")
    genomes = read_table(config_path.parent / config["genomes"])
    ids = [g["genome_id"] for g in genomes]
    if len(ids) != len(set(ids)):
        raise ValueError("duplicate genome ID")
    for genome in genomes:
        if not re.fullmatch(r"GCF_\d{9}\.\d+", genome["genome_id"]):
            raise ValueError("genome IDs must be versioned RefSeq assembly accessions")
        if not re.fullmatch(r"[a-f0-9]{32}", genome["source_md5"]):
            raise ValueError("each reference needs its frozen source MD5")
    composition = read_table(config_path.parent / config["composition"])
    samples = {}
    seen = set()
    for row in composition:
        sample, genome = row["sample_id"], row["genome_id"]
        if not ID.fullmatch(sample) or genome not in ids:
            raise ValueError("invalid sample or unknown composition genome")
        if (sample, genome) in seen:
            raise ValueError("duplicate sample/genome composition entry")
        seen.add((sample, genome))
        pairs = int(row["pairs"])
        if pairs <= 0:
            raise ValueError("composition counts must be positive")
        item = samples.setdefault(sample, {"group": row["community"], "counts": {}})
        if item["group"] != row["community"]:
            raise ValueError("sample has inconsistent metadata groups")
        item["counts"][genome] = pairs
    if any(sum(s["counts"].values()) != config["read_pairs_per_sample"] for s in samples.values()):
        raise ValueError("sample pair counts do not match the configured total")
    roles = {g["genome_id"]: g["role"] for g in genomes}
    if any(roles[g] == "distractor" for s in samples.values() for g in s["counts"]):
        raise ValueError("catalogue-only distractor appears in simulated reads")
    for name, excluded in config["catalogues"].items():
        if not ID.fullmatch(name) or not set(excluded) <= set(ids):
            raise ValueError("invalid catalogue case")
    resources = config["resources"]
    limits = ["threads", "memory_gib", "scratch_gib", "genome_download_limit_mib", "simulation_jobs", "sample_jobs"]
    if any(resources[k] <= 0 for k in limits):
        raise ValueError("resource limits must be positive")
    if resources["simulation_jobs"] > resources["threads"]:
        raise ValueError("simulation jobs exceed thread budget")
    if resources["sample_jobs"] * (resources["simulation_jobs"] + 1) > resources["threads"]:
        raise ValueError("sample preparation concurrency exceeds thread budget")
    return config, genomes, samples


def source_url(genome):
    accession = genome["genome_id"]
    digits = accession.split("_")[1].split(".")[0]
    stem = f'{accession}_{genome["assembly_name"]}'
    return f"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/{digits[:3]}/{digits[3:6]}/{digits[6:]}/{stem}/{stem}_genomic.fna.gz"


def download(url, destination, max_bytes):
    for attempt in range(5):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "meta-co-assembler-pilot/0.1"})
            with urllib.request.urlopen(request, timeout=60) as response, destination.open("wb") as out:
                count = 0
                while block := response.read(1024 * 1024):
                    count += len(block)
                    if count > max_bytes:
                        raise ValueError("genome download byte limit exceeded")
                    out.write(block)
            return count
        except (urllib.error.URLError, TimeoutError) as error:
            if isinstance(error, urllib.error.HTTPError) and error.code not in [429, 500, 502, 503, 504]:
                raise
            if attempt == 4:
                raise
            time.sleep(min(2 ** attempt, 16))


def fasta_records(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as stream:
        name, sequence = None, []
        for line in stream:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(sequence)
                fields = line[1:].split()
                if not fields:
                    raise ValueError(f"{path}: FASTA record without a name")
                name, sequence = fields[0], []
            elif name is None:
                if line.strip():
                    raise ValueError(f"{path}: FASTA sequence before the first header")
            else:
                sequence.append("".join(line.split()).upper())
        if name is not None:
            yield name, "".join(sequence)


def fastq_records(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as stream:
        while header := stream.readline():
            sequence, plus, quality = stream.readline(), stream.readline(), stream.readline()
            if not header.startswith("@") or not plus.startswith("+") or not quality:
                raise ValueError(f"{path}: malformed FASTQ")
            sequence, quality = sequence.rstrip("\r\n"), quality.rstrip("\r\n")
            if not sequence or len(sequence) != len(quality):
                raise ValueError(f"{path}: sequence/quality length mismatch")
            fields = header[1:].split()
            if not fields:
                raise ValueError(f"{path}: FASTQ record without a name")
            yield fields[0], sequence, quality


def pair_name(name):
    return name[:-2] if name.endswith(("/1", "/2")) else name


def paired_records(read1, read2):
    first, second = fastq_records(read1), fastq_records(read2)
    while True:
        one, two = next(first, None), next(second, None)
        if one is None and two is None:
            return
        if one is None or two is None or pair_name(one[0]) != pair_name(two[0]):
            raise ValueError("missing or mismatched FASTQ mate")
        yield one, two


def simulate_part(task):
    command, prefix, count, genome = task
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    with Path(str(prefix) + ".log").open("xb") as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, env=env, check=True)
    actual = sum(1 for _ in paired_records(str(prefix) + "1.fq", str(prefix) + "2.fq"))
    if actual < count:
        raise ValueError(f"{genome}: ART produced {actual} pairs, below target {count}")
    return {"prefix": prefix, "pairs": count, "generated_pairs": actual, "genome": genome, "command": command}


def mix_parts(parts, read1, read2, sample, length, seed, compression_level=3):
    rng = random.Random(seed)
    streams = [paired_records(str(p["prefix"]) + "1.fq", str(p["prefix"]) + "2.fq") for p in parts]
    remaining = [p["pairs"] for p in parts]
    counts = {}
    with ExitStack() as stack:
        writers = [
            stack.enter_context(gzip.GzipFile(filename="", mode="wb", fileobj=stack.enter_context(p.open("xb")), mtime=0, compresslevel=compression_level))
            for p in [read1, read2]
        ]
        total, ordinal = sum(remaining), 0
        while total:
            pick = rng.randrange(total)
            index = 0
            while pick >= remaining[index]:
                pick -= remaining[index]
                index += 1
            pair = next(streams[index], None)
            if pair is None:
                raise ValueError("ART produced fewer pairs than requested")
            genome = parts[index]["genome"]
            ordinal += 1
            for mate, (record, writer) in enumerate(zip(pair, writers), 1):
                _, sequence, quality = record
                if len(sequence) != length:
                    raise ValueError("ART produced an unexpected read length")
                writer.write(f"@{sample}|{genome}|{ordinal:09d}/{mate}\n{sequence}\n+\n{quality}\n".encode())
            counts[genome] = counts.get(genome, 0) + 1
            remaining[index] -= 1
            total -= 1
            if ordinal % 250000 == 0:
                print(f"Writing {sample}: {ordinal:,} synchronized pairs", flush=True)
        for stream, part in zip(streams, parts):
            surplus = sum(1 for _ in stream)
            if surplus != part.get("generated_pairs", part["pairs"]) - part["pairs"]:
                raise ValueError("ART surplus count differs from the verified part count")
    return counts


def simulate_sample(config, destination, source_records, sample, info):
    simulation = config["simulation"]
    print(f"Simulating {sample}: {config['read_pairs_per_sample']:,} read pairs", flush=True)
    with tempfile.TemporaryDirectory(prefix=f"{sample}-", dir=destination) as temp:
        tasks = []
        for genome, pairs in sorted(info["counts"].items()):
            records = source_records[genome]
            for index, ((path, _), count) in enumerate(zip(records, allocate(pairs, [n for _, n in records]))):
                if count == 0:
                    continue
                prefix = Path(temp) / f"{genome}-{index}-"
                # ART may skip fragments overlapping masked bases instead of replacing them.
                # Generate a deterministic surplus; consume exactly the designed count.
                surplus = max(simulation["art_min_surplus_pairs"], int(count * simulation["art_surplus_fraction"]) + 1)
                command = [shutil.which("art_illumina"), "-ss", simulation["profile"], "-p",
                           "-l", str(simulation["read_length"]), "-m", str(simulation["fragment_mean"]),
                           "-s", str(simulation["fragment_sd"]), "-rs", str(derived_seed(config["seed"], sample, genome, index)),
                           "-c", str(count + surplus), "-na", "-i", str(path), "-o", str(prefix)]
                tasks.append((command, prefix, count, genome))
        with ThreadPoolExecutor(max_workers=config["resources"]["simulation_jobs"]) as pool:
            parts = list(pool.map(simulate_part, tasks))
        for part in parts:
            log = Path(str(part["prefix"]) + ".log")
            shutil.copyfile(log, destination / "logs" / f'{sample}-{log.name}')
        read1, read2 = [destination / "reads" / f"{sample}_R{mate}.fastq.gz" for mate in [1, 2]]
        counts = mix_parts(parts, read1, read2, sample, simulation["read_length"],
                           derived_seed(config["seed"], sample, "mixture"), simulation["gzip_level"])
        if counts != info["counts"]:
            raise ValueError("generated source composition differs from design")
        report = {"sample_id": sample, "source_pairs": counts, "pairs": sum(counts.values()),
                  "read1_sha256": digest(read1), "read2_sha256": digest(read2),
                  "commands": [part["command"] for part in parts],
                  "art_generated_pairs": sum(part["generated_pairs"] for part in parts)}
    print(f"Completed {sample}", flush=True)
    return report


def prepare(config_path, destination):
    config_path, destination = Path(config_path).resolve(), Path(destination).resolve()
    config, genomes, samples = load_design(config_path)
    if destination.exists():
        raise ValueError("output exists; choose a fresh destination (existing data is never overwritten)")
    if not shutil.which("art_illumina"):
        raise ValueError("art_illumina is missing; synchronize and use meta-co-assembler-env")
    destination.parent.mkdir(parents=True, exist_ok=True)
    total_pairs = len(samples) * config["read_pairs_per_sample"]
    estimate = total_pairs * (4 * config["simulation"]["read_length"] + 160) + 2 * 1024 ** 3
    if estimate > config["resources"]["scratch_gib"] * 1024 ** 3:
        raise ValueError("estimated preparation size exceeds scratch budget")
    if shutil.disk_usage(destination.parent).free < estimate:
        raise ValueError("insufficient free space for pilot preparation")
    destination.mkdir()
    for name in ["downloads", "references", "records", "reads", "logs", "catalogues"]:
        (destination / name).mkdir()
    (destination / "PREPARING").write_text("Incomplete until dataset.json is written.\n")
    source_records, reference_rows, source_manifest = {}, [], []
    remaining_download = config["resources"]["genome_download_limit_mib"] * 1024 ** 2
    for genome in genomes:
        gid = genome["genome_id"]
        print(f"Downloading frozen reference {gid}", flush=True)
        archive = destination / "downloads" / f"{gid}.fna.gz"
        remaining_download -= download(source_url(genome), archive, remaining_download)
        if digest(archive, "md5") != genome["source_md5"]:
            raise ValueError(f"{gid}: frozen source MD5 mismatch")
        records = list(fasta_records(archive))
        if len(records) != int(genome["records"]) or sum(len(seq) for _, seq in records) != int(genome["bases"]):
            raise ValueError(f"{gid}: reference size/record count mismatch")
        normalized = destination / "references" / f"{gid}.fna"
        per_record = destination / "records" / gid
        per_record.mkdir()
        source_records[gid] = []
        with normalized.open("x") as out:
            for index, (name, sequence) in enumerate(records):
                qualified = f"{gid}__{name}"
                text = f">{qualified}\n" + "\n".join(sequence[i:i+80] for i in range(0, len(sequence), 80)) + "\n"
                out.write(text)
                record_path = per_record / f"{index}.fna"
                record_path.write_text(text)
                source_records[gid].append((record_path, len(sequence)))
                reference_rows.append({"genome_id": gid, "record_id": qualified, "source_record_id": name, "length": len(sequence)})
        source_manifest.append({**genome, "url": source_url(genome), "source_sha256": digest(archive), "fasta_sha256": digest(normalized)})
        time.sleep(1)
    write_table(destination / "record-map.tsv", ["genome_id", "record_id", "source_record_id", "length"], reference_rows)
    fasta_digests = {g["genome_id"]: g["fasta_sha256"] for g in source_manifest}
    for case in config["catalogues"]:
        directory = destination / "catalogues" / case
        directory.mkdir()
        rows = catalogue_rows(config, genomes, fasta_digests, case)
        write_table(directory / "catalogue.tsv",
                    ["reference_id", "fasta", "taxid", "species_taxid", "sha256"], rows)
        (directory / "references.txt").write_text("".join(row["fasta"] + "\n" for row in rows))
    simulation = config["simulation"]
    with ThreadPoolExecutor(max_workers=config["resources"]["sample_jobs"]) as pool:
        reports = list(pool.map(
            lambda item: simulate_sample(config, destination, source_records, *item),
            sorted(samples.items()),
        ))
    write_table(destination / "samples.tsv", ["sample_id", "read1", "read2"],
                [{"sample_id": s, "read1": f"reads/{s}_R1.fastq.gz", "read2": f"reads/{s}_R2.fastq.gz"} for s in sorted(samples)])
    write_table(destination / "metadata.tsv", ["sample_id", "community"],
                [{"sample_id": s, "community": samples[s]["group"]} for s in sorted(samples)])
    inventory = []
    if os.environ.get("CONDA_PREFIX"):
        for record in sorted((Path(os.environ["CONDA_PREFIX"]) / "conda-meta").glob("*.json")):
            package = json.loads(record.read_text())
            inventory.append({key: package[key] for key in ["name", "version", "build"]})
    write_json(destination / "dataset.json", {
        "schema_version": 1, "status": "complete", "config_sha256": digest(config_path),
        "genomes_table_sha256": digest(config_path.parent / config["genomes"]),
        "composition_sha256": digest(config_path.parent / config["composition"]),
        "preparation_script_sha256": digest(Path(__file__)),
        "seed": config["seed"], "simulation": simulation, "references": source_manifest,
        "samples": reports, "packages": inventory,
        "truth": "read headers encode sample|source assembly|pair ordinal/mate; full genomes are assembly truth"
    })
    (destination / "PREPARING").unlink()
    print(f"Prepared {len(samples)} samples, {total_pairs:,} pairs at {destination}", flush=True)


def catalogue_rows(config, genomes, fasta_digests, case):
    """The exact catalogue prepare() writes for one case, derived from the design."""
    excluded = set(config["catalogues"][case])
    return [{"reference_id": g["genome_id"], "fasta": f'../../references/{g["genome_id"]}.fna',
             "taxid": g["taxid"], "species_taxid": g["species_taxid"], "sha256": fasta_digests[g["genome_id"]]}
            for g in genomes if g["genome_id"] not in excluded]


def verify_catalogues(config, genomes, fasta_digests, destination):
    """Holdouts must be physically absent from every index input set."""
    root = destination / "catalogues"
    present = sorted(path.name for path in root.iterdir())
    if present != sorted(config["catalogues"]):
        raise ValueError(f"catalogue cases {present} differ from the design {sorted(config['catalogues'])}")
    for case, excluded in config["catalogues"].items():
        expected = catalogue_rows(config, genomes, fasta_digests, case)
        actual = read_table(root / case / "catalogue.tsv")
        leaked = sorted({row.get("reference_id") for row in actual} & set(excluded))
        if leaked:
            raise ValueError(f"{case}: held-out references {leaked} are present in the catalogue")
        if actual != expected:
            raise ValueError(f"{case}: catalogue.tsv differs from the frozen design")
        listed = (root / case / "references.txt").read_text().splitlines()
        if listed != [row["fasta"] for row in expected]:
            raise ValueError(f"{case}: references.txt differs from catalogue.tsv")


def verify(config_path, destination):
    config, genomes, samples = load_design(config_path)
    destination = Path(destination)
    manifest = json.loads((destination / "dataset.json").read_text())
    if (destination / "PREPARING").exists() or manifest.get("status") != "complete" or manifest.get("config_sha256") != digest(config_path):
        raise ValueError("dataset is incomplete or uses another configuration")
    for key, filename in [("genomes_table_sha256", config["genomes"]), ("composition_sha256", config["composition"])]:
        if manifest.get(key) != digest(Path(config_path).resolve().parent / filename):
            raise ValueError("frozen genome registry or composition has changed")
    # Check against the design, not the manifest's own lists, so omissions are detected.
    reported = sorted(report["sample_id"] for report in manifest["samples"])
    if reported != sorted(samples):
        raise ValueError(f"manifest samples {reported} differ from the design {sorted(samples)}")
    genome_ids = [g["genome_id"] for g in genomes]
    references = {genome["genome_id"]: genome for genome in manifest["references"]}
    if sorted(references) != sorted(genome_ids) or len(references) != len(manifest["references"]):
        raise ValueError("manifest references differ from the frozen genome registry")
    for report in manifest["samples"]:
        sample = report["sample_id"]
        read1, read2 = [destination / "reads" / f"{sample}_R{mate}.fastq.gz" for mate in [1, 2]]
        for path, key in [(read1, "read1_sha256"), (read2, "read2_sha256")]:
            if digest(path) != report[key]:
                raise ValueError(f"{sample}: read file checksum mismatch")
        counts = {}
        for ordinal, pair in enumerate(paired_records(read1, read2), 1):
            for mate, (name, sequence, _) in enumerate(pair, 1):
                head, genome, number = name.split("|")
                if head != sample or number != f"{ordinal:09d}/{mate}" or len(sequence) != config["simulation"]["read_length"]:
                    raise ValueError("invalid sample-aware read identity or length")
            counts[genome] = counts.get(genome, 0) + 1
        if counts != samples[sample]["counts"]:
            raise ValueError(f"{sample}: source pair counts differ from truth")
        print(f"Verified {sample}: {sum(counts.values()):,} synchronized pairs", flush=True)
    fasta_digests = {}
    for genome_id in genome_ids:
        fasta_digests[genome_id] = digest(destination / "references" / f"{genome_id}.fna")
        if fasta_digests[genome_id] != references[genome_id]["fasta_sha256"]:
            raise ValueError("reference checksum mismatch")
    verify_catalogues(config, genomes, fasta_digests, destination)
    print(f"Verified {len(genome_ids)} references and {len(config['catalogues'])} catalogue cases; "
          "held-out references are absent", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=["validate", "prepare", "verify"])
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output-dir", type=Path, default=Path("data/pilot"))
    args = parser.parse_args()
    try:
        if args.action == "prepare":
            prepare(args.config, args.output_dir)
        elif args.action == "verify":
            verify(args.config, args.output_dir)
        else:
            _, genomes, samples = load_design(args.config)
            print(f"Valid design: {len(genomes)} references, {len(samples)} samples")
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
