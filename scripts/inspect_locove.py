#!/usr/bin/env python3
"""Measure LoCovE candidate-window volume; never assign or discard reads."""

import argparse
import ast
import json
import math
from pathlib import Path
import struct
import zipfile


def read_array(archive, name):
    data = archive.read(name + ".npy")
    if data[:8] != b"\x93NUMPY\x01\x00":
        raise ValueError("expected LoCovE's NPY v1 output")
    length = int.from_bytes(data[8:10], "little")
    header = ast.literal_eval(data[10:10 + length].decode().strip())
    if header["fortran_order"] or len(header["shape"]) != 1:
        raise ValueError("expected one-dimensional C-order array")
    fmt = {"<i8": "<q", "<f8": "<d"}.get(header["descr"])
    if fmt is None:
        raise ValueError("unsupported LoCovE array type")
    values = [value[0] for value in struct.iter_unpack(fmt, data[10 + length:])]
    if len(values) != header["shape"][0]:
        raise ValueError("array size mismatch")
    return values


def merged_bases(intervals):
    total, end = 0, -1
    for start, stop in sorted(intervals):
        if start < 0 or stop < start:
            raise ValueError("invalid half-open interval")
        total += max(0, stop - max(start, end))
        end = max(end, stop)
    return total


def inspect(directory):
    directory = Path(directory)
    mapping = json.loads((directory / "record-map.json").read_text())
    with zipfile.ZipFile(directory / "segments.npz") as archive:
        arrays = {name: read_array(archive, name) for name in [
            "seg_start", "seg_end", "seg_record", "n_anchors", "breadth_anchor"
        ]}
    if len({len(values) for values in arrays.values()}) != 1:
        raise ValueError("segment arrays differ in length")
    rows = list(zip(*(arrays[name] for name in [
        "seg_start", "seg_end", "seg_record", "n_anchors", "breadth_anchor"
    ])))
    # LoCovE's record map has no lengths; infer them from segment ends. A record
    # without segments has unknown length and coverage, so it cannot be scored.
    lengths = [0] * len(mapping["records"])
    seen = [False] * len(lengths)
    for start, end, record, _, _ in rows:
        if not 0 <= record < len(lengths) or start < 0 or end < start:
            raise ValueError("invalid record or coordinate")
        lengths[record] = max(lengths[record], end)
        seen[record] = True
    missing = [name for name, present in zip(mapping["records"], seen) if not present]
    if missing:
        raise ValueError(f"{len(missing)} record(s) have no segments, e.g. {missing[:3]}; their coverage is unknown")
    if not sum(lengths):
        raise ValueError("segments cover no reference bases")
    result = {
        "schema_version": 1, "segments": len(rows), "records": len(lengths),
        "reference_bases": sum(lengths),
        "no_anchor_segments": sum(anchors == 0 for _, _, _, anchors, _ in rows),
        "analysis_only": True, "per_read_assignments": False, "experiments": []
    }
    for threshold in [0.0, 0.1, 0.5, 0.8, 0.95]:
        for padding in [0, 500]:
            intervals = [[] for _ in lengths]
            for start, end, record, anchors, breadth in rows:
                # Unknown segments stay in the candidate set. NaN is not zero.
                if anchors == 0 or not math.isfinite(breadth) or breadth > threshold:
                    intervals[record].append((max(0, start - padding), min(lengths[record], end + padding)))
            bases = sum(merged_bases(items) for items in intervals)
            result["experiments"].append({
                "breadth_strictly_greater_than": threshold, "flank_bases": padding,
                "retain_unknown": True, "candidate_bases": bases,
                "fraction_of_reference": bases / sum(lengths)
            })
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("coverage_dir", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    text = json.dumps(inspect(args.coverage_dir), indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("x") as stream:
            stream.write(text)
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
