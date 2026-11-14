//! Sample grouping by k-mer sketch similarity (M8).
//!
//! Each sample's reads (both mates) are streamed once into a FracMinHash
//! sketch: the canonical k-mers whose 64-bit hash is at most 2^64 / scaled.
//! A sketch keeps at most `max_hashes` hashes. When it would hold more, its
//! scale doubles and the larger hashes are dropped, which leaves exactly the
//! subset a sketch made at the larger scale would hold; samples are then
//! compared at the largest scale any of them reached.
//!
//! Hashes seen in only one sample, mostly sequencing errors, are ignored:
//! similarity is the Jaccard index of the hashes each sample shares with at
//! least one other sample. This avoids an abundance filter, which would also
//! remove genomes sequenced below about 1x. Shared hashes are found by
//! sorting (hash, sample) entries. Each sample keeps only its most similar
//! neighbours, so no all-pairs matrix is built. `cluster` turns that sparse
//! graph into groups.

pub mod cluster;

use std::{
    collections::BTreeMap,
    path::Path,
    sync::{
        Mutex, PoisonError,
        atomic::{AtomicUsize, Ordering},
    },
    time::Instant,
};

use serde_json::{Value, json};

use crate::{
    assembly::{AssemblyGroup, inputs::ReadPair},
    cli::{GroupingArgs, SketchPolicy},
    discovery::extract::FastqReader,
    error::Error,
    execution::CancelToken,
};
use cluster::{Caps, Policy};

/// Bytes a sketching worker may hold per allowed hash: a buffer of twice
/// `max_hashes` 64-bit hashes that is sorted and deduplicated in place.
const WORKER_BYTES_PER_HASH: u64 = 16;
/// Bytes per hash once sketching ends: the sketch itself (8), then the sorted
/// (hash, sample) index (16) and the shared-hash lists (8) built from it.
const INDEX_BYTES_PER_HASH: u64 = 32;
const RECORDS_BETWEEN_CANCEL_CHECKS: u64 = 4096;

/// FracMinHash settings.
#[derive(Clone, Copy, Debug)]
pub struct SketchParams {
    pub k: usize,
    pub scaled: u64,
    pub seed: u64,
    pub max_hashes: usize,
}

/// One sample's sketch: sorted, distinct hashes at `scaled`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Sketch {
    pub hashes: Vec<u64>,
    pub scaled: u64,
}

/// MurmurHash3's 64-bit finalizer: a bijection with full avalanche.
fn fmix64(mut value: u64) -> u64 {
    value ^= value >> 33;
    value = value.wrapping_mul(0xff51_afd7_ed55_8ccd);
    value ^= value >> 33;
    value = value.wrapping_mul(0xc4ce_b9fe_1a85_ec53);
    value ^ (value >> 33)
}

/// The largest hash kept at `scaled`.
pub fn max_hash(scaled: u64) -> u64 {
    u64::MAX / scaled
}

const fn base_code(base: u8) -> u8 {
    match base {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 4,
    }
}

/// Accumulates one sample's sketch within a fixed buffer.
pub struct Sketcher {
    k: usize,
    mask: u64,
    shift: u32,
    key: u64,
    scaled: u64,
    max_hash: u64,
    cap: usize,
    buffer: Vec<u64>,
}

impl Sketcher {
    pub fn new(params: &SketchParams) -> Self {
        let bits = 2 * params.k as u32;
        Self {
            k: params.k,
            mask: if bits == 64 {
                u64::MAX
            } else {
                (1 << bits) - 1
            },
            shift: bits - 2,
            key: fmix64(params.seed),
            scaled: params.scaled,
            max_hash: max_hash(params.scaled),
            cap: params.max_hashes,
            buffer: Vec::new(),
        }
    }

    /// Adds the canonical k-mers of one sequence; k-mers with a base other
    /// than A, C, G, or T are skipped.
    pub fn add_sequence(&mut self, sequence: &[u8]) {
        let (mut forward, mut reverse, mut length) = (0u64, 0u64, 0usize);
        for &base in sequence {
            let code = base_code(base);
            if code > 3 {
                length = 0;
                continue;
            }
            forward = ((forward << 2) | u64::from(code)) & self.mask;
            reverse = (reverse >> 2) | (u64::from(3 - code) << self.shift);
            length += 1;
            if length >= self.k {
                let hash = fmix64(forward.min(reverse) ^ self.key);
                if hash <= self.max_hash {
                    self.push(hash);
                }
            }
        }
    }

    fn push(&mut self, hash: u64) {
        if self.buffer.len() == 2 * self.cap {
            self.compact();
        }
        if self.buffer.len() == self.buffer.capacity() {
            let target = (2 * self.cap).min((2 * self.buffer.capacity()).max(1024));
            self.buffer.reserve_exact(target - self.buffer.len());
        }
        self.buffer.push(hash);
    }

    /// Deduplicates, then doubles the scale until at most `cap` hashes remain.
    fn compact(&mut self) {
        self.buffer.sort_unstable();
        self.buffer.dedup();
        while self.buffer.len() > self.cap {
            self.scaled = self.scaled.saturating_mul(2);
            self.max_hash = max_hash(self.scaled);
            let limit = self.max_hash;
            let keep = self.buffer.partition_point(|&hash| hash <= limit);
            self.buffer.truncate(keep);
        }
    }

    pub fn finish(mut self) -> Sketch {
        self.compact();
        self.buffer.shrink_to_fit();
        Sketch {
            hashes: self.buffer,
            scaled: self.scaled,
        }
    }
}

/// Streams both mates of a sample into one sketch.
pub fn sketch_reads(
    read1: &Path,
    read2: &Path,
    params: &SketchParams,
    cancel: &CancelToken,
) -> Result<Sketch, Error> {
    let mut sketcher = Sketcher::new(params);
    for path in [read1, read2] {
        let mut reader = FastqReader::open(path)?;
        let mut records = 0u64;
        while let Some(record) = reader.next_record()? {
            sketcher.add_sequence(record.lines[1].trim_ascii_end());
            records += 1;
            if records % RECORDS_BETWEEN_CANCEL_CHECKS == 0 && cancel.is_cancelled() {
                return Err(Error::Interrupted(
                    "sketching was interrupted; no output was created".into(),
                ));
            }
        }
    }
    Ok(sketcher.finish())
}

/// A link between two samples (indices into the sorted inputs, `a < b`).
#[derive(Clone, Debug, PartialEq)]
pub struct Edge {
    pub a: usize,
    pub b: usize,
    pub shared: u64,
    pub similarity: f64,
}

/// The sparse similarity graph of sketches compared at one scale.
#[derive(Debug)]
pub struct Graph {
    pub scaled: u64,
    pub hashes: Vec<u64>,
    /// Hashes each sample shares with at least one other sample.
    pub shared: Vec<u64>,
    pub edges: Vec<Edge>,
}

/// Compares sketches at their largest scale and keeps, for every sample,
/// the `neighbors` most similar samples (ties go to the earlier sample).
pub fn similarity_graph(sketches: Vec<Sketch>, neighbors: usize) -> Graph {
    let samples = sketches.len();
    let scaled = sketches.iter().map(|s| s.scaled).max().unwrap_or(1);
    let limit = max_hash(scaled);
    let mut hashes = Vec::with_capacity(samples);
    let mut entries: Vec<(u64, u32)> = Vec::new();
    for (index, sketch) in sketches.into_iter().enumerate() {
        let kept = &sketch.hashes[..sketch.hashes.partition_point(|&hash| hash <= limit)];
        hashes.push(kept.len() as u64);
        entries.extend(kept.iter().map(|&hash| (hash, index as u32)));
    }
    entries.sort_unstable();
    // Each run of equal hashes lists the samples holding that hash.
    let mut members: Vec<u32> = Vec::new();
    let mut starts: Vec<usize> = vec![0];
    let mut runs_of: Vec<Vec<usize>> = vec![Vec::new(); samples];
    let mut index = 0;
    while index < entries.len() {
        let end = index + entries[index..].partition_point(|e| e.0 == entries[index].0);
        if end - index >= 2 {
            let run = starts.len() - 1;
            for &(_, sample) in &entries[index..end] {
                members.push(sample);
                runs_of[sample as usize].push(run);
            }
            starts.push(members.len());
        }
        index = end;
    }
    drop(entries);
    let shared: Vec<u64> = runs_of.iter().map(|runs| runs.len() as u64).collect();
    let mut links: BTreeMap<(usize, usize), (u64, f64)> = BTreeMap::new();
    let mut counts = vec![0u64; samples];
    let mut touched = Vec::new();
    for (sample, runs) in runs_of.iter().enumerate() {
        for &run in runs {
            for &other in &members[starts[run]..starts[run + 1]] {
                let other = other as usize;
                if other != sample {
                    if counts[other] == 0 {
                        touched.push(other);
                    }
                    counts[other] += 1;
                }
            }
        }
        let mut candidates: Vec<(f64, usize, u64)> = touched
            .iter()
            .map(|&other| {
                let common = counts[other];
                let union = shared[sample] + shared[other] - common;
                (common as f64 / union as f64, other, common)
            })
            .collect();
        candidates.sort_by(|x, y| y.0.total_cmp(&x.0).then(x.1.cmp(&y.1)));
        for &(similarity, other, common) in candidates.iter().take(neighbors) {
            links.insert((sample.min(other), sample.max(other)), (common, similarity));
        }
        for other in touched.drain(..) {
            counts[other] = 0;
        }
    }
    Graph {
        scaled,
        hashes,
        shared,
        edges: links
            .into_iter()
            .map(|((a, b), (shared, similarity))| Edge {
                a,
                b,
                shared,
                similarity,
            })
            .collect(),
    }
}

/// Sketches samples in parallel, at most `workers` at a time, in input order.
fn sketch_all(
    inputs: &[ReadPair],
    params: &SketchParams,
    workers: usize,
    cancel: &CancelToken,
) -> Result<Vec<Sketch>, Error> {
    let next = AtomicUsize::new(0);
    let results: Vec<Mutex<Option<Result<Sketch, Error>>>> =
        inputs.iter().map(|_| Mutex::new(None)).collect();
    std::thread::scope(|scope| {
        for _ in 0..workers {
            scope.spawn(|| {
                loop {
                    let index = next.fetch_add(1, Ordering::Relaxed);
                    let Some(input) = inputs.get(index) else {
                        break;
                    };
                    let result = sketch_reads(&input.read1, &input.read2, params, cancel);
                    if result.is_err() {
                        cancel.cancel();
                    }
                    *results[index]
                        .lock()
                        .unwrap_or_else(PoisonError::into_inner) = Some(result);
                }
            });
        }
    });
    let mut sketches = Vec::with_capacity(inputs.len());
    let mut first_error = None;
    for slot in results {
        match slot.into_inner().unwrap_or_else(PoisonError::into_inner) {
            Some(Ok(sketch)) => sketches.push(sketch),
            // Prefer the error that caused cancellation over the
            // interruptions it triggered in other workers.
            Some(Err(error))
                if first_error
                    .as_ref()
                    .is_none_or(|e| matches!(e, Error::Interrupted(_))) =>
            {
                first_error = Some(error);
            }
            Some(Err(_)) | None => {}
        }
    }
    match first_error {
        Some(error) => Err(error),
        None if sketches.len() == inputs.len() => Ok(sketches),
        None => Err(Error::Interrupted(
            "sketching was interrupted; no output was created".into(),
        )),
    }
}

/// Concurrent sketches that fit `memory_bytes`, and the resulting bound.
pub fn sketch_workers(
    samples: usize,
    max_hashes: usize,
    threads: usize,
    memory_bytes: u64,
) -> Result<(usize, u64), Error> {
    let per_hash = max_hashes as u64;
    let index = samples as u64 * per_hash * INDEX_BYTES_PER_HASH;
    let worker = per_hash * WORKER_BYTES_PER_HASH;
    let fit = memory_bytes.saturating_sub(index) / worker;
    let workers = (fit as usize).min(threads).min(samples.max(1));
    if workers == 0 {
        let needed = index + worker;
        return Err(Error::InvalidInput(format!(
            "sketching {samples} samples with --sketch-max-hashes {max_hashes} needs up to {:.1} GiB; lower --sketch-max-hashes or raise --memory-gib",
            needed as f64 / f64::from(1u32 << 30)
        )));
    }
    Ok((workers, index + workers as u64 * worker))
}

/// Resources and provenance for sketch grouping.
pub struct SketchContext<'a> {
    pub threads: usize,
    pub memory_mib: u64,
    /// What the sketches stand for: `discovery-pools` or `original-reads`.
    pub represents: &'a str,
    pub cancel: &'a CancelToken,
}

/// Groups sorted inputs by sketch similarity within the caps and returns the
/// groups with a provenance record.
pub fn sketch_groups(
    inputs: Vec<ReadPair>,
    grouping: &GroupingArgs,
    sample_cap: Option<usize>,
    context: &SketchContext,
) -> Result<(Vec<AssemblyGroup>, Value), Error> {
    let options = &grouping.sketch_options;
    let params = SketchParams {
        k: options.sketch_k,
        scaled: options.sketch_scaled.get(),
        seed: options.sketch_seed,
        max_hashes: options.sketch_max_hashes.get(),
    };
    let (workers, memory_bound) = sketch_workers(
        inputs.len(),
        params.max_hashes,
        context.threads,
        context.memory_mib * 1024 * 1024,
    )?;
    let started = Instant::now();
    let sketches = sketch_all(&inputs, &params, workers, context.cancel)?;
    let own_scales: Vec<u64> = sketches.iter().map(|s| s.scaled).collect();
    let graph = similarity_graph(sketches, options.sketch_neighbors.get());
    let policy = match options.sketch_policy {
        SketchPolicy::Merge => Policy::Merge,
        SketchPolicy::Threshold => Policy::Threshold(
            options
                .sketch_min_similarity
                .expect("validated: threshold has a minimum"),
        ),
        SketchPolicy::Seeded => Policy::Seeded,
    };
    let caps = Caps {
        samples: sample_cap,
        bases: grouping.max_group_bases.map(std::num::NonZeroU64::get),
    };
    let bases: Vec<u64> = inputs.iter().map(|i| i.bases).collect();
    let clusters = cluster::cluster(&bases, &graph.shared, &graph.edges, policy, caps);
    let seconds = started.elapsed().as_secs_f64();
    let similarity: BTreeMap<(usize, usize), f64> = graph
        .edges
        .iter()
        .map(|e| ((e.a, e.b), e.similarity))
        .collect();
    let id = |index: usize| inputs[index].sample_id.clone();
    let mut groups = Vec::with_capacity(clusters.len());
    let mut group_records = Vec::with_capacity(clusters.len());
    for (number, found) in clusters.iter().enumerate() {
        let group_id = format!("group-{:06}", number + 1);
        let representative = found.seed.unwrap_or(found.members[0]);
        let pairs = found.members.len() * (found.members.len() - 1) / 2;
        let mut total = 0.0;
        for (position, &a) in found.members.iter().enumerate() {
            for &b in &found.members[position + 1..] {
                total += similarity.get(&(a, b)).copied().unwrap_or(0.0);
            }
        }
        group_records.push(json!({
            "group_id": group_id,
            "representative": id(representative),
            "samples": found.members.iter().map(|&m| id(m)).collect::<Vec<_>>(),
            "mean_within_similarity": (pairs > 0).then(|| total / pairs as f64),
        }));
        groups.push(AssemblyGroup {
            group_id,
            grouped_by: "sketch".into(),
            value: id(representative),
            members: found.members.iter().map(|&m| inputs[m].clone()).collect(),
        });
    }
    let record = json!({
        "represents": context.represents,
        "method": "FracMinHash",
        "k": params.k,
        "scaled": params.scaled,
        "compared_at_scaled": graph.scaled,
        "seed": params.seed,
        "hash": "MurmurHash3 fmix64 of the 2-bit canonical k-mer XOR fmix64(seed); keep hashes <= floor((2^64 - 1) / scaled)",
        "max_hashes_per_sample": params.max_hashes,
        "similarity": "Jaccard index of the sketch hashes each sample shares with at least one other sample",
        "neighbors": options.sketch_neighbors.get(),
        "policy": options.sketch_policy,
        "min_similarity": options.sketch_min_similarity,
        "linkage": match options.sketch_policy {
            SketchPolicy::Seeded => "seed's neighbours by similarity",
            _ => "average over the neighbour graph (missing links count as 0)",
        },
        "workers": workers,
        "memory_bound_bytes": memory_bound,
        "seconds": (seconds * 1000.0).round() / 1000.0,
        "samples": inputs.iter().enumerate().map(|(index, input)| json!({
            "sample_id": input.sample_id,
            "own_scaled": own_scales[index],
            "hashes": graph.hashes[index],
            "shared_hashes": graph.shared[index],
        })).collect::<Vec<_>>(),
        "edges": graph.edges.iter().map(|e| json!({
            "samples": [id(e.a), id(e.b)],
            "shared_hashes": e.shared,
            "similarity": e.similarity,
        })).collect::<Vec<_>>(),
        "groups": group_records,
    });
    Ok((groups, record))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn params(k: usize, scaled: u64, max_hashes: usize) -> SketchParams {
        SketchParams {
            k,
            scaled,
            seed: 42,
            max_hashes,
        }
    }

    fn sketch_of(sequences: &[&[u8]], params: &SketchParams) -> Sketch {
        let mut sketcher = Sketcher::new(params);
        for sequence in sequences {
            sketcher.add_sequence(sequence);
        }
        sketcher.finish()
    }

    /// A deterministic pseudo-random DNA sequence.
    fn genome(seed: u64, length: usize) -> Vec<u8> {
        let mut state = seed;
        (0..length)
            .map(|_| {
                state = fmix64(state.wrapping_add(0x9e37_79b9_7f4a_7c15));
                b"ACGT"[(state >> 62) as usize]
            })
            .collect()
    }

    fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
        sequence
            .iter()
            .rev()
            .map(|base| match base {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                other => *other,
            })
            .collect()
    }

    #[test]
    fn strands_and_case_give_one_sketch_and_k_mers_with_n_are_skipped() {
        let p = params(15, 1, 100_000);
        let forward = genome(1, 500);
        let lower = forward.to_ascii_lowercase();
        let both = sketch_of(&[&forward], &p);
        assert_eq!(both.hashes.len(), 486, "500 - 15 + 1 distinct k-mers");
        assert_eq!(sketch_of(&[&reverse_complement(&forward)], &p), both);
        assert_eq!(sketch_of(&[&lower], &p), both);
        let mut broken = forward.clone();
        broken[100] = b'N';
        assert_eq!(sketch_of(&[&broken], &p).hashes.len(), 486 - 15);
        assert!(sketch_of(&[b"ACGTN"], &p).hashes.is_empty());
    }

    #[test]
    fn the_seed_changes_hashes_and_scaled_keeps_a_consistent_subset() {
        let sequence = genome(2, 20_000);
        let all = sketch_of(&[&sequence], &params(21, 1, 1_000_000));
        let mut other_seed = params(21, 1, 1_000_000);
        other_seed.seed = 7;
        assert_ne!(sketch_of(&[&sequence], &other_seed).hashes, all.hashes);
        let scaled = sketch_of(&[&sequence], &params(21, 10, 1_000_000));
        let expected: Vec<u64> = all
            .hashes
            .iter()
            .copied()
            .filter(|&h| h <= max_hash(10))
            .collect();
        assert_eq!(scaled.hashes, expected);
        assert!((1_000..3_000).contains(&scaled.hashes.len()), "about 1/10");
    }

    #[test]
    fn the_hash_cap_doubles_the_scale_and_bounds_the_buffer() {
        let sequence = genome(3, 50_000);
        let full = sketch_of(&[&sequence], &params(21, 1, 1_000_000));
        let mut sketcher = Sketcher::new(&params(21, 1, 1_000));
        sketcher.add_sequence(&sequence);
        assert!(
            sketcher.buffer.capacity() <= 2_000,
            "the buffer never exceeds twice the cap"
        );
        let capped = sketcher.finish();
        assert!(capped.hashes.len() <= 1_000);
        assert!(
            capped.scaled >= 64,
            "50,000 k-mers need a scale of at least 50"
        );
        let expected: Vec<u64> = full
            .hashes
            .iter()
            .copied()
            .filter(|&h| h <= max_hash(capped.scaled))
            .collect();
        assert_eq!(
            capped.hashes, expected,
            "downsampling equals sketching at the larger scale"
        );
    }

    #[test]
    fn similarity_ignores_unshared_hashes_and_compares_at_the_largest_scale() {
        let sketch = |hashes: &[u64], scaled| Sketch {
            hashes: hashes.to_vec(),
            scaled,
        };
        // Hash 9 is private to sample 0 and never counts.
        let graph = similarity_graph(
            vec![
                sketch(&[1, 2, 3, 9], 1),
                sketch(&[1, 2, 4], 1),
                sketch(&[3, 4], 1),
                sketch(&[], 1),
            ],
            10,
        );
        assert_eq!(graph.shared, [3, 3, 2, 0]);
        let edges: Vec<_> = graph
            .edges
            .iter()
            .map(|e| (e.a, e.b, e.shared, e.similarity))
            .collect();
        assert_eq!(
            edges,
            [
                (0, 1, 2, 2.0 / 4.0),
                (0, 2, 1, 1.0 / 4.0),
                (1, 2, 1, 1.0 / 4.0)
            ]
        );
        let one_neighbor = similarity_graph(
            vec![
                sketch(&[1, 2, 3, 9], 1),
                sketch(&[1, 2, 4], 1),
                sketch(&[3, 4], 1),
            ],
            1,
        );
        assert_eq!(
            one_neighbor
                .edges
                .iter()
                .map(|e| (e.a, e.b))
                .collect::<Vec<_>>(),
            [(0, 1), (0, 2)],
            "0 and 1 pick each other; 2 ties and keeps the earlier sample"
        );
        let limit = max_hash(4);
        let rescaled = similarity_graph(
            vec![sketch(&[1, limit + 1], 1), sketch(&[1, limit + 1], 4)],
            10,
        );
        assert_eq!(rescaled.scaled, 4);
        assert_eq!(
            rescaled.hashes,
            [1, 1],
            "hashes above the common scale are dropped"
        );
    }

    #[test]
    fn workers_fit_the_memory_budget_or_fail_clearly() {
        let mib = 1024 * 1024;
        assert_eq!(sketch_workers(4, 1_000, 8, 1024 * mib).unwrap().0, 4);
        let (workers, bound) = sketch_workers(30, 1_000_000, 64, 1024 * mib).unwrap();
        assert_eq!(workers, 7, "(1 GiB - 30 x 32 MB) / 16 MB");
        assert!(bound <= 1024 * mib);
        let error = sketch_workers(1_000, 1_000_000, 8, 1024 * mib).unwrap_err();
        assert!(error.to_string().contains("--sketch-max-hashes"), "{error}");
    }
}
