//! Accepting assembled contigs into the cumulative contig set.
//!
//! A candidate must be at least the minimum length. Containment is judged
//! from minimap2 `-x asm5 -c` alignments: contig A is contained in B when one
//! alignment covers at least 95% of A at identity (matches / block length)
//! of at least 95%. A new contig contained in an accepted contig, or in a
//! longer new contig (ties: the smaller ID wins), is redundant. An accepted
//! contig contained in a newly accepted one is superseded: it stays in its
//! round's outputs but leaves the cumulative set.

use std::{
    collections::{BTreeMap, HashMap, HashSet},
    fs::{self, File},
    io::{BufRead, BufReader, BufWriter, Write},
    path::{Path, PathBuf},
};

use serde_json::{Value, json};

use crate::{
    assembly::contigs::{FastaReader, file_stats},
    catalogue::sha256_file,
};

pub const CANDIDATES: &str = "candidates.fa";
pub const CANDIDATE_COUNTS: &str = "candidates.json";
pub const CUMULATIVE: &str = "cumulative.fa";
pub const TABLE: &str = "contigs.tsv";
pub const SUMMARY: &str = "accept.json";
pub const NEW_VS_NEW: &str = "new_vs_new.paf";
pub const NEW_VS_OLD: &str = "new_vs_old.paf";
pub const CONTAINMENT_MIN_IDENTITY: f64 = 0.95;
pub const CONTAINMENT_MIN_COVERAGE: f64 = 0.95;
pub const SCHEMA_VERSION: u32 = 1;

/// Contigs offered for acceptance: a FASTA and the prefix that makes its IDs
/// unique in the run (empty when they already are).
#[derive(Clone, Debug)]
pub struct CandidateSource {
    pub path: PathBuf,
    pub prefix: String,
}

/// Writes candidates at least `min_length` long to `candidates.fa` with
/// their final IDs, and their counts to `candidates.json`.
pub fn write_candidates(
    dir: &Path,
    sources: &[CandidateSource],
    min_length: u64,
) -> Result<u64, String> {
    let path = dir.join(CANDIDATES);
    let write_error = |error: std::io::Error| format!("cannot write {}: {error}", path.display());
    let mut out = BufWriter::new(File::create_new(&path).map_err(write_error)?);
    let (mut offered, mut short, mut kept) = (0u64, 0u64, 0u64);
    let mut ids = HashSet::new();
    for source in sources {
        let mut reader = FastaReader::open(&source.path)?;
        while let Some(record) = reader.next_record()? {
            offered += 1;
            if record.length() < min_length {
                short += 1;
                continue;
            }
            let id = format!("{}{}", source.prefix, record.id);
            if !ids.insert(id.clone()) {
                return Err(format!("duplicate candidate contig ID {id:?}"));
            }
            record.write_as(&mut out, &id).map_err(write_error)?;
            kept += 1;
        }
    }
    out.flush().map_err(write_error)?;
    let counts = json!({
        "sources": sources.iter().map(|s| &s.path).collect::<Vec<_>>(),
        "min_length": min_length,
        "offered": offered,
        "short": short,
        "candidates": kept,
    });
    fs::write(dir.join(CANDIDATE_COUNTS), counts.to_string() + "\n")
        .map_err(|error| format!("cannot write {CANDIDATE_COUNTS}: {error}"))?;
    Ok(kept)
}

/// A containment relation between two contigs.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Containment {
    /// A new contig inside an already accepted one.
    NewInOld { new: String, old: String },
    /// An accepted contig inside a new one.
    OldInNew { old: String, new: String },
    /// A new contig inside another new one.
    NewInNew { inner: String, outer: String },
}

fn field<T: std::str::FromStr>(fields: &[&str], index: usize, line: &str) -> Result<T, String> {
    fields
        .get(index)
        .and_then(|value| value.parse().ok())
        .ok_or_else(|| format!("malformed PAF line {line:?}"))
}

/// Containments in a PAF whose queries are new contigs. Targets are accepted
/// contigs when `targets_old`, otherwise new contigs (self hits ignored).
pub fn parse_paf(
    input: impl BufRead,
    targets_old: bool,
    edges: &mut Vec<Containment>,
) -> Result<(), String> {
    for line in input.lines() {
        let line = line.map_err(|error| format!("cannot read PAF: {error}"))?;
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        let (query, target) = (fields[0].to_owned(), field::<String>(&fields, 5, &line)?);
        let qlen: u64 = field(&fields, 1, &line)?;
        let (qstart, qend): (u64, u64) = (field(&fields, 2, &line)?, field(&fields, 3, &line)?);
        let tlen: u64 = field(&fields, 6, &line)?;
        let (tstart, tend): (u64, u64) = (field(&fields, 7, &line)?, field(&fields, 8, &line)?);
        let (matches, block): (u64, u64) = (field(&fields, 9, &line)?, field(&fields, 10, &line)?);
        if query == target || block == 0 || qlen == 0 || tlen == 0 {
            continue;
        }
        #[allow(clippy::cast_precision_loss)]
        let covered =
            |span: u64, length: u64| span as f64 >= CONTAINMENT_MIN_COVERAGE * length as f64;
        #[allow(clippy::cast_precision_loss)]
        if (matches as f64) < CONTAINMENT_MIN_IDENTITY * block as f64 {
            continue;
        }
        if covered(qend.saturating_sub(qstart), qlen) {
            edges.push(if targets_old {
                Containment::NewInOld {
                    new: query.clone(),
                    old: target.clone(),
                }
            } else {
                Containment::NewInNew {
                    inner: query.clone(),
                    outer: target.clone(),
                }
            });
        }
        if covered(tend.saturating_sub(tstart), tlen) {
            edges.push(if targets_old {
                Containment::OldInNew {
                    old: target,
                    new: query,
                }
            } else {
                Containment::NewInNew {
                    inner: target,
                    outer: query,
                }
            });
        }
    }
    Ok(())
}

#[derive(Debug, Default, PartialEq, Eq)]
pub struct Resolution {
    pub accepted: Vec<String>,
    /// (new contig, the contig it is contained in)
    pub redundant: Vec<(String, String)>,
    /// (accepted contig, the new contig that contains it)
    pub superseded: Vec<(String, String)>,
}

/// Decides candidates in order of decreasing length, then ID.
pub fn resolve(candidates: &[(String, u64)], edges: &[Containment]) -> Resolution {
    let mut order: Vec<&(String, u64)> = candidates.iter().collect();
    order.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
    let rank: HashMap<&str, usize> = order
        .iter()
        .enumerate()
        .map(|(index, (id, _))| (id.as_str(), index))
        .collect();
    let mut in_old: BTreeMap<&str, Vec<&str>> = BTreeMap::new();
    let mut in_new: BTreeMap<&str, Vec<&str>> = BTreeMap::new();
    let mut old_in: BTreeMap<&str, Vec<&str>> = BTreeMap::new();
    for edge in edges {
        match edge {
            Containment::NewInOld { new, old } => in_old.entry(new).or_default().push(old),
            Containment::NewInNew { inner, outer } => in_new.entry(inner).or_default().push(outer),
            Containment::OldInNew { old, new } => old_in.entry(old).or_default().push(new),
        }
    }
    let mut resolution = Resolution::default();
    let mut accepted: HashSet<&str> = HashSet::new();
    for (id, _) in &order {
        let own_rank = rank[id.as_str()];
        let old_container = in_old.get(id.as_str()).and_then(|olds| olds.iter().min());
        let new_container = in_new
            .get(id.as_str())
            .into_iter()
            .flatten()
            .filter(|outer| rank.get(*outer).is_some_and(|&r| r < own_rank))
            .min_by_key(|outer| rank[*outer]);
        match old_container.or(new_container) {
            Some(container) => resolution
                .redundant
                .push((id.clone(), (*container).to_owned())),
            None => {
                accepted.insert(id.as_str());
                resolution.accepted.push(id.clone());
            }
        }
    }
    for (old, containers) in old_in {
        if let Some(new) = containers
            .iter()
            .filter(|new| accepted.contains(*new))
            .min_by_key(|new| rank[*new])
        {
            resolution
                .superseded
                .push((old.to_owned(), (*new).to_owned()));
        }
    }
    resolution
}

fn read_counts(dir: &Path) -> Result<Value, String> {
    let text = fs::read_to_string(dir.join(CANDIDATE_COUNTS))
        .map_err(|error| format!("cannot read {CANDIDATE_COUNTS}: {error}"))?;
    serde_json::from_str(&text).map_err(|error| error.to_string())
}

/// Resolves containment, writes the new cumulative set, the contig table,
/// and `accept.json`, and removes the candidate FASTA (a derived copy).
pub fn finish(dir: &Path, previous: Option<&Path>, round: u32) -> Result<Value, String> {
    let counts = read_counts(dir)?;
    let mut candidates = Vec::new();
    let mut reader = FastaReader::open(&dir.join(CANDIDATES))?;
    while let Some(record) = reader.next_record()? {
        let length = record.length();
        candidates.push((record.id, length));
    }
    let mut edges = Vec::new();
    for (name, targets_old) in [(NEW_VS_NEW, false), (NEW_VS_OLD, true)] {
        let path = dir.join(name);
        if path.is_file() {
            let file = File::open(&path).map_err(|error| format!("cannot read {name}: {error}"))?;
            parse_paf(BufReader::new(file), targets_old, &mut edges)?;
        }
    }
    let resolution = resolve(&candidates, &edges);
    let superseded: HashMap<&str, &str> = resolution
        .superseded
        .iter()
        .map(|(old, new)| (old.as_str(), new.as_str()))
        .collect();
    let accepted: HashSet<&str> = resolution.accepted.iter().map(String::as_str).collect();
    let path = dir.join(CUMULATIVE);
    let write_error = |error: std::io::Error| format!("cannot write {}: {error}", path.display());
    let mut out = BufWriter::new(File::create_new(&path).map_err(write_error)?);
    let mut previous_contigs = 0u64;
    let mut table = String::from("contig_id\tlength\tstatus\tcontained_in\n");
    if let Some(previous) = previous {
        let mut reader = FastaReader::open(previous)?;
        while let Some(record) = reader.next_record()? {
            previous_contigs += 1;
            match superseded.get(record.id.as_str()) {
                Some(new) => {
                    table.push_str(&format!(
                        "{}\t{}\tsuperseded\t{new}\n",
                        record.id,
                        record.length()
                    ));
                }
                None => record.write_as(&mut out, &record.id).map_err(write_error)?,
            }
        }
    }
    let redundant: HashMap<&str, &str> = resolution
        .redundant
        .iter()
        .map(|(id, container)| (id.as_str(), container.as_str()))
        .collect();
    let mut reader = FastaReader::open(&dir.join(CANDIDATES))?;
    while let Some(record) = reader.next_record()? {
        if accepted.contains(record.id.as_str()) {
            record.write_as(&mut out, &record.id).map_err(write_error)?;
            table.push_str(&format!("{}\t{}\taccepted\t\n", record.id, record.length()));
        } else {
            let container = redundant
                .get(record.id.as_str())
                .copied()
                .unwrap_or_default();
            table.push_str(&format!(
                "{}\t{}\tredundant\t{container}\n",
                record.id,
                record.length()
            ));
        }
    }
    out.flush().map_err(write_error)?;
    drop(out);
    fs::write(dir.join(TABLE), table).map_err(|error| format!("cannot write {TABLE}: {error}"))?;
    fs::remove_file(dir.join(CANDIDATES))
        .map_err(|error| format!("cannot remove {CANDIDATES}: {error}"))?;
    let stats = file_stats(&path)?;
    let cumulative = stats.contigs;
    let summary = json!({
        "schema_version": SCHEMA_VERSION,
        "producer": crate::execution::run_dir::PRODUCER,
        "round": round,
        "policy": {
            "min_contig_length": counts["min_length"],
            "containment": {
                "tool": "minimap2",
                "preset": "asm5",
                "min_identity": CONTAINMENT_MIN_IDENTITY,
                "min_coverage_of_contained": CONTAINMENT_MIN_COVERAGE,
            },
        },
        "counts": {
            "previous": previous_contigs,
            "offered": counts["offered"],
            "short": counts["short"],
            "candidates": candidates.len(),
            "redundant": resolution.redundant.len(),
            "accepted": resolution.accepted.len(),
            "superseded": resolution.superseded.len(),
            "cumulative": cumulative,
        },
        "changed": !resolution.accepted.is_empty(),
        "cumulative": stats,
        "cumulative_sha256": sha256_file(&path).map_err(|error| error.to_string())?,
    });
    let text = serde_json::to_string_pretty(&summary).map_err(|error| error.to_string())?;
    fs::write(dir.join(SUMMARY), text + "\n").map_err(|error| error.to_string())?;
    Ok(summary)
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use super::*;

    /// (query, qlen, qstart, qend, target, tlen, tstart, tend, matches, block)
    type Hit<'a> = (&'a str, u64, u64, u64, &'a str, u64, u64, u64, u64, u64);

    fn paf(lines: &[Hit<'_>]) -> String {
        lines
            .iter()
            .map(|(q, ql, qs, qe, t, tl, ts, te, m, b)| {
                format!("{q}\t{ql}\t{qs}\t{qe}\t+\t{t}\t{tl}\t{ts}\t{te}\t{m}\t{b}\t60\n")
            })
            .collect()
    }

    fn ids(pairs: &[(String, String)]) -> Vec<(&str, &str)> {
        pairs
            .iter()
            .map(|(a, b)| (a.as_str(), b.as_str()))
            .collect()
    }

    #[test]
    fn paf_containment_needs_identity_and_coverage_of_the_contained_contig() {
        let text = paf(&[
            // new n1 (1000) fully inside old o1 (5000) at 99% identity
            ("n1", 1000, 0, 1000, "o1", 5000, 100, 1100, 990, 1000),
            // n2 covers all of old o2: o2 inside n2
            ("n2", 4000, 500, 1500, "o2", 1000, 0, 1000, 1000, 1000),
            // too divergent
            ("n3", 1000, 0, 1000, "o3", 5000, 0, 1000, 900, 1000),
            // too little of n4 aligned
            ("n4", 1000, 0, 900, "o4", 5000, 0, 900, 900, 900),
        ]);
        let mut edges = Vec::new();
        parse_paf(Cursor::new(text), true, &mut edges).unwrap();
        assert_eq!(
            edges,
            vec![
                Containment::NewInOld {
                    new: "n1".into(),
                    old: "o1".into()
                },
                Containment::OldInNew {
                    old: "o2".into(),
                    new: "n2".into()
                },
            ]
        );
        let mut edges = Vec::new();
        let text = paf(&[
            ("a", 1000, 0, 1000, "b", 1000, 0, 1000, 1000, 1000),
            ("a", 1000, 0, 1000, "a", 1000, 0, 1000, 1000, 1000),
        ]);
        parse_paf(Cursor::new(text), false, &mut edges).unwrap();
        assert_eq!(edges.len(), 2, "mutual containment, self hit ignored");
        assert!(parse_paf(Cursor::new("a\t10\n"), false, &mut Vec::new()).is_err());
    }

    #[test]
    fn redundant_new_contigs_are_dropped_and_contained_old_contigs_superseded() {
        let candidates = vec![
            ("long".to_owned(), 5000),
            ("inside-long".to_owned(), 1200),
            ("inside-old".to_owned(), 1500),
            ("twin-a".to_owned(), 2000),
            ("twin-b".to_owned(), 2000),
            ("novel".to_owned(), 1000),
        ];
        let edges = vec![
            Containment::NewInNew {
                inner: "inside-long".into(),
                outer: "long".into(),
            },
            Containment::NewInOld {
                new: "inside-old".into(),
                old: "old-2".into(),
            },
            Containment::NewInOld {
                new: "inside-old".into(),
                old: "old-1".into(),
            },
            Containment::NewInNew {
                inner: "twin-a".into(),
                outer: "twin-b".into(),
            },
            Containment::NewInNew {
                inner: "twin-b".into(),
                outer: "twin-a".into(),
            },
            Containment::OldInNew {
                old: "old-3".into(),
                new: "long".into(),
            },
            Containment::OldInNew {
                old: "old-4".into(),
                new: "inside-old".into(),
            },
        ];
        let resolution = resolve(&candidates, &edges);
        assert_eq!(resolution.accepted, ["long", "twin-a", "novel"]);
        assert_eq!(
            ids(&resolution.redundant),
            [
                ("twin-b", "twin-a"),
                ("inside-old", "old-1"),
                ("inside-long", "long")
            ]
        );
        // old-4 is inside a contig that was itself redundant, so it stays.
        assert_eq!(ids(&resolution.superseded), [("old-3", "long")]);
    }

    #[test]
    fn finishing_rewrites_the_cumulative_set_and_reports_counts() {
        let dir = std::env::temp_dir().join(format!("mca-accept-{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        let seq = |n: usize| "ACGT".repeat(n);
        fs::write(
            dir.join("previous.fa"),
            format!(">old-1\n{}\n>old-2\n{}\n", seq(300), seq(400)),
        )
        .unwrap();
        fs::write(
            dir.join("round.fa"),
            format!(">k1 len=2000\n{}\n>k2\n{}\n>k3\nACGT\n", seq(500), seq(300)),
        )
        .unwrap();
        let sources = [CandidateSource {
            path: dir.join("round.fa"),
            prefix: "r02.g__".into(),
        }];
        assert_eq!(write_candidates(&dir, &sources, 1000).unwrap(), 2);
        fs::write(
            dir.join(NEW_VS_OLD),
            paf(&[
                (
                    "r02.g__k1",
                    2000,
                    400,
                    1600,
                    "old-1",
                    1200,
                    0,
                    1200,
                    1200,
                    1200,
                ),
                (
                    "r02.g__k2",
                    1200,
                    0,
                    1200,
                    "old-2",
                    1600,
                    0,
                    1200,
                    1200,
                    1200,
                ),
            ]),
        )
        .unwrap();
        let summary = finish(&dir, Some(&dir.join("previous.fa")), 2).unwrap();
        assert_eq!(summary["counts"]["offered"], 3);
        assert_eq!(summary["counts"]["short"], 1);
        assert_eq!(summary["counts"]["accepted"], 1);
        assert_eq!(summary["counts"]["redundant"], 1);
        assert_eq!(summary["counts"]["superseded"], 1);
        assert_eq!(summary["counts"]["cumulative"], 2);
        assert_eq!(summary["changed"], true);
        let cumulative = fs::read_to_string(dir.join(CUMULATIVE)).unwrap();
        assert!(cumulative.starts_with(">old-2\n"));
        assert!(cumulative.contains(">r02.g__k1 len=2000\n"));
        assert!(!cumulative.contains("old-1"));
        let table = fs::read_to_string(dir.join(TABLE)).unwrap();
        assert!(table.contains("old-1\t1200\tsuperseded\tr02.g__k1\n"));
        assert!(table.contains("r02.g__k2\t1200\tredundant\told-2\n"));
        assert!(!dir.join(CANDIDATES).exists());
        fs::remove_dir_all(&dir).unwrap();
    }
}
