//! Pair classification against the cumulative accepted contigs.
//!
//! - `explained` (residual-only): both mates have an alignment, primary or
//!   secondary, with identity and aligned fraction at least the thresholds to
//!   any accepted contig. MAPQ, pairing, and which contig are ignored, since
//!   ambiguity among overlapping contigs is no reason to reassemble.
//! - `certain` (recruitment): both primaries pass the thresholds with known
//!   MAPQ at least the minimum, the mates are properly paired on one contig,
//!   there is no supplementary record, and no qualifying secondary alignment
//!   lies on another contig.
//!
//! Every other pair stays in the next round's pool with both mates. The
//! classifier streams the pool and the alignments in lockstep and copies
//! records verbatim, so pair identity and sample identity are preserved.

use std::{
    collections::BTreeMap,
    io::{BufRead, BufReader, Read, Write},
    path::Path,
};

use serde_json::{Value, json};

use crate::{
    cli::PairRule,
    discovery::{
        decide::{Alignment, PROPER_PAIR, SECONDARY, SUPPLEMENTARY, UNMAPPED, measure},
        extract::{FastqReader, SamGroups, create, finish, open_text, sequence_length},
    },
    error::Error,
};

pub const POOL_R1: &str = "pool_R1.fastq.gz";
pub const POOL_R2: &str = "pool_R2.fastq.gz";
pub const DECISIONS: &str = "decisions.tsv.gz";
pub const SUMMARY: &str = "classification.json";
pub const SCHEMA_VERSION: u32 = 1;
const DECISION_COLUMNS: &str = "pair\tread_name\tresolved\treason\tcontig\n";

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Thresholds {
    pub min_identity: f64,
    pub min_aligned_fraction: f64,
    pub min_mapq: u8,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Classification {
    /// The pair needs no further assembly.
    pub resolved: bool,
    pub reason: &'static str,
    /// Mate 1's primary contig, when mapped.
    pub contig: Option<String>,
}

fn qualifies(record: &Alignment, thresholds: &Thresholds, allow_hard: bool) -> bool {
    record.flag & UNMAPPED == 0
        && measure(record, allow_hard).is_some_and(|m| {
            m.identity >= thresholds.min_identity
                && m.aligned_fraction >= thresholds.min_aligned_fraction
        })
}

/// Classifies one pair from all of its alignment records.
///
/// Returns `Err` for structurally impossible input, such as a mate without
/// exactly one primary record.
pub fn classify(
    records: &[Alignment],
    rule: PairRule,
    thresholds: &Thresholds,
) -> Result<Classification, String> {
    let mut primaries: [Option<&Alignment>; 2] = [None, None];
    for record in records.iter().filter(|r| r.is_primary()) {
        let mate = record
            .mate()
            .ok_or("alignment record is not marked as mate 1 or 2")?;
        if primaries[usize::from(mate - 1)].replace(record).is_some() {
            return Err(format!(
                "mate {mate} has more than one primary alignment; adjacent pairs may share a read name"
            ));
        }
    }
    let [Some(first), Some(second)] = primaries else {
        return Err("a mate has no primary alignment record".into());
    };
    let contig = (first.flag & UNMAPPED == 0).then(|| first.contig.clone());
    let outcome = |resolved, reason| {
        Ok(Classification {
            resolved,
            reason,
            contig: contig.clone(),
        })
    };
    match rule {
        PairRule::Explained => {
            let mate_explained = |mate: u8| {
                records.iter().any(|r| {
                    r.mate() == Some(mate)
                        && r.flag & SUPPLEMENTARY == 0
                        && qualifies(r, thresholds, true)
                })
            };
            if mate_explained(1) && mate_explained(2) {
                outcome(true, "explained")
            } else if first.flag & UNMAPPED != 0 || second.flag & UNMAPPED != 0 {
                outcome(false, "mate-unmapped")
            } else {
                outcome(false, "below-thresholds")
            }
        }
        PairRule::Certain => {
            if first.flag & UNMAPPED != 0 || second.flag & UNMAPPED != 0 {
                return outcome(false, "mate-unmapped");
            }
            for primary in [first, second] {
                if primary.nm.is_none() {
                    return outcome(false, "missing-nm");
                }
                let Some(measured) = measure(primary, false) else {
                    return outcome(false, "unsupported-cigar");
                };
                if measured.identity < thresholds.min_identity {
                    return outcome(false, "low-identity");
                }
                if measured.aligned_fraction < thresholds.min_aligned_fraction {
                    return outcome(false, "low-aligned-fraction");
                }
            }
            if first.mapq == 255 || second.mapq == 255 {
                return outcome(false, "mapq-unavailable");
            }
            if first.mapq < thresholds.min_mapq || second.mapq < thresholds.min_mapq {
                return outcome(false, "low-mapq");
            }
            if records.iter().any(|r| r.flag & SUPPLEMENTARY != 0) {
                return outcome(false, "supplementary-alignment");
            }
            if first.flag & PROPER_PAIR == 0 || second.flag & PROPER_PAIR == 0 {
                return outcome(false, "discordant-pair");
            }
            if first.contig != second.contig {
                return outcome(false, "mates-different-contigs");
            }
            let elsewhere = records.iter().any(|r| {
                r.flag & SECONDARY != 0
                    && r.contig != first.contig
                    && qualifies(r, thresholds, true)
            });
            if elsewhere {
                return outcome(false, "alternative-contig");
            }
            outcome(true, "certain")
        }
    }
}

pub struct Request<'a> {
    pub sample_id: &'a str,
    pub read1: &'a Path,
    pub read2: &'a Path,
    pub rule: PairRule,
    pub thresholds: Thresholds,
    /// Decisions of the previous check over the same pool, row for row.
    pub previous_decisions: Option<&'a Path>,
    /// When false, no alignments are read and every pair stays unresolved.
    pub alignments: bool,
    pub output_dir: &'a Path,
}

fn invalid(message: String) -> Error {
    Error::InvalidInput(message)
}

/// Rows of a previous decisions table: (read name, resolved).
struct PreviousDecisions {
    input: Box<dyn BufRead>,
    path: String,
    line: String,
}

impl PreviousDecisions {
    fn open(path: &Path) -> Result<Self, Error> {
        let mut table = Self {
            input: open_text(path)?,
            path: path.display().to_string(),
            line: String::new(),
        };
        match table.read_line()? {
            Some(header) if format!("{header}\n") == DECISION_COLUMNS => Ok(table),
            _ => Err(invalid(format!("{}: not a decisions table", table.path))),
        }
    }

    fn read_line(&mut self) -> Result<Option<String>, Error> {
        self.line.clear();
        let read = self
            .input
            .read_line(&mut self.line)
            .map_err(|source| Error::Io {
                context: format!("cannot read {}", self.path),
                source,
            })?;
        Ok((read > 0).then(|| self.line.trim_end_matches(['\n', '\r']).to_owned()))
    }

    fn next(&mut self) -> Result<Option<(String, bool)>, Error> {
        let Some(line) = self.read_line()? else {
            return Ok(None);
        };
        let fields: Vec<&str> = line.split('\t').collect();
        match fields.as_slice() {
            [_, name, "1", ..] => Ok(Some(((*name).to_owned(), true))),
            [_, name, "0", ..] => Ok(Some(((*name).to_owned(), false))),
            _ => Err(invalid(format!("{}: malformed row {line:?}", self.path))),
        }
    }
}

/// Classifies every pair of a pool, writes the unresolved pairs as the next
/// pool, and reports counts, progress against the previous round's input,
/// and whether the next input equals the previous one.
pub fn classify_pool(request: &Request<'_>, alignments: impl Read) -> Result<Value, Error> {
    let mut read1 = FastqReader::open(request.read1)?;
    let mut read2 = FastqReader::open(request.read2)?;
    let mut groups = SamGroups::new(BufReader::with_capacity(1 << 20, alignments));
    let mut previous = request
        .previous_decisions
        .map(PreviousDecisions::open)
        .transpose()?;
    let path = |name: &str| request.output_dir.join(name);
    let (mut out1, mut out2, mut table) = (
        create(&path(POOL_R1))?,
        create(&path(POOL_R2))?,
        create(&path(DECISIONS))?,
    );
    let write_error = |source| Error::Io {
        context: format!("cannot write outputs in {}", request.output_dir.display()),
        source,
    };
    table
        .write_all(DECISION_COLUMNS.as_bytes())
        .map_err(write_error)?;
    let (mut pairs, mut resolved, mut bases_in, mut bases_retained) = (0u64, 0u64, 0u64, 0u64);
    let (mut previous_input, mut previous_input_resolved) = (0u64, 0u64);
    let mut unchanged = true;
    let mut reasons: BTreeMap<&str, u64> = BTreeMap::new();
    loop {
        let (first, second) = match (read1.next_record()?, read2.next_record()?) {
            (None, None) => break,
            (Some(first), Some(second)) => (first, second),
            _ => {
                return Err(invalid(format!(
                    "{} and {} have different numbers of records after {pairs} pairs",
                    request.read1.display(),
                    request.read2.display()
                )));
            }
        };
        pairs += 1;
        if first.name != second.name {
            return Err(invalid(format!(
                "pair {pairs}: mate names differ ({:?} and {:?})",
                first.name, second.name
            )));
        }
        let pair_bases = sequence_length(&first) + sequence_length(&second);
        bases_in += pair_bases;
        let decision = if request.alignments {
            let Some((name, records)) = groups.next_group()? else {
                return Err(invalid(format!(
                    "alignments ended before pair {pairs} ({:?})",
                    first.name
                )));
            };
            if name != first.name {
                return Err(invalid(format!(
                    "alignment order differs from the FASTQ at pair {pairs}: expected {:?}, found {name:?}",
                    first.name
                )));
            }
            classify(&records, request.rule, &request.thresholds)
                .map_err(|problem| invalid(format!("pair {pairs} ({name:?}): {problem}")))?
        } else {
            Classification {
                resolved: false,
                reason: "not-mapped",
                contig: None,
            }
        };
        // Without a previous table, the previous round's input is this pool.
        let was_resolved = match previous.as_mut() {
            None => false,
            Some(table) => match table.next()? {
                Some((name, was)) if name == first.name => was,
                Some((name, _)) => {
                    return Err(invalid(format!(
                        "pair {pairs}: previous decisions list {name:?}, the pool has {:?}",
                        first.name
                    )));
                }
                None => {
                    return Err(invalid(format!(
                        "previous decisions end before pair {pairs}"
                    )));
                }
            },
        };
        if !was_resolved {
            previous_input += 1;
            previous_input_resolved += u64::from(decision.resolved);
        }
        unchanged &= was_resolved == decision.resolved;
        *reasons.entry(decision.reason).or_default() += 1;
        if decision.resolved {
            resolved += 1;
        } else {
            bases_retained += pair_bases;
            for line in &first.lines {
                out1.write_all(line).map_err(write_error)?;
            }
            for line in &second.lines {
                out2.write_all(line).map_err(write_error)?;
            }
        }
        writeln!(
            table,
            "{pairs}\t{}\t{}\t{}\t{}",
            first.name,
            u8::from(decision.resolved),
            decision.reason,
            decision.contig.unwrap_or_default()
        )
        .map_err(write_error)?;
    }
    if request.alignments {
        if let Some((name, _)) = groups.next_group()? {
            return Err(invalid(format!(
                "alignments continue after the last FASTQ pair with {name:?}"
            )));
        }
    }
    if let Some(table) = previous.as_mut() {
        if let Some((name, _)) = table.next()? {
            return Err(invalid(format!(
                "previous decisions continue after the last pair with {name:?}"
            )));
        }
    }
    let digest1 = finish(out1, &path(POOL_R1))?;
    let digest2 = finish(out2, &path(POOL_R2))?;
    let digest_table = finish(table, &path(DECISIONS))?;
    let summary = json!({
        "schema_version": SCHEMA_VERSION,
        "producer": crate::execution::run_dir::PRODUCER,
        "sample_id": request.sample_id,
        "rule": request.rule,
        "thresholds": {
            "min_identity": request.thresholds.min_identity,
            "min_aligned_fraction": request.thresholds.min_aligned_fraction,
            "min_mapq": request.thresholds.min_mapq,
        },
        "mode": if request.alignments { "alignments" } else { "no-alignments" },
        "counts": {
            "pairs_in": pairs,
            "pairs_resolved": resolved,
            "pairs_retained": pairs - resolved,
            "bases_in": bases_in,
            "bases_retained": bases_retained,
        },
        "reasons": reasons,
        "previous_input": {"pairs": previous_input, "resolved_now": previous_input_resolved},
        "input_unchanged": unchanged,
        "outputs": {"read1": POOL_R1, "read2": POOL_R2, "decisions": DECISIONS},
        "sha256": {"read1": digest1, "read2": digest2, "decisions": digest_table},
    });
    let summary_path = path(SUMMARY);
    let mut file = std::fs::File::create_new(&summary_path).map_err(|source| Error::Io {
        context: format!("cannot create {}", summary_path.display()),
        source,
    })?;
    serde_json::to_writer_pretty(&mut file, &summary).map_err(Error::Json)?;
    file.write_all(b"\n").map_err(write_error)?;
    Ok(summary)
}

#[cfg(test)]
mod tests {
    use std::{
        fs,
        io::Cursor,
        path::PathBuf,
        sync::atomic::{AtomicU64, Ordering},
    };

    use flate2::read::GzDecoder;

    use super::*;

    const T: Thresholds = Thresholds {
        min_identity: 0.98,
        min_aligned_fraction: 0.9,
        min_mapq: 20,
    };

    fn record(flag: u16, contig: &str, mapq: u8, cigar: &str, nm: u32) -> Alignment {
        Alignment {
            flag,
            contig: contig.into(),
            mapq,
            cigar: cigar.into(),
            nm: Some(nm),
        }
    }

    /// A properly paired, uniquely placed pair on `contig`.
    fn proper(contig: &str) -> Vec<Alignment> {
        vec![
            record(99, contig, 60, "100M", 0),
            record(147, contig, 60, "100M", 0),
        ]
    }

    fn unmapped() -> Vec<Alignment> {
        vec![record(77, "*", 0, "*", 0), record(141, "*", 0, "*", 0)]
    }

    fn reason(records: &[Alignment], rule: PairRule) -> (&'static str, bool) {
        let c = classify(records, rule, &T).unwrap();
        (c.reason, c.resolved)
    }

    #[test]
    fn explained_needs_both_mates_on_any_contig_and_ignores_ambiguity() {
        assert_eq!(
            reason(&proper("c1"), PairRule::Explained),
            ("explained", true)
        );
        // Mates on different contigs with MAPQ 0 still count as explained.
        let split = vec![
            record(65, "c1", 0, "100M", 0),
            record(129, "c2", 0, "100M", 1),
        ];
        assert_eq!(reason(&split, PairRule::Explained), ("explained", true));
        // A qualifying secondary explains a mate whose primary falls short.
        let mut rescued = proper("c1");
        rescued[0].nm = Some(5);
        assert_eq!(
            reason(&rescued, PairRule::Explained),
            ("below-thresholds", false)
        );
        rescued.push(record(256 + 64 + 1, "c2", 0, "100M", 0));
        assert_eq!(reason(&rescued, PairRule::Explained), ("explained", true));
        let mut clipped = proper("c1");
        clipped[1].cigar = "80M20S".into();
        assert_eq!(
            reason(&clipped, PairRule::Explained),
            ("below-thresholds", false)
        );
        assert_eq!(
            reason(&unmapped(), PairRule::Explained),
            ("mate-unmapped", false)
        );
        let mut half = proper("c1");
        half[1] = record(133, "*", 0, "*", 0);
        assert_eq!(reason(&half, PairRule::Explained), ("mate-unmapped", false));
    }

    #[test]
    fn certain_needs_a_unique_proper_pair_on_one_contig() {
        assert_eq!(reason(&proper("c1"), PairRule::Certain), ("certain", true));
        let mut same_contig_repeat = proper("c1");
        same_contig_repeat.push(record(256 + 64 + 1, "c1", 0, "100M", 0));
        assert_eq!(
            reason(&same_contig_repeat, PairRule::Certain),
            ("certain", true)
        );
        let mut elsewhere = proper("c1");
        elsewhere.push(record(256 + 128 + 1, "c2", 0, "100M", 1));
        assert_eq!(
            reason(&elsewhere, PairRule::Certain),
            ("alternative-contig", false)
        );
        let mut weak_elsewhere = proper("c1");
        weak_elsewhere.push(record(256 + 128 + 1, "c2", 0, "100M", 9));
        assert_eq!(
            reason(&weak_elsewhere, PairRule::Certain),
            ("certain", true)
        );
        let split = vec![
            record(97, "c1", 60, "100M", 0),
            record(145, "c2", 60, "100M", 0),
        ];
        assert_eq!(
            reason(&split, PairRule::Certain),
            ("discordant-pair", false)
        );
        let mut split_proper = split.clone();
        split_proper[0].flag |= 2;
        split_proper[1].flag |= 2;
        assert_eq!(
            reason(&split_proper, PairRule::Certain),
            ("mates-different-contigs", false)
        );
        let mut low = proper("c1");
        low[1].mapq = 3;
        assert_eq!(reason(&low, PairRule::Certain), ("low-mapq", false));
        low[1].mapq = 255;
        assert_eq!(reason(&low, PairRule::Certain), ("mapq-unavailable", false));
        let mut divergent = proper("c1");
        divergent[0].nm = Some(3);
        assert_eq!(
            reason(&divergent, PairRule::Certain),
            ("low-identity", false)
        );
        let mut clipped = proper("c1");
        clipped[0].cigar = "50S50M".into();
        assert_eq!(
            reason(&clipped, PairRule::Certain),
            ("low-aligned-fraction", false)
        );
        let mut chimeric = proper("c1");
        chimeric.push(record(2048 + 64 + 1, "c2", 60, "30M70H", 0));
        assert_eq!(
            reason(&chimeric, PairRule::Certain),
            ("supplementary-alignment", false)
        );
        let mut missing = proper("c1");
        missing[0].nm = None;
        assert_eq!(reason(&missing, PairRule::Certain), ("missing-nm", false));
        assert_eq!(
            reason(&unmapped(), PairRule::Certain),
            ("mate-unmapped", false)
        );
        assert!(classify(&proper("c1")[..1], PairRule::Certain, &T).is_err());
        let doubled = [proper("c1"), proper("c1")].concat();
        assert!(classify(&doubled, PairRule::Explained, &T).is_err());
    }

    static NEXT: AtomicU64 = AtomicU64::new(0);

    fn dir() -> PathBuf {
        let path = std::env::temp_dir().join(format!(
            "mca-classify-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        let _ = fs::remove_dir_all(&path);
        fs::create_dir_all(path.join("out")).unwrap();
        path
    }

    fn pool(dir: &Path, names: &[&str]) -> (PathBuf, PathBuf) {
        let write = |mate: u8| {
            let path = dir.join(format!("in_R{mate}.fastq"));
            let text: String = names
                .iter()
                .map(|n| format!("@{n}/{mate}\nACGTACGTAC\n+\nIIIIIIIIII\n"))
                .collect();
            fs::write(&path, text).unwrap();
            path
        };
        (write(1), write(2))
    }

    fn sam(pairs: &[(&str, bool)]) -> String {
        pairs
            .iter()
            .map(|(name, mapped)| {
                if *mapped {
                    format!(
                        "{name}\t99\tc1\t1\t60\t10M\t=\t11\t20\t*\t*\tNM:i:0\n{name}\t147\tc1\t11\t60\t10M\t=\t1\t-20\t*\t*\tNM:i:0\n"
                    )
                } else {
                    format!("{name}\t77\t*\t0\t0\t*\t*\t0\t0\t*\t*\n{name}\t141\t*\t0\t0\t*\t*\t0\t0\t*\t*\n")
                }
            })
            .collect()
    }

    fn gunzip(path: &Path) -> String {
        let mut text = String::new();
        GzDecoder::new(fs::File::open(path).unwrap())
            .read_to_string(&mut text)
            .unwrap();
        text
    }

    fn request<'a>(
        read1: &'a Path,
        read2: &'a Path,
        output_dir: &'a Path,
        previous: Option<&'a Path>,
    ) -> Request<'a> {
        Request {
            sample_id: "S1",
            read1,
            read2,
            rule: PairRule::Certain,
            thresholds: T,
            previous_decisions: previous,
            alignments: true,
            output_dir,
        }
    }

    #[test]
    fn pools_keep_unresolved_pairs_verbatim_and_report_progress_against_the_previous_input() {
        let dir = dir();
        let (r1, r2) = pool(&dir, &["a", "b", "c", "d"]);
        let first = dir.join("out");
        let summary = classify_pool(
            &request(&r1, &r2, &first, None),
            Cursor::new(sam(&[
                ("a", true),
                ("b", false),
                ("c", false),
                ("d", false),
            ])),
        )
        .unwrap();
        assert_eq!(summary["counts"]["pairs_resolved"], 1);
        assert_eq!(summary["counts"]["pairs_retained"], 3);
        assert_eq!(summary["counts"]["bases_retained"], 60);
        assert_eq!(
            summary["previous_input"],
            json!({"pairs": 4, "resolved_now": 1})
        );
        assert_eq!(summary["input_unchanged"], false);
        assert_eq!(
            gunzip(&first.join(POOL_R1)),
            "@b/1\nACGTACGTAC\n+\nIIIIIIIIII\n@c/1\nACGTACGTAC\n+\nIIIIIIIIII\n@d/1\nACGTACGTAC\n+\nIIIIIIIIII\n"
        );
        assert!(gunzip(&first.join(POOL_R2)).starts_with("@b/2\n"));
        let table = gunzip(&first.join(DECISIONS));
        assert!(table.starts_with(DECISION_COLUMNS));
        assert!(table.contains("1\ta\t1\tcertain\tc1\n"));
        assert!(table.contains("2\tb\t0\tmate-unmapped\t\n"));

        // Recruitment re-checks the same pool: c is newly certain, a is not
        // any more. Progress counts only the previous input (b, c, d).
        let second = dir.join("second");
        fs::create_dir(&second).unwrap();
        let decisions = first.join(DECISIONS);
        let summary = classify_pool(
            &request(&r1, &r2, &second, Some(&decisions)),
            Cursor::new(sam(&[
                ("a", false),
                ("b", false),
                ("c", true),
                ("d", false),
            ])),
        )
        .unwrap();
        assert_eq!(
            summary["previous_input"],
            json!({"pairs": 3, "resolved_now": 1})
        );
        assert_eq!(summary["input_unchanged"], false);
        assert_eq!(summary["counts"]["pairs_retained"], 3);

        // The same decisions again: the next input equals the previous one.
        let third = dir.join("third");
        fs::create_dir(&third).unwrap();
        let decisions = second.join(DECISIONS);
        let summary = classify_pool(
            &request(&r1, &r2, &third, Some(&decisions)),
            Cursor::new(sam(&[
                ("a", false),
                ("b", false),
                ("c", true),
                ("d", false),
            ])),
        )
        .unwrap();
        assert_eq!(summary["input_unchanged"], true);
        assert_eq!(
            summary["previous_input"],
            json!({"pairs": 3, "resolved_now": 0})
        );
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn mismatched_alignments_previous_tables_and_empty_pools() {
        let dir = dir();
        let (r1, r2) = pool(&dir, &["a", "b"]);
        let out = |name: &str| {
            let path = dir.join(name);
            fs::create_dir(&path).unwrap();
            path
        };
        let o = out("order");
        let error = classify_pool(
            &request(&r1, &r2, &o, None),
            Cursor::new(sam(&[("b", true), ("a", true)])),
        )
        .unwrap_err();
        assert!(error.to_string().contains("alignment order differs"));
        let o = out("short");
        let error = classify_pool(
            &request(&r1, &r2, &o, None),
            Cursor::new(sam(&[("a", true)])),
        )
        .unwrap_err();
        assert!(error.to_string().contains("alignments ended"));
        let table = dir.join("previous.tsv");
        fs::write(
            &table,
            format!("{DECISION_COLUMNS}1\ta\t1\tcertain\tc1\n2\tz\t0\tx\t\n"),
        )
        .unwrap();
        let o = out("previous");
        let error = classify_pool(
            &request(&r1, &r2, &o, Some(&table)),
            Cursor::new(sam(&[("a", true), ("b", true)])),
        )
        .unwrap_err();
        assert!(error.to_string().contains("previous decisions list \"z\""));
        fs::write(&table, "not a table\n").unwrap();
        let o = out("bad-table");
        assert!(classify_pool(&request(&r1, &r2, &o, Some(&table)), Cursor::new("")).is_err());

        let (e1, e2) = pool(&dir.join("previous"), &[]);
        let o = out("empty");
        let mut empty = request(&e1, &e2, &o, None);
        empty.alignments = false;
        let summary = classify_pool(&empty, std::io::empty()).unwrap();
        assert_eq!(summary["counts"]["pairs_in"], 0);
        assert_eq!(summary["input_unchanged"], true);
        assert_eq!(summary["previous_input"]["pairs"], 0);
        fs::remove_dir_all(dir).unwrap();
    }
}
