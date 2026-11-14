//! Lockstep extraction of discovery pairs.
//!
//! The mapper's SAM output is read from stdin in input order, one group of
//! records per read pair, alongside the original R1/R2 FASTQ files. Each pair's
//! name is checked against its alignment group, so any reordering fails
//! instead of mixing reads. Retained pairs are copied verbatim (header,
//! sequence, separator, and quality lines) to gzip FASTQ. Memory is bounded
//! by one pair's records; nothing is collated or indexed by read name.

use std::{
    collections::BTreeMap,
    fs::File,
    io::{self, BufRead, BufReader, Read, Write},
    path::{Path, PathBuf},
};

use flate2::{Compression, read::MultiGzDecoder, write::GzEncoder};
use serde_json::json;

use super::decide::{Alignment, Contigs, Decision, Outcome, Policy, decide};
use crate::{digest::HashingWriter, error::Error};

pub const DISCOVERY_R1: &str = "discovery_R1.fastq.gz";
pub const DISCOVERY_R2: &str = "discovery_R2.fastq.gz";
pub const DECISIONS: &str = "decisions.tsv.gz";
pub const SUMMARY: &str = "extraction.json";
pub const SCHEMA_VERSION: u32 = 1;
const DECISION_COLUMNS: &str =
    "pair_index\tread_name\toutcome\treason\treference_id\tspecies_taxid\tmapq1\tmapq2\n";

fn invalid(message: String) -> Error {
    Error::InvalidInput(message)
}

/// Opens plain or gzip input, detected from the first two bytes.
pub fn open_text(path: &Path) -> Result<Box<dyn BufRead>, Error> {
    let io_error = |source| Error::Io {
        context: format!("cannot read {}", path.display()),
        source,
    };
    let mut file = BufReader::with_capacity(1 << 20, File::open(path).map_err(io_error)?);
    let gzip = file
        .fill_buf()
        .map_err(io_error)?
        .starts_with(&[0x1f, 0x8b]);
    Ok(if gzip {
        Box::new(BufReader::with_capacity(1 << 20, MultiGzDecoder::new(file)))
    } else {
        Box::new(file)
    })
}

/// One FASTQ record exactly as read, plus its normalized pair name.
pub struct FastqRecord {
    pub lines: [Vec<u8>; 4],
    pub name: String,
}

/// Bases in a record's sequence line, excluding the line ending.
pub fn sequence_length(record: &FastqRecord) -> u64 {
    record.lines[1].trim_ascii_end().len() as u64
}

pub struct FastqReader<'a> {
    input: Box<dyn BufRead + 'a>,
    path: PathBuf,
    records: u64,
}

/// Read name without the leading `@`, comments, or a `/1` or `/2` suffix.
pub fn pair_name(header: &[u8]) -> String {
    let text = String::from_utf8_lossy(header);
    let token = text
        .trim_end_matches(['\n', '\r'])
        .trim_start_matches('@')
        .split_whitespace()
        .next()
        .unwrap_or_default()
        .to_owned();
    token
        .strip_suffix("/1")
        .or_else(|| token.strip_suffix("/2"))
        .unwrap_or(&token)
        .to_owned()
}

impl FastqReader<'static> {
    pub fn open(path: &Path) -> Result<Self, Error> {
        Ok(Self::from_reader(path, open_text(path)?))
    }
}

impl<'a> FastqReader<'a> {
    /// Reads already-decompressed text; `path` labels diagnostics.
    pub fn from_reader(path: &Path, input: Box<dyn BufRead + 'a>) -> Self {
        Self {
            input,
            path: path.to_owned(),
            records: 0,
        }
    }

    pub fn next_record(&mut self) -> Result<Option<FastqRecord>, Error> {
        let mut lines: [Vec<u8>; 4] = Default::default();
        for (index, line) in lines.iter_mut().enumerate() {
            let read = self
                .input
                .read_until(b'\n', line)
                .map_err(|source| Error::Io {
                    context: format!("cannot read {}", self.path.display()),
                    source,
                })?;
            if read == 0 {
                if index == 0 {
                    return Ok(None);
                }
                return Err(self.malformed("truncated record"));
            }
            if !line.ends_with(b"\n") {
                line.push(b'\n');
            }
        }
        let content = |line: &Vec<u8>| line.trim_ascii_end().len();
        if !lines[0].starts_with(b"@") || !lines[2].starts_with(b"+") {
            return Err(self.malformed("missing '@' header or '+' separator"));
        }
        if content(&lines[1]) == 0 || content(&lines[1]) != content(&lines[3]) {
            return Err(self.malformed("sequence and quality lengths differ"));
        }
        let name = pair_name(&lines[0]);
        if name.is_empty() {
            return Err(self.malformed("record has no name"));
        }
        self.records += 1;
        Ok(Some(FastqRecord { lines, name }))
    }

    fn malformed(&self, problem: &str) -> Error {
        invalid(format!(
            "{}: malformed FASTQ after record {}: {problem}",
            self.path.display(),
            self.records
        ))
    }
}

/// Groups consecutive SAM records that share a read name.
pub struct SamGroups<R: BufRead> {
    input: R,
    pending: Option<(String, Alignment)>,
    line: String,
    finished: bool,
}

impl<R: BufRead> SamGroups<R> {
    pub fn new(input: R) -> Self {
        Self {
            input,
            pending: None,
            line: String::new(),
            finished: false,
        }
    }

    fn next_record(&mut self) -> Result<Option<(String, Alignment)>, Error> {
        loop {
            self.line.clear();
            let read = self
                .input
                .read_line(&mut self.line)
                .map_err(|source| Error::Io {
                    context: "cannot read alignments from stdin".into(),
                    source,
                })?;
            if read == 0 {
                return Ok(None);
            }
            let line = self.line.trim_end_matches(['\n', '\r']);
            if line.is_empty() || line.starts_with('@') {
                continue;
            }
            return Alignment::parse(line).map(Some).map_err(invalid);
        }
    }

    pub fn next_group(&mut self) -> Result<Option<(String, Vec<Alignment>)>, Error> {
        if self.finished {
            return Ok(None);
        }
        let Some((name, first)) = (match self.pending.take() {
            Some(record) => Some(record),
            None => self.next_record()?,
        }) else {
            self.finished = true;
            return Ok(None);
        };
        let mut records = vec![first];
        loop {
            match self.next_record()? {
                Some((next, record)) if next == name => records.push(record),
                Some(other) => {
                    self.pending = Some(other);
                    break;
                }
                None => {
                    self.finished = true;
                    break;
                }
            }
        }
        Ok(Some((name, records)))
    }
}

/// Reads `contig, reference_id, species_taxid, ...` with a header row.
pub fn read_contigs(path: &Path) -> Result<Contigs, Error> {
    let text = std::fs::read_to_string(path).map_err(|source| Error::Io {
        context: format!("cannot read {}", path.display()),
        source,
    })?;
    let mut contigs = Contigs::new();
    for (number, line) in text.lines().enumerate().skip(1) {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 {
            return Err(invalid(format!(
                "{}: line {} is malformed",
                path.display(),
                number + 1
            )));
        }
        if contigs
            .insert(
                fields[0].to_owned(),
                (fields[1].to_owned(), fields[2].to_owned()),
            )
            .is_some()
        {
            return Err(invalid(format!(
                "{}: duplicate contig {:?}",
                path.display(),
                fields[0]
            )));
        }
    }
    Ok(contigs)
}

pub struct Request<'a> {
    pub sample_id: &'a str,
    pub read1: &'a Path,
    pub read2: &'a Path,
    /// `None` when the sample has no selected references: every pair is kept.
    pub contigs: Option<&'a Path>,
    pub policy: Policy,
    pub output_dir: &'a Path,
}

pub(crate) type Output = GzEncoder<HashingWriter<io::BufWriter<File>>>;

pub(crate) fn create(path: &Path) -> Result<Output, Error> {
    let file = File::create_new(path).map_err(|source| Error::Io {
        context: format!("cannot create {}", path.display()),
        source,
    })?;
    Ok(GzEncoder::new(
        HashingWriter::new(io::BufWriter::with_capacity(1 << 20, file)),
        Compression::fast(),
    ))
}

/// Completes the gzip stream and returns the SHA-256 of the file's bytes.
pub(crate) fn finish(encoder: Output, path: &Path) -> Result<String, Error> {
    let io_error = |source| Error::Io {
        context: format!("cannot finish {}", path.display()),
        source,
    };
    let (writer, digest) = encoder.finish().map_err(io_error)?.finish();
    writer
        .into_inner()
        .map_err(|error| io_error(error.into_error()))?
        .sync_all()
        .map_err(io_error)?;
    Ok(digest)
}

/// Runs the decision stream over `alignments` (SAM on stdin in practice).
pub fn extract(request: &Request<'_>, alignments: impl Read) -> Result<serde_json::Value, Error> {
    let contigs = request.contigs.map(read_contigs).transpose()?;
    let mut read1 = FastqReader::open(request.read1)?;
    let mut read2 = FastqReader::open(request.read2)?;
    let mut groups = SamGroups::new(BufReader::with_capacity(1 << 20, alignments));
    let path = |name: &str| request.output_dir.join(name);
    let (mut out1, mut out2, mut table) = (
        create(&path(DISCOVERY_R1))?,
        create(&path(DISCOVERY_R2))?,
        create(&path(DECISIONS))?,
    );
    let write_error = |source| Error::Io {
        context: format!("cannot write outputs in {}", request.output_dir.display()),
        source,
    };
    table
        .write_all(DECISION_COLUMNS.as_bytes())
        .map_err(write_error)?;

    let mut pairs = 0u64;
    let (mut bases_in, mut bases_retained) = (0u64, 0u64);
    let mut outcomes: BTreeMap<&str, u64> = BTreeMap::new();
    let mut reasons: BTreeMap<&str, u64> = BTreeMap::new();
    let mut assigned_by_reference: BTreeMap<String, u64> = BTreeMap::new();
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
        let pair_bases = sequence_length(&first) + sequence_length(&second);
        bases_in += pair_bases;
        if first.name != second.name {
            return Err(invalid(format!(
                "pair {pairs}: mate names differ ({:?} and {:?})",
                first.name, second.name
            )));
        }
        let decision: Decision = match &contigs {
            None => Decision {
                outcome: Outcome::DiscoveryUnassigned,
                reason: "no-selected-references",
                reference: None,
                species: None,
                mapq: [None, None],
            },
            Some(contigs) => {
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
                decide(&records, contigs, &request.policy)
                    .map_err(|problem| invalid(format!("pair {pairs} ({name:?}): {problem}")))?
            }
        };
        *outcomes.entry(decision.outcome.label()).or_default() += 1;
        *reasons.entry(decision.reason).or_default() += 1;
        if decision.outcome == Outcome::ReferenceAssigned {
            if let Some(reference) = &decision.reference {
                *assigned_by_reference.entry(reference.clone()).or_default() += 1;
            }
        } else {
            bases_retained += pair_bases;
            for line in &first.lines {
                out1.write_all(line).map_err(write_error)?;
            }
            for line in &second.lines {
                out2.write_all(line).map_err(write_error)?;
            }
        }
        let text = |value: &Option<String>| value.clone().unwrap_or_default();
        let mapq = |value: Option<u8>| value.map(|q| q.to_string()).unwrap_or_default();
        writeln!(
            table,
            "{pairs}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            first.name,
            decision.outcome.label(),
            decision.reason,
            text(&decision.reference),
            text(&decision.species),
            mapq(decision.mapq[0]),
            mapq(decision.mapq[1]),
        )
        .map_err(write_error)?;
    }
    if contigs.is_some() {
        if let Some((name, _)) = groups.next_group()? {
            return Err(invalid(format!(
                "alignments continue after the last FASTQ pair with {name:?}"
            )));
        }
    }
    let digest1 = finish(out1, &path(DISCOVERY_R1))?;
    let digest2 = finish(out2, &path(DISCOVERY_R2))?;
    let digest_table = finish(table, &path(DECISIONS))?;

    let assigned = outcomes.get("reference-assigned").copied().unwrap_or(0);
    let retained = pairs - assigned;
    let summary = json!({
        "schema_version": SCHEMA_VERSION,
        "producer": crate::execution::run_dir::PRODUCER,
        "sample_id": request.sample_id,
        "mode": if contigs.is_some() { "alignments" } else { "no-selected-references" },
        "policy": {
            "min_identity": request.policy.min_identity,
            "min_aligned_fraction": request.policy.min_aligned_fraction,
            "min_mapq": request.policy.min_mapq,
            "secondary_alignments": request.policy.max_secondary,
            "same_species_multimapping_is_assigned": true,
            "retain_both_mates_if_either_is_not_confident": true,
        },
        "counts": {
            "pairs_in": pairs,
            "reads_in": pairs * 2,
            "pairs_assigned": assigned,
            "pairs_retained": retained,
            "reads_retained": retained * 2,
            "bases_in": bases_in,
            "bases_retained": bases_retained,
            "pairs_by_outcome": outcomes,
        },
        "reasons": reasons,
        "assigned_by_reference": assigned_by_reference,
        "outputs": {"read1": DISCOVERY_R1, "read2": DISCOVERY_R2, "decisions": DECISIONS},
        "sha256": {"read1": digest1, "read2": digest2, "decisions": digest_table},
    });
    let summary_path = path(SUMMARY);
    let mut file = File::create_new(&summary_path).map_err(|source| Error::Io {
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
        sync::atomic::{AtomicU64, Ordering},
    };

    use flate2::read::GzDecoder;

    use super::*;

    static NEXT: AtomicU64 = AtomicU64::new(0);

    struct Dir(PathBuf, PathBuf);

    impl Dir {
        fn new() -> Self {
            let path = std::env::temp_dir().join(format!(
                "mca-extract-{}-{}",
                std::process::id(),
                NEXT.fetch_add(1, Ordering::Relaxed)
            ));
            let _ = fs::remove_dir_all(&path);
            fs::create_dir_all(path.join("out")).unwrap();
            let out = path.join("out");
            Self(path, out)
        }

        fn write(&self, name: &str, text: &str) -> PathBuf {
            let path = self.0.join(name);
            fs::write(&path, text).unwrap();
            path
        }

        fn write_gz(&self, name: &str, text: &str) -> PathBuf {
            let path = self.0.join(name);
            let mut encoder = GzEncoder::new(File::create(&path).unwrap(), Compression::default());
            encoder.write_all(text.as_bytes()).unwrap();
            encoder.finish().unwrap();
            path
        }

        fn read_gz(&self, name: &str) -> String {
            let mut text = String::new();
            GzDecoder::new(File::open(self.0.join("out").join(name)).unwrap())
                .read_to_string(&mut text)
                .unwrap();
            text
        }
    }

    impl Drop for Dir {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.0);
        }
    }

    const POLICY: Policy = Policy {
        min_identity: 0.98,
        min_aligned_fraction: 0.9,
        min_mapq: 20,
        max_secondary: 5,
    };

    fn fastq(names: &[&str], mate: u8) -> String {
        names
            .iter()
            .map(|name| format!("@{name}/{mate} comment:{mate}\nACGTACGTAC\n+\nIIIIIIIIII\n"))
            .collect()
    }

    fn sam_pair(name: &str, contig: &str, flags: [u16; 2]) -> String {
        flags
            .iter()
            .map(|flag| format!("{name}\t{flag}\t{contig}\t1\t60\t10=\t=\t1\t10\tACGTACGTAC\tIIIIIIIIII\tNM:i:0\n"))
            .collect()
    }

    fn unmapped(name: &str) -> String {
        [77, 141]
            .iter()
            .map(|flag| format!("{name}\t{flag}\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n"))
            .collect()
    }

    fn setup(dir: &Dir, names: &[&str]) -> (PathBuf, PathBuf, PathBuf) {
        let contigs = dir.write(
            "contigs.tsv",
            "contig\treference_id\tspecies_taxid\tlength\nrefA::chr\trefA\t562\t100\n",
        );
        (
            dir.write_gz("r1.fq.gz", &fastq(names, 1)),
            dir.write("r2.fq", &fastq(names, 2)),
            contigs,
        )
    }

    fn request<'a>(
        dir: &'a Dir,
        files: &'a (PathBuf, PathBuf, PathBuf),
        with_contigs: bool,
    ) -> Request<'a> {
        Request {
            sample_id: "S1",
            read1: &files.0,
            read2: &files.1,
            contigs: with_contigs.then_some(files.2.as_path()),
            policy: POLICY,
            output_dir: &dir.1,
        }
    }

    #[test]
    fn retained_pairs_are_copied_verbatim_and_counts_are_conserved() {
        let dir = Dir::new();
        let files = setup(&dir, &["p1", "p2", "p3"]);
        let sam = format!(
            "@HD\tVN:1.6\n{}{}{}",
            sam_pair("p1", "refA::chr", [99, 147]),
            unmapped("p2"),
            sam_pair("p3", "refA::chr", [65, 129])
        );
        let summary = extract(&request(&dir, &files, true), Cursor::new(sam)).unwrap();
        assert_eq!(summary["counts"]["pairs_in"], 3);
        assert_eq!(summary["counts"]["pairs_assigned"], 1);
        assert_eq!(summary["counts"]["pairs_retained"], 2);
        assert_eq!(summary["counts"]["reads_retained"], 4);
        assert_eq!(summary["reasons"]["mate-unmapped"], 1);
        assert_eq!(summary["reasons"]["discordant-pair"], 1);
        assert_eq!(summary["assigned_by_reference"]["refA"], 1);
        assert_eq!(summary["counts"]["bases_in"], 60);
        assert_eq!(summary["counts"]["bases_retained"], 40);
        let digest = crate::catalogue::sha256_file(&dir.1.join(DISCOVERY_R1)).unwrap();
        assert_eq!(summary["sha256"]["read1"], digest.as_str());
        assert_eq!(
            dir.read_gz(DISCOVERY_R1),
            fastq(&["p2", "p3"], 1),
            "headers with comments are kept"
        );
        assert_eq!(dir.read_gz(DISCOVERY_R2), fastq(&["p2", "p3"], 2));
        let table = dir.read_gz(DECISIONS);
        let rows: Vec<_> = table.lines().collect();
        assert_eq!(rows.len(), 4);
        assert_eq!(
            rows[1],
            "1\tp1\treference-assigned\tunique\trefA\t562\t60\t60"
        );
        assert!(rows[2].starts_with("2\tp2\tdiscovery-unassigned\tmate-unmapped"));
        assert!(dir.0.join("out").join(SUMMARY).is_file());
    }

    #[test]
    fn without_selected_references_every_pair_is_kept() {
        let dir = Dir::new();
        let files = setup(&dir, &["p1", "p2"]);
        let summary = extract(&request(&dir, &files, false), io::empty()).unwrap();
        assert_eq!(summary["mode"], "no-selected-references");
        assert_eq!(summary["counts"]["pairs_retained"], 2);
        assert_eq!(dir.read_gz(DISCOVERY_R1), fastq(&["p1", "p2"], 1));
    }

    #[test]
    fn empty_inputs_produce_empty_valid_pools() {
        let dir = Dir::new();
        let files = setup(&dir, &[]);
        let summary = extract(&request(&dir, &files, true), Cursor::new("@HD\tVN:1.6\n")).unwrap();
        assert_eq!(summary["counts"]["pairs_in"], 0);
        assert_eq!(dir.read_gz(DISCOVERY_R1), "");
    }

    #[test]
    fn order_count_and_name_problems_fail_instead_of_mixing_reads() {
        let swapped = format!(
            "{}{}",
            sam_pair("p2", "refA::chr", [99, 147]),
            sam_pair("p1", "refA::chr", [99, 147])
        );
        let short = sam_pair("p1", "refA::chr", [99, 147]);
        let long = format!(
            "{short}{}{}",
            sam_pair("p2", "refA::chr", [99, 147]),
            sam_pair("p3", "refA::chr", [99, 147])
        );
        for (name, sam, problem) in [
            ("swapped", swapped, "alignment order differs"),
            ("short", short, "alignments ended before pair 2"),
            (
                "long",
                long,
                "alignments continue after the last FASTQ pair",
            ),
        ] {
            let dir = Dir::new();
            let files = setup(&dir, &["p1", "p2"]);
            let error = extract(&request(&dir, &files, true), Cursor::new(sam)).unwrap_err();
            assert!(error.to_string().contains(problem), "{name}: {error}");
        }
        let dir = Dir::new();
        let mut files = setup(&dir, &["p1", "p2"]);
        files.1 = dir.write("r2-short.fq", &fastq(&["p1"], 2));
        let error = extract(&request(&dir, &files, false), io::empty()).unwrap_err();
        assert!(error.to_string().contains("different numbers of records"));
        files.1 = dir.write("r2-other.fq", &fastq(&["p1", "zz"], 2));
        let dir2 = Dir::new();
        let error = extract(&request(&dir2, &files, false), io::empty()).unwrap_err();
        assert!(error.to_string().contains("mate names differ"));
    }

    #[test]
    fn adjacent_pairs_sharing_a_name_are_rejected() {
        let dir = Dir::new();
        let files = setup(&dir, &["dup", "dup"]);
        let sam = format!(
            "{}{}",
            sam_pair("dup", "refA::chr", [99, 147]),
            sam_pair("dup", "refA::chr", [99, 147])
        );
        let error = extract(&request(&dir, &files, true), Cursor::new(sam)).unwrap_err();
        assert!(
            error.to_string().contains("more than one primary"),
            "{error}"
        );
    }

    #[test]
    fn malformed_fastq_and_existing_outputs_are_errors() {
        let dir = Dir::new();
        let mut files = setup(&dir, &["p1"]);
        files.0 = dir.write("bad.fq", "@p1/1\nACGT\n+\nII\n");
        assert!(
            extract(&request(&dir, &files, false), io::empty())
                .unwrap_err()
                .to_string()
                .contains("lengths differ")
        );
        let dir = Dir::new();
        let files = setup(&dir, &["p1"]);
        fs::write(dir.0.join("out").join(DISCOVERY_R1), "old").unwrap();
        assert!(extract(&request(&dir, &files, false), io::empty()).is_err());
        assert_eq!(
            fs::read_to_string(dir.0.join("out").join(DISCOVERY_R1)).unwrap(),
            "old"
        );
    }

    #[test]
    fn pair_names_drop_comments_and_mate_suffixes_only() {
        assert_eq!(pair_name(b"@A01|g|1/1 x:y\n"), "A01|g|1");
        assert_eq!(pair_name(b"@read/2\n"), "read");
        assert_eq!(pair_name(b"@read/3\n"), "read/3");
        assert_eq!(pair_name(b"@read 1:N:0\r\n"), "read");
    }
}
