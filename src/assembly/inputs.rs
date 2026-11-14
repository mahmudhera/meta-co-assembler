//! One streaming pass over a read pair that checks mates, counts pairs and
//! bases, and fingerprints both files.

use std::{
    fs::File,
    io::{BufRead, BufReader},
    path::{Path, PathBuf},
};

use flate2::read::MultiGzDecoder;
use serde::Serialize;

use crate::{
    digest::HashingReader,
    discovery::extract::{FastqReader, sequence_length},
    error::Error,
};

/// One sample's reads as an assembly input.
#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub struct ReadPair {
    pub sample_id: String,
    pub read1: PathBuf,
    pub read2: PathBuf,
    pub pairs: u64,
    pub bases: u64,
    pub read1_sha256: String,
    pub read2_sha256: String,
}

fn open<'a>(path: &Path, hashing: &'a mut HashingReader<File>) -> Result<FastqReader<'a>, Error> {
    let mut buffered = BufReader::with_capacity(1 << 20, hashing);
    let gzip = buffered
        .fill_buf()
        .map_err(|source| Error::Io {
            context: format!("cannot read {}", path.display()),
            source,
        })?
        .starts_with(&[0x1f, 0x8b]);
    let input: Box<dyn BufRead + 'a> = if gzip {
        Box::new(BufReader::with_capacity(
            1 << 20,
            MultiGzDecoder::new(buffered),
        ))
    } else {
        Box::new(buffered)
    };
    Ok(FastqReader::from_reader(path, input))
}

fn hashing(path: &Path) -> Result<HashingReader<File>, Error> {
    File::open(path)
        .map(HashingReader::new)
        .map_err(|source| Error::Io {
            context: format!("cannot open {}", path.display()),
            source,
        })
}

/// Reads both mates in lockstep. Mismatched names, different record counts,
/// or malformed FASTQ are input errors.
pub fn scan(sample_id: &str, read1: &Path, read2: &Path) -> Result<ReadPair, Error> {
    let (mut hash1, mut hash2) = (hashing(read1)?, hashing(read2)?);
    let (mut pairs, mut bases) = (0u64, 0u64);
    {
        let mut first = open(read1, &mut hash1)?;
        let mut second = open(read2, &mut hash2)?;
        loop {
            match (first.next_record()?, second.next_record()?) {
                (None, None) => break,
                (Some(one), Some(two)) => {
                    pairs += 1;
                    if one.name != two.name {
                        return Err(Error::InvalidInput(format!(
                            "{sample_id}: pair {pairs} has mate names {:?} and {:?}",
                            one.name, two.name
                        )));
                    }
                    bases += sequence_length(&one) + sequence_length(&two);
                }
                _ => {
                    return Err(Error::InvalidInput(format!(
                        "{sample_id}: {} and {} have different numbers of records",
                        read1.display(),
                        read2.display()
                    )));
                }
            }
        }
    }
    Ok(ReadPair {
        sample_id: sample_id.to_owned(),
        read1: read1.to_owned(),
        read2: read2.to_owned(),
        pairs,
        bases,
        read1_sha256: hash1.digest(),
        read2_sha256: hash2.digest(),
    })
}

/// Scans samples in parallel, at most `threads` at a time, in input order.
pub fn scan_all(
    samples: &[(String, PathBuf, PathBuf)],
    threads: usize,
) -> Result<Vec<ReadPair>, Error> {
    let next = std::sync::atomic::AtomicUsize::new(0);
    let results: Vec<std::sync::Mutex<Option<Result<ReadPair, Error>>>> = samples
        .iter()
        .map(|_| std::sync::Mutex::new(None))
        .collect();
    std::thread::scope(|scope| {
        for _ in 0..threads.clamp(1, samples.len().max(1)) {
            scope.spawn(|| {
                loop {
                    let index = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    let Some((id, read1, read2)) = samples.get(index) else {
                        break;
                    };
                    let result = scan(id, read1, read2);
                    *results[index]
                        .lock()
                        .unwrap_or_else(std::sync::PoisonError::into_inner) = Some(result);
                }
            });
        }
    });
    results
        .into_iter()
        .map(|slot| {
            slot.into_inner()
                .unwrap_or_else(std::sync::PoisonError::into_inner)
                .expect("every sample is scanned")
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use std::{fs, io::Write};

    use flate2::{Compression, write::GzEncoder};

    use super::*;
    use crate::catalogue::sha256_file;

    fn dir(name: &str) -> PathBuf {
        let path = std::env::temp_dir().join(format!("mca-scan-{}-{name}", std::process::id()));
        let _ = fs::remove_dir_all(&path);
        fs::create_dir_all(&path).unwrap();
        path
    }

    fn fastq(mate: u8, names: &[&str]) -> String {
        names
            .iter()
            .map(|n| format!("@{n}/{mate} x\nACGTA\n+\nIIIII\n"))
            .collect()
    }

    #[test]
    fn counts_pairs_and_bases_and_fingerprints_raw_file_bytes() {
        let root = dir("ok");
        let (r1, r2) = (root.join("a_R1.fq.gz"), root.join("a_R2.fq"));
        let mut encoder = GzEncoder::new(fs::File::create(&r1).unwrap(), Compression::default());
        encoder
            .write_all(fastq(1, &["p1", "p2", "p3"]).as_bytes())
            .unwrap();
        encoder.finish().unwrap();
        fs::write(&r2, fastq(2, &["p1", "p2", "p3"])).unwrap();
        let scan = scan("A", &r1, &r2).unwrap();
        assert_eq!((scan.pairs, scan.bases), (3, 30));
        assert_eq!(scan.read1_sha256, sha256_file(&r1).unwrap());
        assert_eq!(scan.read2_sha256, sha256_file(&r2).unwrap());
        let both = scan_all(
            &[("A".into(), r1.clone(), r2.clone()), ("B".into(), r1, r2)],
            4,
        )
        .unwrap();
        assert_eq!(both.len(), 2);
        assert_eq!(both[1].sample_id, "B");
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn mismatched_or_uneven_mates_are_input_errors() {
        let root = dir("bad");
        let (r1, r2) = (root.join("a_R1.fq"), root.join("a_R2.fq"));
        fs::write(&r1, fastq(1, &["p1", "p2"])).unwrap();
        fs::write(&r2, fastq(2, &["p1", "zz"])).unwrap();
        assert!(
            scan("A", &r1, &r2)
                .unwrap_err()
                .to_string()
                .contains("mate names")
        );
        fs::write(&r2, fastq(2, &["p1"])).unwrap();
        assert!(
            scan_all(&[("A".into(), r1, r2)], 2)
                .unwrap_err()
                .to_string()
                .contains("different numbers")
        );
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn empty_pools_are_valid_and_count_zero() {
        let root = dir("empty");
        let (r1, r2) = (root.join("a_R1.fq"), root.join("a_R2.fq"));
        fs::write(&r1, "").unwrap();
        fs::write(&r2, "").unwrap();
        let scan = scan("A", &r1, &r2).unwrap();
        assert_eq!((scan.pairs, scan.bases), (0, 0));
        fs::remove_dir_all(root).unwrap();
    }
}
