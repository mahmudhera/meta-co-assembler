//! Contig FASTA normalization and basic statistics.

use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    path::Path,
};

use serde::Serialize;

/// Contiguity and composition of one assembly. Lengths of all contigs are
/// held in memory (proportional to contigs, not reads).
#[derive(Clone, Debug, Default, PartialEq, Serialize)]
pub struct ContigStats {
    pub contigs: u64,
    pub total_bases: u64,
    pub max_length: u64,
    pub n50: u64,
    pub l50: u64,
    pub contigs_at_least_1000: u64,
    pub bases_in_contigs_at_least_1000: u64,
    /// GC share of A/C/G/T bases; absent for an empty assembly.
    pub gc_percent: Option<f64>,
}

fn stats(mut lengths: Vec<u64>, gc: u64, acgt: u64) -> ContigStats {
    lengths.sort_unstable_by(|a, b| b.cmp(a));
    let total: u64 = lengths.iter().sum();
    let (mut running, mut n50, mut l50) = (0u64, 0u64, 0u64);
    for (index, &length) in lengths.iter().enumerate() {
        running += length;
        if 2 * running >= total {
            n50 = length;
            l50 = index as u64 + 1;
            break;
        }
    }
    let long: Vec<u64> = lengths.iter().copied().filter(|&l| l >= 1000).collect();
    #[allow(clippy::cast_precision_loss)]
    let gc_percent = (acgt > 0).then(|| (gc as f64 / acgt as f64 * 1e4).round() / 100.0);
    ContigStats {
        contigs: lengths.len() as u64,
        total_bases: total,
        max_length: lengths.first().copied().unwrap_or(0),
        n50,
        l50,
        contigs_at_least_1000: long.len() as u64,
        bases_in_contigs_at_least_1000: long.iter().sum(),
        gc_percent,
    }
}

/// Copies an assembler's contigs to `output`, renaming each to
/// `{prefix}__{original id}` (the original header text is kept after it) so
/// IDs stay unique across groups. Fails on malformed FASTA instead of
/// silently dropping records.
pub fn rewrite(input: &Path, output: &Path, prefix: &str) -> Result<ContigStats, String> {
    let read_error = |error: std::io::Error| format!("cannot read {}: {error}", input.display());
    let write_error = |error: std::io::Error| format!("cannot write {}: {error}", output.display());
    let reader = BufReader::with_capacity(1 << 20, File::open(input).map_err(read_error)?);
    let mut writer =
        BufWriter::with_capacity(1 << 20, File::create_new(output).map_err(write_error)?);
    let (mut lengths, mut gc, mut acgt) = (Vec::new(), 0u64, 0u64);
    let mut current: Option<u64> = None;
    let finish = |current: &mut Option<u64>, lengths: &mut Vec<u64>| -> Result<(), String> {
        match current.take() {
            Some(0) => Err(format!("{}: a contig has no sequence", input.display())),
            Some(length) => {
                lengths.push(length);
                Ok(())
            }
            None => Ok(()),
        }
    };
    for (number, line) in reader.lines().enumerate() {
        let line = line.map_err(read_error)?;
        let line = line.trim_end();
        if let Some(header) = line.strip_prefix('>') {
            finish(&mut current, &mut lengths)?;
            let mut parts = header.splitn(2, char::is_whitespace);
            let id = parts.next().unwrap_or_default();
            if id.is_empty() {
                return Err(format!(
                    "{}: line {} has an empty contig name",
                    input.display(),
                    number + 1
                ));
            }
            let rest = parts.next().map(str::trim).filter(|rest| !rest.is_empty());
            match rest {
                Some(rest) => writeln!(writer, ">{prefix}__{id} {rest}"),
                None => writeln!(writer, ">{prefix}__{id}"),
            }
            .map_err(write_error)?;
            current = Some(0);
        } else if line.is_empty() {
            continue;
        } else {
            let Some(length) = current.as_mut() else {
                return Err(format!(
                    "{}: sequence before the first header",
                    input.display()
                ));
            };
            for byte in line.bytes() {
                match byte.to_ascii_uppercase() {
                    b'G' | b'C' => {
                        gc += 1;
                        acgt += 1;
                    }
                    b'A' | b'T' => acgt += 1,
                    b'N' => {}
                    _ => {
                        return Err(format!(
                            "{}: line {} contains {:?}, which is not a nucleotide",
                            input.display(),
                            number + 1,
                            char::from(byte)
                        ));
                    }
                }
            }
            *length += line.len() as u64;
            writeln!(writer, "{line}").map_err(write_error)?;
        }
    }
    finish(&mut current, &mut lengths)?;
    writer.flush().map_err(write_error)?;
    Ok(stats(lengths, gc, acgt))
}

#[cfg(test)]
mod tests {
    use std::{fs, path::PathBuf};

    use super::*;

    fn dir(name: &str) -> PathBuf {
        let path = std::env::temp_dir().join(format!("mca-contigs-{}-{name}", std::process::id()));
        let _ = fs::remove_dir_all(&path);
        fs::create_dir_all(&path).unwrap();
        path
    }

    #[test]
    fn contigs_are_prefixed_and_measured() {
        let root = dir("ok");
        let input = root.join("final.contigs.fa");
        let long = "ACGT".repeat(300);
        fs::write(
            &input,
            format!(">k141_1 flag=1 multi=2.0 len=1200\n{long}\n>k141_2\nGGCC\nNA\n"),
        )
        .unwrap();
        let stats = rewrite(&input, &root.join("contigs.fa"), "group-000001").unwrap();
        assert_eq!(
            fs::read_to_string(root.join("contigs.fa")).unwrap(),
            format!(
                ">group-000001__k141_1 flag=1 multi=2.0 len=1200\n{long}\n>group-000001__k141_2\nGGCC\nNA\n"
            )
        );
        assert_eq!(
            (stats.contigs, stats.total_bases, stats.max_length),
            (2, 1206, 1200)
        );
        assert_eq!((stats.n50, stats.l50), (1200, 1));
        assert_eq!(
            (
                stats.contigs_at_least_1000,
                stats.bases_in_contigs_at_least_1000
            ),
            (1, 1200)
        );
        // 604 G/C of 1205 A/C/G/T bases.
        assert_eq!(stats.gc_percent, Some(50.12));
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn an_empty_assembly_is_valid_with_empty_statistics() {
        let root = dir("empty");
        fs::write(root.join("final.contigs.fa"), "").unwrap();
        let stats = rewrite(
            &root.join("final.contigs.fa"),
            &root.join("contigs.fa"),
            "g",
        )
        .unwrap();
        assert_eq!(stats, ContigStats::default());
        assert_eq!(fs::read(root.join("contigs.fa")).unwrap(), b"");
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn n50_uses_half_of_total_length() {
        let s = stats(vec![2, 3, 4, 5, 6], 0, 0);
        // Total 20: 6 + 5 = 11 >= 10, so N50 is 5 and L50 is 2.
        assert_eq!((s.n50, s.l50), (5, 2));
    }

    #[test]
    fn malformed_fasta_fails() {
        let root = dir("bad");
        for (name, text) in [
            ("headerless", "ACGT\n>c\nAC\n"),
            ("empty contig", ">a\n>b\nAC\n"),
            ("nameless", ">\nAC\n"),
            ("not dna", ">a\nAC-T\n"),
        ] {
            let input = root.join(format!("{name}.fa"));
            fs::write(&input, text).unwrap();
            assert!(
                rewrite(&input, &root.join(format!("{name}.out.fa")), "g").is_err(),
                "{name}"
            );
        }
        fs::remove_dir_all(root).unwrap();
    }
}
