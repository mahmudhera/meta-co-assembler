//! Per-contig abundance as strobealign `--aemb` computes it, and the
//! per-sample coverage matrix the binners read.
//!
//! Each read adds its full length, split evenly among the distinct contigs of
//! its best-scoring alignments (primary plus secondaries with the same `AS`);
//! a contig's abundance is that total divided by its length. Supplementary
//! and unmapped records add nothing. Strobealign computes this for read
//! pairs; from minimap2 SAM it is computed per mate, which gives the same
//! result when both mates align to the same contigs.

use std::{
    collections::HashMap,
    fs,
    io::{BufRead, Write},
    path::Path,
};

const UNMAPPED: u16 = 0x4;
const FIRST: u16 = 0x40;
const SECOND: u16 = 0x80;
const SUPPLEMENTARY: u16 = 0x800;

/// The fields of one SAM record that abundance uses.
#[derive(Debug, PartialEq)]
struct Record {
    name: String,
    flag: u16,
    contig: Option<usize>,
    score: Option<i64>,
    query_length: u64,
}

/// Read length from a CIGAR, counting clipped bases; `None` for `*`.
fn cigar_query_length(cigar: &str) -> Result<Option<u64>, String> {
    if cigar == "*" {
        return Ok(None);
    }
    let (mut total, mut number) = (0u64, 0u64);
    for byte in cigar.bytes() {
        if byte.is_ascii_digit() {
            number = number * 10 + u64::from(byte - b'0');
            continue;
        }
        if matches!(byte, b'M' | b'I' | b'S' | b'=' | b'X' | b'H') {
            total += number;
        } else if !matches!(byte, b'D' | b'N' | b'P') {
            return Err(format!("CIGAR {cigar:?} has an unknown operation"));
        }
        number = 0;
    }
    Ok(Some(total))
}

fn parse_record(line: &str, contigs: &HashMap<String, usize>) -> Result<Record, String> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 11 {
        return Err(format!(
            "SAM record has {} fields; expected at least 11",
            fields.len()
        ));
    }
    let flag: u16 = fields[1]
        .parse()
        .map_err(|_| format!("invalid FLAG {:?}", fields[1]))?;
    let contig = if flag & UNMAPPED != 0 || fields[2] == "*" {
        None
    } else {
        Some(
            *contigs
                .get(fields[2])
                .ok_or_else(|| format!("contig {:?} has no @SQ header", fields[2]))?,
        )
    };
    let sequence = (fields[9] != "*").then(|| fields[9].len() as u64);
    let query_length = cigar_query_length(fields[5])?.or(sequence).unwrap_or(0);
    let score = fields[11..]
        .iter()
        .find_map(|tag| tag.strip_prefix("AS:i:"))
        .map(|value| {
            value
                .parse::<i64>()
                .map_err(|_| format!("invalid AS tag {value:?}"))
        })
        .transpose()?;
    Ok(Record {
        name: fields[0].to_owned(),
        flag,
        contig,
        score,
        query_length,
    })
}

/// Adds one fragment's records (every record with one read name).
fn add_fragment(records: &[Record], bases: &mut [f64]) {
    for mate in [FIRST, SECOND, 0] {
        let mine: Vec<&Record> = records
            .iter()
            .filter(|r| r.flag & (FIRST | SECOND) == mate)
            .collect();
        let length = mine.iter().map(|r| r.query_length).max().unwrap_or(0);
        let aligned: Vec<&&Record> = mine
            .iter()
            .filter(|r| r.contig.is_some() && r.flag & SUPPLEMENTARY == 0)
            .collect();
        let Some(best) = aligned.iter().map(|r| r.score.unwrap_or(i64::MIN)).max() else {
            continue;
        };
        let mut hits: Vec<usize> = aligned
            .iter()
            .filter(|r| r.score.unwrap_or(i64::MIN) == best)
            .filter_map(|r| r.contig)
            .collect();
        hits.sort_unstable();
        hits.dedup();
        #[allow(clippy::cast_precision_loss)]
        let share = length as f64 / hits.len() as f64;
        for contig in hits {
            bases[contig] += share;
        }
    }
}

/// Abundance of every `@SQ` contig, in header order, from SAM whose records
/// of one read name are consecutive (as mappers write them).
pub fn from_sam(input: impl BufRead) -> Result<Vec<(String, f64)>, String> {
    let mut names: Vec<String> = Vec::new();
    let mut lengths: Vec<u64> = Vec::new();
    let mut index: HashMap<String, usize> = HashMap::new();
    let mut bases: Vec<f64> = Vec::new();
    let mut fragment: Vec<Record> = Vec::new();
    for (number, line) in input.lines().enumerate() {
        let line = line.map_err(|error| format!("cannot read SAM: {error}"))?;
        if let Some(header) = line.strip_prefix("@SQ\t") {
            let field = |key: &str| {
                header
                    .split('\t')
                    .find_map(|f| f.strip_prefix(key))
                    .map(str::to_owned)
            };
            let (Some(name), Some(length)) = (field("SN:"), field("LN:")) else {
                return Err(format!("SAM line {}: @SQ lacks SN or LN", number + 1));
            };
            let length: u64 = length
                .parse()
                .map_err(|_| format!("SAM line {}: invalid LN {length:?}", number + 1))?;
            if length == 0 || index.insert(name.clone(), names.len()).is_some() {
                return Err(format!(
                    "SAM line {}: contig {name:?} is empty or repeated",
                    number + 1
                ));
            }
            names.push(name);
            lengths.push(length);
            bases.push(0.0);
            continue;
        }
        if line.starts_with('@') || line.is_empty() {
            continue;
        }
        let record = parse_record(&line, &index)
            .map_err(|error| format!("SAM line {}: {error}", number + 1))?;
        if fragment
            .first()
            .is_some_and(|first| first.name != record.name)
        {
            add_fragment(&fragment, &mut bases);
            fragment.clear();
        }
        fragment.push(record);
    }
    add_fragment(&fragment, &mut bases);
    if names.is_empty() {
        return Err("SAM has no @SQ headers".into());
    }
    #[allow(clippy::cast_precision_loss)]
    Ok(names
        .into_iter()
        .zip(bases.iter().zip(&lengths))
        .map(|(name, (&total, &length))| (name, total / length as f64))
        .collect())
}

/// Writes `contig<TAB>abundance` lines, as strobealign `--aemb` does.
pub fn write_tsv(out: &mut impl Write, rows: &[(String, f64)]) -> std::io::Result<()> {
    for (name, value) in rows {
        writeln!(out, "{name}\t{value:.6}")?;
    }
    Ok(())
}

/// Reads an abundance file into contig -> abundance.
pub fn read_tsv(path: &Path) -> Result<HashMap<String, f64>, String> {
    let text = fs::read_to_string(path)
        .map_err(|error| format!("cannot read {}: {error}", path.display()))?;
    let mut values = HashMap::new();
    for (number, line) in text.lines().enumerate() {
        let mut fields = line.split('\t');
        let (Some(name), Some(value), None) = (fields.next(), fields.next(), fields.next()) else {
            return Err(format!(
                "{} line {}: expected contig and abundance",
                path.display(),
                number + 1
            ));
        };
        let value: f64 = value
            .parse()
            .ok()
            .filter(|v: &f64| v.is_finite() && *v >= 0.0)
            .ok_or_else(|| format!("{} line {}: invalid abundance", path.display(), number + 1))?;
        if values.insert(name.to_owned(), value).is_some() {
            return Err(format!("{}: contig {name:?} is repeated", path.display()));
        }
    }
    if values.is_empty() {
        return Err(format!("{} is empty", path.display()));
    }
    Ok(values)
}

/// Whole-contig abundance per sample, in `contigs` order. With `split`, the
/// files hold SemiBin2's contig halves (`ID_1`, `ID_2`), averaged here.
pub fn matrix(
    contigs: &[String],
    files: &[(String, &Path)],
    split: bool,
) -> Result<Vec<Vec<f64>>, String> {
    let mut columns = Vec::with_capacity(files.len());
    for (sample, path) in files {
        let values = read_tsv(path)?;
        let mut column = Vec::with_capacity(contigs.len());
        for contig in contigs {
            let value = if split {
                let halves = [format!("{contig}_1"), format!("{contig}_2")];
                let found: Vec<f64> = halves
                    .iter()
                    .filter_map(|h| values.get(h))
                    .copied()
                    .collect();
                (!found.is_empty()).then(|| found.iter().sum::<f64>() / found.len() as f64)
            } else {
                values.get(contig).copied()
            };
            column.push(value.ok_or_else(|| {
                format!(
                    "{}: no abundance for contig {contig:?} of {sample}",
                    path.display()
                )
            })?);
        }
        columns.push(column);
    }
    Ok((0..contigs.len())
        .map(|row| columns.iter().map(|column| column[row]).collect())
        .collect())
}

/// Writes the matrix with a `contigName` header, the format MetaBAT 2 reads
/// with `--cvExt`.
pub fn write_matrix(
    out: &mut impl Write,
    samples: &[String],
    contigs: &[String],
    rows: &[Vec<f64>],
) -> std::io::Result<()> {
    write!(out, "contigName")?;
    for sample in samples {
        write!(out, "\t{sample}")?;
    }
    writeln!(out)?;
    for (contig, row) in contigs.iter().zip(rows) {
        write!(out, "{contig}")?;
        for value in row {
            write!(out, "\t{value:.6}")?;
        }
        writeln!(out)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    const HEADER: &str = "@HD\tVN:1.6\n@SQ\tSN:a\tLN:100\n@SQ\tSN:b\tLN:100\n@SQ\tSN:c\tLN:200\n";

    fn sam(records: &[&str]) -> String {
        let mut text = HEADER.to_owned();
        for record in records {
            text.push_str(record);
            text.push('\n');
        }
        text
    }

    fn abundance(text: &str) -> Vec<f64> {
        from_sam(text.as_bytes())
            .unwrap()
            .into_iter()
            .map(|(_, v)| v)
            .collect()
    }

    #[test]
    fn ties_split_a_read_and_lower_scores_do_not() {
        let text = sam(&[
            // r1 mate 1: tied between a and b; a weaker hit on c is ignored.
            "r1\t65\ta\t1\t60\t50M\t=\t1\t0\t*\t*\tAS:i:100",
            "r1\t321\tb\t1\t0\t50M\t=\t1\t0\t*\t*\tAS:i:100",
            "r1\t321\tc\t1\t0\t50M\t=\t1\t0\t*\t*\tAS:i:80",
            // r1 mate 2: soft clips count toward the read length.
            "r1\t129\tc\t1\t60\t10S40M\t=\t1\t0\t*\t*\tAS:i:80",
            // r2: unmapped mate adds nothing; supplementary records add nothing.
            "r2\t77\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII",
            "r2\t2177\ta\t1\t0\t20H30M\t*\t0\t0\t*\t*\tAS:i:60",
        ]);
        let values = abundance(&text);
        assert_eq!(
            values,
            [0.25, 0.25, 0.25],
            "a, b: 25 bases each; c: 50 of 200"
        );
    }

    #[test]
    fn records_without_scores_count_once_and_malformed_sam_is_an_error() {
        let text = sam(&["r1\t0\tc\t1\t60\t100M\t*\t0\t0\t*\t*"]);
        assert_eq!(abundance(&text), [0.0, 0.0, 0.5]);
        for bad in [
            sam(&["r1\t0\tz\t1\t60\t10M\t*\t0\t0\t*\t*"]),
            sam(&["r1\t0\ta\t1\t60\t10Q\t*\t0\t0\t*\t*"]),
            sam(&["r1\tx\ta\t1\t60\t10M\t*\t0\t0\t*\t*"]),
            "r1\t4\t*\t0\t0\t*\t*\t0\t0\tA\tI\n".to_owned(),
            "@SQ\tSN:a\tLN:0\n".to_owned(),
        ] {
            assert!(from_sam(bad.as_bytes()).is_err(), "{bad:?}");
        }
    }

    #[test]
    fn matrices_follow_contig_order_and_average_split_halves() {
        let dir = std::env::temp_dir().join(format!("mca-abundance-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let (whole, split) = (dir.join("whole.tsv"), dir.join("split.tsv"));
        fs::write(&whole, "b\t2.000000\na\t1.000000\n").unwrap();
        fs::write(&split, "a_1\t1.0\na_2\t3.0\nb_1\t4.0\nb_2\t4.0\n").unwrap();
        let contigs = vec!["a".to_owned(), "b".to_owned()];
        let rows = matrix(
            &contigs,
            &[
                ("S1".into(), whole.as_path()),
                ("S2".into(), split.as_path()),
            ],
            false,
        );
        assert!(rows.is_err(), "split names are not whole contigs");
        let rows = matrix(&contigs, &[("S1".into(), whole.as_path())], false).unwrap();
        assert_eq!(rows, [vec![1.0], vec![2.0]]);
        let rows = matrix(&contigs, &[("S2".into(), split.as_path())], true).unwrap();
        assert_eq!(rows, [vec![2.0], vec![4.0]]);
        let mut out = Vec::new();
        write_matrix(&mut out, &["S2".into()], &contigs, &rows).unwrap();
        assert_eq!(
            String::from_utf8(out).unwrap(),
            "contigName\tS2\na\t2.000000\nb\t4.000000\n"
        );
        fs::write(&whole, "a\t1\na\t2\n").unwrap();
        assert!(read_tsv(&whole).is_err(), "repeated contigs are errors");
        fs::remove_dir_all(&dir).unwrap();
    }
}
