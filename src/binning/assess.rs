//! Bin tables: normalized bins, CheckM2 quality, ANI dereplication, and the
//! catalogue comparison. Pure parsing and decisions; no tools run here.

use std::{
    cmp::Ordering,
    collections::HashMap,
    fs::File,
    io::{BufWriter, Write},
    path::{Path, PathBuf},
};

use serde::{Deserialize, Serialize};

use crate::assembly::contigs::{FastaReader, stats};

/// Quality tiers by CheckM2 completeness and contamination only, after the
/// MIMAG thresholds (Bowers et al. 2017). rRNA and tRNA genes are not
/// assessed, so `high` is not a MIMAG high-quality draft.
pub fn tier(completeness: f64, contamination: f64) -> &'static str {
    if completeness >= 90.0 && contamination < 5.0 {
        "high"
    } else if completeness >= 50.0 && contamination < 10.0 {
        "medium"
    } else {
        "low"
    }
}

/// dRep-style score used to choose a cluster's representative.
pub fn score(completeness: f64, contamination: f64) -> f64 {
    completeness - 5.0 * contamination
}

/// One normalized bin of one binning unit.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Bin {
    pub bin_id: String,
    pub unit: String,
    /// The binner's own file name.
    pub source: String,
    pub contigs: u64,
    pub bases: u64,
    pub n50: u64,
    pub gc_percent: Option<f64>,
}

fn number_in(name: &str) -> Option<u64> {
    let digits: String = name
        .chars()
        .rev()
        .skip_while(|c| !c.is_ascii_digit())
        .take_while(char::is_ascii_digit)
        .collect();
    digits.chars().rev().collect::<String>().parse().ok()
}

/// Each contig's bin: (contig ID, bin ID).
pub type Membership = Vec<(String, String)>;

/// Copies the binner's bins (`raw`, any order) to `out_dir/BIN_ID.fa` with
/// IDs `UNIT__binNNN` in the binner's numeric order and headers reduced to
/// contig IDs. Returns the bins and each contig's bin.
pub fn normalize(
    raw: &[PathBuf],
    unit: &str,
    out_dir: &Path,
) -> Result<(Vec<Bin>, Membership), String> {
    let mut ordered: Vec<&PathBuf> = raw.iter().collect();
    let name = |p: &PathBuf| {
        p.file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_default()
    };
    ordered.sort_by_key(|p| (number_in(&name(p)), name(p)));
    let (mut bins, mut members) = (Vec::new(), Vec::new());
    for (index, path) in ordered.into_iter().enumerate() {
        let bin_id = format!("{unit}__bin{:03}", index + 1);
        let target = out_dir.join(format!("{bin_id}.fa"));
        let write_error =
            |error: std::io::Error| format!("cannot write {}: {error}", target.display());
        let mut out = BufWriter::new(File::create(&target).map_err(write_error)?);
        let mut reader = FastaReader::open(path)?;
        let (mut lengths, mut gc, mut acgt) = (Vec::new(), 0u64, 0u64);
        while let Some(record) = reader.next_record()? {
            writeln!(out, ">{}", record.id).map_err(write_error)?;
            out.write_all(&record.sequence).map_err(write_error)?;
            out.write_all(b"\n").map_err(write_error)?;
            lengths.push(record.length());
            for byte in &record.sequence {
                match byte.to_ascii_uppercase() {
                    b'G' | b'C' => {
                        gc += 1;
                        acgt += 1;
                    }
                    b'A' | b'T' => acgt += 1,
                    _ => {}
                }
            }
            members.push((record.id, bin_id.clone()));
        }
        out.flush().map_err(write_error)?;
        if lengths.is_empty() {
            return Err(format!("{} holds no contigs", path.display()));
        }
        let stats = stats(lengths, gc, acgt);
        bins.push(Bin {
            bin_id,
            unit: unit.to_owned(),
            source: name(path),
            contigs: stats.contigs,
            bases: stats.total_bases,
            n50: stats.n50,
            gc_percent: stats.gc_percent,
        });
    }
    Ok((bins, members))
}

/// CheckM2's estimate for one bin.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Quality {
    pub completeness: f64,
    pub contamination: f64,
    pub model: String,
}

fn columns(header: &str, wanted: &[&str], what: &str) -> Result<Vec<usize>, String> {
    let names: Vec<&str> = header.split('\t').collect();
    wanted
        .iter()
        .map(|w| {
            names
                .iter()
                .position(|n| n == w)
                .ok_or_else(|| format!("{what} lacks the {w} column"))
        })
        .collect()
}

/// Parses CheckM2's `quality_report.tsv` into bin name -> estimate.
pub fn parse_checkm2(text: &str) -> Result<HashMap<String, Quality>, String> {
    let mut lines = text.lines();
    let header = lines.next().ok_or("the CheckM2 report is empty")?;
    let index = columns(
        header,
        &[
            "Name",
            "Completeness",
            "Contamination",
            "Completeness_Model_Used",
        ],
        "the CheckM2 report",
    )?;
    let mut result = HashMap::new();
    for line in lines.filter(|l| !l.trim().is_empty()) {
        let fields: Vec<&str> = line.split('\t').collect();
        let get = |i: usize| {
            fields
                .get(index[i])
                .copied()
                .ok_or("a CheckM2 row is short")
        };
        let number = |i: usize| -> Result<f64, String> {
            get(i)?
                .parse::<f64>()
                .ok()
                .filter(|v| v.is_finite() && *v >= 0.0)
                .ok_or_else(|| format!("invalid CheckM2 value in {line:?}"))
        };
        result.insert(
            get(0)?.to_owned(),
            Quality {
                completeness: number(1)?,
                contamination: number(2)?,
                model: get(3)?.to_owned(),
            },
        );
    }
    Ok(result)
}

/// One skani comparison; aligned fractions are percentages.
#[derive(Clone, Debug, PartialEq)]
pub struct Pair {
    pub reference: PathBuf,
    pub query: PathBuf,
    pub ani: f64,
    pub af_reference: f64,
    pub af_query: f64,
}

/// Parses skani `dist` or sparse `triangle` output. Only the first five
/// columns are read: record names may contain tabs.
pub fn parse_skani(text: &str) -> Result<Vec<Pair>, String> {
    let mut lines = text.lines();
    let Some(header) = lines.next() else {
        return Ok(Vec::new());
    };
    let expected = [
        "Ref_file",
        "Query_file",
        "ANI",
        "Align_fraction_ref",
        "Align_fraction_query",
    ];
    if header.split('\t').take(5).ne(expected) {
        return Err(format!("unexpected skani header {header:?}"));
    }
    lines
        .filter(|l| !l.trim().is_empty())
        .map(|line| {
            let fields: Vec<&str> = line.splitn(6, '\t').collect();
            let number = |i: usize| -> Result<f64, String> {
                fields
                    .get(i)
                    .and_then(|v| v.parse::<f64>().ok())
                    .filter(|v| v.is_finite())
                    .ok_or_else(|| format!("invalid skani row {line:?}"))
            };
            Ok(Pair {
                reference: PathBuf::from(fields[0]),
                query: PathBuf::from(fields.get(1).ok_or("short skani row")?),
                ani: number(2)?,
                af_reference: number(3)?,
                af_query: number(4)?,
            })
        })
        .collect()
}

/// A bin's dereplication outcome.
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct Cluster {
    pub cluster: usize,
    pub representative: String,
    pub is_representative: bool,
    /// ANI to the representative (absent for the representative itself).
    pub ani: Option<f64>,
    pub aligned_fraction: Option<f64>,
}

/// Greedy dereplication: bins in order of score (then size, then ID); each
/// unassigned bin founds a cluster and takes every unassigned bin linked to
/// it at `min_ani` or more over at least `min_af` of the smaller bin.
/// `pairs` holds (bin index, bin index, ANI, AF of the first, AF of the second).
pub fn dereplicate(
    ids: &[String],
    bases: &[u64],
    scores: &[f64],
    pairs: &[(usize, usize, f64, f64, f64)],
    min_ani: f64,
    min_af: f64,
) -> Vec<Cluster> {
    let mut links: HashMap<(usize, usize), (f64, f64)> = HashMap::new();
    for &(a, b, ani, af_a, af_b) in pairs {
        if a == b {
            continue;
        }
        let smaller = if bases[a] <= bases[b] { af_a } else { af_b };
        if ani >= min_ani && smaller >= min_af * 100.0 {
            let key = (a.min(b), a.max(b));
            let entry = links.entry(key).or_insert((ani, smaller));
            if ani > entry.0 {
                *entry = (ani, smaller);
            }
        }
    }
    let mut order: Vec<usize> = (0..ids.len()).collect();
    order.sort_by(|&a, &b| {
        scores[b]
            .partial_cmp(&scores[a])
            .unwrap_or(Ordering::Equal)
            .then(bases[b].cmp(&bases[a]))
            .then(ids[a].cmp(&ids[b]))
    });
    let mut result: Vec<Option<Cluster>> = vec![None; ids.len()];
    let mut clusters = 0;
    for &rep in &order {
        if result[rep].is_some() {
            continue;
        }
        clusters += 1;
        result[rep] = Some(Cluster {
            cluster: clusters,
            representative: ids[rep].clone(),
            is_representative: true,
            ani: None,
            aligned_fraction: None,
        });
        for &other in &order {
            if result[other].is_some() {
                continue;
            }
            if let Some(&(ani, af)) = links.get(&(rep.min(other), rep.max(other))) {
                result[other] = Some(Cluster {
                    cluster: clusters,
                    representative: ids[rep].clone(),
                    is_representative: false,
                    ani: Some(ani),
                    aligned_fraction: Some(af / 100.0),
                });
            }
        }
    }
    result
        .into_iter()
        .map(|c| c.expect("every bin is assigned"))
        .collect()
}

/// A bin's closest catalogue reference.
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct CatalogueMatch {
    /// `catalogue-species-match`, `below-species-threshold`, or
    /// `no-catalogue-hit` (skani reported no comparison, roughly below 80% ANI).
    pub status: &'static str,
    pub reference_id: Option<String>,
    pub species_taxid: Option<String>,
    pub ani: Option<f64>,
    /// Fraction of the bin aligned to the reference.
    pub aligned_fraction_bin: Option<f64>,
    pub aligned_fraction_reference: Option<f64>,
}

/// The best catalogue hit of each bin: qualifying hits first, then highest
/// ANI, then the larger aligned fraction of the bin. `pairs` has catalogue
/// FASTAs as references and bins as queries.
pub fn catalogue_matches(
    bins: &[(String, PathBuf)],
    references: &HashMap<PathBuf, (String, String)>,
    pairs: &[Pair],
    min_ani: f64,
    min_af: f64,
) -> Result<Vec<CatalogueMatch>, String> {
    let mut best: HashMap<&Path, &Pair> = HashMap::new();
    let qualifies = |p: &Pair| p.ani >= min_ani && p.af_query >= min_af * 100.0;
    for pair in pairs {
        if !references.contains_key(&pair.reference) {
            return Err(format!(
                "skani reported unknown reference {}",
                pair.reference.display()
            ));
        }
        let key = (qualifies(pair), pair.ani, pair.af_query);
        let better = best.get(pair.query.as_path()).is_none_or(|current| {
            let other = (qualifies(current), current.ani, current.af_query);
            key.partial_cmp(&other) == Some(Ordering::Greater)
        });
        if better {
            best.insert(pair.query.as_path(), pair);
        }
    }
    Ok(bins
        .iter()
        .map(|(_, path)| match best.get(path.as_path()) {
            None => CatalogueMatch {
                status: "no-catalogue-hit",
                reference_id: None,
                species_taxid: None,
                ani: None,
                aligned_fraction_bin: None,
                aligned_fraction_reference: None,
            },
            Some(pair) => {
                let (reference_id, species) = &references[&pair.reference];
                CatalogueMatch {
                    status: if qualifies(pair) {
                        "catalogue-species-match"
                    } else {
                        "below-species-threshold"
                    },
                    reference_id: Some(reference_id.clone()),
                    species_taxid: Some(species.clone()),
                    ani: Some(pair.ani),
                    aligned_fraction_bin: Some(pair.af_query / 100.0),
                    aligned_fraction_reference: Some(pair.af_reference / 100.0),
                }
            }
        })
        .collect())
}

#[cfg(test)]
mod tests {
    use std::fs;

    use super::*;

    #[test]
    fn tiers_follow_the_mimag_completeness_and_contamination_thresholds() {
        assert_eq!(tier(90.0, 4.99), "high");
        assert_eq!(tier(89.99, 1.0), "medium");
        assert_eq!(tier(95.0, 5.0), "medium");
        assert_eq!(tier(50.0, 9.99), "medium");
        assert_eq!(tier(99.0, 10.0), "low");
        assert_eq!(tier(49.9, 0.0), "low");
        assert!((score(90.0, 2.0) - 80.0).abs() < 1e-9);
    }

    #[test]
    fn bins_are_renamed_in_numeric_order_with_bare_contig_headers() {
        let dir = std::env::temp_dir().join(format!("mca-bins-{}", std::process::id()));
        let (raw, out) = (dir.join("raw"), dir.join("out"));
        fs::create_dir_all(&raw).unwrap();
        fs::create_dir_all(&out).unwrap();
        fs::write(raw.join("bin.10.fa"), ">c3 depth=1\tx\nGGGG\n").unwrap();
        fs::write(raw.join("bin.2.fa"), ">c1 flag=1\nACGT\nAC\n>c2\nTT\n").unwrap();
        let (bins, members) = normalize(
            &[raw.join("bin.10.fa"), raw.join("bin.2.fa")],
            "group-000001",
            &out,
        )
        .unwrap();
        assert_eq!(
            bins.iter()
                .map(|b| (b.bin_id.as_str(), b.source.as_str(), b.contigs, b.bases))
                .collect::<Vec<_>>(),
            [
                ("group-000001__bin001", "bin.2.fa", 2, 8),
                ("group-000001__bin002", "bin.10.fa", 1, 4)
            ]
        );
        assert_eq!(
            fs::read_to_string(out.join("group-000001__bin001.fa")).unwrap(),
            ">c1\nACGTAC\n>c2\nTT\n"
        );
        assert_eq!(
            members[2],
            ("c3".to_owned(), "group-000001__bin002".to_owned())
        );
        fs::write(raw.join("empty.fa"), "").unwrap();
        assert!(normalize(&[raw.join("empty.fa")], "u", &out).is_err());
        fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn checkm2_and_skani_reports_parse_by_column() {
        let report = "Name\tCompleteness\tContamination\tCompleteness_Model_Used\tExtra\n\
                      b1\t99.5\t0.2\tNeural Network (Specific Model)\tx\n";
        let quality = parse_checkm2(report).unwrap();
        assert_eq!(quality["b1"].completeness, 99.5);
        assert!(parse_checkm2("Name\tCompleteness\n").is_err());
        assert!(
            parse_checkm2(
                "Name\tCompleteness\tContamination\tCompleteness_Model_Used\nb\tx\t1\tm\n"
            )
            .is_err()
        );
        let skani = "Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\tRef_name\tQuery_name\n\
                     r.fa\tq.fa\t99.10\t80.00\t60.50\tc1\tdepth=1\tq1\n";
        let pairs = parse_skani(skani).unwrap();
        assert_eq!(pairs[0].query, PathBuf::from("q.fa"));
        assert_eq!((pairs[0].ani, pairs[0].af_query), (99.1, 60.5));
        assert!(parse_skani("").unwrap().is_empty());
        assert!(parse_skani("Ref\tQuery\n").is_err());
    }

    #[test]
    fn dereplication_keeps_the_best_bin_and_needs_ani_and_aligned_fraction() {
        let ids: Vec<String> = ["a", "b", "c", "d"]
            .iter()
            .map(|s| (*s).to_owned())
            .collect();
        let bases = [1000, 2000, 1500, 500];
        let scores = [50.0, 90.0, 90.0, 10.0];
        let pairs = [
            (0, 1, 99.5, 80.0, 40.0), // a (smaller) aligns 80% to b
            (1, 2, 99.9, 30.0, 40.0), // c (smaller) aligns only 40%
            (2, 3, 98.0, 90.0, 95.0), // below 99% ANI
        ];
        let clusters = dereplicate(&ids, &bases, &scores, &pairs, 99.0, 0.5);
        let view: Vec<_> = clusters
            .iter()
            .map(|c| (c.cluster, c.representative.as_str(), c.is_representative))
            .collect();
        assert_eq!(
            view,
            [
                (1, "b", false),
                (1, "b", true),
                (2, "c", true),
                (3, "d", true)
            ],
            "b outscores c by size at equal score"
        );
        assert_eq!(clusters[0].aligned_fraction, Some(0.8));
    }

    #[test]
    fn catalogue_matches_need_species_level_ani_over_half_the_bin() {
        let refs: HashMap<PathBuf, (String, String)> = [
            (PathBuf::from("r1.fa"), ("R1".to_owned(), "562".to_owned())),
            (PathBuf::from("r2.fa"), ("R2".to_owned(), "573".to_owned())),
        ]
        .into();
        let pair = |r: &str, q: &str, ani, af_q| Pair {
            reference: PathBuf::from(r),
            query: PathBuf::from(q),
            ani,
            af_reference: 50.0,
            af_query: af_q,
        };
        let bins = vec![
            ("b1".to_owned(), PathBuf::from("b1.fa")),
            ("b2".to_owned(), PathBuf::from("b2.fa")),
            ("b3".to_owned(), PathBuf::from("b3.fa")),
        ];
        let pairs = [
            pair("r1.fa", "b1.fa", 99.0, 30.0),
            pair("r2.fa", "b1.fa", 96.0, 70.0),
            pair("r1.fa", "b2.fa", 94.0, 90.0),
        ];
        let matches = catalogue_matches(&bins, &refs, &pairs, 95.0, 0.5).unwrap();
        assert_eq!(
            matches
                .iter()
                .map(|m| (m.status, m.reference_id.as_deref()))
                .collect::<Vec<_>>(),
            [
                ("catalogue-species-match", Some("R2")),
                ("below-species-threshold", Some("R1")),
                ("no-catalogue-hit", None)
            ],
            "a qualifying hit beats a higher ANI over too little of the bin"
        );
        assert!(
            catalogue_matches(
                &bins,
                &refs,
                &[pair("x.fa", "b1.fa", 99.0, 99.0)],
                95.0,
                0.5
            )
            .is_err()
        );
    }
}
