//! Pair-level assignment decisions.
//!
//! A pair leaves the discovery pool only when both mates are confidently
//! explained by selected references of one species:
//!
//! - each mate's primary alignment is mapped, has NM and a CIGAR without
//!   hard clips or skips, identity `1 - NM/(M+=+X+I+D)` and aligned fraction
//!   `(M+=+X+I)/(M+=+X+I+S)` at least the thresholds, and a known MAPQ;
//! - the mates are properly paired and on references of the same species,
//!   with no supplementary (chimeric) records;
//! - every qualifying secondary alignment is in that species, so multi-mapping
//!   between strains of one species or repeats within one genome still counts
//!   as explained (owner decision, M4);
//! - a mate with MAPQ below the minimum whose secondaries reached the
//!   reporting cap may have unreported alternatives, so the pair is retained.
//!
//! Otherwise both mates are retained, as unassigned or ambiguous.

use std::collections::HashMap;

/// Contig name in the alignment reference -> (reference ID, species taxid).
pub type Contigs = HashMap<String, (String, String)>;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Policy {
    pub min_identity: f64,
    pub min_aligned_fraction: f64,
    pub min_mapq: u8,
    /// Secondary alignments requested from the mapper per read.
    pub max_secondary: usize,
}

pub const UNMAPPED: u16 = 0x4;
pub const PROPER_PAIR: u16 = 0x2;
pub const FIRST: u16 = 0x40;
pub const SECOND: u16 = 0x80;
pub const SECONDARY: u16 = 0x100;
pub const SUPPLEMENTARY: u16 = 0x800;

/// The fields of one SAM record that decisions use.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Alignment {
    pub flag: u16,
    pub contig: String,
    pub mapq: u8,
    pub cigar: String,
    pub nm: Option<u32>,
}

impl Alignment {
    /// Parses a SAM line's mandatory fields and the NM tag.
    pub fn parse(line: &str) -> Result<(String, Self), String> {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 11 {
            return Err(format!(
                "SAM record has {} fields; expected at least 11",
                fields.len()
            ));
        }
        let flag = fields[1]
            .parse()
            .map_err(|_| format!("invalid SAM FLAG {:?}", fields[1]))?;
        let mapq = fields[4]
            .parse()
            .map_err(|_| format!("invalid SAM MAPQ {:?}", fields[4]))?;
        let nm = fields[11..]
            .iter()
            .find_map(|tag| tag.strip_prefix("NM:i:"))
            .map(|value| {
                value
                    .parse()
                    .map_err(|_| format!("invalid NM tag {value:?}"))
            })
            .transpose()?;
        Ok((
            fields[0].to_owned(),
            Self {
                flag,
                contig: fields[2].to_owned(),
                mapq,
                cigar: fields[5].to_owned(),
                nm,
            },
        ))
    }

    pub fn mate(&self) -> Option<u8> {
        match (self.flag & FIRST != 0, self.flag & SECOND != 0) {
            (true, false) => Some(1),
            (false, true) => Some(2),
            _ => None,
        }
    }

    pub fn is_primary(&self) -> bool {
        self.flag & (SECONDARY | SUPPLEMENTARY) == 0
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Outcome {
    ReferenceAssigned,
    DiscoveryUnassigned,
    DiscoveryAmbiguous,
}

impl Outcome {
    pub fn label(self) -> &'static str {
        match self {
            Self::ReferenceAssigned => "reference-assigned",
            Self::DiscoveryUnassigned => "discovery-unassigned",
            Self::DiscoveryAmbiguous => "discovery-ambiguous",
        }
    }

    pub fn retained(self) -> bool {
        self != Self::ReferenceAssigned
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Decision {
    pub outcome: Outcome,
    pub reason: &'static str,
    /// The reference of mate 1's primary alignment, when mapped.
    pub reference: Option<String>,
    pub species: Option<String>,
    pub mapq: [Option<u8>; 2],
}

pub struct Measured {
    pub identity: f64,
    pub aligned_fraction: f64,
}

/// Identity and aligned fraction from CIGAR and NM. `None` when they cannot
/// be established. Hard clips count as unaligned query only if `allow_hard`.
pub fn measure(alignment: &Alignment, allow_hard: bool) -> Option<Measured> {
    let nm = f64::from(alignment.nm?);
    let (mut aligned, mut deleted, mut clipped) = (0u64, 0u64, 0u64);
    let mut number = 0u64;
    for byte in alignment.cigar.bytes() {
        if byte.is_ascii_digit() {
            number = number
                .checked_mul(10)?
                .checked_add(u64::from(byte - b'0'))?;
            continue;
        }
        match byte {
            b'M' | b'=' | b'X' | b'I' => aligned += number,
            b'D' => deleted += number,
            b'S' => clipped += number,
            b'H' if allow_hard => clipped += number,
            _ => return None, // N, P, unexpected H, or malformed
        }
        number = 0;
    }
    if aligned == 0 || number != 0 {
        return None;
    }
    Some(Measured {
        identity: 1.0 - nm / (aligned + deleted) as f64,
        aligned_fraction: aligned as f64 / (aligned + clipped) as f64,
    })
}

fn qualifies(measured: &Measured, policy: &Policy) -> bool {
    measured.identity >= policy.min_identity
        && measured.aligned_fraction >= policy.min_aligned_fraction
}

fn decision(outcome: Outcome, reason: &'static str) -> Decision {
    Decision {
        outcome,
        reason,
        reference: None,
        species: None,
        mapq: [None, None],
    }
}

/// Decides one pair from all of its alignment records.
///
/// Returns `Err` for structurally impossible input, such as a mate without
/// exactly one primary record or a contig outside the selected references.
pub fn decide(
    records: &[Alignment],
    contigs: &Contigs,
    policy: &Policy,
) -> Result<Decision, String> {
    let mut primaries: [Option<&Alignment>; 2] = [None, None];
    for record in records.iter().filter(|r| r.is_primary()) {
        let mate = record
            .mate()
            .ok_or("alignment record is not marked as mate 1 or 2")?;
        let slot = &mut primaries[usize::from(mate - 1)];
        if slot.replace(record).is_some() {
            return Err(format!(
                "mate {mate} has more than one primary alignment; adjacent pairs may share a read name"
            ));
        }
    }
    let [Some(first), Some(second)] = primaries else {
        return Err("a mate has no primary alignment record".into());
    };
    let mapq = [Some(first.mapq), Some(second.mapq)];
    let with_mapq = |mut d: Decision| {
        d.mapq = mapq;
        d
    };
    if first.flag & UNMAPPED != 0 || second.flag & UNMAPPED != 0 {
        return Ok(with_mapq(decision(
            Outcome::DiscoveryUnassigned,
            "mate-unmapped",
        )));
    }
    let lookup = |alignment: &Alignment| {
        contigs.get(&alignment.contig).ok_or_else(|| {
            format!(
                "alignment to {:?}, which is not a selected reference contig",
                alignment.contig
            )
        })
    };
    let (reference, species) = lookup(first)?.clone();
    let (_, second_species) = lookup(second)?;
    let mut located = with_mapq(decision(Outcome::DiscoveryUnassigned, ""));
    located.reference = Some(reference);
    located.species = Some(species.clone());
    let retain = |outcome, reason| {
        let mut d = located.clone();
        d.outcome = outcome;
        d.reason = reason;
        Ok(d)
    };

    for primary in [first, second] {
        if primary.nm.is_none() {
            return retain(Outcome::DiscoveryUnassigned, "missing-nm");
        }
        let Some(measured) = measure(primary, false) else {
            return retain(Outcome::DiscoveryUnassigned, "unsupported-cigar");
        };
        if measured.identity < policy.min_identity {
            return retain(Outcome::DiscoveryUnassigned, "low-identity");
        }
        if measured.aligned_fraction < policy.min_aligned_fraction {
            return retain(Outcome::DiscoveryUnassigned, "low-aligned-fraction");
        }
    }
    if first.mapq == 255 || second.mapq == 255 {
        return retain(Outcome::DiscoveryAmbiguous, "mapq-unavailable");
    }
    if records.iter().any(|r| r.flag & SUPPLEMENTARY != 0) {
        return retain(Outcome::DiscoveryAmbiguous, "supplementary-alignment");
    }
    if first.flag & PROPER_PAIR == 0 || second.flag & PROPER_PAIR == 0 {
        return retain(Outcome::DiscoveryAmbiguous, "discordant-pair");
    }
    if *second_species != species {
        return retain(Outcome::DiscoveryAmbiguous, "mates-different-species");
    }
    let mut alternatives = false;
    let mut secondaries = [0usize; 2];
    for record in records.iter().filter(|r| r.flag & SECONDARY != 0) {
        if let Some(mate) = record.mate() {
            secondaries[usize::from(mate - 1)] += 1;
        }
        if record.flag & UNMAPPED != 0 {
            continue;
        }
        if !measure(record, true).is_some_and(|m| qualifies(&m, policy)) {
            continue;
        }
        alternatives = true;
        if lookup(record)?.1 != species {
            return retain(Outcome::DiscoveryAmbiguous, "cross-species-alternative");
        }
    }
    for (index, primary) in [first, second].into_iter().enumerate() {
        if primary.mapq < policy.min_mapq && secondaries[index] >= policy.max_secondary {
            return retain(Outcome::DiscoveryAmbiguous, "alternatives-truncated");
        }
    }
    let low_mapq = first.mapq < policy.min_mapq || second.mapq < policy.min_mapq;
    let reason = match (alternatives, low_mapq) {
        (true, _) => "same-species-alternatives",
        (false, true) => "low-mapq-no-reported-alternative",
        (false, false) => "unique",
    };
    retain(Outcome::ReferenceAssigned, reason)
}

#[cfg(test)]
mod tests {
    use super::*;

    const POLICY: Policy = Policy {
        min_identity: 0.98,
        min_aligned_fraction: 0.90,
        min_mapq: 20,
        max_secondary: 3,
    };

    fn contigs() -> Contigs {
        [
            ("ecoliA::chr", "ecoliA", "562"),
            ("ecoliB::chr", "ecoliB", "562"),
            ("ecoliA::plasmid", "ecoliA", "562"),
            ("bsub::chr", "bsub", "1423"),
        ]
        .into_iter()
        .map(|(contig, reference, species)| {
            (
                contig.to_owned(),
                (reference.to_owned(), species.to_owned()),
            )
        })
        .collect()
    }

    fn record(flag: u16, contig: &str, mapq: u8, cigar: &str, nm: Option<u32>) -> Alignment {
        Alignment {
            flag,
            contig: contig.into(),
            mapq,
            cigar: cigar.into(),
            nm,
        }
    }

    /// A properly paired, perfect primary pair on `contig`.
    fn pair(contig: &str, mapq: u8) -> Vec<Alignment> {
        vec![
            record(0x1 | 0x2 | 0x40, contig, mapq, "150=", Some(0)),
            record(0x1 | 0x2 | 0x80, contig, mapq, "150=", Some(0)),
        ]
    }

    fn secondary(mate: u16, contig: &str, nm: u32) -> Alignment {
        record(0x1 | SECONDARY | mate, contig, 0, "150M", Some(nm))
    }

    fn decide_ok(records: &[Alignment]) -> Decision {
        decide(records, &contigs(), &POLICY).unwrap()
    }

    #[test]
    fn unique_confident_pairs_are_assigned_to_their_reference() {
        let d = decide_ok(&pair("ecoliA::chr", 60));
        assert_eq!(
            (d.outcome, d.reason),
            (Outcome::ReferenceAssigned, "unique")
        );
        assert_eq!(d.reference.as_deref(), Some("ecoliA"));
        assert_eq!(d.species.as_deref(), Some("562"));
        assert_eq!(d.mapq, [Some(60), Some(60)]);
    }

    #[test]
    fn same_species_multimapping_and_repeats_are_assigned() {
        let mut strains = pair("ecoliA::chr", 0);
        strains.extend([
            secondary(0x40, "ecoliB::chr", 0),
            secondary(0x80, "ecoliB::chr", 1),
        ]);
        assert_eq!(decide_ok(&strains).reason, "same-species-alternatives");
        let mut repeat = pair("ecoliA::chr", 3);
        repeat.push(secondary(0x40, "ecoliA::plasmid", 0));
        assert_eq!(decide_ok(&repeat).outcome, Outcome::ReferenceAssigned);
        let d = decide_ok(&pair("bsub::chr", 1));
        assert_eq!(
            (d.outcome, d.reason),
            (
                Outcome::ReferenceAssigned,
                "low-mapq-no-reported-alternative"
            )
        );
    }

    #[test]
    fn cross_species_alternatives_make_the_pair_ambiguous_unless_they_fail_thresholds() {
        let mut close = pair("ecoliA::chr", 0);
        close.push(secondary(0x80, "bsub::chr", 2));
        let d = decide_ok(&close);
        assert_eq!(
            (d.outcome, d.reason),
            (Outcome::DiscoveryAmbiguous, "cross-species-alternative")
        );
        let mut weak = pair("ecoliA::chr", 40);
        weak.push(secondary(0x80, "bsub::chr", 10)); // 93% identity: not a qualifying alternative
        assert_eq!(decide_ok(&weak).reason, "unique");
    }

    #[test]
    fn a_reached_secondary_cap_with_low_mapq_is_ambiguous() {
        let mut capped = pair("ecoliA::chr", 0);
        capped.extend((0..3).map(|_| secondary(0x40, "ecoliB::chr", 0)));
        assert_eq!(decide_ok(&capped).reason, "alternatives-truncated");
        let mut confident = pair("ecoliA::chr", 30);
        confident.extend((0..3).map(|_| secondary(0x40, "ecoliB::chr", 0)));
        assert_eq!(decide_ok(&confident).outcome, Outcome::ReferenceAssigned);
    }

    #[test]
    fn either_mate_failing_retains_both() {
        let cases: Vec<(Vec<Alignment>, Outcome, &str)> = vec![
            (
                vec![
                    record(0x1 | 0x40, "ecoliA::chr", 60, "150=", Some(0)),
                    record(0x1 | 0x4 | 0x80, "*", 0, "*", None),
                ],
                Outcome::DiscoveryUnassigned,
                "mate-unmapped",
            ),
            (
                vec![
                    record(0x3 | 0x40, "ecoliA::chr", 60, "150M", Some(4)),
                    record(0x3 | 0x80, "ecoliA::chr", 60, "150M", Some(0)),
                ],
                Outcome::DiscoveryUnassigned,
                "low-identity",
            ),
            (
                vec![
                    record(0x3 | 0x40, "ecoliA::chr", 60, "120M30S", Some(0)),
                    record(0x3 | 0x80, "ecoliA::chr", 60, "150M", Some(0)),
                ],
                Outcome::DiscoveryUnassigned,
                "low-aligned-fraction",
            ),
            (
                vec![
                    record(0x3 | 0x40, "ecoliA::chr", 60, "150M", None),
                    record(0x3 | 0x80, "ecoliA::chr", 60, "150M", Some(0)),
                ],
                Outcome::DiscoveryUnassigned,
                "missing-nm",
            ),
            (
                vec![
                    record(0x3 | 0x40, "ecoliA::chr", 60, "140M10H", Some(0)),
                    record(0x3 | 0x80, "ecoliA::chr", 60, "150M", Some(0)),
                ],
                Outcome::DiscoveryUnassigned,
                "unsupported-cigar",
            ),
            (
                vec![
                    record(0x1 | 0x40, "ecoliA::chr", 60, "150=", Some(0)),
                    record(0x1 | 0x80, "ecoliA::chr", 60, "150=", Some(0)),
                ],
                Outcome::DiscoveryAmbiguous,
                "discordant-pair",
            ),
            (
                vec![
                    record(0x3 | 0x40, "ecoliA::chr", 60, "150=", Some(0)),
                    record(0x3 | 0x80, "bsub::chr", 60, "150=", Some(0)),
                ],
                Outcome::DiscoveryAmbiguous,
                "mates-different-species",
            ),
            (
                vec![
                    record(0x3 | 0x40, "ecoliA::chr", 255, "150=", Some(0)),
                    record(0x3 | 0x80, "ecoliA::chr", 60, "150=", Some(0)),
                ],
                Outcome::DiscoveryAmbiguous,
                "mapq-unavailable",
            ),
        ];
        for (records, outcome, reason) in cases {
            let d = decide_ok(&records);
            assert_eq!((d.outcome, d.reason), (outcome, reason));
            assert!(d.outcome.retained());
        }
        let mut chimeric = pair("ecoliA::chr", 60);
        chimeric.push(record(
            0x1 | SUPPLEMENTARY | 0x40,
            "bsub::chr",
            60,
            "70H80M",
            Some(0),
        ));
        assert_eq!(decide_ok(&chimeric).reason, "supplementary-alignment");
    }

    #[test]
    fn identity_counts_deletions_and_aligned_fraction_counts_soft_clips() {
        // 140 aligned + 5 deleted, NM 2: identity 1 - 2/145 > 0.98; aligned 140/150 > 0.90.
        let records = vec![
            record(0x3 | 0x40, "ecoliA::chr", 60, "10S70=5D70=", Some(2)),
            record(0x3 | 0x80, "ecoliA::chr", 60, "150=", Some(0)),
        ];
        assert_eq!(decide_ok(&records).outcome, Outcome::ReferenceAssigned);
    }

    #[test]
    fn impossible_structures_are_errors_not_decisions() {
        let mut two_primaries = pair("ecoliA::chr", 60);
        two_primaries.push(record(0x3 | 0x40, "ecoliA::chr", 60, "150=", Some(0)));
        assert!(
            decide(&two_primaries, &contigs(), &POLICY)
                .unwrap_err()
                .contains("more than one primary")
        );
        assert!(decide(&pair("ecoliA::chr", 60)[..1], &contigs(), &POLICY).is_err());
        assert!(
            decide(&pair("unknown::chr", 60), &contigs(), &POLICY)
                .unwrap_err()
                .contains("not a selected reference")
        );
    }

    #[test]
    fn sam_lines_parse_flags_mapq_cigar_and_nm() {
        let (name, a) = Alignment::parse(
            "r1\t99\tecoliA::chr\t5\t60\t150=\t=\t100\t250\tACGT\tIIII\tAS:i:300\tNM:i:2",
        )
        .unwrap();
        assert_eq!(name, "r1");
        assert_eq!((a.flag, a.mapq, a.nm, a.mate()), (99, 60, Some(2), Some(1)));
        assert!(Alignment::parse("r1\t99\tchr").is_err());
        assert!(Alignment::parse("r1\tx\tc\t1\t60\t1M\t*\t0\t0\tA\tI").is_err());
    }
}
