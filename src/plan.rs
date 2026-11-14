use std::collections::BTreeMap;

use serde::Serialize;

use crate::{
    cli::{IterationStrategy, PlanArgs},
    error::Error,
    manifest::{
        Metadata, Sample, ValidationReport, read_metadata, read_samples, validation_report,
    },
    tools::{MappingMode, Profiler, ToolSelection},
};

#[derive(Debug, Serialize)]
pub struct ResearchPlan {
    pub schema_version: u32,
    pub executable: bool,
    pub summary: ValidationReport,
    pub samples: Vec<Sample>,
    pub metadata: BTreeMap<String, BTreeMap<String, String>>,
    pub group_by: String,
    pub groups: Vec<Group>,
    pub resources: Resources,
    pub tools: ToolSelection,
    pub read_extractor: &'static str,
    pub iteration: Iteration,
    pub read_retention: &'static str,
    pub stages: Vec<&'static str>,
    pub unresolved: Vec<&'static str>,
}

#[derive(Debug, Serialize)]
pub struct Group {
    pub group_id: String,
    pub metadata_value: String,
    pub sample_ids: Vec<String>,
}

#[derive(Debug, Serialize)]
pub struct Resources {
    pub threads: usize,
    pub memory_gib: u32,
    pub jobs: usize,
    pub max_group_samples: usize,
}

#[derive(Debug, Serialize)]
pub struct Iteration {
    pub strategy: IterationStrategy,
    pub max_rounds: u32,
}

pub fn make_plan(args: &PlanArgs) -> Result<ResearchPlan, Error> {
    args.validate()?;
    let samples = read_samples(&args.samples)?;
    let metadata = read_metadata(&args.metadata, &samples)?;
    let groups = group_samples(
        &samples,
        &metadata,
        &args.group_by,
        args.max_group_samples.get(),
    )?;
    let mut unresolved = vec![
        "reference-catalogue-and-taxonomy-crosswalk",
        "profiler-adapter",
        "reference-mapper-adapter",
        "read-extractor-adapter",
        "coassembler-adapter",
        "assignment-and-recruitment-thresholds",
        "convergence-criteria",
    ];
    if args.tools.profiler == Profiler::Locove
        || args.tools.mapping_mode == MappingMode::LocoveGuided
    {
        unresolved.push("locove-evidence-policy-and-source-revision");
    }
    if args.tools.mapping_mode == MappingMode::LocoveGuided {
        unresolved.extend([
            "region-padding-and-coordinate-liftover",
            "full-reference-confidence-validation",
            "guided-mapping-correctness-and-speed-evaluation",
        ]);
    }
    Ok(ResearchPlan {
        schema_version: 2,
        executable: false,
        summary: validation_report(&samples, Some(&metadata)),
        samples,
        metadata: metadata.rows,
        group_by: args.group_by.clone(),
        groups,
        tools: args.tools.clone(),
        read_extractor: "samtools-with-pair-level-decision",
        resources: Resources {
            threads: args.threads.get(),
            memory_gib: args.memory_gib.get(),
            jobs: args.jobs.get(),
            max_group_samples: args.max_group_samples.get(),
        },
        iteration: Iteration {
            strategy: args.iteration_strategy,
            max_rounds: args.max_rounds.get(),
        },
        read_retention: "retain-both-if-either-unassigned",
        stages: vec![
            "profile-references-per-sample",
            "map-selected-references-per-sample",
            "retain-unassigned-read-pairs-per-sample",
            "group-discovery-pools-by-metadata",
            "coassemble-per-group",
            "map-cumulative-contigs-and-evaluate-next-round",
        ],
        unresolved,
    })
}

fn group_samples(
    samples: &[Sample],
    metadata: &Metadata,
    group_by: &str,
    limit: usize,
) -> Result<Vec<Group>, Error> {
    if !metadata.columns.iter().any(|column| column == group_by) {
        return Err(Error::InvalidInput(format!(
            "metadata has no grouping column {group_by:?}"
        )));
    }
    let mut by_value: BTreeMap<&str, Vec<String>> = BTreeMap::new();
    for sample in samples {
        let value = &metadata.rows[&sample.sample_id][group_by];
        if value.trim().is_empty() {
            return Err(Error::InvalidInput(format!(
                "{}: grouping column {group_by:?} is empty",
                sample.sample_id
            )));
        }
        by_value
            .entry(value)
            .or_default()
            .push(sample.sample_id.clone());
    }
    let mut groups = Vec::new();
    for (value, mut ids) in by_value {
        ids.sort();
        for chunk in ids.chunks(limit) {
            groups.push(Group {
                group_id: format!("group-{:06}", groups.len() + 1),
                metadata_value: value.to_owned(),
                sample_ids: chunk.to_vec(),
            });
        }
    }
    Ok(groups)
}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;

    use super::*;

    fn inputs(rows: &[(&str, &str)]) -> (Vec<Sample>, Metadata) {
        let samples = rows
            .iter()
            .map(|(id, _)| Sample {
                sample_id: (*id).into(),
                read1: PathBuf::new(),
                read2: PathBuf::new(),
            })
            .collect();
        let rows = rows
            .iter()
            .map(|(id, value)| {
                (
                    (*id).into(),
                    BTreeMap::from([("habitat".into(), (*value).into())]),
                )
            })
            .collect();
        (
            samples,
            Metadata {
                columns: vec!["habitat".into()],
                rows,
            },
        )
    }

    #[test]
    fn batching_preserves_every_sample_once_and_keeps_metadata_groups_separate() {
        let (samples, metadata) = inputs(&[
            ("d", "soil"),
            ("b", "soil"),
            ("c", "soil"),
            ("a", "water"),
            ("e", "soil"),
        ]);
        let groups = group_samples(&samples, &metadata, "habitat", 3).unwrap();
        assert_eq!(groups.len(), 3);
        assert_eq!(groups[0].sample_ids, ["b", "c", "d"]);
        assert_eq!(groups[1].sample_ids, ["e"]);
        assert_eq!(groups[2].sample_ids, ["a"]);
        assert_eq!(groups[0].metadata_value, "soil");
        assert_eq!(groups[2].metadata_value, "water");
        let mut ids: Vec<_> = groups
            .iter()
            .flat_map(|g| g.sample_ids.iter().cloned())
            .collect();
        ids.sort();
        assert_eq!(ids, ["a", "b", "c", "d", "e"]);
        assert!(groups.iter().all(|g| g.sample_ids.len() <= 3));
    }

    #[test]
    fn group_membership_and_ids_do_not_depend_on_sample_row_order() {
        let (mut samples, metadata) = inputs(&[("z", "water"), ("b", "soil"), ("a", "soil")]);
        let first = group_samples(&samples, &metadata, "habitat", 1).unwrap();
        samples.reverse();
        let second = group_samples(&samples, &metadata, "habitat", 1).unwrap();
        assert_eq!(
            serde_json::to_value(&first).unwrap(),
            serde_json::to_value(&second).unwrap()
        );
        assert_eq!(first[0].group_id, "group-000001");
    }

    #[test]
    fn missing_column_or_empty_group_values_cannot_create_an_implicit_group() {
        let (samples, metadata) = inputs(&[("a", "soil")]);
        assert!(group_samples(&samples, &metadata, "absent", 1).is_err());
        for value in ["", "  "] {
            let (samples, metadata) = inputs(&[("a", value)]);
            assert!(group_samples(&samples, &metadata, "habitat", 1).is_err());
        }
    }

    #[test]
    fn raw_metadata_values_never_become_output_path_components() {
        let (samples, metadata) = inputs(&[("a", "../../soil"), ("b", "SOIL"), ("c", "soil")]);
        let groups = group_samples(&samples, &metadata, "habitat", 32).unwrap();
        assert_eq!(groups.len(), 3);
        assert!(
            groups
                .iter()
                .all(|group| group.group_id.starts_with("group-"))
        );
        assert_eq!(groups[0].metadata_value, "../../soil");
    }
}
