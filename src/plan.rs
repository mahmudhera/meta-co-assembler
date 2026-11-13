use std::collections::BTreeMap;

use serde::Serialize;

use crate::{
    cli::{IterationStrategy, PlanArgs},
    error::Error,
    manifest::{
        Metadata, Sample, ValidationReport, read_metadata, read_samples, validation_report,
    },
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
    Ok(ResearchPlan {
        schema_version: 1,
        executable: false,
        summary: validation_report(&samples, Some(&metadata)),
        samples,
        metadata: metadata.rows,
        group_by: args.group_by.clone(),
        groups,
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
        unresolved: vec![
            "reference-catalogue-and-taxonomy-crosswalk",
            "profiler",
            "reference-mapper",
            "read-extractor",
            "coassembler",
            "assignment-and-recruitment-thresholds",
            "convergence-criteria",
        ],
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
