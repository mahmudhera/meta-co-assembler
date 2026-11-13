use std::{
    collections::{BTreeMap, BTreeSet},
    fs::{self, File},
    path::{Path, PathBuf},
};

use csv::{Reader, ReaderBuilder, StringRecord};
use serde::Serialize;

use crate::error::Error;

#[derive(Debug, Serialize)]
pub struct Sample {
    pub sample_id: String,
    pub read1: PathBuf,
    pub read2: PathBuf,
}

#[derive(Debug)]
pub struct Metadata {
    pub columns: Vec<String>,
    pub rows: BTreeMap<String, BTreeMap<String, String>>,
}

#[derive(Debug, Serialize)]
pub struct ValidationReport {
    pub schema_version: u32,
    pub sample_count: usize,
    pub metadata_columns: Vec<String>,
    pub missing_metadata_values: BTreeMap<String, usize>,
}

fn open_tsv(path: &Path) -> Result<Reader<File>, Error> {
    ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .map_err(|source| Error::Csv {
            context: format!("cannot open TSV {}", path.display()),
            source,
        })
}

fn csv_error(path: &Path, source: csv::Error) -> Error {
    Error::Csv {
        context: format!("invalid TSV {}", path.display()),
        source,
    }
}

fn headers(reader: &mut Reader<File>, path: &Path) -> Result<StringRecord, Error> {
    let headers = reader.headers().map_err(|e| csv_error(path, e))?.clone();
    let mut seen = BTreeSet::new();
    for header in &headers {
        if header.trim().is_empty() || header.trim() != header || !seen.insert(header) {
            return Err(Error::InvalidInput(format!(
                "{}: column names must be nonempty, unique, and have no surrounding whitespace",
                path.display()
            )));
        }
    }
    Ok(headers)
}

fn column(headers: &StringRecord, name: &str, path: &Path) -> Result<usize, Error> {
    headers
        .iter()
        .position(|header| header == name)
        .ok_or_else(|| {
            Error::InvalidInput(format!(
                "{}: missing required column {name}",
                path.display()
            ))
        })
}

fn validate_id(id: &str, path: &Path) -> Result<(), Error> {
    let mut bytes = id.bytes();
    if !bytes.next().is_some_and(|b| b.is_ascii_alphanumeric())
        || !bytes.all(|b| b.is_ascii_alphanumeric() || b"._-".contains(&b))
    {
        return Err(Error::InvalidInput(format!(
            "{}: invalid sample_id {id:?}; use an ASCII letter/digit followed by letters, digits, ., _, or -",
            path.display()
        )));
    }
    Ok(())
}

fn resolve_read(base: &Path, raw: &str, sample: &str) -> Result<PathBuf, Error> {
    if raw.is_empty() {
        return Err(Error::InvalidInput(format!("{sample}: read path is empty")));
    }
    let path = base.join(raw);
    let supported = [".fastq", ".fq", ".fastq.gz", ".fq.gz"]
        .iter()
        .any(|extension| raw.ends_with(extension));
    if !supported {
        return Err(Error::InvalidInput(format!(
            "{sample}: {} must have a .fastq, .fq, .fastq.gz, or .fq.gz extension",
            path.display()
        )));
    }
    let resolved = fs::canonicalize(&path).map_err(|source| Error::Io {
        context: format!("{sample}: cannot resolve read file {}", path.display()),
        source,
    })?;
    let metadata = fs::metadata(&resolved).map_err(|source| Error::Io {
        context: format!("{sample}: cannot inspect read file {}", resolved.display()),
        source,
    })?;
    if !metadata.is_file() || metadata.len() == 0 {
        return Err(Error::InvalidInput(format!(
            "{sample}: {} must be a nonempty regular file",
            resolved.display()
        )));
    }
    // Opening checks readability without inspecting or buffering the reads.
    File::open(&resolved).map_err(|source| Error::Io {
        context: format!("{sample}: cannot open read file {}", resolved.display()),
        source,
    })?;
    Ok(resolved)
}

pub fn read_samples(path: &Path) -> Result<Vec<Sample>, Error> {
    let mut reader = open_tsv(path)?;
    let headers = headers(&mut reader, path)?;
    let id_col = column(&headers, "sample_id", path)?;
    let read1_col = column(&headers, "read1", path)?;
    let read2_col = column(&headers, "read2", path)?;
    if headers.len() != 3 {
        return Err(Error::InvalidInput(format!(
            "{}: sample table must contain exactly sample_id, read1, read2; put other columns in metadata",
            path.display()
        )));
    }
    let base = path.parent().unwrap_or_else(|| Path::new("."));
    let mut samples = Vec::new();
    let mut ids = BTreeSet::new();
    let mut files = BTreeSet::new();
    for record in reader.records() {
        let record = record.map_err(|e| csv_error(path, e))?;
        let sample_id = &record[id_col];
        validate_id(sample_id, path)?;
        if !ids.insert(sample_id.to_owned()) {
            return Err(Error::InvalidInput(format!(
                "{}: duplicate sample_id {sample_id:?}",
                path.display()
            )));
        }
        let read1 = resolve_read(base, &record[read1_col], sample_id)?;
        let read2 = resolve_read(base, &record[read2_col], sample_id)?;
        for read in [&read1, &read2] {
            if !files.insert(read.clone()) {
                return Err(Error::InvalidInput(format!(
                    "{sample_id}: read file {} is reused as another mate or sample",
                    read.display()
                )));
            }
        }
        samples.push(Sample {
            sample_id: sample_id.to_owned(),
            read1,
            read2,
        });
    }
    if samples.is_empty() {
        return Err(Error::InvalidInput(format!(
            "{}: no samples found",
            path.display()
        )));
    }
    samples.sort_by(|a, b| a.sample_id.cmp(&b.sample_id));
    Ok(samples)
}

pub fn read_metadata(path: &Path, samples: &[Sample]) -> Result<Metadata, Error> {
    let mut reader = open_tsv(path)?;
    let headers = headers(&mut reader, path)?;
    let id_col = column(&headers, "sample_id", path)?;
    if headers.len() < 2 {
        return Err(Error::InvalidInput(format!(
            "{}: metadata needs at least one column besides sample_id",
            path.display()
        )));
    }
    let expected: BTreeSet<_> = samples.iter().map(|s| s.sample_id.as_str()).collect();
    let mut rows = BTreeMap::new();
    for record in reader.records() {
        let record = record.map_err(|e| csv_error(path, e))?;
        let id = &record[id_col];
        if !expected.contains(id) {
            return Err(Error::InvalidInput(format!(
                "{}: metadata contains unknown sample_id {id:?}",
                path.display()
            )));
        }
        let mut values = BTreeMap::new();
        for (i, header) in headers.iter().enumerate() {
            if i != id_col {
                values.insert(header.to_owned(), record[i].to_owned());
            }
        }
        if rows.insert(id.to_owned(), values).is_some() {
            return Err(Error::InvalidInput(format!(
                "{}: duplicate metadata sample_id {id:?}",
                path.display()
            )));
        }
    }
    for id in expected {
        if !rows.contains_key(id) {
            return Err(Error::InvalidInput(format!(
                "{}: missing metadata for sample_id {id:?}",
                path.display()
            )));
        }
    }
    let mut columns: Vec<_> = headers
        .iter()
        .filter(|h| *h != "sample_id")
        .map(str::to_owned)
        .collect();
    columns.sort();
    Ok(Metadata { columns, rows })
}

pub fn validation_report(samples: &[Sample], metadata: Option<&Metadata>) -> ValidationReport {
    let metadata_columns = metadata.map(|m| m.columns.clone()).unwrap_or_default();
    let mut missing_metadata_values = BTreeMap::new();
    if let Some(metadata) = metadata {
        for column in &metadata.columns {
            let missing = metadata
                .rows
                .values()
                .filter(|row| row[column].trim().is_empty())
                .count();
            missing_metadata_values.insert(column.clone(), missing);
        }
    }
    ValidationReport {
        schema_version: 1,
        sample_count: samples.len(),
        metadata_columns,
        missing_metadata_values,
    }
}
