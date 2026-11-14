//! Frozen reference catalogues: exact reference IDs, FASTA files, taxonomy,
//! and content digests. Discovery is relative to exactly these sequences.
//!
//! Format (UTF-8 TSV, any column order):
//! `reference_id  fasta  taxid  species_taxid  sha256`.
//! FASTA paths resolve against the catalogue's directory.

use std::{
    collections::{BTreeMap, BTreeSet},
    fs::{self, File},
    io::{self, Read},
    path::{Path, PathBuf},
};

use csv::ReaderBuilder;
use serde::Serialize;
use sha2::{Digest, Sha256};

use crate::error::Error;

pub const COLUMNS: [&str; 5] = ["reference_id", "fasta", "taxid", "species_taxid", "sha256"];
const FASTA_SUFFIXES: [&str; 3] = [".fna", ".fa", ".fasta"];

#[derive(Clone, Debug, Serialize)]
pub struct Reference {
    pub reference_id: String,
    /// Canonical absolute path.
    pub fasta: PathBuf,
    pub taxid: String,
    pub species_taxid: String,
    pub sha256: String,
    pub fasta_bytes: u64,
}

#[derive(Clone, Debug, Serialize)]
pub struct Catalogue {
    pub path: PathBuf,
    /// Digest of the catalogue table itself, which lists every FASTA digest.
    pub sha256: String,
    pub references: Vec<Reference>,
}

impl Catalogue {
    pub fn total_fasta_bytes(&self) -> u64 {
        self.references.iter().fold(0, |sum, reference| {
            sum.saturating_add(reference.fasta_bytes)
        })
    }

    /// Reference indices by species taxonomy ID, in catalogue order.
    pub fn species_members(&self) -> BTreeMap<&str, Vec<usize>> {
        let mut members: BTreeMap<&str, Vec<usize>> = BTreeMap::new();
        for (index, reference) in self.references.iter().enumerate() {
            members
                .entry(reference.species_taxid.as_str())
                .or_default()
                .push(index);
        }
        members
    }
}

/// Streams a file through SHA-256 and returns lowercase hex.
pub fn sha256_file(path: &Path) -> io::Result<String> {
    let mut hasher = Sha256::new();
    let mut file = File::open(path)?;
    let mut buffer = vec![0; 1024 * 1024];
    loop {
        match file.read(&mut buffer) {
            Ok(0) => break,
            Ok(read) => hasher.update(&buffer[..read]),
            Err(error) if error.kind() == io::ErrorKind::Interrupted => {}
            Err(error) => return Err(error),
        }
    }
    Ok(hasher
        .finalize()
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect())
}

/// LoCovE's genome name: the file name with `.gz`, `.fna`, `.fa`, and
/// `.fasta` removed in that order (LoCovE `refindex.rs`, pinned revision).
pub fn locove_genome_name(path: &Path) -> String {
    let mut name = path
        .file_name()
        .map_or_else(|| path.to_string_lossy(), |name| name.to_string_lossy())
        .into_owned();
    for suffix in [".gz", ".fna", ".fa", ".fasta"] {
        if let Some(stripped) = name.strip_suffix(suffix) {
            name = stripped.to_owned();
        }
    }
    name
}

fn invalid(path: &Path, message: impl std::fmt::Display) -> Error {
    Error::InvalidInput(format!("{}: {message}", path.display()))
}

fn valid_id(id: &str) -> bool {
    let mut bytes = id.bytes();
    bytes.next().is_some_and(|b| b.is_ascii_alphanumeric())
        && bytes.all(|b| b.is_ascii_alphanumeric() || b"._-".contains(&b))
}

fn valid_taxid(value: &str) -> bool {
    !value.is_empty() && value.bytes().all(|b| b.is_ascii_digit()) && !value.starts_with('0')
}

fn valid_sha256(value: &str) -> bool {
    value.len() == 64
        && value
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}

/// Parses and checks a catalogue's structure and files. Digests are checked
/// separately by [`verify_digests`] because that reads every FASTA.
pub fn read_catalogue(path: &Path) -> Result<Catalogue, Error> {
    let sha256 = sha256_file(path).map_err(|source| Error::Io {
        context: format!("cannot read catalogue {}", path.display()),
        source,
    })?;
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .map_err(|source| Error::Csv {
            context: format!("cannot open catalogue {}", path.display()),
            source,
        })?;
    let headers = reader
        .headers()
        .map_err(|source| Error::Csv {
            context: format!("invalid catalogue {}", path.display()),
            source,
        })?
        .clone();
    let names: Vec<&str> = headers.iter().collect();
    let mut sorted = names.clone();
    sorted.sort_unstable();
    let mut expected = COLUMNS;
    expected.sort_unstable();
    if sorted != expected {
        return Err(invalid(
            path,
            format!("columns must be exactly {}", COLUMNS.join(", ")),
        ));
    }
    let column = |name: &str| names.iter().position(|header| *header == name).unwrap_or(0);
    let (id_col, fasta_col, taxid_col, species_col, sha_col) = (
        column("reference_id"),
        column("fasta"),
        column("taxid"),
        column("species_taxid"),
        column("sha256"),
    );
    let base = path.parent().unwrap_or_else(|| Path::new("."));
    let mut references = Vec::new();
    let (mut ids, mut files, mut digests, mut locove_names) = (
        BTreeSet::new(),
        BTreeSet::new(),
        BTreeSet::new(),
        BTreeSet::new(),
    );
    for record in reader.records() {
        let record = record.map_err(|source| Error::Csv {
            context: format!("invalid catalogue {}", path.display()),
            source,
        })?;
        let id = &record[id_col];
        if !valid_id(id) {
            return Err(invalid(path, format!("invalid reference_id {id:?}")));
        }
        if !ids.insert(id.to_owned()) {
            return Err(invalid(path, format!("duplicate reference_id {id:?}")));
        }
        let raw = &record[fasta_col];
        if raw.is_empty() || raw.chars().any(char::is_control) {
            return Err(invalid(path, format!("{id}: invalid FASTA path {raw:?}")));
        }
        let plain = raw.strip_suffix(".gz").unwrap_or(raw);
        if !FASTA_SUFFIXES.iter().any(|suffix| plain.ends_with(suffix)) {
            return Err(invalid(
                path,
                format!("{id}: FASTA {raw:?} must end in .fna, .fa, or .fasta, optionally .gz"),
            ));
        }
        let fasta = fs::canonicalize(base.join(raw)).map_err(|source| Error::Io {
            context: format!("{id}: cannot resolve FASTA {}", base.join(raw).display()),
            source,
        })?;
        let metadata = fs::metadata(&fasta).map_err(|source| Error::Io {
            context: format!("{id}: cannot inspect FASTA {}", fasta.display()),
            source,
        })?;
        if !metadata.is_file() || metadata.len() == 0 {
            return Err(invalid(
                path,
                format!("{id}: {} must be a nonempty regular file", fasta.display()),
            ));
        }
        if fasta.to_str().is_none() {
            return Err(invalid(
                path,
                format!("{id}: FASTA path must be valid UTF-8"),
            ));
        }
        if !files.insert(fasta.clone()) {
            return Err(invalid(
                path,
                format!("{id}: FASTA {} is listed twice", fasta.display()),
            ));
        }
        if !locove_names.insert(locove_genome_name(&fasta)) {
            return Err(invalid(
                path,
                format!("{id}: FASTA file names must be unique after removing extensions"),
            ));
        }
        let (taxid, species_taxid, digest) =
            (&record[taxid_col], &record[species_col], &record[sha_col]);
        if !valid_taxid(taxid) || !valid_taxid(species_taxid) {
            return Err(invalid(
                path,
                format!("{id}: taxid and species_taxid must be positive integers"),
            ));
        }
        if !valid_sha256(digest) {
            return Err(invalid(
                path,
                format!("{id}: sha256 must be 64 lowercase hexadecimal digits"),
            ));
        }
        if !digests.insert(digest.to_owned()) {
            return Err(invalid(
                path,
                format!(
                    "{id}: identical FASTA content under two reference IDs makes selection ambiguous"
                ),
            ));
        }
        references.push(Reference {
            reference_id: id.to_owned(),
            fasta,
            taxid: taxid.to_owned(),
            species_taxid: species_taxid.to_owned(),
            sha256: digest.to_owned(),
            fasta_bytes: metadata.len(),
        });
    }
    if references.is_empty() {
        return Err(invalid(path, "no references found"));
    }
    Ok(Catalogue {
        path: fs::canonicalize(path).map_err(|source| Error::Io {
            context: format!("cannot resolve catalogue {}", path.display()),
            source,
        })?,
        sha256,
        references,
    })
}

/// Recomputes every FASTA digest. A mismatch means the catalogue no longer
/// describes the sequences on disk, so selection results would be unbound.
pub fn verify_digests(catalogue: &Catalogue) -> Result<(), Error> {
    for reference in &catalogue.references {
        let actual = sha256_file(&reference.fasta).map_err(|source| Error::Io {
            context: format!(
                "{}: cannot read FASTA {}",
                reference.reference_id,
                reference.fasta.display()
            ),
            source,
        })?;
        if actual != reference.sha256 {
            return Err(Error::InvalidInput(format!(
                "{}: FASTA {} has SHA-256 {actual}, but the catalogue records {}",
                reference.reference_id,
                reference.fasta.display(),
                reference.sha256
            )));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::sync::atomic::{AtomicU64, Ordering};

    use super::*;

    static NEXT: AtomicU64 = AtomicU64::new(0);

    struct Dir(PathBuf);

    impl Dir {
        fn new() -> Self {
            let path = std::env::temp_dir().join(format!(
                "mca-catalogue-{}-{}",
                std::process::id(),
                NEXT.fetch_add(1, Ordering::Relaxed)
            ));
            let _ = fs::remove_dir_all(&path);
            fs::create_dir_all(path.join("refs")).unwrap();
            Self(path)
        }

        fn fasta(&self, name: &str, body: &str) -> String {
            fs::write(self.0.join("refs").join(name), body).unwrap();
            sha256_file(&self.0.join("refs").join(name)).unwrap()
        }

        fn catalogue(&self, rows: &[String]) -> PathBuf {
            let path = self.0.join("catalogue.tsv");
            let mut text = "reference_id\tfasta\ttaxid\tspecies_taxid\tsha256\n".to_owned();
            for row in rows {
                text.push_str(row);
                text.push('\n');
            }
            fs::write(&path, text).unwrap();
            path
        }
    }

    impl Drop for Dir {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.0);
        }
    }

    fn row(id: &str, file: &str, species: &str, digest: &str) -> String {
        format!("{id}\trefs/{file}\t{species}1\t{species}\t{digest}")
    }

    #[test]
    fn sha256_matches_the_standard_test_vector() {
        let dir = Dir::new();
        let path = dir.0.join("abc");
        fs::write(&path, "abc").unwrap();
        assert_eq!(
            sha256_file(&path).unwrap(),
            "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
        );
    }

    #[test]
    fn locove_names_strip_known_extensions_in_order() {
        for (file, name) in [
            ("GCF_000005845.2.fna", "GCF_000005845.2"),
            ("a.fna.gz", "a"),
            ("b.fasta.gz", "b"),
            ("c.fa", "c"),
            ("d.txt", "d.txt"),
        ] {
            assert_eq!(locove_genome_name(Path::new(file)), name, "{file}");
        }
    }

    #[test]
    fn valid_catalogues_resolve_paths_and_group_species() {
        let dir = Dir::new();
        let a = dir.fasta("a.fna", ">a\nACGT\n");
        let b = dir.fasta("b.fa.gz", ">b\nACGA\n");
        let c = dir.fasta("c.fasta", ">c\nTTTT\n");
        let path = dir.catalogue(&[
            row("ref.A", "a.fna", "562", &a),
            row("ref-B", "b.fa.gz", "562", &b),
            row("C", "c.fasta", "1423", &c),
        ]);
        let catalogue = read_catalogue(&path).unwrap();
        verify_digests(&catalogue).unwrap();
        assert_eq!(catalogue.references.len(), 3);
        assert!(catalogue.references[0].fasta.is_absolute());
        assert_eq!(catalogue.species_members()["562"], [0, 1]);
        // Each test FASTA is ">x\nACGT\n": 8 bytes.
        assert_eq!(catalogue.total_fasta_bytes(), 3 * 8);
        assert_eq!(catalogue.sha256, sha256_file(&path).unwrap());
    }

    #[test]
    fn malformed_catalogues_are_input_errors() {
        let dir = Dir::new();
        let a = dir.fasta("a.fna", ">a\nACGT\n");
        let b = dir.fasta("b.fna", ">b\nACGA\n");
        fs::write(dir.0.join("refs/a.txt"), ">a\nACGT\n").unwrap();
        fs::write(dir.0.join("refs/empty.fna"), "").unwrap();
        fs::create_dir_all(dir.0.join("refs/sub")).unwrap();
        fs::write(dir.0.join("refs/sub/a.fna"), ">x\nGGGG\n").unwrap();
        let other = sha256_file(&dir.0.join("refs/sub/a.fna")).unwrap();
        for (name, rows) in [
            ("no rows", vec![]),
            ("bad id", vec![row("../a", "a.fna", "562", &a)]),
            (
                "duplicate id",
                vec![row("a", "a.fna", "562", &a), row("a", "b.fna", "562", &b)],
            ),
            (
                "same file twice",
                vec![row("a", "a.fna", "562", &a), row("b", "./a.fna", "562", &b)],
            ),
            (
                "same content twice",
                vec![row("a", "a.fna", "562", &a), row("b", "b.fna", "562", &a)],
            ),
            (
                "same LoCovE name",
                vec![
                    row("a", "a.fna", "562", &a),
                    row("x", "sub/a.fna", "562", &other),
                ],
            ),
            ("bad extension", vec![row("a", "a.txt", "562", &a)]),
            ("empty file", vec![row("a", "empty.fna", "562", &a)]),
            ("bad taxid", vec![row("a", "a.fna", "0562", &a)]),
            ("short digest", vec![row("a", "a.fna", "562", "abc")]),
            (
                "uppercase digest",
                vec![row("a", "a.fna", "562", &a.to_uppercase())],
            ),
        ] {
            let path = dir.catalogue(&rows.iter().map(String::from).collect::<Vec<_>>());
            let error = read_catalogue(&path).unwrap_err();
            assert_eq!(error.exit_code(), 2, "{name}: {error}");
        }
        fs::write(
            dir.0.join("catalogue.tsv"),
            "reference_id\tfasta\ttaxid\tsha256\n",
        )
        .unwrap();
        assert_eq!(
            read_catalogue(&dir.0.join("catalogue.tsv"))
                .unwrap_err()
                .exit_code(),
            2
        );
    }

    #[test]
    fn missing_fasta_is_an_operational_error_and_changed_content_is_detected() {
        let dir = Dir::new();
        let a = dir.fasta("a.fna", ">a\nACGT\n");
        let path = dir.catalogue(&[row("a", "missing.fna", "562", &a)]);
        assert_eq!(read_catalogue(&path).unwrap_err().exit_code(), 1);
        let path = dir.catalogue(&[row("a", "a.fna", "562", &a)]);
        let catalogue = read_catalogue(&path).unwrap();
        fs::write(dir.0.join("refs/a.fna"), ">a\nACGG\n").unwrap();
        let error = verify_digests(&catalogue).unwrap_err();
        assert_eq!(error.exit_code(), 2);
        assert!(error.to_string().contains("but the catalogue records"));
    }
}
