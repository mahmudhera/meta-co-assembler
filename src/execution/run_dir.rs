//! Output-directory ownership, layout, and atomic publication.
//!
//! A run owns its output directory through a lock. An existing, non-empty
//! directory is refused unless `force` is set; then only a directory that a
//! previous run marked as its own is cleared. Resume is not supported.

use std::{
    fs::{self, File, OpenOptions},
    io::{self, Write},
    path::{Component, Path, PathBuf},
};

use serde::Serialize;

use crate::error::Error;

pub const MANIFEST: &str = "run.json";
pub const PRODUCER: &str = "meta-co-assembler";
const LOCK_FILE: &str = ".meta-co-assembler.lock";
const SCRATCH_MARKER: &str = ".meta-co-assembler-scratch";

fn io_error(context: String) -> impl FnOnce(io::Error) -> Error {
    move |source| Error::Io { context, source }
}

/// Absolute locations used by a run. Paths are resolved through existing
/// symbolic links so that containment checks compare real locations.
#[derive(Clone, Debug)]
pub struct RunLayout {
    output_dir: PathBuf,
    scratch_dir: PathBuf,
}

impl RunLayout {
    /// `scratch_dir` defaults to `OUTPUT/scratch`. Nothing is created.
    pub fn new(output_dir: &Path, scratch_dir: Option<&Path>) -> Result<Self, Error> {
        let output_dir = resolve(output_dir)?;
        let scratch_dir = match scratch_dir {
            Some(path) => resolve(path)?,
            None => output_dir.join("scratch"),
        };
        Ok(Self {
            output_dir,
            scratch_dir,
        })
    }

    pub fn output_dir(&self) -> &Path {
        &self.output_dir
    }

    pub fn scratch_dir(&self) -> &Path {
        &self.scratch_dir
    }

    pub fn manifest(&self) -> PathBuf {
        self.output_dir.join(MANIFEST)
    }

    /// Published, validated stage outputs. Exists only for complete stages.
    pub fn stage_dir(&self, id: &str) -> PathBuf {
        self.output_dir.join("stages").join(id)
    }

    /// Where a stage runs. Remains after failure or interruption.
    pub fn partial_dir(&self, id: &str) -> PathBuf {
        self.output_dir.join("partial").join(id)
    }

    /// Temporary space for one task; also its `TMPDIR`.
    pub fn task_scratch_dir(&self, id: &str) -> PathBuf {
        self.scratch_dir.join(id)
    }

    fn scratch_is_internal(&self) -> bool {
        self.scratch_dir.starts_with(&self.output_dir)
    }
}

/// Makes a path absolute, canonicalizing its longest existing ancestor.
/// Components below that ancestor must be plain names.
fn resolve(path: &Path) -> Result<PathBuf, Error> {
    let absolute = std::path::absolute(path)
        .map_err(io_error(format!("cannot resolve {}", path.display())))?;
    let mut existing = absolute.as_path();
    let mut missing = Vec::new();
    while !existing.exists() {
        match (existing.parent(), existing.file_name()) {
            (Some(parent), Some(name)) => {
                missing.push(name.to_owned());
                existing = parent;
            }
            _ => break,
        }
    }
    let mut resolved = fs::canonicalize(existing)
        .map_err(io_error(format!("cannot resolve {}", existing.display())))?;
    for name in missing.iter().rev() {
        if !matches!(
            Path::new(name).components().next(),
            Some(Component::Normal(_))
        ) {
            return Err(Error::InvalidInput(format!(
                "{}: use a path without '.' or '..' below its existing parent",
                path.display()
            )));
        }
        resolved.push(name);
    }
    Ok(resolved)
}

/// Exclusive ownership of an output directory for the lifetime of a run.
pub(crate) struct RunLock {
    #[cfg(unix)]
    _lock: nix::fcntl::Flock<File>,
    #[cfg(not(unix))]
    path: PathBuf,
}

#[cfg(not(unix))]
impl Drop for RunLock {
    fn drop(&mut self) {
        let _ = fs::remove_file(&self.path);
    }
}

fn lock(output_dir: &Path) -> Result<RunLock, Error> {
    let path = output_dir.join(LOCK_FILE);
    let busy = || {
        Error::Execution(format!(
            "another run is using {}; wait for it to finish",
            output_dir.display()
        ))
    };
    #[cfg(unix)]
    {
        use nix::{
            errno::Errno,
            fcntl::{Flock, FlockArg},
        };
        let file = OpenOptions::new()
            .create(true)
            .truncate(false)
            .write(true)
            .open(&path)
            .map_err(io_error(format!("cannot open lock {}", path.display())))?;
        match Flock::lock(file, FlockArg::LockExclusiveNonblock) {
            Ok(lock) => Ok(RunLock { _lock: lock }),
            Err((_, Errno::EWOULDBLOCK)) => Err(busy()),
            Err((_, errno)) => Err(Error::Execution(format!(
                "cannot lock {}: {errno}",
                path.display()
            ))),
        }
    }
    #[cfg(not(unix))]
    {
        // Best effort: a crash leaves this file behind for manual removal.
        match OpenOptions::new().write(true).create_new(true).open(&path) {
            Ok(_) => Ok(RunLock { path }),
            Err(error) if error.kind() == io::ErrorKind::AlreadyExists => Err(busy()),
            Err(error) => Err(Error::Io {
                context: format!("cannot create lock {}", path.display()),
                source: error,
            }),
        }
    }
}

fn entries_except(dir: &Path, keep: &str) -> Result<Vec<PathBuf>, Error> {
    let mut entries = Vec::new();
    for entry in fs::read_dir(dir).map_err(io_error(format!("cannot list {}", dir.display())))? {
        let entry = entry.map_err(io_error(format!("cannot list {}", dir.display())))?;
        if entry.file_name() != keep {
            entries.push(entry.path());
        }
    }
    entries.sort();
    Ok(entries)
}

fn remove_entries(entries: &[PathBuf]) -> Result<(), Error> {
    for entry in entries {
        let metadata = match entry.symlink_metadata() {
            Ok(metadata) => metadata,
            // Already removed, for example when old and new scratch coincide.
            Err(error) if error.kind() == io::ErrorKind::NotFound => continue,
            Err(error) => {
                return Err(Error::Io {
                    context: format!("cannot inspect {}", entry.display()),
                    source: error,
                });
            }
        };
        // remove_dir_all does not follow symbolic links.
        let result = if metadata.is_dir() {
            fs::remove_dir_all(entry)
        } else {
            fs::remove_file(entry)
        };
        result.map_err(io_error(format!("cannot remove {}", entry.display())))?;
    }
    Ok(())
}

fn previous_run(output_dir: &Path) -> Option<serde_json::Value> {
    let text = fs::read_to_string(output_dir.join(MANIFEST)).ok()?;
    let value: serde_json::Value = serde_json::from_str(&text).ok()?;
    (value["producer"] == PRODUCER).then_some(value)
}

/// `Ok(None)` for an empty directory, `Ok(Some(previous manifest))` for a
/// previous run that `force` allows deleting, and an error otherwise.
fn reusable(output: &Path, force: bool) -> Result<Option<serde_json::Value>, Error> {
    if entries_except(output, LOCK_FILE)?.is_empty() {
        return Ok(None);
    }
    if !force {
        return Err(Error::InvalidInput(format!(
            "{} is not empty; existing outputs are never overwritten. Pass --force to delete the previous run and rerun everything",
            output.display()
        )));
    }
    previous_run(output).map(Some).ok_or_else(|| {
        Error::InvalidInput(format!(
            "refusing to delete {}: it has no meta-co-assembler {MANIFEST}",
            output.display()
        ))
    })
}

/// Entries of an external scratch directory that `force` may clear. Internal
/// scratch is cleared with the output directory.
fn scratch_leftovers(layout: &RunLayout, force: bool) -> Result<Vec<PathBuf>, Error> {
    let scratch = layout.scratch_dir();
    if layout.scratch_is_internal() {
        return Ok(Vec::new());
    }
    let leftovers = entries_except(scratch, SCRATCH_MARKER)?;
    let marked = scratch.join(SCRATCH_MARKER).is_file();
    if !leftovers.is_empty() && !(force && marked) {
        return Err(Error::InvalidInput(format!(
            "scratch directory {} is not empty{}",
            scratch.display(),
            if marked {
                "; pass --force to clear the previous run's scratch"
            } else {
                "; it was not created by meta-co-assembler"
            }
        )));
    }
    Ok(leftovers)
}

pub(crate) struct Prepared {
    pub lock: RunLock,
    pub removed: Vec<String>,
}

/// Claims the output and scratch directories. Refuses to touch anything that
/// contains a protected input or that was not created by a previous run.
pub(crate) fn prepare(
    layout: &RunLayout,
    force: bool,
    protected: &[PathBuf],
) -> Result<Prepared, Error> {
    let output = layout.output_dir();
    let scratch = layout.scratch_dir();
    if scratch == output || output.starts_with(scratch) {
        return Err(Error::InvalidInput(format!(
            "scratch directory {} must not contain the output directory",
            scratch.display()
        )));
    }
    let mut inputs = Vec::with_capacity(protected.len());
    for input in protected {
        let input = fs::canonicalize(input).map_err(io_error(format!(
            "cannot resolve protected input {}",
            input.display()
        )))?;
        for owned in [output, scratch] {
            if input.starts_with(owned) {
                return Err(Error::InvalidInput(format!(
                    "{} contains input {}; choose a separate output location so inputs cannot be deleted",
                    owned.display(),
                    input.display()
                )));
            }
        }
        inputs.push(input);
    }
    if output.exists() && !output.is_dir() {
        return Err(Error::InvalidInput(format!(
            "{} exists and is not a directory",
            output.display()
        )));
    }
    // Check before creating the lock file so a refused, unrelated directory
    // is left exactly as it was; check again once the lock is held.
    if output.exists() {
        reusable(output, force)?;
    }
    if scratch.exists() {
        scratch_leftovers(layout, force)?;
    }
    fs::create_dir_all(output).map_err(io_error(format!("cannot create {}", output.display())))?;
    let lock = lock(output)?;

    let mut removed = Vec::new();
    let existing = entries_except(output, LOCK_FILE)?;
    let previous = reusable(output, force)?;
    let leftovers = if scratch.exists() {
        scratch_leftovers(layout, force)?
    } else {
        Vec::new()
    };
    if let Some(previous) = previous {
        if let Some(old_scratch) = previous["scratch_dir"].as_str().map(PathBuf::from) {
            let holds_input = inputs.iter().any(|input| input.starts_with(&old_scratch));
            if !old_scratch.starts_with(output)
                && !holds_input
                && old_scratch.join(SCRATCH_MARKER).is_file()
            {
                let entries = entries_except(&old_scratch, SCRATCH_MARKER)?;
                remove_entries(&entries)?;
                removed.extend(
                    entries
                        .iter()
                        .map(|path| path.to_string_lossy().into_owned()),
                );
            }
        }
        remove_entries(&existing)?;
        removed.extend(
            existing
                .iter()
                .map(|path| path.to_string_lossy().into_owned()),
        );
    }

    remove_entries(&leftovers)?;
    removed.extend(
        leftovers
            .iter()
            .map(|path| path.to_string_lossy().into_owned()),
    );
    fs::create_dir_all(scratch)
        .map_err(io_error(format!("cannot create {}", scratch.display())))?;
    if !layout.scratch_is_internal() {
        File::create(scratch.join(SCRATCH_MARKER))
            .map_err(io_error(format!("cannot mark {}", scratch.display())))?;
    }
    for name in ["stages", "partial"] {
        let dir = output.join(name);
        fs::create_dir(&dir).map_err(io_error(format!("cannot create {}", dir.display())))?;
    }
    Ok(Prepared { lock, removed })
}

/// Replaces `path` atomically: readers see the old or the new file, never a
/// partial write.
pub(crate) fn write_json_atomic(path: &Path, value: &impl Serialize) -> Result<(), Error> {
    let name = path
        .file_name()
        .map_or_else(Default::default, |n| n.to_string_lossy());
    let temporary = path.with_file_name(format!(".{name}.tmp-{}", std::process::id()));
    let write = || -> io::Result<()> {
        let mut file = File::create(&temporary)?;
        serde_json::to_writer_pretty(&mut file, value)?;
        file.write_all(b"\n")?;
        file.sync_all()?;
        fs::rename(&temporary, path)?;
        sync_parent(path)
    };
    write().map_err(io_error(format!("cannot write {}", path.display())))
}

#[cfg(unix)]
fn sync_dir(dir: &Path) -> io::Result<()> {
    File::open(dir)?.sync_all()
}

#[cfg(not(unix))]
fn sync_dir(_dir: &Path) -> io::Result<()> {
    Ok(())
}

fn sync_parent(path: &Path) -> io::Result<()> {
    path.parent().map_or(Ok(()), sync_dir)
}

/// Flushes every regular file and directory below `dir` to stable storage.
#[cfg(unix)]
fn sync_tree(dir: &Path) -> io::Result<()> {
    for entry in fs::read_dir(dir)? {
        let path = entry?.path();
        let metadata = path.symlink_metadata()?;
        if metadata.is_dir() {
            sync_tree(&path)?;
        } else if metadata.is_file() {
            File::open(&path)?.sync_all()?;
        }
    }
    sync_dir(dir)
}

#[cfg(not(unix))]
fn sync_tree(_dir: &Path) -> io::Result<()> {
    Ok(())
}

/// Moves a validated partial directory to its final location in one rename.
/// Never replaces an existing final directory.
pub(crate) fn publish(partial: &Path, published: &Path) -> io::Result<()> {
    if published.symlink_metadata().is_ok() {
        return Err(io::Error::new(
            io::ErrorKind::AlreadyExists,
            format!("{} already exists", published.display()),
        ));
    }
    sync_tree(partial)?;
    fs::rename(partial, published)?;
    sync_parent(published)?;
    sync_parent(partial)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn temp(name: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("mca-rundir-{}-{name}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir(&dir).unwrap();
        fs::canonicalize(dir).unwrap()
    }

    #[test]
    fn layout_resolves_missing_components_and_rejects_dot_dot_below_them() {
        let root = temp("resolve");
        let layout = RunLayout::new(&root.join("a/b"), None).unwrap();
        assert_eq!(layout.output_dir(), root.join("a/b"));
        assert_eq!(layout.scratch_dir(), root.join("a/b/scratch"));
        assert_eq!(layout.stage_dir("s1"), root.join("a/b/stages/s1"));
        assert_eq!(layout.partial_dir("s1"), root.join("a/b/partial/s1"));
        assert!(RunLayout::new(&root.join("missing/../x"), None).is_err());
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn atomic_json_replaces_the_whole_file() {
        let root = temp("json");
        let path = root.join("run.json");
        write_json_atomic(&path, &serde_json::json!({"n": 1})).unwrap();
        write_json_atomic(&path, &serde_json::json!({"n": 2})).unwrap();
        let value: serde_json::Value = serde_json::from_slice(&fs::read(&path).unwrap()).unwrap();
        assert_eq!(value["n"], 2);
        assert_eq!(
            fs::read_dir(&root).unwrap().count(),
            1,
            "temporary file left behind"
        );
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn publication_never_replaces_an_existing_directory() {
        let root = temp("publish");
        fs::create_dir_all(root.join("partial/s")).unwrap();
        fs::write(root.join("partial/s/out"), "new").unwrap();
        fs::create_dir_all(root.join("stages/s")).unwrap();
        assert!(publish(&root.join("partial/s"), &root.join("stages/s")).is_err());
        assert!(root.join("partial/s/out").exists());
        fs::remove_dir(root.join("stages/s")).unwrap();
        publish(&root.join("partial/s"), &root.join("stages/s")).unwrap();
        assert_eq!(
            fs::read_to_string(root.join("stages/s/out")).unwrap(),
            "new"
        );
        assert!(!root.join("partial/s").exists());
        fs::remove_dir_all(root).unwrap();
    }
}
