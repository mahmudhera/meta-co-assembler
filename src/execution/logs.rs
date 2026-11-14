//! Drain child-process streams to capped log files.
//!
//! A drainer always reads its pipe to EOF so a chatty child can never block,
//! but it keeps at most `cap` leading bytes and a `tail`-byte ring buffer.

use std::{
    collections::VecDeque,
    fs::File,
    io::{self, Read, Write},
    path::{Path, PathBuf},
    sync::mpsc,
    thread,
    time::Duration,
};

use serde::Serialize;

const CHUNK: usize = 64 * 1024;

/// What happened to one captured stream.
#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub struct LogRecord {
    pub path: String,
    pub bytes_seen: u64,
    pub bytes_written: u64,
    pub truncated: bool,
    /// False when the stream was still open after the child ended (an
    /// escaped descendant kept it) and the drainer was abandoned.
    pub complete: bool,
    pub error: Option<String>,
}

/// Copies `source` into a new file at `path`, keeping the first `cap` bytes
/// and the last `tail` bytes. Memory use is bounded by `tail` plus one chunk.
pub fn drain(mut source: impl Read, path: &Path, cap: u64, tail: usize) -> LogRecord {
    let mut record = LogRecord {
        path: path.to_string_lossy().into_owned(),
        bytes_seen: 0,
        bytes_written: 0,
        truncated: false,
        complete: true,
        error: None,
    };
    let mut file = match File::create_new(path) {
        Ok(file) => Some(file),
        Err(error) => {
            record.error = Some(format!("cannot create log: {error}"));
            None
        }
    };
    let mut ring: VecDeque<u8> = VecDeque::with_capacity(tail.min(CHUNK));
    let mut buffer = vec![0; CHUNK];
    loop {
        let read = match source.read(&mut buffer) {
            Ok(0) => break,
            Ok(read) => read,
            Err(error) if error.kind() == io::ErrorKind::Interrupted => continue,
            Err(error) => {
                record
                    .error
                    .get_or_insert(format!("cannot read stream: {error}"));
                break;
            }
        };
        record.bytes_seen += read as u64;
        let chunk = &buffer[..read];
        let room = cap.saturating_sub(record.bytes_written).min(read as u64) as usize;
        if room > 0 && write_or_disable(&mut file, &chunk[..room], &mut record) {
            record.bytes_written += room as u64;
        }
        if room < read {
            record.truncated = true;
            ring.extend(&chunk[room..]);
            let excess = ring.len().saturating_sub(tail);
            ring.drain(..excess);
        }
    }
    if record.truncated {
        let kept = ring.len() as u64;
        let omitted = record.bytes_seen - record.bytes_written - kept;
        let note = format!(
            "\n[meta-co-assembler: log capped at {cap} bytes; omitted {omitted} bytes; last {kept} bytes follow]\n"
        );
        write_or_disable(&mut file, note.as_bytes(), &mut record);
        let (front, back) = ring.as_slices();
        write_or_disable(&mut file, front, &mut record);
        write_or_disable(&mut file, back, &mut record);
    }
    if let Some(file) = file.as_mut() {
        if let Err(error) = file.flush() {
            record
                .error
                .get_or_insert(format!("cannot flush log: {error}"));
        }
    }
    record
}

/// Returns whether the bytes were written. After a write error the file is
/// dropped, but the caller keeps draining so the child never blocks.
fn write_or_disable(file: &mut Option<File>, bytes: &[u8], record: &mut LogRecord) -> bool {
    let Some(handle) = file.as_mut() else {
        return false;
    };
    match handle.write_all(bytes) {
        Ok(()) => true,
        Err(error) => {
            record
                .error
                .get_or_insert(format!("cannot write log: {error}"));
            *file = None;
            false
        }
    }
}

/// A drainer running on its own thread.
pub struct Drainer {
    path: PathBuf,
    receiver: mpsc::Receiver<LogRecord>,
}

impl Drainer {
    pub fn spawn(source: impl Read + Send + 'static, path: PathBuf, cap: u64, tail: usize) -> Self {
        let (sender, receiver) = mpsc::channel();
        let target = path.clone();
        thread::spawn(move || {
            let _ = sender.send(drain(source, &target, cap, tail));
        });
        Self { path, receiver }
    }

    /// Waits for EOF. If a surviving descendant holds the stream open, the
    /// drainer is abandoned after `wait` and the record is marked incomplete.
    pub fn finish(self, wait: Duration) -> LogRecord {
        self.receiver
            .recv_timeout(wait)
            .unwrap_or_else(|_| LogRecord {
                path: self.path.to_string_lossy().into_owned(),
                bytes_seen: 0,
                bytes_written: 0,
                truncated: false,
                complete: false,
                error: Some("stream still open after the process ended; capture abandoned".into()),
            })
    }
}

#[cfg(test)]
mod tests {
    use std::{fs, io::Cursor};

    use super::*;

    fn temp_path(name: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("mca-logs-{}-{name}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir(&dir).unwrap();
        dir.join("stream.log")
    }

    #[test]
    fn short_streams_are_copied_exactly() {
        let path = temp_path("short");
        let record = drain(Cursor::new(b"hello\n".to_vec()), &path, 100, 10);
        assert_eq!(fs::read(&path).unwrap(), b"hello\n");
        assert_eq!((record.bytes_seen, record.bytes_written), (6, 6));
        assert!(!record.truncated && record.complete && record.error.is_none());
        fs::remove_dir_all(path.parent().unwrap()).unwrap();
    }

    #[test]
    fn long_streams_keep_head_and_tail_within_bounds() {
        let path = temp_path("long");
        let data: Vec<u8> = (0..300_000u32).map(|i| (i % 251) as u8).collect();
        let record = drain(Cursor::new(data.clone()), &path, 1000, 500);
        let log = fs::read(&path).unwrap();
        assert!(record.truncated);
        assert_eq!(record.bytes_seen, 300_000);
        assert_eq!(record.bytes_written, 1000);
        assert!(log.starts_with(&data[..1000]));
        assert!(log.ends_with(&data[data.len() - 500..]));
        assert!(log.len() < 1000 + 500 + 200);
        let text = String::from_utf8_lossy(&log);
        assert!(text.contains(&format!("omitted {} bytes", 300_000 - 1000 - 500)));
        fs::remove_dir_all(path.parent().unwrap()).unwrap();
    }

    #[test]
    fn an_existing_log_is_never_overwritten_but_the_stream_is_still_drained() {
        let path = temp_path("exists");
        fs::write(&path, "previous").unwrap();
        let record = drain(Cursor::new(vec![7; 10_000]), &path, 100, 10);
        assert_eq!(record.bytes_seen, 10_000);
        assert_eq!(record.bytes_written, 0);
        assert!(record.error.unwrap().contains("cannot create log"));
        assert_eq!(fs::read(&path).unwrap(), b"previous");
        fs::remove_dir_all(path.parent().unwrap()).unwrap();
    }

    #[test]
    fn a_zero_cap_keeps_only_the_tail() {
        let path = temp_path("zero");
        let record = drain(Cursor::new(b"abcdefghij".to_vec()), &path, 0, 3);
        assert_eq!(record.bytes_written, 0);
        assert!(fs::read_to_string(&path).unwrap().ends_with("hij"));
        fs::remove_dir_all(path.parent().unwrap()).unwrap();
    }
}
