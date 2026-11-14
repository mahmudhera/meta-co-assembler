//! SHA-256 over streams that are read or written anyway, so large files are
//! fingerprinted without an extra pass.

use std::io::{self, Read, Write};

use sha2::{Digest, Sha256};

pub fn hex(bytes: &[u8]) -> String {
    bytes.iter().map(|byte| format!("{byte:02x}")).collect()
}

/// Hashes every byte written through it.
pub struct HashingWriter<W> {
    inner: W,
    hasher: Sha256,
}

impl<W> HashingWriter<W> {
    pub fn new(inner: W) -> Self {
        Self {
            inner,
            hasher: Sha256::new(),
        }
    }

    /// The writer and the lowercase hex digest of everything written.
    pub fn finish(self) -> (W, String) {
        (self.inner, hex(&self.hasher.finalize()))
    }
}

impl<W: Write> Write for HashingWriter<W> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let written = self.inner.write(buf)?;
        self.hasher.update(&buf[..written]);
        Ok(written)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.inner.flush()
    }
}

/// Hashes every byte read through it.
pub struct HashingReader<R> {
    inner: R,
    hasher: Sha256,
}

impl<R> HashingReader<R> {
    pub fn new(inner: R) -> Self {
        Self {
            inner,
            hasher: Sha256::new(),
        }
    }

    pub fn digest(&self) -> String {
        hex(&self.hasher.clone().finalize())
    }
}

impl<R: Read> Read for HashingReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let read = self.inner.read(buf)?;
        self.hasher.update(&buf[..read]);
        Ok(read)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn readers_and_writers_match_the_standard_vector() {
        let expected = "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad";
        let mut writer = HashingWriter::new(Vec::new());
        writer.write_all(b"ab").unwrap();
        writer.write_all(b"c").unwrap();
        let (bytes, digest) = writer.finish();
        assert_eq!((bytes.as_slice(), digest.as_str()), (&b"abc"[..], expected));
        let mut reader = HashingReader::new(&b"abc"[..]);
        let mut sink = Vec::new();
        reader.read_to_end(&mut sink).unwrap();
        assert_eq!(reader.digest(), expected);
    }
}
