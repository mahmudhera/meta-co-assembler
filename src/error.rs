use std::{fmt, io};

#[derive(Debug)]
pub enum Error {
    InvalidInput(String),
    Io { context: String, source: io::Error },
    Csv { context: String, source: csv::Error },
    Json(serde_json::Error),
    Unavailable(&'static str),
}

impl Error {
    pub fn exit_code(&self) -> u8 {
        match self {
            Self::InvalidInput(_) => 2,
            Self::Csv { source, .. } if !source.is_io_error() => 2,
            Self::Unavailable(_) => 3,
            _ => 1,
        }
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidInput(message) => f.write_str(message),
            Self::Unavailable(message) => f.write_str(message),
            Self::Io { context, source } => write!(f, "{context}: {source}"),
            Self::Csv { context, source } => write!(f, "{context}: {source}"),
            Self::Json(source) => write!(f, "cannot write JSON: {source}"),
        }
    }
}

impl std::error::Error for Error {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Io { source, .. } => Some(source),
            Self::Csv { source, .. } => Some(source),
            Self::Json(source) => Some(source),
            _ => None,
        }
    }
}
