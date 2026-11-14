use std::{
    io::{self, Write},
    process::ExitCode,
};

use clap::Parser;
use meta_co_assembler::{
    cli::{Cli, Command},
    error::Error,
    execution::{CancelToken, RunStatus},
    manifest::{read_metadata, read_samples, validation_report},
    plan::make_plan,
    profiling::run_profile,
};
use serde::Serialize;

fn write_json(value: &impl Serialize) -> Result<(), Error> {
    let mut stdout = io::stdout().lock();
    serde_json::to_writer_pretty(&mut stdout, value).map_err(Error::Json)?;
    writeln!(stdout).map_err(|source| Error::Io {
        context: "cannot write stdout".into(),
        source,
    })
}

fn execute(cli: Cli) -> Result<(), Error> {
    match cli.command {
        Command::Validate(args) => {
            let samples = read_samples(&args.samples)?;
            let metadata = args
                .metadata
                .as_deref()
                .map(|path| read_metadata(path, &samples))
                .transpose()?;
            write_json(&validation_report(&samples, metadata.as_ref()))
        }
        Command::Plan(args) => write_json(&make_plan(&args)?),
        Command::Profile(args) => {
            let cancel = CancelToken::new();
            cancel.cancel_on_termination_signals()?;
            let outcome = run_profile(&args, &cancel)?;
            write_json(&outcome.summary)?;
            let manifest = outcome.summary["manifest"].as_str().unwrap_or_default();
            let reason = outcome.report.stop_reason.unwrap_or_default();
            match outcome.report.status {
                RunStatus::Interrupted => Err(Error::Interrupted(format!(
                    "profiling was interrupted; partial work is recorded in {manifest}"
                ))),
                RunStatus::Succeeded => Ok(()),
                _ => Err(Error::RunFailed(format!(
                    "profiling failed: {reason}; see {manifest}"
                ))),
            }
        }
        Command::Run(args) => {
            args.plan.validate()?;
            Err(Error::Unavailable(
                "pipeline execution is unavailable: execution adapters are not implemented; selected tools are recorded for planning only; no outputs were created"
                    .into(),
            ))
        }
    }
}

fn main() -> ExitCode {
    match execute(Cli::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("error: {error}");
            ExitCode::from(error.exit_code())
        }
    }
}
