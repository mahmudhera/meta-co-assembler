use std::{
    io::{self, Write},
    process::ExitCode,
};

use clap::Parser;
use meta_co_assembler::{
    assembly::{run_assemble_reads, run_coassemble},
    cli::{Cli, Command},
    discovery::{run_discover, run_extract_pairs},
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

/// Prints an executed run's summary and maps its status to an exit category.
fn finish_run(
    name: &str,
    summary: &serde_json::Value,
    status: RunStatus,
    stop_reason: Option<String>,
) -> Result<(), Error> {
    write_json(summary)?;
    let manifest = summary["manifest"].as_str().unwrap_or_default();
    match status {
        RunStatus::Succeeded => Ok(()),
        RunStatus::Interrupted => Err(Error::Interrupted(format!(
            "{name} was interrupted; partial work is recorded in {manifest}"
        ))),
        _ => Err(Error::RunFailed(format!(
            "{name} failed: {}; see {manifest}",
            stop_reason.unwrap_or_default()
        ))),
    }
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
            finish_run(
                "profiling",
                &outcome.summary,
                outcome.report.status,
                outcome.report.stop_reason,
            )
        }
        Command::Discover(args) => {
            let cancel = CancelToken::new();
            cancel.cancel_on_termination_signals()?;
            let outcome = run_discover(&args, &cancel)?;
            finish_run(
                "discovery",
                &outcome.summary,
                outcome.report.status,
                outcome.report.stop_reason,
            )
        }
        Command::ExtractPairs(args) => write_json(&run_extract_pairs(&args)?),
        Command::Coassemble(args) => {
            let cancel = CancelToken::new();
            cancel.cancel_on_termination_signals()?;
            let outcome = run_coassemble(&args, &cancel)?;
            finish_run(
                "coassembly",
                &outcome.summary,
                outcome.report.status,
                outcome.report.stop_reason,
            )
        }
        Command::AssembleReads(args) => {
            let cancel = CancelToken::new();
            cancel.cancel_on_termination_signals()?;
            let outcome = run_assemble_reads(&args, &cancel)?;
            finish_run(
                "assembly",
                &outcome.summary,
                outcome.report.status,
                outcome.report.stop_reason,
            )
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
