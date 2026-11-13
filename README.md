# meta-co-assembler

An exploratory Rust project for coassembling reference-unassigned metagenomic
read pairs across related samples. The first research endpoint is contigs.
MAG recovery and evaluation follow only if the assembly approach proves useful.

AI agents help develop the software. The pipeline itself will run without an
LLM, on a single node, with explicit CPU and memory budgets.

## Current state

The CLI validates sample and metadata tables and prints deterministic metadata
grouping plans. Profiling, mapping, read extraction, assembly, and iteration are
not implemented. No external biological tool or reference catalogue has been
selected. Select these with the project owner before implementing adapters.

## Build

Install Rust 1.85 or newer, then:

```sh
cargo build --locked --release
./target/release/meta-co-assembler --help
```

On Windows, the executable is `target/release/meta-co-assembler.exe`.
The Rust CLI targets Linux, macOS, and Windows; initial biological execution
will be validated on Linux.

## Try the input contracts

Provide already processed paired-end Illumina reads in FASTQ or gzip FASTQ.
Copy [the sample table](examples/samples.tsv) and
[the metadata table](examples/metadata.tsv) and replace the placeholder paths.
Relative read paths resolve against the sample table's directory.

```sh
meta-co-assembler validate --samples samples.tsv --metadata metadata.tsv

meta-co-assembler plan \
  --samples samples.tsv \
  --metadata metadata.tsv \
  --group-by habitat \
  --memory-gib 64 \
  --threads 8 \
  --jobs 2 \
  --max-group-samples 16 \
  --iteration-strategy residual-only \
  --max-rounds 1 > plan.json
```

Use `--iteration-strategy recruitment` to describe the alternative iterative
strategy. Neither strategy executes yet. `plan` reads tables and file metadata,
not FASTQ contents; it prints JSON to stdout and writes no project files.

Unassigned means not confidently explained by the selected references under an
explicit mapping policy. It does not mean novel. Both mates will be retained
when either mate is unassigned or ambiguous. Samples remain traceable through
grouping and assembly.

See the [CLI contract](docs/cli-contract.md) for input formats, defaults, exit
codes, and future execution boundaries.

## Development checks

```sh
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo test --locked --all-targets
```

The tests use tiny, temporary fixtures and require no reference database,
biological tool, or real sequencing data. They exercise CLI contracts,
resource validation, sample/metadata integrity, deterministic grouping,
iteration selection, and the unavailable execution boundary.

[GitHub Actions](.github/workflows/ci.yml) builds release executables and runs
tests on Linux, macOS, and Windows. Separate jobs test Rust 1.85 and check
formatting/Clippy. Executables are retained as build artifacts for seven days.
The workflow runs on pushes, pull requests, and manual dispatch.

Start future development with the [milestone roadmap](ROADMAP.md) and
[intended architecture](docs/architecture.md). Development agents enter through
[AGENTS.md](AGENTS.md) and the shared [agent guide](agents/AGENT.md).
The owner creates each feature branch manually; later milestones are
conditional on the contig pilot demonstrating value.
