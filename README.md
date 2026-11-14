# meta-co-assembler

An exploratory Rust project for coassembling reference-unassigned metagenomic
read pairs across related samples. The first research endpoint is contigs.
MAG recovery and evaluation follow only if the assembly approach proves useful.

AI agents help develop the software. The pipeline itself will run without an
LLM, on a single node, with explicit CPU and memory budgets.

## Current state

The CLI validates sample and metadata tables and prints deterministic metadata
grouping plans. Profiling, mapping, read extraction, assembly, and iteration are
not implemented. M1 defines a controlled synthetic pilot and selectable tools:
Sylph/Kraken2/experimental LoCovE, Strobealign/Minimap2, and
MEGAHIT/metaSPAdes, with Samtools for alignment I/O. The initial pilot uses
Sylph, Strobealign, and MEGAHIT. See the [pilot protocol](experiments/pilot/README.md)
and [LoCovE investigation](docs/research/locove.md).

## Environment and build

Create the project environment from the repository root, then activate it:

```sh
conda env create --file environment.yml
conda activate meta-co-assembler-env
```

If the environment already exists, synchronize it with the definition:

```sh
conda env update --name meta-co-assembler-env --file environment.yml --prune
```

The [environment definition](environment.yml) pins Rust 1.99.1, including Cargo,
rustfmt, and Clippy, the native C compiler/linker, Python, and the released
pilot tools. Conda manages these tools; `Cargo.toml` and `Cargo.lock` manage
Rust library dependencies.
The compiler metapackage's version is pinned; platform-specific transitive
packages are resolved by Conda, so this file is not an exact platform lockfile.
LoCovE is source-built into this environment using its
[pinned revision and explicit setup](experiments/pilot/README.md#locove-setup-and-investigation).
New tool choices still require the project owner's input.

Run builds and development checks inside `meta-co-assembler-env`:

```sh
cargo build --locked --release
./target/release/meta-co-assembler --help
```

On Windows, the executable is `target/release/meta-co-assembler.exe`.
The Rust CLI targets Linux, macOS, and Windows; initial biological execution
will be validated on Linux. The Conda environment is initially verified on
Linux. Native-toolchain CI also tests the minimum supported Rust version, 1.85.

For scripts and agent sessions that do not retain shell activation, use
`conda run`:

```sh
conda run --no-capture-output --name meta-co-assembler-env cargo build --locked --release
conda run --no-capture-output --name meta-co-assembler-env cargo test --locked --all-targets
```

Confirm that `cargo --version` and `rustc --version` report the pinned Rust
version before recording results. A rustup or system toolchain that appears
earlier on `PATH`, for example from `~/.cargo/env`, can shadow the environment's
toolchain even under `conda run`. CI fails if the resolved versions differ.

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

Tool choices are independent plan options: `--profiler sylph|kraken2|locove`,
`--mapper strobealign|minimap2`, and `--assembler megahit|metaspades`.
Use `--mapping-mode locove-guided` to record the experimental guidance design.
Every plan remains non-executable; selecting a backend does not run it.

Unassigned means not confidently explained by the selected references under an
explicit mapping policy. It does not mean novel. Both mates will be retained
when either mate is unassigned or ambiguous. Samples remain traceable through
grouping and assembly.

See the [CLI contract](docs/cli-contract.md) for input formats, defaults, exit
codes, and future execution boundaries.

## Development checks

Activate `meta-co-assembler-env` first, or prefix each command with
`conda run --no-capture-output --name meta-co-assembler-env`.

```sh
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo test --locked --all-targets
python -m unittest discover -s tests/pilot -v
```

The tests use tiny, temporary fixtures and require no reference database,
biological tool, or real sequencing data. They exercise CLI contracts,
resource validation, sample/metadata integrity, deterministic grouping,
iteration selection, and the unavailable execution boundary.

[GitHub Actions](.github/workflows/ci.yml) builds release executables and runs
tests on Linux, macOS, and Windows. Separate jobs test Rust 1.85 and check
formatting/Clippy. A Linux job creates `meta-co-assembler-env` from
`environment.yml` and runs the development checks and a release build inside it.
Executables are retained as build artifacts for seven days.
The workflow runs on pushes, pull requests, and manual dispatch.

Start future development with the [milestone roadmap](ROADMAP.md) and
[intended architecture](docs/architecture.md). Development agents enter through
[AGENTS.md](AGENTS.md) and the shared [agent guide](agents/AGENT.md).
The owner creates each feature branch manually; later milestones are
conditional on the contig pilot demonstrating value.
