# Agent entry point

Read [agents/AGENT.md](agents/AGENT.md) before changing this repository. It is
the shared development guide. Consult the relevant role guide when needed:

- [Architecture](agents/software_architect.md)
- [Implementation](agents/implementation.md)
- [Testing](agents/testing.md)
- [Research](agents/research.md)

Use the Conda environment `meta-co-assembler-env` for local builds, tests,
and research commands. See [environment.yml](environment.yml) and the
[setup instructions](README.md#environment-and-build). For sessions without
persistent activation, prefix commands with
`conda run --no-capture-output --name meta-co-assembler-env`.

These are role guidelines, not autonomous runtime agents. A single development
agent may use all of them; a separate architect agent is optional.

Follow the [CLI contract](docs/cli-contract.md),
[architecture](docs/architecture.md), and [roadmap](ROADMAP.md). The project
owner's current instructions take precedence over repository guidance.

The [M1 decision](docs/decisions/0001-pilot.md) records approved pilot tools,
selectable alternatives, and the small RefSeq catalogue. Reuse those choices.
Ask the owner about new tools, production catalogues, or unselected scientific
policies. LoCovE guidance is experimental; preserve the confidence requirements
in [the interface contract](docs/pilot-interfaces.md).

The initial bootstrap was requested on main. Subsequent feature branches are
created manually by the owner; work on the provided branch and requested
milestone. Do not automatically start the next milestone.
