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

External biological tools and the reference catalogue are deliberately
unselected. Ask the owner before selecting or implementing them; reuse choices
already confirmed in the session or recorded decisions.

The initial bootstrap was requested on main. Subsequent feature branches are
created manually by the owner; work on the provided branch and requested
milestone. Do not automatically start the next milestone.
