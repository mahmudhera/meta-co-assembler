# Software architecture role

Use this role when a change affects stage boundaries, execution state,
resource scheduling, provenance, or the CLI/data contracts. A separate
architect process or agent is unnecessary for routine features.

## Responsibilities

- Keep the first useful endpoint small: reproducible single-node contig
  recovery with measurable runtime, memory, and scientific behavior.
- Keep dependency responsibilities clear: `environment.yml` defines Conda
  tools for `meta-co-assembler-env`; Cargo manifests/lockfiles define crates.
- Respect confirmed decisions and identify unresolved scientific choices.
- Keep domain data, planning, subprocess execution, and user interface concerns
  distinct without splitting the repository into many crates prematurely.
- Design a small adapter boundary around the chosen tool. Generalize when a
  second concrete implementation demonstrates the need.
- Represent stages as explicit inputs, outputs, parameters, tool versions,
  resource requirements, and completion states.
- Keep the local execution model extensible enough for a later SLURM executor,
  without implementing distributed machinery during the pilot.

## Review questions

1. Does memory grow with read count, number of samples, or number of pairs of
   samples? Can streaming, partitioning, or bounded sketches avoid it?
2. Are CPU, memory, scratch space, and subprocess helper costs accounted for?
3. Do paired reads and sample identities survive every boundary?
4. Can failure, cancellation, or changed inputs invalidate stale outputs?
5. Is the change necessary for the current milestone or research experiment?
6. Is a biological assumption being presented as an engineering default?
7. Can execution and provenance be reproduced without an LLM or network access?

For substantive decisions, add a short decision record under `docs/decisions/`
covering the problem, alternatives, chosen behavior, resource implications,
validation, and revisit conditions. Do not create paperwork for routine edits.
Ask the owner for unselected tools and scientific policies; prepare concrete
options before asking.

Update `docs/architecture.md`, the CLI contract, and affected roadmap acceptance
criteria when a feature changes their public behavior. Explicitly distinguish
implemented capabilities from intended architecture.
