# Research role

Keep the project focused on whether reference subtraction and sample grouping
improve discovery at acceptable resource cost.

Start with the approved `experiments/pilot/` protocol and settings. Separate
LoCovE coverage/region measurements from claims of mapping acceleration.
Keep full-reference controls and measure false confident assignment caused
by removing competing reference windows. Unknown regions remain uncertain.

## Scientific boundaries

- A profiler's detected taxa are candidates for reference selection, not proof
  that all corresponding reads should be removed.
- Unmapped/ambiguous reads are not automatically novel, microbial, or useful.
- Unassigned to selected references is different from absent from the full
  frozen catalogue. Use separate terms and measurements.
- Recovered contigs are not automatically MAGs. Binning, quality assessment,
  dereplication, and catalogue comparison are later work.
- Preserve both mates when either is unassigned. Evaluate any loss of shared
  conserved sequence caused by reference subtraction.

## Research protocol

1. State the hypothesis and the smallest experiment that could disprove it.
2. Define datasets, reference holdouts, comparators, resource limits, metrics,
   and success/failure criteria with the owner before expensive experiments.
3. Compare tools and database choices using primary papers and official
   documentation, including memory, throughput, pairing, and licensing needs.
4. Ask the owner to choose unselected tools/algorithms/thresholds. Record the
   decision and experiment configuration so later agents reuse it.
5. Run reproducible experiments with fixed data versions, parameters, and seeds.
6. Report negative results and failure cases as clearly as successful ones.

Run experiments in `meta-co-assembler-env`. Add confirmed external tools to
`environment.yml` and record the actual resolved package/tool versions with
each run; the portable environment specification is not a full platform lockfile.

Test single-sample/all-read and grouped/all-read assembly baselines where
feasible, using comparable resource budgets. Separate the effect of subtraction
from grouping and iteration. Include rare populations, related strains, and
references deliberately omitted from the catalogue; holdouts simulate catalogue
absence and do not demonstrate real-world taxonomic novelty.

Prefer a decision to continue, revise, or stop over premature polishing.
If the pilot fails, summarize what failed before adding more pipeline features.
