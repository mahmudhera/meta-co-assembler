#!/bin/sh
# M7 factorial: input {all reads, discovery} x grouping {per sample, metadata, pooled},
# plus pooled iteration (residual-only, recruitment) on grouped discovery.
# Datasets: pilot (species holdout), split-genome scenario, hard benchmark.
# Every command at 16 threads / 64 GiB, serially, under /usr/bin/time -v. BIN defaults to the
# release build (M7 used commit 67d293c). Outputs go to work/m7 (gitignored).
set -u
cd "$(dirname "$0")/../.."
M=work/m7
B="${BIN:-target/release/meta-co-assembler}"
R="--threads 16 --memory-gib 64"
step() { d=$1; name=$2; shift 2; start=$(date +%s); echo "== $d/$name started $(date -Is)" >> $M/log.txt
  /usr/bin/time -v -o $M/$d/$name.time "$@" > $M/$d/$name.json 2> $M/$d/$name.err; code=$?
  echo "== $d/$name exit=$code seconds=$(( $(date +%s) - start ))" >> $M/log.txt; }
run_dataset() { d=$1 samples=$2 metadata=$3 catalogue=$4; mkdir -p $M/$d
  step $d discover $B discover --samples $samples --catalogue $catalogue --output-dir $M/$d/discover $R --jobs 4
  step $d discovery-grouped $B coassemble --discovery-dir $M/$d/discover --metadata $metadata --group-by community --output-dir $M/$d/discovery-grouped $R
  step $d discovery-per-sample $B coassemble --discovery-dir $M/$d/discover --per-sample --output-dir $M/$d/discovery-per-sample $R
  step $d discovery-pooled $B coassemble --discovery-dir $M/$d/discover --all-samples --output-dir $M/$d/discovery-pooled $R
  for strategy in residual-only recruitment; do
    step $d iterate-$strategy $B iterate --coassembly-dir $M/$d/discovery-grouped --iteration-strategy $strategy --all-samples --max-rounds 3 --output-dir $M/$d/iterate-$strategy $R
  done
  step $d allreads-grouped $B assemble-reads --samples $samples --metadata $metadata --group-by community --output-dir $M/$d/allreads-grouped $R
  step $d allreads-pooled $B assemble-reads --samples $samples --all-samples --output-dir $M/$d/allreads-pooled $R
  step $d allreads-per-sample $B assemble-reads --samples $samples --per-sample --output-dir $M/$d/allreads-per-sample $R
}
for d in ${DATASETS:-hard pilot split}; do
  case $d in
    pilot) run_dataset pilot data/pilot/samples.tsv data/pilot/metadata.tsv data/pilot/catalogues/species_holdout/catalogue.tsv ;;
    split) run_dataset split data/pilot-split/samples.tsv data/pilot-split/metadata.tsv data/pilot-split/catalogues/split_holdout/catalogue.tsv ;;
    hard) run_dataset hard data/m7-hard/samples.tsv data/m7-hard/metadata.tsv data/m7-hard/catalogues/hard_holdout/catalogue.tsv ;;
  esac
done
echo "== all done $(date -Is)" >> $M/log.txt
