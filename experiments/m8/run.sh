#!/bin/sh
# M8: sketch grouping against metadata and random grouping (design agreed with the owner).
# Part 1, grouping only (`group`, no assembly) on the three M7 datasets, at a 4-sample cap:
#   discovery pools and original reads with the merge and seeded policies, threshold
#   sensitivity on the pools (PARTS grouping), and the memory cap (--sketch-max-hashes) on the
#   hard reads (PARTS memory).
# Part 2, hard benchmark assemblies at a 4-sample cap: sketch grouping and three seeded random
#   groupings (random-groupings.tsv) of the M7 discovery pools. Community grouping and the
#   pooled assembly reuse the M7 runs.
# Needs the M7 discover runs (experiments/m7/run.sh). 16 threads / 64 GiB, serially, under
# /usr/bin/time -v. BIN defaults to the release build (M8 used the feature/sketch-grouping
# implementation, built before it was committed, on 0a744e9). Outputs go to work/m8 (gitignored).
set -u
cd "$(dirname "$0")/../.."
M=work/m8
B="${BIN:-target/release/meta-co-assembler}"
R="--threads 16 --memory-gib 64"
step() { d=$1; name=$2; shift 2; start=$(date +%s); echo "== $d/$name started $(date -Is)" >> $M/log.txt
  /usr/bin/time -v -o $M/$d/$name.time "$@" > $M/$d/$name.json 2> $M/$d/$name.err; code=$?
  echo "== $d/$name exit=$code seconds=$(( $(date +%s) - start ))" >> $M/log.txt; }
samples_of() { case $1 in pilot) echo data/pilot/samples.tsv ;; split) echo data/pilot-split/samples.tsv ;;
  hard) echo data/m7-hard/samples.tsv ;; esac; }

PARTS="${PARTS:-grouping memory assemblies}"
has() { case " $PARTS " in *" $1 "*) return 0 ;; *) return 1 ;; esac; }
has grouping && for d in pilot split hard; do
  mkdir -p $M/$d
  pools="--discovery-dir work/m7/$d/discover"
  reads="--samples $(samples_of $d)"
  for policy in merge seeded; do
    step $d group-pools-$policy $B group $pools --sketch --sketch-policy $policy --max-group-samples 4 $R
    step $d group-reads-$policy $B group $reads --sketch --sketch-policy $policy --max-group-samples 4 $R
  done
  for t in 0.05 0.1 0.5 0.8; do
    step $d group-pools-threshold-$t $B group $pools --sketch --sketch-policy threshold \
      --sketch-min-similarity $t --max-group-samples 4 $R
  done
done
mkdir -p $M/hard
# Memory cap: at the default scale, and at scale 1 (every k-mer kept), where only the cap bounds
# the sketch. Runs after the assemblies (PARTS=memory) so that it does not share the CPUs with them.
has memory && for cap in 1000 10000 100000; do
  step hard group-reads-merge-cap$cap $B group --samples data/m7-hard/samples.tsv --sketch \
    --sketch-max-hashes $cap --max-group-samples 4 $R
done
has memory && for cap in 100000 1000000 10000000; do
  step hard group-reads-merge-scaled1-cap$cap $B group --samples data/m7-hard/samples.tsv --sketch \
    --sketch-scaled 1 --sketch-max-hashes $cap --max-group-samples 4 $R
done
if has assemblies; then
  d=hard
  step $d sketch-merge $B coassemble --discovery-dir work/m7/hard/discover --sketch \
    --max-group-samples 4 --output-dir $M/$d/sketch-merge $R
  for s in 1 2 3; do
    step $d random-$s $B coassemble --discovery-dir work/m7/hard/discover \
      --metadata experiments/m8/random-groupings.tsv --group-by random_$s --max-group-samples 4 \
      --output-dir $M/$d/random-$s $R
  done
fi
echo "== all done $(date -Is)" >> $M/log.txt
