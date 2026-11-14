#!/bin/sh
# Evaluate the M8 hard-benchmark arms with the M7 evaluator: MetaQUAST per community with the gate
# against grouped all-read assembly, and the global lens (LENSES selects; default both).
# Community grouping and pooling are the M7 runs.
set -u
cd "$(dirname "$0")/../.."
M=work/m8/hard
M7=work/m7/hard
LENSES="${LENSES:-community global}"
case " $LENSES " in *" community "*)
python scripts/evaluate_assemblies.py \
  --arm sketch-merge=$M/sketch-merge \
  --arm random-1=$M/random-1 --arm random-2=$M/random-2 --arm random-3=$M/random-3 \
  --arm community=$M7/discovery-grouped \
  --arm pooled=$M7/discovery-pooled \
  --arm allreads-grouped=$M7/allreads-grouped \
  --discovery-arm sketch-merge --discovery-arm random-1 --discovery-arm random-2 \
  --discovery-arm random-3 --discovery-arm community --discovery-arm pooled \
  --baseline-arm allreads-grouped \
  --metadata data/m7-hard/metadata.tsv --group-by community \
  --composition experiments/m7-hard/composition.tsv --references data/m7-hard/references \
  --catalogue data/m7-hard/catalogues/hard_holdout/catalogue.tsv \
  --genomes experiments/m7-hard/genomes.tsv --workdir $M/quast --threads 16 --output $M/evaluation.json
echo "== evaluate-assemblies hard exit=$?" >> work/m8/log.txt
;; esac
# Global lens: each arm's whole contig output against every source genome (one community "all").
# The per-community lens credits a community with every group that holds any of its samples, so
# groupings that spread a community over several groups get several assemblies of a genome shared
# by all communities; here every arm is credited with all of its groups.
case " $LENSES " in *" global "*)
awk -F'\t' 'BEGIN{OFS="\t"} NR==1{print "sample_id","scope"; next} {print $1,"all"}' \
  data/m7-hard/metadata.tsv > $M/metadata-all.tsv
python scripts/evaluate_assemblies.py \
  --arm sketch-merge=$M/sketch-merge \
  --arm random-1=$M/random-1 --arm random-2=$M/random-2 --arm random-3=$M/random-3 \
  --arm community=$M7/discovery-grouped \
  --arm pooled=$M7/discovery-pooled \
  --arm allreads-grouped=$M7/allreads-grouped \
  --discovery-arm sketch-merge --discovery-arm random-1 --discovery-arm random-2 \
  --discovery-arm random-3 --discovery-arm community --discovery-arm pooled \
  --baseline-arm allreads-grouped \
  --metadata $M/metadata-all.tsv --group-by scope \
  --composition experiments/m7-hard/composition.tsv --references data/m7-hard/references \
  --catalogue data/m7-hard/catalogues/hard_holdout/catalogue.tsv \
  --genomes experiments/m7-hard/genomes.tsv --workdir $M/quast-global --threads 16 \
  --output $M/evaluation-global.json
echo "== evaluate-global hard exit=$?" >> work/m8/log.txt
;; esac
