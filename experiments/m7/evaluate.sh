#!/bin/sh
# Evaluate one M7 dataset ($1 = pilot | split | hard): MetaQUAST per community over the eight
# arms plus iterate round 1 (grouped all-read assembly is the gate baseline; iterate arms are
# shared sets, the evaluator default since M7), and discovery-pool composition.
set -u
cd "$(dirname "$0")/../.."
M=work/m7
d=$1
case $d in
  pilot) meta=data/pilot/metadata.tsv comp=experiments/pilot/composition.tsv refs=data/pilot/references
         cat=data/pilot/catalogues/species_holdout/catalogue.tsv genomes=experiments/pilot/genomes.tsv ;;
  split) meta=data/pilot-split/metadata.tsv comp=data/pilot-split/composition.tsv refs=data/pilot/references
         cat=data/pilot-split/catalogues/split_holdout/catalogue.tsv genomes=experiments/pilot/genomes.tsv ;;
  hard)  meta=data/m7-hard/metadata.tsv comp=experiments/m7-hard/composition.tsv refs=data/m7-hard/references
         cat=data/m7-hard/catalogues/hard_holdout/catalogue.tsv genomes=experiments/m7-hard/genomes.tsv ;;
esac
python scripts/evaluate_discovery.py $M/$d/discover --genomes $genomes --output $M/$d/discovery-eval.json > /dev/null
echo "== evaluate-discovery $d exit=$?" >> $M/log.txt
python scripts/evaluate_assemblies.py \
  --arm discovery-grouped=$M/$d/discovery-grouped \
  --arm discovery-per-sample=$M/$d/discovery-per-sample \
  --arm discovery-pooled=$M/$d/discovery-pooled \
  --arm iterate-round1=$M/$d/iterate-residual-only@1 \
  --arm iterate-residual-only=$M/$d/iterate-residual-only \
  --arm iterate-recruitment=$M/$d/iterate-recruitment \
  --arm allreads-per-sample=$M/$d/allreads-per-sample \
  --arm allreads-grouped=$M/$d/allreads-grouped \
  --arm allreads-pooled=$M/$d/allreads-pooled \
  --discovery-arm discovery-grouped --discovery-arm discovery-per-sample --discovery-arm discovery-pooled \
  --discovery-arm iterate-round1 --discovery-arm iterate-residual-only --discovery-arm iterate-recruitment \
  --baseline-arm allreads-grouped \
  --metadata $meta --group-by community --composition $comp --references $refs \
  --catalogue $cat --genomes $genomes --workdir $M/$d/quast --threads 16 --output $M/$d/evaluation.json
echo "== evaluate-assemblies $d exit=$?" >> $M/log.txt
