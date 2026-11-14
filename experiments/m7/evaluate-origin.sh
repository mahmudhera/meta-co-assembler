#!/bin/sh
# Sensitivity evaluation ($1 = pilot | split | hard): iterate contigs attributed to the group that
# assembled them (--origin-arm, the M6 lens), against the same grouped all-read baseline.
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
python scripts/evaluate_assemblies.py \
  --arm iterate-round1=$M/$d/iterate-residual-only@1 \
  --arm iterate-residual-only=$M/$d/iterate-residual-only \
  --arm iterate-recruitment=$M/$d/iterate-recruitment \
  --arm discovery-pooled=$M/$d/discovery-pooled \
  --arm allreads-grouped=$M/$d/allreads-grouped \
  --origin-arm iterate-round1 --origin-arm iterate-residual-only --origin-arm iterate-recruitment \
  --discovery-arm iterate-round1 --discovery-arm iterate-residual-only --discovery-arm iterate-recruitment \
  --discovery-arm discovery-pooled \
  --baseline-arm allreads-grouped \
  --metadata $meta --group-by community --composition $comp --references $refs \
  --catalogue $cat --genomes $genomes --workdir $M/$d/quast-origin --threads 16 --output $M/$d/evaluation-origin.json
echo "== evaluate-origin $d exit=$?" >> $M/log.txt
