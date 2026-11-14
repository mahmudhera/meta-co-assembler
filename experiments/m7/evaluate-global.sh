#!/bin/sh
# Contig correctness ($1 = pilot | split | hard): each arm's whole contig output against every
# source genome of the dataset (one community "all"), so contigs of a genome confined to one
# community are not scored against another community's congeners.
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
awk -F'\t' 'BEGIN{OFS="\t"} NR==1{print "sample_id","scope"; next} {print $1,"all"}' $meta > $M/$d/metadata-all.tsv
python scripts/evaluate_assemblies.py \
  --arm discovery-grouped=$M/$d/discovery-grouped \
  --arm discovery-per-sample=$M/$d/discovery-per-sample \
  --arm discovery-pooled=$M/$d/discovery-pooled \
  --arm iterate-residual-only=$M/$d/iterate-residual-only \
  --arm iterate-recruitment=$M/$d/iterate-recruitment \
  --arm allreads-per-sample=$M/$d/allreads-per-sample \
  --arm allreads-grouped=$M/$d/allreads-grouped \
  --arm allreads-pooled=$M/$d/allreads-pooled \
  --discovery-arm discovery-pooled --baseline-arm allreads-grouped \
  --metadata $M/$d/metadata-all.tsv --group-by scope --composition $comp --references $refs \
  --catalogue $cat --genomes $genomes --workdir $M/$d/quast-global --threads 16 --output $M/$d/evaluation-global.json
echo "== evaluate-global $d exit=$?" >> $M/log.txt
