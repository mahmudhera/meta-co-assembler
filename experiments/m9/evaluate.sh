#!/bin/sh
# Truth-based evaluation of every M9 `bin` run of a dataset ($1 = pilot | hard).
set -u
cd "$(dirname "$0")/../.."
M=work/m9
d=$1
case $d in
  pilot) comp=experiments/pilot/composition.tsv refs=data/pilot/references
         genomes=experiments/pilot/genomes.tsv cat=data/pilot/catalogues/species_holdout/catalogue.tsv ;;
  hard)  comp=experiments/m7-hard/composition.tsv refs=data/m7-hard/references
         genomes=experiments/m7-hard/genomes.tsv cat=data/m7-hard/catalogues/hard_holdout/catalogue.tsv ;;
esac
runs=""
for dir in $M/$d/*/; do
  name=$(basename "$dir")
  [ -f "$dir/stages/summarize/bins.tsv" ] && runs="$runs --run $name=$dir"
done
python scripts/evaluate_bins.py $runs --composition $comp --references $refs --genomes $genomes \
  --catalogue $cat --threads 16 --output $M/$d/evaluation.json
echo "== evaluate $d exit=$?" >> $M/log.txt
