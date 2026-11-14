#!/bin/sh
# M9: bin the M7 assemblies of the hard benchmark and the pilot (design agreed with the owner).
# Assemblies: pooled discovery coassembly, and pooled and community-grouped all-read assembly
# (the M7 runs in work/m7). Arms per assembly: MetaBAT 2 with strobealign --aemb and with
# minimap2 coverage, and SemiBin2 with strobealign; SemiBin2 with minimap2 on the pooled
# discovery coassembly only (SemiBin2 takes minutes per run on CPU).
# Every `bin` run uses 16 threads / 64 GiB under /usr/bin/time -v; within a dataset runs are
# serial, and the two datasets run concurrently on a 64-CPU node. BIN defaults to the release
# build; CHECKM2_DB to the database downloaded with `checkm2 database --download --path
# data/checkm2`. Outputs go to work/m9 (gitignored).
set -u
cd "$(dirname "$0")/../.."
M=work/m9
B="${BIN:-target/release/meta-co-assembler}"
DB="${CHECKM2_DB:-data/checkm2/CheckM2_database/uniref100.KO.1.dmnd}"
R="--threads 16 --memory-gib 64"
mkdir -p $M
step() { d=$1; name=$2; shift 2; start=$(date +%s); echo "== $d/$name started $(date -Is)" >> $M/log.txt
  /usr/bin/time -v -o $M/$d/$name.time "$@" > $M/$d/$name.json 2> $M/$d/$name.err; code=$?
  echo "== $d/$name exit=$code seconds=$(( $(date +%s) - start ))" >> $M/log.txt; }
catalogue_of() { case $1 in pilot) echo data/pilot/catalogues/species_holdout/catalogue.tsv ;;
  hard) echo data/m7-hard/catalogues/hard_holdout/catalogue.tsv ;; esac; }
binning() { d=$1 arm=$2 binner=$3 mapper=$4
  step $d $arm-$binner-$mapper $B bin --assembly-dir work/m7/$d/$arm --binner $binner --mapper $mapper \
    --catalogue "$(catalogue_of $d)" --checkm2-database "$DB" --output-dir $M/$d/$arm-$binner-$mapper $R; }
run_dataset() { d=$1; mkdir -p $M/$d
  for arm in discovery-pooled allreads-pooled allreads-grouped; do
    for mapper in strobealign minimap2; do binning $d $arm metabat2 $mapper; done
  done
  for arm in discovery-pooled allreads-pooled allreads-grouped; do binning $d $arm semibin2 strobealign; done
  binning $d discovery-pooled semibin2 minimap2
}
for d in ${DATASETS:-hard pilot}; do run_dataset $d & done
wait
echo "== all done $(date -Is)" >> $M/log.txt
