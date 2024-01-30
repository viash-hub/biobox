#!/bin/bash

set -e

dir_in="$meta_resources_dir/test_data"

echo "> Run arriba with blacklist"
"$meta_executable" \
  --bam "$dir_in/A.bam" \
  --genome "$dir_in/genome.fasta" \
  --gene_annotation "$dir_in/annotation.gtf" \
  --blacklist "$dir_in/blacklist.tsv" \
  --fusions "fusions.tsv" \
  --fusions_discarded "fusions_discarded.tsv" \
  --interesting_contigs "1,2"

echo ">> Checking output"
[ ! -f "fusions.tsv" ] && echo "Output file fusions.tsv does not exist" && exit 1
[ ! -f "fusions_discarded.tsv" ] && echo "Output file fusions_discarded.tsv does not exist" && exit 1

echo ">> Check if output is empty"
[ ! -s "fusions.tsv" ] && echo "Output file fusions.tsv is empty" && exit 1
[ ! -s "fusions_discarded.tsv" ] && echo "Output file fusions_discarded.tsv is empty" && exit 1

rm fusions.tsv fusions_discarded.tsv

echo "> Run arriba without blacklist"
"$meta_executable" \
  --bam "$dir_in/A.bam" \
  --genome "$dir_in/genome.fasta" \
  --gene_annotation "$dir_in/annotation.gtf" \
  --fusions "fusions.tsv" \
  --fusions_discarded "fusions_discarded.tsv" \
  --interesting_contigs "1,2" \
  --disable_filters blacklist

echo ">> Checking output"
[ ! -f "fusions.tsv" ] && echo "Output file fusions.tsv does not exist" && exit 1
[ ! -f "fusions_discarded.tsv" ] && echo "Output file fusions_discarded.tsv does not exist" && exit 1

echo ">> Check if output is empty"
[ ! -s "fusions.tsv" ] && echo "Output file fusions.tsv is empty" && exit 1
[ ! -s "fusions_discarded.tsv" ] && echo "Output file fusions_discarded.tsv is empty" && exit 1

echo "> Test successful"