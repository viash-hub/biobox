#!/bin/bash

set -e

dir_in="$meta_resources_dir/test_data"

echo "> Run featureCounts"
"$meta_executable" \
  --input "$dir_in/a.bam" \
  --annotation "$dir_in/annotation.gtf" \
  --output_counts "features.tsv" \
  --output_summary "summary.tsv" \
  --output_junctions "junction_counts.txt" \
  --ref_fasta "$dir_in/genome.fasta" \
  --overlapping \
  --frac_overlap 0.2 \
  --paired \
  --strand 0

echo ">> Checking output"
[ ! -f "features.tsv" ] && echo "Output file features.tsv does not exist" && exit 1
[ ! -f "summary.tsv" ] && echo "Output file summary.tsv does not exist" && exit 1
[ ! -f "junction_counts.txt" ] && echo "Output file junction_counts.txt does not exist" && exit 1

echo ">> Check if output is empty"
[ ! -s "features.tsv" ] && echo "Output file features.tsv is empty" && exit 1
[ ! -s "summary.tsv" ] && echo "Output file summary.tsv is empty" && exit 1
[ ! -s "junction_counts.txt" ] && echo "Output file junction_counts.txt is empty" && exit 1

echo "> Test successful"