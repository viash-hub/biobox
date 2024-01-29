#!/bin/bash

set -e

dir_in="$meta_resources_dir/test_data"

echo "> Run featurCounts"
"$meta_executable" \
  --input "$dir_in/a.bam" \
  --annotation "$dir_in/annotation.gtf" \
  --output "features.txt" \
  --output_summary "features.txt.summary" \
  --output_junctions "features.txt.jcounts" \
  --junctions \
  --ref_fasta "$dir_in/genome.fasta" \
  --overlapping \
  --frac_overlap 0.2 \
  --paired \
  --strand 0

echo ">> Checking output"
[ ! -f "features.txt" ] && echo "Output file features.txt does not exist" && exit 1
[ ! -f "features.txt.summary" ] && echo "Output file features.txt.summary does not exist" && exit 1
[ ! -f "features.txt.jcounts" ] && echo "Output file features.txt.jcounts does not exist" && exit 1

echo ">> Check if output is empty"
[ ! -s "features.txt" ] && echo "Output file features.txt is empty" && exit 1
[ ! -s "features.txt.summary" ] && echo "Output file features.txt.summary is empty" && exit 1
[ ! -s "features.txt.jcounts" ] && echo "Output file features.txt.jcounts is empty" && exit 1

echo "> Test successful"