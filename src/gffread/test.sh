#!/bin/bash

## VIASH START
## VIASH END

set -e

test_output_dir="${meta_resources_dir}/test_data/test_output"
test_dir="${meta_resources_dir}/test_data"
expected_output_dir="${meta_resources_dir}/test_data/output"

mkdir -p "$test_output_dir"


################################################################################

echo "> Test 1 - Read annotation file, output GFF" 

"$meta_executable" \
  --expose_dups \
  --outfile "$test_output_dir/ann_simple.gff" \
  --input "$test_dir/sequence.gff3" 
 

echo ">> Check if output exists"
[ ! -f "$test_output_dir/ann_simple.gff" ] \
    && echo "Output file test_output/ann_simple.gff does not exist" && exit 1

echo ">> Check if output is empty"
[ ! -s "$test_output_dir/ann_simple.gff" ] \
    && echo "Output file test_output/ann_simple.gff is empty" && exit 1

echo ">> Compare output to expected output"

# compare file expect lines starting with "#"
diff <(grep -v "^#" "$expected_output_dir/ann_simple.gff") \
    <(grep -v "^#" "$test_output_dir/ann_simple.gff") || \
    (echo "Output file ann_simple.gff does not match expected output" && exit 1)

################################################################################

echo "> Test 2 - Read annotation file, output GTF"

"$meta_executable" \
  --gtf_output \
  --outfile "$test_output_dir/annotation.gtf" \
  --input "$test_dir/sequence.gff3"

echo ">> Check if output exists"
[ ! -f "$test_output_dir/annotation.gtf" ] \
    && echo "Output file test_output/annotation.gtf does not exist" && exit 1

echo ">> Check if output is empty"
[ ! -s "$test_output_dir/annotation.gtf" ] \
    && echo "Output file test_output/annotation.gtf is empty" && exit 1

echo ">> Compare output to expected output"
diff "$expected_output_dir/annotation.gtf" "$test_output_dir/annotation.gtf" || \
    (echo "Output file annotation.gtf does not match expected output" && exit 1)

################################################################################

echo "> Test 3 - Generate fasta file from annotation file"


"$meta_executable" \
  --genome "$test_dir/sequence.fasta" \
  --spliced_exons "$test_output_dir/transcripts.fa" \
  --outfile "$test_output_dir/output.gff" \
  --input "$test_dir/sequence.gff3"

echo ">> Check if output exists"
[ ! -f "$test_output_dir/transcripts.fa" ] \
    && echo "Output file transcripts.fa does not exist" && exit 1

echo ">> Check if output is empty"
[ ! -s "$test_output_dir/transcripts.fa" ] \
    && echo "Output file transcripts.fa is empty" && exit 1

echo ">> Compare output to expected output"
diff "$expected_output_dir/transcripts.fa" "$test_output_dir/transcripts.fa" || \
    (echo "Output file transcripts.fa does not match expected output" && exit 1)

################################################################################

echo "> Test 4 - Generate table from GFF annotation file"

"$meta_executable" \
  --table "@id;@chr;@start;@end;@strand;@exons;Name;gene;product" \
  --outfile "$test_output_dir/annotation.tbl" \
  --input "$test_dir/sequence.gff3"

echo ">> Check if output exists"
[ ! -f "$test_output_dir/annotation.tbl" ] \
    && echo "Output file test_output/annotation.tbl does not exist" && exit 1

echo ">> Check if output is empty"
[ ! -s "$test_output_dir/annotation.tbl" ] \
    && echo "Output file test_output/annotation.tbl is empty" && exit 1

echo ">> Compare output to expected output"
diff "$expected_output_dir/annotation.tbl" "$test_output_dir/annotation.tbl" || \
    (echo "Output file annotation.tbl does not match expected output" && exit 1)

################################################################################

rm -r "$test_output_dir"

echo "> All tests successful"

exit 0
