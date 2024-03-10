#!/bin/bash

## VIASH START
## VIASH END

set -e

test_output_dir="$meta_resources_dir/test_data/test_output/"
test_dir = "$meta_resources_dir/test_data/"
expected_output_dir = "$meta_resources_dir/test_data/output/"
echo "> Run gffread with a single genome"

################################################################################

echo "> Test 1 - Read annotation file, output GFF" 

"$meta_executable" \
    -E annotation.gff \
    -o "$test_output_dir/ann_simple.gff"

echo ">> Check if output exists"
[ ! -f "$test_output_dir/ann_simple.gff" ] && echo "Output file test_output/ann_simple.gff does not exist" && exit 1

echo ">> Check if output is empty"
[ ! -s "$test_output_dir/ann_simple.gff" ] && echo "Output file test_output/ann_simple.gff is empty" && exit 1

echo ">> Compare output to expected output"
diff "$expected_output_dir/ann_simple.gff" "$test_output_dir/ann_simple.gff" || (echo "Output file ann_simple.gff does not match expected output" && exit 1)

################################################################################

echo "> Test 2 - Read annotation file, output GTF"

"$meta_executable" \
    annotation.gff \
    -T \
    -o "$test_output_dir/annotation.gtf"

echo ">> Check if output exists"
[ ! -f "$test_output_dir/annotation.gtf" ] && echo "Output file test_output/annotation.gtf does not exist" && exit 1

echo ">> Check if output is empty"
[ ! -s "$test_output_dir/annotation.gtf" ] && echo "Output file test_output/annotation.gtf is empty" && exit 1

echo ">> Compare output to expected output"
diff "$expected_output_dir/annotation.gtf" "$test_output_dir/annotation.gtf" || (echo "Output file annotation.gtf does not match expected output" && exit 1)

################################################################################

echo "> Test 3 - Generate fasta file from annotation file"

"$meta_executable" \
    -w "$test_output_dir/transcripts.fa" \
    -g genome.fa \
    annotation.gff

echo ">> Check if output exists"
[ ! -f "$test_output_dir/transcripts.fa" ] && echo "Output file test_output/transcripts.fa does not exist" && exit 1

echo ">> Check if output is empty"
[ ! -s "$test_output_dir/transcripts.fa" ] && echo "Output file test_output/transcripts.fa is empty" && exit 1

echo ">> Compare output to expected output"
diff "$expected_output_dir/transcripts.fa" "$test_output_dir/transcripts.fa" || (echo "Output file transcripts.fa does not match expected output" && exit 1)

################################################################################

echo "> Test 4 - Generate table from GFF annotation file"

"$meta_executable" \
    --table @id,@chr,@start,@end,@strand,@exons,Name,gene,product \
    -o "$test_output_dir/annotation.tbl" \
    annotation.gff

echo ">> Check if output exists"
[ ! -f "$test_output_dir/annotation.tbl" ] && echo "Output file test_output/annotation.tbl does not exist" && exit 1

echo ">> Check if output is empty"
[ ! -s "$test_output_dir/annotation.tbl" ] && echo "Output file test_output/annotation.tbl is empty" && exit 1

echo ">> Compare output to expected output"
diff "$expected_output_dir/annotation.tbl" "$test_output_dir/annotation.tbl" || (echo "Output file annotation.tbl does not match expected output" && exit 1)

################################################################################

echo "> All tests successful"