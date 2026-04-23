#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

source "$meta_resources_dir/test_helpers.sh"

setup_test_env

log "Starting tests for $meta_name"

gtf="$meta_resources_dir/test_data/test.gtf"
sam="$meta_resources_dir/test_data/test.sam"
bam="$meta_resources_dir/test_data/test.bam"
cram="$meta_resources_dir/test_data/test.cram"

# Test 1: Basic run
log "Starting TEST 1: Basic run"
mkdir -p "$meta_temp_dir/outdir1"
"$meta_executable" \
  --input "$sam" \
  --gtf "$gtf" \
  --outdir "$meta_temp_dir/outdir1" \
  --sample_name "test_sample"

check_dir_exists "$meta_temp_dir/outdir1" "output directory"
[ -z "$(ls -A "$meta_temp_dir/outdir1")" ] && echo "Output directory is empty" && exit 1
log "✅ TEST 1 completed successfully"

# Test 2: Flat output layout
log "Starting TEST 2: Flat output layout"
mkdir -p "$meta_temp_dir/outdir2"
"$meta_executable" \
  --input "$sam" \
  --gtf "$gtf" \
  --outdir "$meta_temp_dir/outdir2" \
  --flat_output \
  --sample_name "test_sample"

check_dir_exists "$meta_temp_dir/outdir2" "flat output directory"
[ -z "$(ls -A "$meta_temp_dir/outdir2")" ] && echo "Flat output directory is empty" && exit 1
log "✅ TEST 2 completed successfully"

# Test 3: BAM input
log "Starting TEST 3: BAM input"
mkdir -p "$meta_temp_dir/outdir3"
"$meta_executable" \
  --input "$bam" \
  --gtf "$gtf" \
  --outdir "$meta_temp_dir/outdir3" \
  --sample_name "test_sample"

check_dir_exists "$meta_temp_dir/outdir3" "BAM output directory"
[ -z "$(ls -A "$meta_temp_dir/outdir3")" ] && echo "BAM output directory is empty" && exit 1
log "✅ TEST 3 completed successfully"

# Test 4: CRAM input
log "Starting TEST 4: CRAM input"
mkdir -p "$meta_temp_dir/outdir4"
"$meta_executable" \
  --input "$cram" \
  --gtf "$gtf" \
  --reference "$meta_resources_dir/test_data/reference.fasta" \
  --outdir "$meta_temp_dir/outdir4" \
  --sample_name "test_sample"

check_dir_exists "$meta_temp_dir/outdir4" "CRAM output directory"
[ -z "$(ls -A "$meta_temp_dir/outdir4")" ] && echo "CRAM output directory is empty" && exit 1
log "✅ TEST 4 completed successfully"

print_test_summary "All tests completed successfully"
