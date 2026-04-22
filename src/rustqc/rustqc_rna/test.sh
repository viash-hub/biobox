#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

source "$meta_resources_dir/test_helpers.sh"

setup_test_env

log "Starting tests for $meta_name"

# Create a minimal GTF annotation
log "Creating test data..."
{
  printf 'chr1\ttest\tgene\t1000\t1999\t.\t+\t.\tgene_id "gene1"; gene_name "GENE1";\n'
  printf 'chr1\ttest\ttranscript\t1000\t1999\t.\t+\t.\tgene_id "gene1"; transcript_id "tx1"; gene_name "GENE1";\n'
  printf 'chr1\ttest\texon\t1000\t1499\t.\t+\t.\tgene_id "gene1"; transcript_id "tx1"; exon_number "1";\n'
  printf 'chr1\ttest\texon\t1500\t1999\t.\t+\t.\tgene_id "gene1"; transcript_id "tx1"; exon_number "2";\n'
} > "$meta_temp_dir/test.gtf"

# Create a minimal SAM file (rustqc supports SAM directly)
{
  printf '@HD\tVN:1.6\tSO:coordinate\n'
  printf '@SQ\tSN:chr1\tLN:5000\n'
  printf '@PG\tID:test\tPN:test\tVN:0.1\n'
  for i in $(seq 1 30); do
    pos=$((1000 + (i - 1) * 5))
    printf "read%d\t0\tchr1\t%d\t60\t50M\t*\t0\t0\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\tNH:i:1\n" \
      "$i" "$pos"
  done
} > "$meta_temp_dir/test.sam"

# Test 1: Basic run
log "Starting TEST 1: Basic run"
mkdir -p "$meta_temp_dir/outdir1"
"$meta_executable" \
  --input "$meta_temp_dir/test.sam" \
  --gtf "$meta_temp_dir/test.gtf" \
  --outdir "$meta_temp_dir/outdir1" \
  --skip_dup_check \
  --skip_tin \
  --skip_preseq \
  --sample_name "test_sample"

check_dir_exists "$meta_temp_dir/outdir1" "output directory"
[ -z "$(ls -A "$meta_temp_dir/outdir1")" ] && echo "Output directory is empty" && exit 1
log "✅ TEST 1 completed successfully"

# Test 2: Flat output layout
log "Starting TEST 2: Flat output layout"
mkdir -p "$meta_temp_dir/outdir2"
"$meta_executable" \
  --input "$meta_temp_dir/test.sam" \
  --gtf "$meta_temp_dir/test.gtf" \
  --outdir "$meta_temp_dir/outdir2" \
  --skip_dup_check \
  --skip_tin \
  --skip_preseq \
  --flat_output \
  --sample_name "test_sample"

check_dir_exists "$meta_temp_dir/outdir2" "flat output directory"
[ -z "$(ls -A "$meta_temp_dir/outdir2")" ] && echo "Flat output directory is empty" && exit 1
log "✅ TEST 2 completed successfully"

print_test_summary "All tests completed successfully"
