#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

source "$meta_resources_dir/test_helpers.sh"

setup_test_env

log "Starting tests for $meta_name"

test_data="$meta_resources_dir/test_data"
bam="$test_data/test.paired_end.sorted.bam"
bai="$test_data/test.paired_end.sorted.bam.bai"
cram="$test_data/test.cram"
crai="$test_data/test.cram.crai"
fasta="$test_data/reference.fasta"
fasta_fai="$test_data/reference.fasta.fai"

##############################################################
log "Starting TEST 1: Basic run"
##############################################################
"$meta_executable" \
  --input "$bam" \
  --input_index "$bai" \
  --output_summary "$meta_temp_dir/test1_summary.txt" \
  --output_global_dist "$meta_temp_dir/test1_global_dist.txt" \
  --output_per_base "$meta_temp_dir/test1_per_base.bed.gz" \
  --output_per_base_index "$meta_temp_dir/test1_per_base.bed.gz.csi"

check_file_exists "$meta_temp_dir/test1_summary.txt" "Summary output"
check_file_not_empty "$meta_temp_dir/test1_summary.txt" "Summary output"
check_file_contains "$meta_temp_dir/test1_summary.txt" "^total" "Summary output"
check_file_exists "$meta_temp_dir/test1_global_dist.txt" "Global distribution output"
check_file_not_empty "$meta_temp_dir/test1_global_dist.txt" "Global distribution output"
check_file_exists "$meta_temp_dir/test1_per_base.bed.gz" "Per-base output"
check_file_exists "$meta_temp_dir/test1_per_base.bed.gz.csi" "Per-base output index"
log "✅ TEST 1 completed successfully"

##############################################################
log "Starting TEST 2: Windowed coverage with --by"
##############################################################
"$meta_executable" \
  --input "$bam" \
  --input_index "$bai" \
  --by 100 \
  --fragment_mode \
  --output_regions "$meta_temp_dir/test2_regions.bed.gz" \
  --output_regions_index "$meta_temp_dir/test2_regions.bed.gz.csi" \
  --output_region_dist "$meta_temp_dir/test2_region_dist.txt"

check_file_exists "$meta_temp_dir/test2_regions.bed.gz" "Regions output"
check_file_not_empty "$meta_temp_dir/test2_regions.bed.gz" "Regions output"
check_file_exists "$meta_temp_dir/test2_regions.bed.gz.csi" "Regions output index"
check_file_exists "$meta_temp_dir/test2_region_dist.txt" "Region distribution output"
log "✅ TEST 2 completed successfully"

##############################################################
log "Starting TEST 3: --no_per_base skips per-base output"
##############################################################
"$meta_executable" \
  --input "$bam" \
  --input_index "$bai" \
  --no_per_base \
  --output_summary "$meta_temp_dir/test3_summary.txt"

check_file_exists "$meta_temp_dir/test3_summary.txt" "Summary output (no per-base)"
log "✅ TEST 3 completed successfully"

##############################################################
log "Starting TEST 4: --quantize with custom --quantize_labels"
##############################################################
"$meta_executable" \
  --input "$bam" \
  --input_index "$bai" \
  --quantize "0:1:5:100:500" \
  --quantize_labels "NO_COVERAGE,LOW_COVERAGE,CALLABLE,HIGH_COVERAGE" \
  --output_quantized "$meta_temp_dir/test4_quantized.bed.gz" \
  --output_quantized_index "$meta_temp_dir/test4_quantized.bed.gz.csi"

check_file_exists "$meta_temp_dir/test4_quantized.bed.gz" "Quantized output"
check_file_exists "$meta_temp_dir/test4_quantized.bed.gz.csi" "Quantized output index"
gunzip -c "$meta_temp_dir/test4_quantized.bed.gz" >"$meta_temp_dir/test4_quantized.bed"
check_file_contains "$meta_temp_dir/test4_quantized.bed" "CALLABLE" "Quantized output"
log "✅ TEST 4 completed successfully"

##############################################################
log "Starting TEST 5: --thresholds"
##############################################################
"$meta_executable" \
  --input "$bam" \
  --input_index "$bai" \
  --by 100 \
  --thresholds "1,5,10" \
  --output_thresholds "$meta_temp_dir/test5_thresholds.bed.gz" \
  --output_thresholds_index "$meta_temp_dir/test5_thresholds.bed.gz.csi"

check_file_exists "$meta_temp_dir/test5_thresholds.bed.gz" "Thresholds output"
check_file_exists "$meta_temp_dir/test5_thresholds.bed.gz.csi" "Thresholds output index"
gunzip -c "$meta_temp_dir/test5_thresholds.bed.gz" >"$meta_temp_dir/test5_thresholds.bed"
check_file_contains "$meta_temp_dir/test5_thresholds.bed" "1X" "Thresholds output"
log "✅ TEST 5 completed successfully"

##############################################################
log "Starting TEST 6: CRAM input with --fasta"
##############################################################
"$meta_executable" \
  --input "$cram" \
  --input_index "$crai" \
  --fasta "$fasta" \
  --fasta_index "$fasta_fai" \
  --output_summary "$meta_temp_dir/test6_summary.txt"

check_file_exists "$meta_temp_dir/test6_summary.txt" "Summary output (CRAM input)"
check_file_not_empty "$meta_temp_dir/test6_summary.txt" "Summary output (CRAM input)"
log "✅ TEST 6 completed successfully"

##############################################################
log "Starting TEST 7: Read filtering arguments"
##############################################################
"$meta_executable" \
  --input "$bam" \
  --input_index "$bai" \
  --mapq 20 \
  --flag 1796 \
  --include_flag 2 \
  --min_frag_len 50 \
  --max_frag_len 500 \
  --chrom "chr1" \
  --fast_mode \
  --use_median \
  --dist_precision 4 \
  --output_summary "$meta_temp_dir/test7_summary.txt" \
  --output_global_dist "$meta_temp_dir/test7_global_dist.txt"

check_file_exists "$meta_temp_dir/test7_summary.txt" "Summary output (read filtering)"
check_file_not_empty "$meta_temp_dir/test7_summary.txt" "Summary output (read filtering)"
check_file_matches_regex "$meta_temp_dir/test7_global_dist.txt" "\.[0-9]{4}$" "Global distribution output (custom precision)"
log "✅ TEST 7 completed successfully"

##############################################################
log "Starting TEST 8: Error handling for a missing input file"
##############################################################
set +e
"$meta_executable" \
  --input "$meta_temp_dir/does_not_exist.bam" \
  2>"$meta_temp_dir/test8_stderr.txt"
exit_code=$?
set -e

if [ "$exit_code" -eq 0 ]; then
  log_error "Expected a non-zero exit code for a missing input file, got 0"
  exit 1
fi
log "✅ TEST 8 completed successfully"

##############################################################
log "Starting TEST 9: Index auto-discovery when --input_index is omitted"
##############################################################
"$meta_executable" \
  --input "$bam" \
  --output_summary "$meta_temp_dir/test9_summary.txt"

check_file_exists "$meta_temp_dir/test9_summary.txt" "Summary output (auto-discovered index)"
check_file_not_empty "$meta_temp_dir/test9_summary.txt" "Summary output (auto-discovered index)"
log "✅ TEST 9 completed successfully"

print_test_summary "All tests completed successfully"
