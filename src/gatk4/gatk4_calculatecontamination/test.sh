#!/bin/bash

## VIASH START
## VIASH END

# Source the centralized test helpers
source "$meta_resources_dir/test_helpers.sh"
source "$meta_resources_dir/gatk4/test_helpers.sh"

# Set up the test environment and the test data directory
gatk4_test_init

# --- Build shared pileup summary fixtures ---
# Hand-crafted pileup summary tables in the GetPileupSummaries output format:
# a "#<METADATA>SAMPLE=<name>" line, then contig/position/ref_count/
# alt_count/other_alt_count/allele_frequency columns. The tumor sites
# alternate between mostly-reference sites and mostly-alternate sites, with a
# few reads of the other allele that stand in for cross-sample contamination.
# The normal table has the same sites but is purely homozygous.
log "Writing tumor and normal pileup summary tables..."
{
  echo "#<METADATA>SAMPLE=tumor"
  printf 'contig\tposition\tref_count\talt_count\tother_alt_count\tallele_frequency\n'
} > "$test_data_dir/tumor_pileups.table"
{
  echo "#<METADATA>SAMPLE=normal"
  printf 'contig\tposition\tref_count\talt_count\tother_alt_count\tallele_frequency\n'
} > "$test_data_dir/normal_pileups.table"

for i in $(seq 0 19); do
  position=$((1000 + i * 100))
  if ((i % 2 == 0)); then
    # Mostly reference allele
    printf 'chr1\t%d\t%d\t%d\t0\t0.05\n' "$position" $((38 + i % 5)) $((i % 3 == 0 ? 0 : 1)) >> "$test_data_dir/tumor_pileups.table"
    printf 'chr1\t%d\t%d\t0\t0\t0.05\n' "$position" $((40 + i % 5)) >> "$test_data_dir/normal_pileups.table"
  else
    # Mostly alternate allele
    printf 'chr1\t%d\t%d\t%d\t0\t0.95\n' "$position" $((1 + i % 3)) $((37 + i % 4)) >> "$test_data_dir/tumor_pileups.table"
    printf 'chr1\t%d\t0\t%d\t0\t0.95\n' "$position" $((39 + i % 4)) >> "$test_data_dir/normal_pileups.table"
  fi
done
check_file_exists "$test_data_dir/tumor_pileups.table" "tumor pileup summary table"
check_file_exists "$test_data_dir/normal_pileups.table" "normal pileup summary table"

# --- Test Case 1: Basic run (tumor only) ---
log "Starting TEST 1: Basic run (tumor only)"

# The --tumor_segmentation path is set without --create_tumor_segmentation
# so the file must not be written
log "Executing $meta_name with the tumor pileups only..."
"$meta_executable" \
  --input "$test_data_dir/tumor_pileups.table" \
  --output "$meta_temp_dir/contamination.table" \
  --tumor_segmentation "$meta_temp_dir/segmentation_not_created.table"

log "Validating TEST 1 outputs..."
check_file_exists "$meta_temp_dir/contamination.table" "contamination table"
check_file_contains "$meta_temp_dir/contamination.table" $'sample\tcontamination\terror' "contamination table header"
check_file_matches_regex "$meta_temp_dir/contamination.table" $'^tumor\t' "contamination table tumor row"
check_file_not_exists "$meta_temp_dir/segmentation_not_created.table" "tumor segmentation table (not created)"

log "✅ TEST 1 completed successfully"

# --- Test Case 2: Matched normal and options ---
log "Starting TEST 2: Matched normal, tumor segmentation and coverage thresholds"

log "Executing $meta_name with --matched_normal, --tumor_segmentation and the coverage ratio thresholds..."
"$meta_executable" \
  --input "$test_data_dir/tumor_pileups.table" \
  --matched_normal "$test_data_dir/normal_pileups.table" \
  --output "$meta_temp_dir/contamination_normal.table" \
  --create_tumor_segmentation \
  --tumor_segmentation "$meta_temp_dir/segmentation.table" \
  --high_coverage_ratio_threshold 2.0 \
  --low_coverage_ratio_threshold 0.25 2>&1 | tee "$meta_temp_dir/test2.log"

log "Validating TEST 2 outputs..."
check_file_exists "$meta_temp_dir/contamination_normal.table" "contamination table (matched normal)"
check_file_matches_regex "$meta_temp_dir/contamination_normal.table" $'^tumor\t' "contamination table tumor row (matched normal)"
check_file_exists "$meta_temp_dir/segmentation.table" "tumor segmentation table"
check_file_contains "$meta_temp_dir/segmentation.table" $'contig\tstart\tend\tminor_allele_fraction' "tumor segmentation table header"
check_file_contains "$meta_temp_dir/test2.log" "matched-normal .*normal_pileups.table" "run log command line (--matched_normal)"
check_file_contains "$meta_temp_dir/test2.log" "high-coverage-ratio-threshold 2.0" "run log command line (--high_coverage_ratio_threshold)"
check_file_contains "$meta_temp_dir/test2.log" "low-coverage-ratio-threshold 0.25" "run log command line (--low_coverage_ratio_threshold)"

log "✅ TEST 2 completed successfully"

# --- Test Case 3: Argument errors ---
log "Starting TEST 3: Argument errors"

check_fails_with "Error: --create_tumor_segmentation requires --tumor_segmentation" "$meta_temp_dir/test3_segmentation.log" \
  "$meta_executable" \
  --input "$test_data_dir/tumor_pileups.table" \
  --output "$meta_temp_dir/contamination_bad.table" \
  --create_tumor_segmentation

log "✅ TEST 3 completed successfully"

print_test_summary "All tests"
