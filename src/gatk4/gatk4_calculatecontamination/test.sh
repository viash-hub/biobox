#!/bin/bash

## VIASH START
## VIASH END

# Source the centralized test helpers
source "$meta_resources_dir/test_helpers.sh"

# Initialize test environment with strict error handling
setup_test_env

#############################################
# Test execution with centralized functions
#############################################

log "Starting tests for $meta_name"

# Create test data directory
test_data_dir="$meta_temp_dir/test_data"
mkdir -p "$test_data_dir"

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

# --- Test Case 1: Tumor only ---
log "Starting TEST 1: Tumor only"

log "Executing $meta_name with the tumor pileups only..."
"$meta_executable" \
  --input "$test_data_dir/tumor_pileups.table" \
  --output "$meta_temp_dir/contamination.table"

log "Validating TEST 1 outputs..."
check_file_exists "$meta_temp_dir/contamination.table" "contamination table"
check_file_contains "$meta_temp_dir/contamination.table" $'sample\tcontamination\terror' "contamination table header"
check_file_matches_regex "$meta_temp_dir/contamination.table" $'^tumor\t' "contamination table tumor row"

log "✅ TEST 1 completed successfully"

# --- Test Case 2: Matched normal and tumor segmentation ---
log "Starting TEST 2: Matched normal and tumor segmentation"

log "Executing $meta_name with --matched_normal and --tumor_segmentation..."
"$meta_executable" \
  --input "$test_data_dir/tumor_pileups.table" \
  --matched_normal "$test_data_dir/normal_pileups.table" \
  --output "$meta_temp_dir/contamination_normal.table" \
  --tumor_segmentation "$meta_temp_dir/segmentation.table"

log "Validating TEST 2 outputs..."
check_file_exists "$meta_temp_dir/contamination_normal.table" "contamination table (matched normal)"
check_file_matches_regex "$meta_temp_dir/contamination_normal.table" $'^tumor\t' "contamination table tumor row (matched normal)"
check_file_exists "$meta_temp_dir/segmentation.table" "tumor segmentation table"
check_file_contains "$meta_temp_dir/segmentation.table" $'contig\tstart\tend\tminor_allele_fraction' "tumor segmentation table header"

log "✅ TEST 2 completed successfully"

# --- Test Case 3: Coverage thresholds, no tumor segmentation ---
log "Starting TEST 3: Coverage thresholds without tumor segmentation"

log "Executing $meta_name with --high_coverage_ratio_threshold and --low_coverage_ratio_threshold..."
"$meta_executable" \
  --input "$test_data_dir/tumor_pileups.table" \
  --output "$meta_temp_dir/contamination_thresholds.table" \
  --high_coverage_ratio_threshold 2.0 \
  --low_coverage_ratio_threshold 0.25

log "Validating TEST 3 outputs..."
check_file_exists "$meta_temp_dir/contamination_thresholds.table" "contamination table (thresholds)"
check_file_matches_regex "$meta_temp_dir/contamination_thresholds.table" $'^tumor\t' "contamination table tumor row (thresholds)"
# Only TEST 2 asked for a tumor segmentation table
segmentation_count=$(find "$meta_temp_dir" -maxdepth 1 -name '*segment*' | wc -l)
if [[ "$segmentation_count" -eq 1 ]]; then
  log "✓ No tumor segmentation table was written without --tumor_segmentation"
else
  log_error "✗ Expected 1 tumor segmentation table (from TEST 2), found $segmentation_count"
  exit 1
fi

log "✅ TEST 3 completed successfully"

print_test_summary "All tests"
