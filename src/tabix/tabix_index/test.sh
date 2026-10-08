#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

source "$meta_resources_dir/test_helpers.sh"
setup_test_env

log "Starting tests for $meta_name"

test_dir="$meta_temp_dir/test_data"
mkdir -p "$test_dir"

create_test_bed "$test_dir/regions.bed" 3
bgzip -c "$test_dir/regions.bed" >"$test_dir/regions.bed.gz"
check_file_exists "$test_dir/regions.bed.gz" "bgzipped test BED"

# --- TEST 1: index a bgzipped BED file using the bed preset ---
log "TEST 1: index a bgzipped BED file using --preset bed"
"$meta_executable" \
  --input "$test_dir/regions.bed.gz" \
  --output_index "$meta_temp_dir/regions.bed.gz.tbi" \
  --preset bed

check_file_exists "$meta_temp_dir/regions.bed.gz.tbi" "TBI index"
check_file_not_empty "$meta_temp_dir/regions.bed.gz.tbi" "TBI index"
log "✅ TEST 1 passed"

# --- TEST 2: create a CSI index instead of TBI ---
log "TEST 2: index with --preset bed --csi"
"$meta_executable" \
  --input "$test_dir/regions.bed.gz" \
  --output_index "$meta_temp_dir/regions.bed.gz.csi" \
  --preset bed \
  --csi

check_file_exists "$meta_temp_dir/regions.bed.gz.csi" "CSI index"
check_file_not_empty "$meta_temp_dir/regions.bed.gz.csi" "CSI index"
log "✅ TEST 2 passed"

# --- TEST 3: explicit sequence/begin/end columns instead of a preset ---
log "TEST 3: index using explicit --sequence/--begin/--end columns"
"$meta_executable" \
  --input "$test_dir/regions.bed.gz" \
  --output_index "$meta_temp_dir/regions_explicit.bed.gz.tbi" \
  --sequence 1 \
  --begin 2 \
  --end 3 \
  --zero_based

check_file_exists "$meta_temp_dir/regions_explicit.bed.gz.tbi" "TBI index (explicit columns)"
check_file_not_empty "$meta_temp_dir/regions_explicit.bed.gz.tbi" "TBI index (explicit columns)"
log "✅ TEST 3 passed"

# --- TEST 4: --force accepts rebuilding an index ---
log "TEST 4: --force re-creates an index"
"$meta_executable" \
  --input "$test_dir/regions.bed.gz" \
  --output_index "$meta_temp_dir/regions.bed.gz.tbi" \
  --preset bed \
  --force

check_file_exists "$meta_temp_dir/regions.bed.gz.tbi" "TBI index (forced rebuild)"
check_file_not_empty "$meta_temp_dir/regions.bed.gz.tbi" "TBI index (forced rebuild)"
log "✅ TEST 4 passed"

print_test_summary "$meta_name tests passed"
