#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

source "$meta_resources_dir/test_helpers.sh"

setup_test_env

log "Starting tests for $meta_name"

###########################################################################

log "Starting TEST 1: basic download"

mkdir -p "$meta_temp_dir/test1"

"$meta_executable" \
  --genome_version AY278488.2 \
  --output "$meta_temp_dir/test1"

check_dir_exists "$meta_temp_dir/test1/AY278488.2" "downloaded database directory"
check_file_exists "$meta_temp_dir/test1/AY278488.2/snpEffectPredictor.bin" "downloaded database"
check_file_not_empty "$meta_temp_dir/test1/AY278488.2/snpEffectPredictor.bin" "downloaded database"

log "✅ TEST 1 completed successfully"

###########################################################################

log "Starting TEST 2: error handling for an unknown genome"

mkdir -p "$meta_temp_dir/test2"

if "$meta_executable" \
  --genome_version NoSuchGenomeXYZ123 \
  --output "$meta_temp_dir/test2" 2>/dev/null; then
  log_error "Component should have failed for an unknown genome version"
  exit 1
else
  log "✓ Component properly failed for an unknown genome version"
fi

log "✅ TEST 2 completed successfully"

###########################################################################

print_test_summary "All snpeff_download tests completed successfully"
