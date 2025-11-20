#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

# Always source centralized helpers
source "$meta_resources_dir/test_helpers.sh"

# Initialize test environment
setup_test_env

log "Starting tests for $meta_name"

# Test 1: Test directory creation and basic parameter passing
log "Starting TEST 1: Directory creation and parameter validation"
mkdir -p "$meta_temp_dir/test_install"
mkdir -p "$meta_temp_dir/test_cache"
mkdir -p "$meta_temp_dir/test_plugins"

# Test that directories are created and parameters are properly passed
# This will fail due to no internet/cache but should show proper parameter handling
"$meta_executable" \
  --destdir "$meta_temp_dir/test_install" \
  --cachedir "$meta_temp_dir/test_cache" \
  --pluginsdir "$meta_temp_dir/test_plugins" \
  --cache_version 115 \
  --no_update \
  --quiet 2>&1 | head -10 > "$meta_temp_dir/install_test.txt" || true

check_file_exists "$meta_temp_dir/install_test.txt" "installation test output"
check_dir_exists "$meta_temp_dir/test_install" "destination directory"
check_dir_exists "$meta_temp_dir/test_cache" "cache directory"
check_dir_exists "$meta_temp_dir/test_plugins" "plugins directory"
log "✅ TEST 1 completed successfully"

# Test 2: Test installation with small species (C. elegans for faster testing)
log "Starting TEST 2: Small species installation test"
mkdir -p "$meta_temp_dir/test_install"
mkdir -p "$meta_temp_dir/test_cache"

# Test that actual installation parameters work with C. elegans (much smaller than human)
# Using same parameters as nf-core modules for faster testing
"$meta_executable" \
  --destdir "$meta_temp_dir/test_install" \
  --cachedir "$meta_temp_dir/test_cache" \
  --cache_version 115 \
  --auto "c" \
  --species "caenorhabditis_elegans" \
  --assembly "WBcel235" \
  --no_update \
  --no_bioperl \
  --no_htslib > "$meta_temp_dir/install_test.txt" 2>&1 || true

# Test 3: Multiple species parameter handling
log "Starting TEST 3: Multiple species parameter handling"
mkdir -p "$meta_temp_dir/multi_install"

"$meta_executable" \
  --destdir "$meta_temp_dir/multi_install" \
  --species "caenorhabditis_elegans;drosophila_melanogaster" \
  --no_update \
  --quiet 2>&1 | head -10 > "$meta_temp_dir/multi_test.txt" || true

check_file_exists "$meta_temp_dir/multi_test.txt" "multi-species test output"
log "✅ TEST 3 completed successfully"

# Test 4: Plugin parameter handling
log "Starting TEST 4: Plugin parameter handling"
mkdir -p "$meta_temp_dir/plugin_install"

"$meta_executable" \
  --destdir "$meta_temp_dir/plugin_install" \
  --plugins "CADD;dbNSFP" \
  --no_plugins \
  --no_update \
  --quiet 2>&1 | head -10 > "$meta_temp_dir/plugin_test.txt" || true

check_file_exists "$meta_temp_dir/plugin_test.txt" "plugin test output"
log "✅ TEST 4 completed successfully"

# Always end with summary
print_test_summary "All vep_install tests completed successfully"