#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

source "$meta_resources_dir/test_helpers.sh"

setup_test_env

log "Starting smoke tests for $meta_name"

# NOTE: This is a smoke test. A full functional test of `seekarctools arc run`
# requires a multi-GB SeekArc reference genome (GRCh38/GRCm38) and FASTQ pairs with
# barcodes matching the selected chemistry, neither of which can be synthesised with
# the test_helpers or downloaded within a self-contained test. These checks instead
# verify that the SeekArc Tools install in the container is functional and that the
# command-line interface the component wraps is wired up correctly.

# Test 1: the bundled seekarctools binary is installed and reports its version
log "Starting TEST 1: seekarctools is installed"
seekarctools --version > "$meta_temp_dir/version.txt" 2>&1
check_file_not_empty "$meta_temp_dir/version.txt" "seekarctools version output"
log "✅ TEST 1 completed successfully"

# Test 2: the 'arc run' subcommand the component wraps is available
log "Starting TEST 2: 'arc run' help is available"
seekarctools arc run --help > "$meta_temp_dir/arc_run_help.txt" 2>&1
check_file_not_empty "$meta_temp_dir/arc_run_help.txt" "arc run help output"
check_file_contains "$meta_temp_dir/arc_run_help.txt" "rnafq1" "arc run help mentions --rnafq1"
log "✅ TEST 2 completed successfully"

print_test_summary "All smoke tests completed successfully"
