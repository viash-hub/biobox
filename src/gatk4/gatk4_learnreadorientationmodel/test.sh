#!/bin/bash

## VIASH START
## VIASH END

# Source the centralized test helpers
source "$meta_resources_dir/test_helpers.sh"
source "$meta_resources_dir/gatk4/test_helpers.sh"

# Set up the test environment and the test data directory
gatk4_test_init

# --- Build shared fixtures: a real F1R2 tarball from Mutect2 ---
create_test_mutect2_outputs "$test_data_dir"
check_file_exists "$test_data_dir/f1r2.tar.gz" "F1R2 tarball from the tumor/normal Mutect2 run"

# A copy of the F1R2 tarball under another name, to stand in for the output of
# a second shard of a scattered Mutect2 run
mkdir -p "$test_data_dir/shard2"
cp "$test_data_dir/f1r2.tar.gz" "$test_data_dir/shard2/f1r2_shard2.tar.gz"

# Check that a read orientation model tarball contains the prior tables
check_orientation_model() {
  local model_path="$1"
  check_file_exists "$model_path" "read orientation model"
  check_file_not_empty "$model_path" "read orientation model"
  if tar tzf "$model_path" | grep -q "orientation_priors"; then
    log "✓ Read orientation model contains an orientation_priors table"
  else
    log_error "✗ Read orientation model does not contain an orientation_priors table: $(tar tzf "$model_path")"
    exit 1
  fi
}

# --- Test Case 1: Basic run (defaults) ---
log "Starting TEST 1: Basic run (defaults)"

log "Executing $meta_name with a single F1R2 tarball..."
"$meta_executable" \
  --input "$test_data_dir/f1r2.tar.gz" \
  --output "$meta_temp_dir/model.tar.gz"

log "Validating TEST 1 outputs..."
check_orientation_model "$meta_temp_dir/model.tar.gz"

log "✅ TEST 1 completed successfully"

# --- Test Case 2: Inputs and tool options ---
log "Starting TEST 2: Two F1R2 tarballs and EM options"

log "Executing $meta_name with two F1R2 tarballs, --convergence_threshold, --max_depth and --num_em_iterations..."
"$meta_executable" \
  --input "$test_data_dir/f1r2.tar.gz" \
  --input "$test_data_dir/shard2/f1r2_shard2.tar.gz" \
  --output "$meta_temp_dir/model_options.tar.gz" \
  --convergence_threshold 0.001 \
  --max_depth 100 \
  --num_em_iterations 10 2>&1 | tee "$meta_temp_dir/test2.log"

log "Validating TEST 2 outputs..."
check_orientation_model "$meta_temp_dir/model_options.tar.gz"
check_file_contains "$meta_temp_dir/test2.log" "input .*f1r2_shard2.tar.gz" "run log command line (second --input)"
check_file_contains "$meta_temp_dir/test2.log" "convergence-threshold 0.001" "run log command line (--convergence_threshold)"
check_file_contains "$meta_temp_dir/test2.log" "max-depth 100" "run log command line (--max_depth)"
check_file_contains "$meta_temp_dir/test2.log" "num-em-iterations 10" "run log command line (--num_em_iterations)"

log "✅ TEST 2 completed successfully"

print_test_summary "All tests"
