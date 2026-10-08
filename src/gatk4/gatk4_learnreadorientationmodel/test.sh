#!/bin/bash

## VIASH START
## VIASH END

# Source the centralized test helpers
source "$meta_resources_dir/test_helpers.sh"
source "$meta_resources_dir/gatk4/test_helpers.sh"

# Initialize test environment with strict error handling
setup_test_env

#############################################
# Test execution with centralized functions
#############################################

log "Starting tests for $meta_name"

# Create test data directory
test_data_dir="$meta_temp_dir/test_data"
mkdir -p "$test_data_dir"

# --- Build shared fixtures: a real F1R2 tarball from Mutect2 ---
create_test_mutect2_outputs "$test_data_dir"
check_file_exists "$test_data_dir/f1r2.tar.gz" "F1R2 tarball from the tumor/normal Mutect2 run"

# A second F1R2 tarball from a tumor-only Mutect2 run, to stand in for the
# output of a second shard of a scattered Mutect2 run
log "Running Mutect2 in tumor-only mode for a second F1R2 tarball..."
gatk Mutect2 \
  --reference "$test_data_dir/reference.fasta" \
  --input "$test_data_dir/tumor.bam" \
  --output "$test_data_dir/tumor_only.vcf" \
  --f1r2-tar-gz "$test_data_dir/f1r2_tumor_only.tar.gz" \
  --verbosity ERROR
check_file_exists "$test_data_dir/f1r2_tumor_only.tar.gz" "F1R2 tarball from the tumor-only Mutect2 run"

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

# --- Test Case 1: Default parameters ---
log "Starting TEST 1: Default parameters"

log "Executing $meta_name with a single F1R2 tarball..."
"$meta_executable" \
  --input "$test_data_dir/f1r2.tar.gz" \
  --output "$meta_temp_dir/model.tar.gz"

log "Validating TEST 1 outputs..."
check_orientation_model "$meta_temp_dir/model.tar.gz"

log "✅ TEST 1 completed successfully"

# --- Test Case 2: Multiple F1R2 tarballs ---
log "Starting TEST 2: Multiple F1R2 tarballs"

log "Executing $meta_name with two F1R2 tarballs..."
"$meta_executable" \
  --input "$test_data_dir/f1r2.tar.gz" \
  --input "$test_data_dir/f1r2_tumor_only.tar.gz" \
  --output "$meta_temp_dir/model_multiple.tar.gz"

log "Validating TEST 2 outputs..."
check_orientation_model "$meta_temp_dir/model_multiple.tar.gz"

log "✅ TEST 2 completed successfully"

# --- Test Case 3: EM options ---
log "Starting TEST 3: EM options"

log "Executing $meta_name with --convergence_threshold, --max_depth and --num_em_iterations..."
"$meta_executable" \
  --input "$test_data_dir/f1r2.tar.gz" \
  --output "$meta_temp_dir/model_options.tar.gz" \
  --convergence_threshold 0.001 \
  --max_depth 100 \
  --num_em_iterations 10

log "Validating TEST 3 outputs..."
check_orientation_model "$meta_temp_dir/model_options.tar.gz"

log "✅ TEST 3 completed successfully"

print_test_summary "All tests"
