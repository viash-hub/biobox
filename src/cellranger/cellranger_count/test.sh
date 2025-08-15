#!/bin/bash

## VIASH START
## VIASH END

# Include test helpers for enhanced logging and utilities
source "${meta_resources_dir}/test_helpers.sh"

# Set up test environment with enhanced error handling 
setup_test_env

# create temporary directory
tmp_dir=$(mktemp -d "${meta_temp_dir}/${meta_name}-XXXXXXXX")
function clean_up {
    log "Cleaning up temporary directory: $tmp_dir"
    rm -rf "${tmp_dir}"
}
trap clean_up EXIT

log "Created temporary directory: $tmp_dir"

log "Copying test data from Cell Ranger installation directory"
test_data_dir="${tmp_dir}/test_data"
mkdir -p "$test_data_dir"

# Define source and destination paths
cellranger_external_dir="/opt/cellranger-9.0.1/external"
source_fastq_dir="$cellranger_external_dir/cellranger_tiny_fastq"
source_ref_dir="$cellranger_external_dir/cellranger_tiny_ref"
dest_fastq_dir="$test_data_dir/cellranger_tiny_fastq"
dest_ref_dir="$test_data_dir/cellranger_tiny_ref"

# Validate source directories exist
check_dir_exists "$source_fastq_dir" "source FASTQ directory"
check_dir_exists "$source_ref_dir" "source reference directory"

# Copy test data
cp -r "$source_fastq_dir" "$test_data_dir"
cp -r "$source_ref_dir" "$test_data_dir"

# Validate copied directories
check_dir_exists "$dest_fastq_dir" "destination FASTQ directory"
check_dir_exists "$dest_ref_dir" "destination reference directory"

log "Test data copied successfully"


## TEST 1: run with folder input
log "Starting TEST 1: Running ${meta_name} with folder input"
test1_output_dir="${tmp_dir}/test1"
mkdir -p "$test1_output_dir"

log "Executing cellranger count with folder input..."
"${meta_executable}" \
  --fastqs "$dest_fastq_dir" \
  --transcriptome "$dest_ref_dir" \
  --output "$test1_output_dir" \
  --lanes 1

log "Validating TEST 1 outputs..."
check_file_exists "$test1_output_dir/filtered_feature_bc_matrix.h5" "filtered feature matrix"
check_file_not_exists "$test1_output_dir/possorted_genome_bam.bam" "BAM file (should be disabled by default)"
check_dir_exists "$test1_output_dir/analysis" "analysis output directory"

log "✅ TEST 1 completed successfully"

## TEST 2: run with individual file input
log "Starting TEST 2: Running ${meta_name} with individual file input"
test2_output_dir="${tmp_dir}/test2"
mkdir -p "$test2_output_dir"

# Define individual FASTQ files
fastq_r1="$dest_fastq_dir/tinygex_S1_L001_R1_001.fastq.gz"
fastq_r2="$dest_fastq_dir/tinygex_S1_L001_R2_001.fastq.gz"

# Validate FASTQ files exist
check_file_exists "$fastq_r1" "R1 FASTQ file"
check_file_exists "$fastq_r2" "R2 FASTQ file"

log "Executing cellranger count with individual file input..."
"${meta_executable}" \
  --fastqs "$fastq_r1" \
  --fastqs "$fastq_r2" \
  --transcriptome "$dest_ref_dir" \
  --output "$test2_output_dir" \
  --no_secondary \
  --create_bam

log "Validating TEST 2 outputs..."
check_file_exists "$test2_output_dir/filtered_feature_bc_matrix.h5" "filtered feature matrix"
check_file_exists "$test2_output_dir/possorted_genome_bam.bam" "BAM file (should exist with --create_bam)"
check_dir_not_exists "$test2_output_dir/analysis" "analysis output directory (should be disabled with --no_secondary)"

log "✅ TEST 2 completed successfully"

print_test_summary "All cellranger_count tests"
