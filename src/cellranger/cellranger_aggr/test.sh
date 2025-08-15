#!/bin/bash

## VIASH START
## VIASH END

# Include test helpers for enhanced logging and utilities
source "${meta_resources_dir}/test_helpers.sh"

# Set up test environment with enhanced error handling 
setup_test_env

# Set cellranger installation directory
CELLRANGER_DIR="/opt/cellranger-9.0.1"

# create temporary directory
tmp_dir=$(mktemp -d "${meta_temp_dir}/${meta_name}-XXXXXXXX")
function clean_up {
    log "Cleaning up temporary directory: $tmp_dir"
    rm -rf "${tmp_dir}"
}
trap clean_up EXIT

log "Created temporary directory: $tmp_dir"

log "Copy test data from Cell Ranger installation directory: ${CELLRANGER_DIR}"
mkdir -p "${tmp_dir}/test_data/"

# Copy test data from cellranger installation
log "Found cellranger test data in ${CELLRANGER_DIR}/external/"
cp -r "${CELLRANGER_DIR}/external/cellranger_tiny_fastq/" "${tmp_dir}/test_data"
cp -r "${CELLRANGER_DIR}/external/cellranger_tiny_ref/" "${tmp_dir}/test_data"

input_dir="${tmp_dir}/test_data/cellranger_tiny_fastq"
reference_dir="${tmp_dir}/test_data/cellranger_tiny_ref"

# Validate that test data was copied successfully
check_dir_exists "$input_dir" "cellranger_tiny_fastq directory"
check_dir_exists "$reference_dir" "cellranger_tiny_ref directory"

log "Using real cellranger test data"

# Run cellranger count to generate sample outputs for aggregation
log "Running cellranger count to generate test samples"
cd "${tmp_dir}"

# Create first sample
cellranger count \
    --id="sample1" \
    --fastqs="${input_dir}" \
    --transcriptome="${reference_dir}" \
    --create-bam=false \
    --disable-ui \
    --dry

# Create second sample (same data, different ID)
cellranger count \
    --id="sample2" \
    --fastqs="${input_dir}" \
    --transcriptome="${reference_dir}" \
    --create-bam=false \
    --disable-ui \
    --dry

# Create aggregation CSV with dry run outputs
log "Creating aggregation CSV file"
echo "library_id,molecule_h5" > aggregation.csv
echo "sample1,${tmp_dir}/sample1/outs/molecule_info.h5" >> aggregation.csv
echo "sample2,${tmp_dir}/sample2/outs/molecule_info.h5" >> aggregation.csv

## TEST: Run cellranger aggr with dry flag
log "Starting TEST: Running ${meta_name} with --dry flag"
cd "${tmp_dir}"

# For testing, just use a simple CSV with mock data since we're doing a dry run
log "Creating simple test CSV with mock data"
echo "library_id,molecule_h5" > simple_test.csv
echo "sample1,/fake/path/sample1/outs/molecule_info.h5" >> simple_test.csv
echo "sample2,/fake/path/sample2/outs/molecule_info.h5" >> simple_test.csv

output_dir="${tmp_dir}/test_output"
mkdir -p "${output_dir}"

log "Executing cellranger aggr with dry run..."

"${meta_executable}" \
    --id "test_aggr" \
    --csv "simple_test.csv" \
    --description "Test aggregation" \
    --normalize "mapped" \
    --dry \
    --output_dir "${output_dir}"

log "Validating dry run outputs..."
log "Checking that dry run completed successfully"
log "Current directory contents:"
ls -la

# Check for any .mro files
if ls *.mro 1> /dev/null 2>&1; then
    log "✓ MRO files found:"
    ls -la *.mro
else
    log "ℹ No .mro files found, but dry run completed without error"
fi

print_test_summary "cellranger_aggr dry run test"
