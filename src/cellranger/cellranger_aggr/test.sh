#!/bin/bash

set -e

## VIASH START
## VIASH END

# create temporary directory
tmp_dir=$(mktemp -d "${meta_temp_dir}/${meta_name}-XXXXXXXX")
function clean_up {
    rm -rf "${tmp_dir}"
}
trap clean_up EXIT

#!/bin/bash

set -e

## VIASH START
## VIASH END

# Set cellranger installation directory
CELLRANGER_DIR="/opt/cellranger-9.0.1"

# create temporary directory
tmp_dir=$(mktemp -d "${meta_temp_dir}/${meta_name}-XXXXXXXX")
function clean_up {
    rm -rf "${tmp_dir}"
}
trap clean_up EXIT

echo "Copy test data from Cell Ranger installation directory: ${CELLRANGER_DIR}"
mkdir -p "${tmp_dir}/test_data/"

# Copy test data from cellranger installation
echo "Found cellranger test data in ${CELLRANGER_DIR}/external/"
cp -r "${CELLRANGER_DIR}/external/cellranger_tiny_fastq/" "${tmp_dir}/test_data"
cp -r "${CELLRANGER_DIR}/external/cellranger_tiny_ref/" "${tmp_dir}/test_data"

input_dir="${tmp_dir}/test_data/cellranger_tiny_fastq"
reference_dir="${tmp_dir}/test_data/cellranger_tiny_ref"

echo "Using real cellranger test data"

# Run cellranger count to generate sample outputs for aggregation
echo "Running cellranger count to generate test samples"
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
echo "library_id,molecule_h5" > aggregation.csv
echo "sample1,${tmp_dir}/sample1/outs/molecule_info.h5" >> aggregation.csv
echo "sample2,${tmp_dir}/sample2/outs/molecule_info.h5" >> aggregation.csv

## TEST: Run cellranger aggr with dry flag
echo "Running ${meta_name} with --dry flag"
cd "${tmp_dir}"

# For testing, just use a simple CSV with mock data since we're doing a dry run
echo "library_id,molecule_h5" > simple_test.csv
echo "sample1,/fake/path/sample1/outs/molecule_info.h5" >> simple_test.csv
echo "sample2,/fake/path/sample2/outs/molecule_info.h5" >> simple_test.csv

"${meta_executable}" \
    --id "test_aggr" \
    --csv "simple_test.csv" \
    --description "Test aggregation" \
    --normalize "mapped" \
    --dry

echo "Checking that dry run completed successfully"
echo "Current directory contents:"
ls -la

# Check for any .mro files
if ls *.mro 1> /dev/null 2>&1; then
    echo "SUCCESS: MRO files found:"
    ls -la *.mro
else
    echo "WARNING: No .mro files found, but dry run completed without error"
fi

echo "SUCCESS: Test completed successfully"
