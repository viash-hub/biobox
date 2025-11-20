#!/bin/bash

set -eo pipefail

## VIASH START
meta_executable="target/native/spaceranger/spaceranger_count/spaceranger_count"
meta_resources_dir="src/spaceranger/spaceranger_count"
## VIASH END

test_data="$meta_resources_dir/test_data"

echo "> Default test run"
"$meta_executable" \
    --id test_spaceranger \
    --transcriptome "$test_data/GRCh38" \
    --fastqs "$test_data/subsampled_fastqs" \
    --probe_set "$test_data/Visium_HD_Human_Pancreas_probe_set.csv" \
    --cytaimage "$test_data/Visium_HD_Human_Pancreas_image.tif" \
    --image "$test_data/Visium_HD_Human_Pancreas_tissue_image.btf" \
    --create_bam false \
    --lanes 1 \
    --slide H1-HBNMBMC \
    --area D1 \
    --nosecondary \
    ---cpus 4 \
    ---memory 18GB


echo "> Checking outputs..."

# Define output directory
OUT_DIR="test_spaceranger/outs"

# Function to check if file exists and is non-empty
check_file() {
    local file=$1
    local description=$2
    echo -n "Checking $description... "
    if [ ! -f "$file" ]; then
        echo "FAIL (file not found)"
        exit 1
    elif [ ! -s "$file" ]; then
        echo "FAIL (file is empty)"
        exit 1
    else
        echo "OK"
    fi
}

# Check essential files
check_file "$OUT_DIR/web_summary.html" "web summary"
check_file "$OUT_DIR/metrics_summary.csv" "metrics summary"

# Check per resolution outputs
for res in "002" "008" "016"; do
    bin_dir="$OUT_DIR/binned_outputs/square_${res}um"
    check_file "$bin_dir/filtered_feature_bc_matrix/matrix.mtx.gz" "${res}um matrix"
done

echo "> All tests passed successfully!"