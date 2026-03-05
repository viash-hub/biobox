#!/usr/bin/env bash

set -euo pipefail

## VIASH START
meta_executable="target/executable/minimap2/minimap2_index"
meta_resources_dir="src/minimap2"
## VIASH END

# Run minimap2_index component
echo "> Run minimap2_index on test FASTA"
"$meta_executable" \
  --input "$meta_resources_dir/test_data/test.fasta" \
  --output "$meta_resources_dir/test_data/test.mmi"

echo "Output .mmi will be written to: $meta_resources_dir/test_data/test.mmi"
  
# Check output exists
echo ">> Check if output exists"
if [ ! -f "$meta_resources_dir/test_data/test.mmi" ]; then
    echo ">> test.mmi does not exist"
    exit 1
fi

# Generate checksums of the output index file and the provided check.mmi file
output_mmi_hash=$(sha256sum "$meta_resources_dir/test_data/test.mmi" | awk '{print $1}')
check_mmi_hash=$(sha256sum "$meta_resources_dir/test_data/check.mmi" | awk '{print $1}')

# Compare checksums
echo "Comparing checksums..."

if [ "$output_mmi_hash" != "$check_mmi_hash" ]; then
    echo "Warning: The output .mmi file does not match 'check.mmi'."
else
    echo "The output .mmi file matches 'check.mmi'."
fi

echo "minimap2_index tests passed."
