#!/usr/bin/env bash

set -euo pipefail

## VIASH START
meta_executable="target/executable/minimap2/minimap2_align"
meta_resources_dir="src/minimap2/minimap2_align"
## VIASH END

# Test PAF output
echo ">> Testing PAF output..."
TEST_DATA="$meta_resources_dir/test_data"
"$meta_executable" \
  --reference "$TEST_DATA/ref.fasta" \
  --query "$TEST_DATA/reads.fastq" \
  --output "$TEST_DATA/test.paf"

if [ ! -s "$TEST_DATA/test.paf" ]; then echo "PAF file is empty"; exit 1; fi
grep -q "seq1" "$TEST_DATA/test.paf" || { echo "PAF alignment failed"; exit 1; }

# Test BAM output
echo ">> Testing BAM output..."
"$meta_executable" \
  --reference "$TEST_DATA/ref.fasta" \
  --query "$TEST_DATA/reads.fastq" \
  --bam true \
  --output "$TEST_DATA/test.bam"

if [ ! -s "$TEST_DATA/test.bam" ]; then echo "BAM file is empty"; exit 1; fi
# Check if BAM index was created
if [ ! -f "$TEST_DATA/test.bam.bai" ]; then echo "BAM index missing"; exit 1; fi

# Check BAM content using samtools
samtools view "$TEST_DATA/test.bam" | grep -q "seq1" || { echo "BAM alignment failed"; exit 1; }

echo "minimap2_align tests passed."
