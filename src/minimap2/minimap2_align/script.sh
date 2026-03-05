#!/usr/bin/env bash

set -euo pipefail

## VIASH START
## VIASH END

echo "Running minimap2_align..."
echo "DEBUG: Output path is $par_output"
echo "DEBUG: BAM mode is $par_bam"

# Set number of threads
threads=${VIASH_META_CPUS:-2}

# Build minimap2 options
opts=()
[ ! -z "${par_preset:-}" ] && opts+=("-x" "$par_preset")
[ "${par_cigar_paf:-false}" = "true" ] && opts+=("-c")
[ "${par_cigar_bam:-false}" = "true" ] && opts+=("-L")

if [ "${par_bam:-false}" = "true" ]; then
  echo "Output format: Sorted BAM"
  # -a is required for SAM/BAM output
  minimap2 \
    -t "$threads" \
    -a \
    "${opts[@]}" \
    "$par_reference" \
    "$par_query" | \
    samtools sort -@ "$threads" -o "$par_output" -
  
  # Create index
  echo "Indexing BAM..."
  samtools index "$par_output"
else
  echo "Output format: PAF"
  minimap2 \
    -t "$threads" \
    "${opts[@]}" \
    "$par_reference" \
    "$par_query" \
    > "$par_output"
fi

echo "Alignment finished successfully."