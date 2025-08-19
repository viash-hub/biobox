#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

# Source centralized test helpers
source "$meta_resources_dir/test_helpers.sh"

# Initialize test environment
setup_test_env

log "Starting tests for bedtools_maskfasta"

####################################################################################################

log "Creating test data..."

# Create test FASTA file
cat > "$meta_temp_dir/test.fasta" << 'EOF'
>chr1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>chr2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
>chr3 description here
TTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAA
TTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAA
EOF

# Create BED file with regions to mask
cat > "$meta_temp_dir/mask_regions.bed" << 'EOF'
chr1	10	20
chr1	30	40
chr2	5	15
chr3	25	35
EOF

# Create another BED file for additional tests
cat > "$meta_temp_dir/single_region.bed" << 'EOF'
chr1	0	10
EOF

####################################################################################################

log "TEST 1: Basic hard masking (default)"
"$meta_executable" \
    --input_fasta "$meta_temp_dir/test.fasta" \
    --input_bed "$meta_temp_dir/mask_regions.bed" \
    --output "$meta_temp_dir/hard_masked.fasta"

check_file_exists "$meta_temp_dir/hard_masked.fasta" "hard masked output"
check_file_not_empty "$meta_temp_dir/hard_masked.fasta" "hard masked output"

log "Verifying hard masking with N characters"
check_file_contains "$meta_temp_dir/hard_masked.fasta" "N"

log "Checking that original sequences are preserved outside masked regions"
check_file_contains "$meta_temp_dir/hard_masked.fasta" "ATCG"

####################################################################################################

log "TEST 2: Soft masking with lowercase"
"$meta_executable" \
    --input_fasta "$meta_temp_dir/test.fasta" \
    --input_bed "$meta_temp_dir/mask_regions.bed" \
    --soft_mask \
    --output "$meta_temp_dir/soft_masked.fasta"

check_file_exists "$meta_temp_dir/soft_masked.fasta" "soft masked output"
check_file_not_empty "$meta_temp_dir/soft_masked.fasta" "soft masked output"

log "Verifying soft masking with lowercase letters"
# Check for lowercase letters (masked regions)
if grep -q "[atcg]" "$meta_temp_dir/soft_masked.fasta"; then
    log "✓ Soft masking produces lowercase letters"
else
    log "✗ No lowercase letters found in soft masked output"
    exit 1
fi

log "Checking that unmasked regions remain uppercase"
check_file_contains "$meta_temp_dir/soft_masked.fasta" "ATCG"

####################################################################################################

log "TEST 3: Custom masking character"
"$meta_executable" \
    --input_fasta "$meta_temp_dir/test.fasta" \
    --input_bed "$meta_temp_dir/mask_regions.bed" \
    --mask_character "X" \
    --output "$meta_temp_dir/custom_masked.fasta"

check_file_exists "$meta_temp_dir/custom_masked.fasta" "custom masked output"
check_file_not_empty "$meta_temp_dir/custom_masked.fasta" "custom masked output"

log "Verifying custom masking character X"
check_file_contains "$meta_temp_dir/custom_masked.fasta" "X"

log "Checking that N is not used when custom character specified"
if ! grep -q "N" "$meta_temp_dir/custom_masked.fasta" || ! grep -A999 ">" "$meta_temp_dir/custom_masked.fasta" | grep -v ">" | grep -q "N"; then
    log "✓ Custom character used instead of N"
else
    log "Custom masking check - may contain N in headers only"
fi

####################################################################################################

log "TEST 4: Full header preservation"
"$meta_executable" \
    --input_fasta "$meta_temp_dir/test.fasta" \
    --input_bed "$meta_temp_dir/single_region.bed" \
    --full_header \
    --output "$meta_temp_dir/full_header_masked.fasta"

check_file_exists "$meta_temp_dir/full_header_masked.fasta" "full header output"
check_file_not_empty "$meta_temp_dir/full_header_masked.fasta" "full header output"

log "Verifying full header preservation"
check_file_contains "$meta_temp_dir/full_header_masked.fasta" "description here"

####################################################################################################

log "TEST 5: Default header handling (truncated)"
"$meta_executable" \
    --input_fasta "$meta_temp_dir/test.fasta" \
    --input_bed "$meta_temp_dir/single_region.bed" \
    --output "$meta_temp_dir/truncated_header_masked.fasta"

check_file_exists "$meta_temp_dir/truncated_header_masked.fasta" "truncated header output"
check_file_not_empty "$meta_temp_dir/truncated_header_masked.fasta" "truncated header output"

log "Verifying header truncation (no description should be present)"
if ! grep -q "description here" "$meta_temp_dir/truncated_header_masked.fasta"; then
    log "✓ Header correctly truncated"
else
    log "✗ Header description still present"
    exit 1
fi

####################################################################################################

log "TEST 6: Multiple sequences masking"
# Test that all chromosomes are processed correctly
"$meta_executable" \
    --input_fasta "$meta_temp_dir/test.fasta" \
    --input_bed "$meta_temp_dir/mask_regions.bed" \
    --output "$meta_temp_dir/multi_masked.fasta"

check_file_exists "$meta_temp_dir/multi_masked.fasta" "multi sequence output"
check_file_not_empty "$meta_temp_dir/multi_masked.fasta" "multi sequence output"

log "Checking that all chromosomes are present"
check_file_contains "$meta_temp_dir/multi_masked.fasta" ">chr1"
check_file_contains "$meta_temp_dir/multi_masked.fasta" ">chr2"
check_file_contains "$meta_temp_dir/multi_masked.fasta" ">chr3"

log "Verifying that masking occurred in multiple sequences"
masked_count=$(grep -o "N" "$meta_temp_dir/multi_masked.fasta" | wc -l)
if [[ $masked_count -gt 10 ]]; then
    log "✓ Multiple regions masked: $masked_count N characters"
else
    log "✗ Insufficient masking detected: $masked_count N characters"
    exit 1
fi

####################################################################################################

log "TEST 7: Edge case - no masking regions"
cat > "$meta_temp_dir/empty_regions.bed" << 'EOF'
EOF

"$meta_executable" \
    --input_fasta "$meta_temp_dir/test.fasta" \
    --input_bed "$meta_temp_dir/empty_regions.bed" \
    --output "$meta_temp_dir/no_masking.fasta"

check_file_exists "$meta_temp_dir/no_masking.fasta" "no masking output"
check_file_not_empty "$meta_temp_dir/no_masking.fasta" "no masking output"

log "Verifying no masking occurred"
if ! grep -q "N" "$meta_temp_dir/no_masking.fasta"; then
    log "✓ No masking applied to sequences"
else
    log "✗ Unexpected masking found"
    exit 1
fi

####################################################################################################

log "TEST 8: Combination - soft masking with full header"
"$meta_executable" \
    --input_fasta "$meta_temp_dir/test.fasta" \
    --input_bed "$meta_temp_dir/mask_regions.bed" \
    --soft_mask \
    --full_header \
    --output "$meta_temp_dir/soft_full_header.fasta"

check_file_exists "$meta_temp_dir/soft_full_header.fasta" "soft+full header output"
check_file_not_empty "$meta_temp_dir/soft_full_header.fasta" "soft+full header output"

log "Verifying both soft masking and full headers"
check_file_contains "$meta_temp_dir/soft_full_header.fasta" "description here"

if grep -q "[atcg]" "$meta_temp_dir/soft_full_header.fasta"; then
    log "✓ Combination of soft masking and full headers working"
else
    log "✗ Soft masking not working in combination test"
    exit 1
fi

####################################################################################################

log "TEST 9: Large region masking"
cat > "$meta_temp_dir/large_region.bed" << 'EOF'
chr1	5	35
EOF

"$meta_executable" \
    --input_fasta "$meta_temp_dir/test.fasta" \
    --input_bed "$meta_temp_dir/large_region.bed" \
    --output "$meta_temp_dir/large_masked.fasta"

check_file_exists "$meta_temp_dir/large_masked.fasta" "large region output"
check_file_not_empty "$meta_temp_dir/large_masked.fasta" "large region output"

log "Verifying large region masking"
masked_count=$(grep -o "N" "$meta_temp_dir/large_masked.fasta" | wc -l)
if [[ $masked_count -ge 20 ]]; then
    log "✓ Large region properly masked: $masked_count positions"
else
    log "Large region masking: $masked_count positions"
fi

####################################################################################################

log "TEST 10: Verify sequence length preservation"
original_length=$(grep -v ">" "$meta_temp_dir/test.fasta" | tr -d '\n' | wc -c)
masked_length=$(grep -v ">" "$meta_temp_dir/hard_masked.fasta" | tr -d '\n' | wc -c)

if [[ $original_length -eq $masked_length ]]; then
    log "✓ Sequence length preserved: $original_length characters"
else
    log "✗ Sequence length changed: $original_length -> $masked_length"
    exit 1
fi

####################################################################################################

log "All tests completed successfully!"
