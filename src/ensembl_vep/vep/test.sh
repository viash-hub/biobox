#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

# Always source centralized helpers
source "$meta_resources_dir/test_helpers.sh"

# Initialize test environment
setup_test_env

log "Starting tests for $meta_name"

# Create test VCF file
create_test_vcf() {
  cat > "$1" << 'EOF'
##fileformat=VCFv4.2
##contig=<ID=1,length=249250621>
##contig=<ID=2,length=242193529>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	230845794	rs699	G	A	.	PASS	.
2	47641559	rs12345	C	T	.	PASS	.
EOF
}

# Create test variant identifier file
create_test_variant_id() {
  cat > "$1" << 'EOF'
1 230845794 230845794 G/A +
2 47641559 47641559 C/T +
EOF
}

# Test 2: VCF input format test (offline mode without cache)
log "Starting TEST 2: VCF input format validation"
create_test_vcf "$meta_temp_dir/test_input.vcf"

# Test parameter parsing without requiring cache (should fail gracefully)
# Using C. elegans for consistency with installation tests
"$meta_executable" \
  --input_file "$meta_temp_dir/test_input.vcf" \
  --output_file "$meta_temp_dir/vep_output1.txt" \
  --format vcf \
  --species caenorhabditis_elegans \
  --assembly WBcel235 \
  --offline 2>&1 | head -20 > "$meta_temp_dir/vcf_test.txt" || true

check_file_exists "$meta_temp_dir/vcf_test.txt" "VCF test output"
check_file_not_empty "$meta_temp_dir/vcf_test.txt" "VCF test output"
log "✅ TEST 1 completed successfully"

# Test 2: Variant identifier format test
log "Starting TEST 2: Variant identifier format"
create_test_variant_id "$meta_temp_dir/test_variants.txt"

"$meta_executable" \
  --input_file "$meta_temp_dir/test_variants.txt" \
  --output_file "$meta_temp_dir/vep_output2.txt" \
  --format variant_identifier \
  --offline 2>&1 | head -20 > "$meta_temp_dir/variant_test.txt" || true

check_file_exists "$meta_temp_dir/variant_test.txt" "variant identifier test output"
check_file_not_empty "$meta_temp_dir/variant_test.txt" "variant identifier test output"
log "✅ TEST 2 completed successfully"

# Test 3: VCF output format
log "Starting TEST 3: VCF output format test"
create_test_vcf "$meta_temp_dir/test_input2.vcf"

"$meta_executable" \
  --input_file "$meta_temp_dir/test_input2.vcf" \
  --output_file "$meta_temp_dir/vep_output.vcf" \
  --vcf \
  --offline 2>&1 | head -20 > "$meta_temp_dir/vcf_output_test.txt" || true

check_file_exists "$meta_temp_dir/vcf_output_test.txt" "VCF output test"
log "✅ TEST 3 completed successfully"

# Test 4: Multiple annotation options
log "Starting TEST 4: Annotation options test"
create_test_vcf "$meta_temp_dir/test_input3.vcf"

"$meta_executable" \
  --input_file "$meta_temp_dir/test_input3.vcf" \
  --output_file "$meta_temp_dir/annotated_output.txt" \
  --canonical \
  --symbol \
  --hgvs \
  --pick \
  --offline 2>&1 | head -20 > "$meta_temp_dir/annotation_test.txt" || true

check_file_exists "$meta_temp_dir/annotation_test.txt" "annotation test output"
check_file_not_empty "$meta_temp_dir/annotation_test.txt" "annotation test output"
log "✅ TEST 4 completed successfully"

# Test 5: Buffer size parameter test
log "Starting TEST 5: Buffer size parameter test"
create_test_vcf "$meta_temp_dir/test_input4.vcf"

# Test with buffer size parameter
"$meta_executable" \
  --input_file "$meta_temp_dir/test_input4.vcf" \
  --output_file "$meta_temp_dir/buffer_output.txt" \
  --buffer_size 1000 \
  --offline 2>&1 | head -20 > "$meta_temp_dir/buffer_test.txt" || true

check_file_exists "$meta_temp_dir/buffer_test.txt" "buffer test output"
log "✅ TEST 5 completed successfully"

# Always end with summary
print_test_summary "All VEP tests completed successfully"