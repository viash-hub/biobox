#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# Always source centralized helpers
source "$meta_resources_dir/test_helpers.sh"

# Initialize test environment
setup_test_env

log "Starting tests for $meta_name"

# Create test VCF files using helpers
create_test_vcf() {
  local output_file="$1"
  local chrom="$2"
  local start_pos="$3"
  local end_pos="$4"
  
  cat > "$output_file" <<EOF
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=$chrom,length=249250621,assembly=b37>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
$chrom	$start_pos	.	G	C	15	PASS	.	GT	0/1
$chrom	$end_pos	.	A	T	20	PASS	.	GT	1/1
EOF
}

# Create test data
log "Creating test VCF files"
create_test_vcf "$meta_temp_dir/input1.vcf" "chr1" 1000 2000
create_test_vcf "$meta_temp_dir/input2.vcf" "chr1" 3000 4000

# Compress and index VCF files for bcftools
bgzip -c "$meta_temp_dir/input1.vcf" > "$meta_temp_dir/input1.vcf.gz"
bgzip -c "$meta_temp_dir/input2.vcf" > "$meta_temp_dir/input2.vcf.gz"
tabix -p vcf "$meta_temp_dir/input1.vcf.gz"
tabix -p vcf "$meta_temp_dir/input2.vcf.gz"

# Create file list
echo "$meta_temp_dir/input1.vcf.gz" > "$meta_temp_dir/file_list.txt"
echo "$meta_temp_dir/input2.vcf.gz" >> "$meta_temp_dir/file_list.txt"

# Test 1: Basic concatenation
log "Starting TEST 1: Basic concatenation"
"$meta_executable" \
  --input "$meta_temp_dir/input1.vcf.gz" \
  --input "$meta_temp_dir/input2.vcf.gz" \
  --output "$meta_temp_dir/output1.vcf"

check_file_exists "$meta_temp_dir/output1.vcf" "basic concatenation output"
check_file_not_empty "$meta_temp_dir/output1.vcf" "basic concatenation output"
check_file_contains "$meta_temp_dir/output1.vcf" "chr1	1000" "first variant"
check_file_contains "$meta_temp_dir/output1.vcf" "chr1	4000" "last variant"
log "✅ TEST 1 completed successfully"

# Test 2: File list input
log "Starting TEST 2: File list input"
"$meta_executable" \
  --file_list "$meta_temp_dir/file_list.txt" \
  --output "$meta_temp_dir/output2.vcf"

check_file_exists "$meta_temp_dir/output2.vcf" "file list output"
check_file_not_empty "$meta_temp_dir/output2.vcf" "file list output"
check_file_contains "$meta_temp_dir/output2.vcf" "chr1	1000" "first variant from file list"
log "✅ TEST 2 completed successfully"

# Test 3: Allow overlaps and output type
log "Starting TEST 3: Allow overlaps with compressed output"
"$meta_executable" \
  --input "$meta_temp_dir/input1.vcf.gz" \
  --input "$meta_temp_dir/input2.vcf.gz" \
  --output "$meta_temp_dir/output3.vcf.gz" \
  --allow_overlaps \
  --output_type "z"

check_file_exists "$meta_temp_dir/output3.vcf.gz" "compressed output"
check_file_not_empty "$meta_temp_dir/output3.vcf.gz" "compressed output"
log "✅ TEST 3 completed successfully"

# Test 4: Remove duplicates
log "Starting TEST 4: Remove duplicates"
"$meta_executable" \
  --input "$meta_temp_dir/input1.vcf.gz" \
  --input "$meta_temp_dir/input1.vcf.gz" \
  --output "$meta_temp_dir/output4.vcf" \
  --allow_overlaps \
  --rm_dups "exact"

check_file_exists "$meta_temp_dir/output4.vcf" "deduplicated output"
check_file_not_empty "$meta_temp_dir/output4.vcf" "deduplicated output"
log "✅ TEST 4 completed successfully"

# Test 5: Naive concatenation
log "Starting TEST 5: Naive concatenation"
"$meta_executable" \
  --input "$meta_temp_dir/input1.vcf.gz" \
  --input "$meta_temp_dir/input2.vcf.gz" \
  --output "$meta_temp_dir/output5.vcf" \
  --naive

check_file_exists "$meta_temp_dir/output5.vcf" "naive concatenation output"
check_file_not_empty "$meta_temp_dir/output5.vcf" "naive concatenation output"
log "✅ TEST 5 completed successfully"

# Test 6: Drop genotypes
log "Starting TEST 6: Drop genotypes"
"$meta_executable" \
  --input "$meta_temp_dir/input1.vcf.gz" \
  --input "$meta_temp_dir/input2.vcf.gz" \
  --output "$meta_temp_dir/output6.vcf" \
  --drop_genotypes

check_file_exists "$meta_temp_dir/output6.vcf" "genotype-free output"
check_file_not_empty "$meta_temp_dir/output6.vcf" "genotype-free output"
log "✅ TEST 6 completed successfully"

# Test 7: Regions filtering
log "Starting TEST 7: Regions filtering"
"$meta_executable" \
  --input "$meta_temp_dir/input1.vcf.gz" \
  --input "$meta_temp_dir/input2.vcf.gz" \
  --output "$meta_temp_dir/output7.vcf" \
  --allow_overlaps \
  --regions "chr1:1000-2000"

check_file_exists "$meta_temp_dir/output7.vcf" "regions filtered output"
check_file_not_empty "$meta_temp_dir/output7.vcf" "regions filtered output"
check_file_contains "$meta_temp_dir/output7.vcf" "chr1	1000" "variant in region"
log "✅ TEST 7 completed successfully"

# Always end with summary
print_test_summary "All tests completed successfully"

