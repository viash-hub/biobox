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
  local has_multiallelics="$2"
  
  if [[ "$has_multiallelics" == "true" ]]; then
    cat > "$output_file" <<EOF
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr1,length=249250621,assembly=b37>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
chr1	1000	.	G	C,A	15	PASS	.	GT	0/1
chr1	2000	.	ATG	A,AT	20	PASS	.	GT	1/2
EOF
  else
    cat > "$output_file" <<EOF
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr1,length=249250621,assembly=b37>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
chr1	1000	.	G	C	15	PASS	.	GT	0/1
chr1	2000	.	A	T	20	PASS	.	GT	1/1
EOF
  fi
}

create_test_reference() {
  local output_file="$1"
  cat > "$output_file" <<EOF
>chr1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
EOF
}

# Create test data
log "Creating test VCF files and reference"
create_test_vcf "$meta_temp_dir/input.vcf" "false"
create_test_vcf "$meta_temp_dir/multiallelic.vcf" "true"
create_test_reference "$meta_temp_dir/reference.fa"

# Compress and index VCF files for bcftools
bgzip -c "$meta_temp_dir/input.vcf" > "$meta_temp_dir/input.vcf.gz"
bgzip -c "$meta_temp_dir/multiallelic.vcf" > "$meta_temp_dir/multiallelic.vcf.gz"
tabix -p vcf "$meta_temp_dir/input.vcf.gz"
tabix -p vcf "$meta_temp_dir/multiallelic.vcf.gz"

# Test 1: Basic normalization with atomize
log "Starting TEST 1: Basic normalization with atomize"
"$meta_executable" \
  --input "$meta_temp_dir/input.vcf.gz" \
  --output "$meta_temp_dir/output1.vcf" \
  --atomize

check_file_exists "$meta_temp_dir/output1.vcf" "basic normalization output"
check_file_not_empty "$meta_temp_dir/output1.vcf" "basic normalization output"
check_file_contains "$meta_temp_dir/output1.vcf" "chr1	1000" "first variant"
log "✅ TEST 1 completed successfully"

# Test 2: Atomize complex variants
log "Starting TEST 2: Atomize complex variants"
"$meta_executable" \
  --input "$meta_temp_dir/multiallelic.vcf.gz" \
  --output "$meta_temp_dir/output2.vcf" \
  --atomize

check_file_exists "$meta_temp_dir/output2.vcf" "atomized output"
check_file_not_empty "$meta_temp_dir/output2.vcf" "atomized output"
log "✅ TEST 2 completed successfully"

# Test 3: Split multiallelics
log "Starting TEST 3: Split multiallelics"
"$meta_executable" \
  --input "$meta_temp_dir/multiallelic.vcf.gz" \
  --output "$meta_temp_dir/output3.vcf" \
  --multiallelics "-both"

check_file_exists "$meta_temp_dir/output3.vcf" "split multiallelic output"
check_file_not_empty "$meta_temp_dir/output3.vcf" "split multiallelic output"
log "✅ TEST 3 completed successfully"

# Test 4: Check reference with fasta
log "Starting TEST 4: Check reference with fasta"
"$meta_executable" \
  --input "$meta_temp_dir/input.vcf.gz" \
  --output "$meta_temp_dir/output4.vcf" \
  --fasta_ref "$meta_temp_dir/reference.fa" \
  --check_ref "w"

check_file_exists "$meta_temp_dir/output4.vcf" "reference checked output"
check_file_not_empty "$meta_temp_dir/output4.vcf" "reference checked output"
log "✅ TEST 4 completed successfully"

# Test 5: Remove duplicates
log "Starting TEST 5: Remove duplicates"
"$meta_executable" \
  --input "$meta_temp_dir/input.vcf.gz" \
  --output "$meta_temp_dir/output5.vcf" \
  --rm_dup "exact"

check_file_exists "$meta_temp_dir/output5.vcf" "deduplicated output"
check_file_not_empty "$meta_temp_dir/output5.vcf" "deduplicated output"
log "✅ TEST 5 completed successfully"

# Test 6: Output format and compression
log "Starting TEST 6: Output format and compression"
"$meta_executable" \
  --input "$meta_temp_dir/input.vcf.gz" \
  --output "$meta_temp_dir/output6.vcf.gz" \
  --output_type "z" \
  --atomize

check_file_exists "$meta_temp_dir/output6.vcf.gz" "compressed output"
check_file_not_empty "$meta_temp_dir/output6.vcf.gz" "compressed output"
log "✅ TEST 6 completed successfully"

# Test 7: Regions filtering
log "Starting TEST 7: Regions filtering"
"$meta_executable" \
  --input "$meta_temp_dir/input.vcf.gz" \
  --output "$meta_temp_dir/output7.vcf" \
  --regions "chr1:1000-1500" \
  --atomize

check_file_exists "$meta_temp_dir/output7.vcf" "regions filtered output"
check_file_not_empty "$meta_temp_dir/output7.vcf" "regions filtered output"
check_file_contains "$meta_temp_dir/output7.vcf" "chr1	1000" "variant in region"
log "✅ TEST 7 completed successfully"

# Always end with summary
print_test_summary "All tests completed successfully"
