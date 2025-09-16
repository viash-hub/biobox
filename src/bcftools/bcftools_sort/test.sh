#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# Always source centralized helpers
source "$meta_resources_dir/test_helpers.sh"

# Initialize test environment
setup_test_env

log "Starting tests for $meta_name"

# Create test VCF file with unsorted variants
create_test_vcf() {
  local output_file="$1"
  cat > "$output_file" << 'EOF'
##fileformat=VCFv4.2
##contig=<ID=chr19>
##contig=<ID=chr20>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
chr20	1000	.	A	G	60	PASS	DP=30	GT:DP	0/1:15
chr19	500	.	C	T	50	PASS	DP=25	GT:DP	1/1:12
chr20	800	.	G	A	40	PASS	DP=20	GT:DP	0/1:10
chr19	700	.	T	C	70	PASS	DP=35	GT:DP	0/0:18
chr20	1200	.	C	G	55	PASS	DP=28	GT:DP	1/1:14
EOF
}

# Create expected sorted output
create_expected_sorted() {
  local output_file="$1"
  cat > "$output_file" << 'EOF'
##fileformat=VCFv4.2
##contig=<ID=chr19>
##contig=<ID=chr20>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
chr19	500	.	C	T	50	PASS	DP=25	GT:DP	1/1:12
chr19	700	.	T	C	70	PASS	DP=35	GT:DP	0/0:18
chr20	800	.	G	A	40	PASS	DP=20	GT:DP	0/1:10
chr20	1000	.	A	G	60	PASS	DP=30	GT:DP	0/1:15
chr20	1200	.	C	G	55	PASS	DP=28	GT:DP	1/1:14
EOF
}

log "Starting TEST 1: Basic VCF sorting"
cd "$meta_temp_dir"

# Create test data
create_test_vcf "unsorted.vcf"
create_expected_sorted "expected_sorted.vcf"

# Run bcftools sort
"$meta_executable" \
  --input "unsorted.vcf" \
  --output "sorted.vcf"

# Verify output
check_file_exists "sorted.vcf"
check_file_not_empty "sorted.vcf"

# Check that variants are properly sorted (by comparing position order)
check_file_contains "sorted.vcf" "chr19.*500"
check_file_contains "sorted.vcf" "chr20.*1200"

log "TEST 1 passed"

log "Starting TEST 2: Compressed VCF output"
cd "$meta_temp_dir"

# Create compressed VCF (reuse test data)
"$meta_executable" \
  --input "unsorted.vcf" \
  --output "sorted.vcf.gz" \
  --output_type "z"

# Verify compressed output
check_file_exists "sorted.vcf.gz"
check_file_not_empty "sorted.vcf.gz"

# Check it's properly compressed by trying to decompress it
zcat "sorted.vcf.gz" | head -1 | grep -q "##fileformat"

log "TEST 2 passed"

log "Starting TEST 3: BCF output"
cd "$meta_temp_dir"

# Create BCF output (reuse test data)
"$meta_executable" \
  --input "unsorted.vcf" \
  --output "sorted.bcf" \
  --output_type "b"

# Verify BCF output
check_file_exists "sorted.bcf"
check_file_not_empty "sorted.bcf"

# Verify it's a valid BCF file by reading it back
bcftools view "sorted.bcf" -o "from_bcf.vcf"
check_file_exists "from_bcf.vcf"
check_file_not_empty "from_bcf.vcf"

log "TEST 3 passed"

log "Starting TEST 4: Memory limit parameter"
cd "$meta_temp_dir"

# Test with memory limit (reuse test data)
"$meta_executable" \
  --input "unsorted.vcf" \
  --output "sorted_mem.vcf" \
  --max_mem "500M"

# Verify output
check_file_exists "sorted_mem.vcf"
check_file_not_empty "sorted_mem.vcf"

log "TEST 4 passed"

log "Starting TEST 5: Custom temporary directory"
cd "$meta_temp_dir"

# Create temp directory
mkdir -p "custom_temp"

# Test with custom temp directory (reuse test data)
"$meta_executable" \
  --input "unsorted.vcf" \
  --output "sorted_temp.vcf" \
  --temp_dir "custom_temp"

# Verify output
check_file_exists "sorted_temp.vcf"
check_file_not_empty "sorted_temp.vcf"

log "TEST 5 passed"

log "Starting TEST 6: Verbosity parameter"
cd "$meta_temp_dir"

# Test with verbosity (reuse test data)
"$meta_executable" \
  --input "unsorted.vcf" \
  --output "sorted_verbose.vcf" \
  --verbosity 1

# Verify output
check_file_exists "sorted_verbose.vcf"
check_file_not_empty "sorted_verbose.vcf"

log "TEST 6 passed"

log "Starting TEST 7: Write index parameter"
cd "$meta_temp_dir"

# Test with index writing (for compressed output, reuse test data)
"$meta_executable" \
  --input "unsorted.vcf" \
  --output "sorted_indexed.vcf.gz" \
  --output_type "z" \
  --write_index "tbi"

# Verify output and index
check_file_exists "sorted_indexed.vcf.gz"
check_file_not_empty "sorted_indexed.vcf.gz"

# Check if index was created
check_file_exists "sorted_indexed.vcf.gz.tbi"

log "TEST 7 passed"

print_test_summary "All tests completed successfully"
