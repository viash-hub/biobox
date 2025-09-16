#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

# Simple helper functions for testing
check_file_exists() {
  [[ -f "$1" ]] || (echo "Error: File '$1' does not exist" && exit 1)
}

check_file_not_empty() {
  [[ -s "$1" ]] || (echo "Error: File '$1' is empty" && exit 1)
}

check_file_contains() {
  grep -q "$2" "$1" || (echo "Error: File '$1' does not contain '$2'" && exit 1)
}

echo "Starting tests for $meta_name"

# Create test VCF file with some basic data
cat > "$meta_temp_dir/test.vcf" << 'EOF'
##fileformat=VCFv4.3
##contig=<ID=chr1,length=248956422>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1	sample2	sample3
chr1	1000	.	A	T	60	PASS	DP=100;AF=0.5	GT:DP	0/1:30	1/1:35	0/0:35
chr1	2000	rs123	G	C	45	PASS	DP=80;AF=0.3	GT:DP	0/1:25	0/0:30	0/1:25
chr1	3000	.	C	G	50	PASS	DP=90;AF=0.4	GT:DP	1/1:30	0/1:30	0/0:30
chr1	4000	.	T	A	30	PASS	DP=60;AF=0.2	GT:DP	0/1:20	0/0:20	0/0:20
chr1	5000	rs456	A	G	40	PASS	DP=70;AF=0.6	GT:DP	1/1:25	0/1:25	1/1:20
EOF

# Compress and index the VCF
bgzip -c "$meta_temp_dir/test.vcf" > "$meta_temp_dir/test.vcf.gz"
tabix -p vcf "$meta_temp_dir/test.vcf.gz"

# Test 1: Basic functionality
echo "TEST 1: Basic functionality"
"$meta_executable" \
  --input "$meta_temp_dir/test.vcf" \
  --output "$meta_temp_dir/output1.txt"

check_file_exists "$meta_temp_dir/output1.txt"
check_file_not_empty "$meta_temp_dir/output1.txt"
check_file_contains "$meta_temp_dir/output1.txt" "SN.*number of samples"
echo "✅ TEST 1 completed successfully"

# Test 2: Compressed input with verbose
echo "TEST 2: Compressed input with verbose"
"$meta_executable" \
  --input "$meta_temp_dir/test.vcf.gz" \
  --output "$meta_temp_dir/output2.txt" \
  --verbose

check_file_exists "$meta_temp_dir/output2.txt"
check_file_not_empty "$meta_temp_dir/output2.txt"
check_file_contains "$meta_temp_dir/output2.txt" "SN.*number of samples"
echo "✅ TEST 2 completed successfully"

# Test 3: With regions and split by ID
echo "TEST 3: Regions and split by ID"
"$meta_executable" \
  --input "$meta_temp_dir/test.vcf.gz" \
  --output "$meta_temp_dir/output3.txt" \
  --regions "chr1:1000-3000" \
  --split_by_id

check_file_exists "$meta_temp_dir/output3.txt"
check_file_not_empty "$meta_temp_dir/output3.txt"
check_file_contains "$meta_temp_dir/output3.txt" "SN.*number of records"
echo "✅ TEST 3 completed successfully"

# Test 4: Advanced options
echo "TEST 4: Advanced options"
"$meta_executable" \
  --input "$meta_temp_dir/test.vcf" \
  --output "$meta_temp_dir/output4.txt" \
  --af_bins "0.1,0.3,0.5,0.7,0.9" \
  --depth "0,100,10" \
  --collapse "snps"

check_file_exists "$meta_temp_dir/output4.txt"
check_file_not_empty "$meta_temp_dir/output4.txt"
check_file_contains "$meta_temp_dir/output4.txt" "SN.*number of SNPs"
echo "✅ TEST 4 completed successfully"

# Test 5: Sample filtering
echo "TEST 5: Sample filtering"
"$meta_executable" \
  --input "$meta_temp_dir/test.vcf.gz" \
  --output "$meta_temp_dir/output5.txt" \
  --samples "sample1,sample2" \
  --include "QUAL >= 40"

check_file_exists "$meta_temp_dir/output5.txt"
check_file_not_empty "$meta_temp_dir/output5.txt"
check_file_contains "$meta_temp_dir/output5.txt" "SN.*number of samples"
echo "✅ TEST 5 completed successfully"

# Test 6: User-defined Ts/Tv
echo "TEST 6: User-defined Ts/Tv"
"$meta_executable" \
  --input "$meta_temp_dir/test.vcf" \
  --output "$meta_temp_dir/output6.txt" \
  --user_tstv "DP:0:100:10" \
  --first_allele_only

check_file_exists "$meta_temp_dir/output6.txt"
check_file_not_empty "$meta_temp_dir/output6.txt"
check_file_contains "$meta_temp_dir/output6.txt" "SN.*number of records"
echo "✅ TEST 6 completed successfully"

# Test 7: Targets vs regions
echo "TEST 7: Targets functionality"
"$meta_executable" \
  --input "$meta_temp_dir/test.vcf.gz" \
  --output "$meta_temp_dir/output7.txt" \
  --targets "chr1:2000-4000" \
  --targets_overlap "pos"

check_file_exists "$meta_temp_dir/output7.txt"
check_file_not_empty "$meta_temp_dir/output7.txt"
check_file_contains "$meta_temp_dir/output7.txt" "SN.*number of records"
echo "✅ TEST 7 completed successfully"

echo "All bcftools_stats tests completed successfully!"


