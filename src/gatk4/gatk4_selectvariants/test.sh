#!/bin/bash

## VIASH START
## VIASH END

# Source the centralized test helpers
source "$meta_resources_dir/test_helpers.sh"

# Initialize test environment with strict error handling
setup_test_env

#############################################
# Test execution with centralized functions
#############################################

log "Starting tests for $meta_name"

# Create test data directory
test_data_dir="$meta_temp_dir/test_data"
mkdir -p "$test_data_dir"

# --- Build a small VCF with 2 SNPs and 2 indels ---
log "Writing hand-crafted test VCF with SNP and indel records..."
cat > "$test_data_dir/variants.vcf" << 'EOF'
##fileformat=VCFv4.2
##contig=<ID=seq1,length=1000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1
seq1	100	snp1	A	G	50	PASS	DP=20	GT:DP	0/1:20
seq1	200	snp2	C	T	60	PASS	DP=25	GT:DP	1/1:25
seq1	300	indel1	A	ATG	55	PASS	DP=22	GT:DP	0/1:22
seq1	400	indel2	ATTG	A	65	PASS	DP=30	GT:DP	1/1:30
EOF
check_file_exists "$test_data_dir/variants.vcf" "hand-crafted test VCF"

# --- Test Case 1: Select only SNPs ---
log "Starting TEST 1: Select only SNP records"

log "Executing $meta_name with --select_type_to_include SNP..."
"$meta_executable" \
  --variant "$test_data_dir/variants.vcf" \
  --output "$meta_temp_dir/snps_only.vcf" \
  --select_type_to_include SNP

log "Validating TEST 1 outputs..."
check_file_exists "$meta_temp_dir/snps_only.vcf" "SNP-only output VCF"
check_file_not_empty "$meta_temp_dir/snps_only.vcf" "SNP-only output VCF"
check_file_contains "$meta_temp_dir/snps_only.vcf" "^##fileformat=VCF" "SNP-only output VCF header"

record_count=$(grep -c "^seq1" "$meta_temp_dir/snps_only.vcf")
if [[ "$record_count" -eq 2 ]]; then
  log "✓ SNP-only output VCF contains expected number of records (2): $meta_temp_dir/snps_only.vcf"
else
  log_error "✗ SNP-only output VCF has $record_count records, expected 2: $meta_temp_dir/snps_only.vcf"
  exit 1
fi

check_file_contains "$meta_temp_dir/snps_only.vcf" "snp1" "SNP-only output VCF (snp1 record)"
check_file_contains "$meta_temp_dir/snps_only.vcf" "snp2" "SNP-only output VCF (snp2 record)"
check_file_not_contains "$meta_temp_dir/snps_only.vcf" "indel1" "SNP-only output VCF (indel1 record should be absent)"
check_file_not_contains "$meta_temp_dir/snps_only.vcf" "indel2" "SNP-only output VCF (indel2 record should be absent)"

log "✅ TEST 1 completed successfully"

# --- Test Case 2: Select only indels ---
log "Starting TEST 2: Select only INDEL records"

log "Executing $meta_name with --select_type_to_include INDEL..."
"$meta_executable" \
  --variant "$test_data_dir/variants.vcf" \
  --output "$meta_temp_dir/indels_only.vcf" \
  --select_type_to_include INDEL

log "Validating TEST 2 outputs..."
check_file_exists "$meta_temp_dir/indels_only.vcf" "INDEL-only output VCF"
check_file_not_empty "$meta_temp_dir/indels_only.vcf" "INDEL-only output VCF"
check_file_contains "$meta_temp_dir/indels_only.vcf" "^##fileformat=VCF" "INDEL-only output VCF header"

record_count=$(grep -c "^seq1" "$meta_temp_dir/indels_only.vcf")
if [[ "$record_count" -eq 2 ]]; then
  log "✓ INDEL-only output VCF contains expected number of records (2): $meta_temp_dir/indels_only.vcf"
else
  log_error "✗ INDEL-only output VCF has $record_count records, expected 2: $meta_temp_dir/indels_only.vcf"
  exit 1
fi

check_file_contains "$meta_temp_dir/indels_only.vcf" "indel1" "INDEL-only output VCF (indel1 record)"
check_file_contains "$meta_temp_dir/indels_only.vcf" "indel2" "INDEL-only output VCF (indel2 record)"
check_file_not_contains "$meta_temp_dir/indels_only.vcf" "snp1" "INDEL-only output VCF (snp1 record should be absent)"
check_file_not_contains "$meta_temp_dir/indels_only.vcf" "snp2" "INDEL-only output VCF (snp2 record should be absent)"

log "✅ TEST 2 completed successfully"

# --- Test Case 3: Select variants by rsID (--keep_ids) ---
log "Starting TEST 3: Select variants by rsID"

log "Executing $meta_name with --keep_ids snp1;indel2..."
"$meta_executable" \
  --variant "$test_data_dir/variants.vcf" \
  --output "$meta_temp_dir/keep_ids.vcf" \
  --keep_ids "snp1;indel2"

log "Validating TEST 3 outputs..."
check_file_exists "$meta_temp_dir/keep_ids.vcf" "keep-ids output VCF"
check_file_not_empty "$meta_temp_dir/keep_ids.vcf" "keep-ids output VCF"

record_count=$(grep -c "^seq1" "$meta_temp_dir/keep_ids.vcf")
if [[ "$record_count" -eq 2 ]]; then
  log "✓ keep-ids output VCF contains expected number of records (2): $meta_temp_dir/keep_ids.vcf"
else
  log_error "✗ keep-ids output VCF has $record_count records, expected 2: $meta_temp_dir/keep_ids.vcf"
  exit 1
fi

check_file_contains "$meta_temp_dir/keep_ids.vcf" "snp1" "keep-ids output VCF (snp1 record)"
check_file_contains "$meta_temp_dir/keep_ids.vcf" "indel2" "keep-ids output VCF (indel2 record)"
check_file_not_contains "$meta_temp_dir/keep_ids.vcf" "snp2" "keep-ids output VCF (snp2 record should be absent)"
check_file_not_contains "$meta_temp_dir/keep_ids.vcf" "indel1" "keep-ids output VCF (indel1 record should be absent)"

log "✅ TEST 3 completed successfully"

# --- Test Case 4: Restrict indel size (--max_indel_size / --min_indel_size) ---
# indel1 (A -> ATG) is a 2bp insertion, indel2 (ATTG -> A) is a 3bp deletion.
log "Starting TEST 4: Restrict indel size with --max_indel_size"

log "Executing $meta_name with --select_type_to_include INDEL --max_indel_size 2..."
"$meta_executable" \
  --variant "$test_data_dir/variants.vcf" \
  --output "$meta_temp_dir/small_indels.vcf" \
  --select_type_to_include INDEL \
  --max_indel_size 2

log "Validating TEST 4 outputs..."
check_file_exists "$meta_temp_dir/small_indels.vcf" "small-indels output VCF"
check_file_not_empty "$meta_temp_dir/small_indels.vcf" "small-indels output VCF"

record_count=$(grep -c "^seq1" "$meta_temp_dir/small_indels.vcf")
if [[ "$record_count" -eq 1 ]]; then
  log "✓ small-indels output VCF contains expected number of records (1): $meta_temp_dir/small_indels.vcf"
else
  log_error "✗ small-indels output VCF has $record_count records, expected 1: $meta_temp_dir/small_indels.vcf"
  exit 1
fi

check_file_contains "$meta_temp_dir/small_indels.vcf" "indel1" "small-indels output VCF (indel1 record)"
check_file_not_contains "$meta_temp_dir/small_indels.vcf" "indel2" "small-indels output VCF (indel2 record should be absent, 3bp > max)"

log "✅ TEST 4 completed successfully"

# --- Test Case 5: Shared GATK engine options ---
log "Starting TEST 5: Shared GATK engine options"

log "Executing $meta_name with --interval_padding, --sites_only_vcf_output, --create_output_variant_index false..."
"$meta_executable" \
  --variant "$test_data_dir/variants.vcf" \
  --output "$meta_temp_dir/output_engine.vcf" \
  --interval_padding 10 \
  --sites_only_vcf_output \
  --create_output_variant_index false

log "Validating TEST 5 outputs..."
check_file_exists "$meta_temp_dir/output_engine.vcf" "output VCF file (engine options)"
check_file_not_empty "$meta_temp_dir/output_engine.vcf" "output VCF file (engine options)"
check_file_not_exists "$meta_temp_dir/output_engine.vcf.idx" "output VCF index (should not exist with --create_output_variant_index false)"

log "Checking that --sites_only_vcf_output dropped the FORMAT/genotype columns..."
chrom_line=$(grep "^#CHROM" "$meta_temp_dir/output_engine.vcf")
if [[ "$chrom_line" == *"FORMAT"* || "$chrom_line" == *"sample1"* ]]; then
  log_error "✗ #CHROM header still contains FORMAT/sample columns despite --sites_only_vcf_output: $chrom_line"
  exit 1
else
  log "✓ #CHROM header has no FORMAT/sample columns (as expected): $chrom_line"
fi

log "✅ TEST 5 completed successfully"

print_test_summary "All tests completed successfully"
