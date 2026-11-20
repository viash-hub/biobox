#!/bin/bash

## VIASH START
## VIASH END

# Source the centralized test helpers and the GATK4-specific helpers
source "$meta_resources_dir/test_helpers.sh"
source "$meta_resources_dir/gatk4/test_helpers.sh"

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

mkdir -p "$meta_temp_dir/indexes"

log "Executing $meta_name with --select_type_to_include SNP..."
"$meta_executable" \
  --variant "$test_data_dir/variants.vcf" \
  --output "$meta_temp_dir/snps_only.vcf" \
  --output_index "$meta_temp_dir/indexes/snps_only.idx" \
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
check_file_exists "$meta_temp_dir/indexes/snps_only.idx" "output index file"
check_file_not_empty "$meta_temp_dir/indexes/snps_only.idx" "output index file"
check_file_not_exists "$meta_temp_dir/snps_only.vcf.idx" "index next to the output VCF (moved to --output_index)"

log "✅ TEST 1 completed successfully"

# --- Test Case 2: Select only indels ---
log "Starting TEST 2: Select only INDEL records"

log "Executing $meta_name with --select_type_to_include INDEL..."
"$meta_executable" \
  --variant "$test_data_dir/variants.vcf" \
  --output "$meta_temp_dir/indels_only.vcf" \
  --output_index "$meta_temp_dir/indels_only.vcf.idx" \
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
  --output_index "$meta_temp_dir/keep_ids.vcf.idx" \
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
  --output_index "$meta_temp_dir/small_indels.vcf.idx" \
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

# --- Test Case 6: Optional --reference trio ---
log "Starting TEST 6: Optional --reference"

log "Creating synthetic reference matching the test VCF's seq1 contig..."
reference_fasta="$test_data_dir/reference.fasta"
create_test_reference "$reference_fasta" 1 1000

log "Executing $meta_name with --reference..."
"$meta_executable" \
  --variant "$test_data_dir/variants.vcf" \
  --output "$meta_temp_dir/with_reference.vcf" \
  --output_index "$meta_temp_dir/with_reference.vcf.idx" \
  --reference "$reference_fasta" \
  --reference_fai "${reference_fasta}.fai" \
  --reference_dict "${reference_fasta%.fasta}.dict" \
  --select_type_to_include SNP

log "Validating TEST 6 outputs..."
check_file_exists "$meta_temp_dir/with_reference.vcf" "output VCF file (with reference)"
check_file_not_empty "$meta_temp_dir/with_reference.vcf" "output VCF file (with reference)"

record_count=$(grep -c "^seq1" "$meta_temp_dir/with_reference.vcf")
if [[ "$record_count" -eq 2 ]]; then
  log "✓ output VCF (with reference) contains expected number of records (2): $meta_temp_dir/with_reference.vcf"
else
  log_error "✗ output VCF (with reference) has $record_count records, expected 2: $meta_temp_dir/with_reference.vcf"
  exit 1
fi

log "✅ TEST 6 completed successfully"

# --- Test Case 7: Bgzipped --variant with a separate index ---
log "Starting TEST 7: Bgzipped --variant with its index in a different directory"

create_test_vcf_gz_with_separate_index "$test_data_dir/variants.vcf" \
  "$test_data_dir/variants.vcf.gz" "$test_data_dir/indexes/variants.vcf.gz.tbi"

log "Executing $meta_name with --variant .vcf.gz and --variant_index..."
"$meta_executable" \
  --variant "$test_data_dir/variants.vcf.gz" \
  --variant_index "$test_data_dir/indexes/variants.vcf.gz.tbi" \
  --output "$meta_temp_dir/from_vcf_gz.vcf" \
  --output_index "$meta_temp_dir/from_vcf_gz.vcf.idx" \
  --select_type_to_include SNP 2>&1 | tee "$meta_temp_dir/test7.log"

log "Validating TEST 7 outputs..."
check_file_exists "$meta_temp_dir/from_vcf_gz.vcf" "output VCF file (from .vcf.gz)"
check_file_contains "$meta_temp_dir/from_vcf_gz.vcf" "^##fileformat=VCF" "output VCF file (from .vcf.gz) header"
check_file_contains "$meta_temp_dir/from_vcf_gz.vcf" "snp1" "output VCF file (from .vcf.gz) (snp1 record)"
check_file_contains "$meta_temp_dir/from_vcf_gz.vcf" "snp2" "output VCF file (from .vcf.gz) (snp2 record)"
check_file_not_contains "$meta_temp_dir/from_vcf_gz.vcf" "indel1" "output VCF file (from .vcf.gz) (indel1 record should be absent)"
check_file_not_contains "$meta_temp_dir/test7.log" "no index was provided" "run log (index was provided)"
check_file_not_exists "$test_data_dir/variants.vcf.gz.tbi" "index next to the input VCF"

log "✅ TEST 7 completed successfully"

# --- Test Case 8: Bgzipped --variant without an index ---
log "Starting TEST 8: Bgzipped --variant without an index"

bgzip -c "$test_data_dir/variants.vcf" > "$test_data_dir/variants_no_index.vcf.gz"

log "Executing $meta_name with --variant .vcf.gz and no index..."
"$meta_executable" \
  --variant "$test_data_dir/variants_no_index.vcf.gz" \
  --output "$meta_temp_dir/from_vcf_gz_no_index.vcf" \
  --output_index "$meta_temp_dir/from_vcf_gz_no_index.vcf.idx" \
  --select_type_to_include SNP 2>&1 | tee "$meta_temp_dir/test8.log"

log "Validating TEST 8 outputs..."
check_file_exists "$meta_temp_dir/from_vcf_gz_no_index.vcf" "output VCF file (from .vcf.gz without index)"
check_file_contains "$meta_temp_dir/from_vcf_gz_no_index.vcf" "snp1" "output VCF file (from .vcf.gz without index) (snp1 record)"
check_file_contains "$meta_temp_dir/test8.log" "Warning: no index was provided for '.*variants_no_index.vcf.gz'" "run log (copy of --variant indexed)"
check_file_not_exists "$test_data_dir/variants_no_index.vcf.gz.tbi" "index next to the input VCF (a copy is indexed instead)"

log "✅ TEST 8 completed successfully"

# --- Test Case 9: --concordance with a separate index ---
log "Starting TEST 9: Bgzipped --concordance with its index in a different directory"

# The comparison VCF has only the snp1 and indel2 records
grep -v -e "snp2" -e "indel1" "$test_data_dir/variants.vcf" > "$test_data_dir/comparison.vcf"
create_test_vcf_gz_with_separate_index "$test_data_dir/comparison.vcf" \
  "$test_data_dir/comparison.vcf.gz" "$test_data_dir/indexes/comparison.vcf.gz.tbi"

log "Executing $meta_name with --concordance and --concordance_index..."
"$meta_executable" \
  --variant "$test_data_dir/variants.vcf" \
  --concordance "$test_data_dir/comparison.vcf.gz" \
  --concordance_index "$test_data_dir/indexes/comparison.vcf.gz.tbi" \
  --output "$meta_temp_dir/concordance.vcf" \
  --output_index "$meta_temp_dir/concordance.vcf.idx"

log "Validating TEST 9 outputs..."
check_file_exists "$meta_temp_dir/concordance.vcf" "concordance output VCF"
record_count=$(grep -c "^seq1" "$meta_temp_dir/concordance.vcf")
if [[ "$record_count" -eq 2 ]]; then
  log "✓ concordance output VCF contains expected number of records (2): $meta_temp_dir/concordance.vcf"
else
  log_error "✗ concordance output VCF has $record_count records, expected 2: $meta_temp_dir/concordance.vcf"
  exit 1
fi
check_file_contains "$meta_temp_dir/concordance.vcf" "snp1" "concordance output VCF (snp1 record)"
check_file_contains "$meta_temp_dir/concordance.vcf" "indel2" "concordance output VCF (indel2 record)"
check_file_not_contains "$meta_temp_dir/concordance.vcf" "snp2" "concordance output VCF (snp2 record should be absent)"
check_file_not_exists "$test_data_dir/comparison.vcf.gz.tbi" "index next to the --concordance VCF"

log "✅ TEST 9 completed successfully"

# --- Test Case 10: An index without its input fails ---
log "Starting TEST 10: --concordance_index without --concordance"

if "$meta_executable" \
  --variant "$test_data_dir/variants.vcf" \
  --concordance_index "$test_data_dir/indexes/comparison.vcf.gz.tbi" \
  --output "$meta_temp_dir/index_only.vcf" > "$meta_temp_dir/test10.log" 2>&1; then
  log_error "✗ $meta_name did not fail with --concordance_index but no --concordance"
  exit 1
fi
check_file_contains "$meta_temp_dir/test10.log" "Error: --concordance_index requires --concordance" "error message"
check_file_not_exists "$meta_temp_dir/index_only.vcf" "output VCF (GATK did not run)"

log "✅ TEST 10 completed successfully"

# --- Test Case 11: Plain --variant without an index ---
log "Starting TEST 11: Plain --variant without an index, with and without --intervals"

mkdir -p "$test_data_dir/plain"
cp "$test_data_dir/variants.vcf" "$test_data_dir/plain/variants.vcf"
echo "seq1:50-250" > "$test_data_dir/plain/region.intervals"

log "Executing $meta_name with a plain --variant and no --intervals..."
"$meta_executable" \
  --variant "$test_data_dir/plain/variants.vcf" \
  --output "$meta_temp_dir/plain_all.vcf" \
  --output_index "$meta_temp_dir/plain_all.vcf.idx" 2> "$meta_temp_dir/test11_all.log"

log "Validating TEST 11 outputs without --intervals..."
check_file_contains "$meta_temp_dir/plain_all.vcf" "indel2" "output VCF (all records)"
check_file_not_contains "$meta_temp_dir/test11_all.log" "no index was provided" "run log (plain VCF used without an index)"
check_file_not_exists "$test_data_dir/plain/variants.vcf.idx" "index next to the plain input"

log "Executing $meta_name with a plain --variant and --intervals..."
"$meta_executable" \
  --variant "$test_data_dir/plain/variants.vcf" \
  --intervals "$test_data_dir/plain/region.intervals" \
  --output "$meta_temp_dir/plain_region.vcf" \
  --output_index "$meta_temp_dir/plain_region.vcf.idx" 2> "$meta_temp_dir/test11_region.log"

log "Validating TEST 11 outputs with --intervals..."
check_file_contains "$meta_temp_dir/plain_region.vcf" "snp2" "output VCF (record in the interval)"
check_file_not_contains "$meta_temp_dir/plain_region.vcf" "indel2" "output VCF (record outside the interval)"
check_file_contains "$meta_temp_dir/test11_region.log" "no index was provided" "run log (indexing warning with --intervals)"
check_file_not_exists "$test_data_dir/plain/variants.vcf.idx" "index next to the plain input"

log "✅ TEST 11 completed successfully"

# --- Test Case 12: Bgzipped output with an output index ---
log "Starting TEST 12: Bgzipped output with its index in a different directory"

log "Executing $meta_name with a .vcf.gz --output and --output_index..."
"$meta_executable" \
  --variant "$test_data_dir/variants.vcf" \
  --output "$meta_temp_dir/snps_only.vcf.gz" \
  --output_index "$meta_temp_dir/indexes/snps_only.tbi" \
  --select_type_to_include SNP

log "Validating TEST 12 outputs..."
check_file_exists "$meta_temp_dir/snps_only.vcf.gz" "bgzipped output VCF"
record_count=$(zcat "$meta_temp_dir/snps_only.vcf.gz" | grep -c "^seq1")
if [[ "$record_count" -eq 2 ]]; then
  log "✓ bgzipped output VCF contains expected number of records (2)"
else
  log_error "✗ bgzipped output VCF has $record_count records, expected 2"
  exit 1
fi
check_file_exists "$meta_temp_dir/indexes/snps_only.tbi" "output index file"
check_file_not_empty "$meta_temp_dir/indexes/snps_only.tbi" "output index file"
check_file_not_exists "$meta_temp_dir/snps_only.vcf.gz.tbi" "index next to the output VCF (moved to --output_index)"

log "✅ TEST 12 completed successfully"

# --- Test Case 13: Missing --output_index ---
log "Starting TEST 13: A missing --output_index fails"

if "$meta_executable" \
  --variant "$test_data_dir/variants.vcf" \
  --output "$meta_temp_dir/no_output_index.vcf" > "$meta_temp_dir/test13.log" 2>&1; then
  log_error "✗ $meta_name did not fail without --output_index"
  exit 1
fi
check_file_contains "$meta_temp_dir/test13.log" "Error: --output_index is required unless --create_output_variant_index is false." "error message"
check_file_not_exists "$meta_temp_dir/no_output_index.vcf" "output VCF (GATK did not run)"

log "✅ TEST 13 completed successfully"

print_test_summary "All tests completed successfully"
