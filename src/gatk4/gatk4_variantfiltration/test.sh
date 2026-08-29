#!/bin/bash

## VIASH START
## VIASH END

# Source the centralized test helpers
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

# --- Prepare shared reference files ---
ref_fasta="$test_data_dir/reference.fasta"
ref_fai="$test_data_dir/reference.fasta.fai"
ref_dict="$test_data_dir/reference.dict"

create_test_reference "$ref_fasta" 1 1000
check_file_exists "$ref_fasta" "reference FASTA file"
check_file_exists "$ref_fai" "reference FASTA index"
check_file_exists "$ref_dict" "reference sequence dictionary"

# create_test_fasta writes a single contig named "seq1" whose sequence is
# the repeating pattern "ATCG" truncated to the requested length. For any
# position that is a multiple of 4 (e.g. 100, 200, 300, 400), the base at
# that position is therefore always "G" -- used as the REF allele below.

# --- Test Case 1: --filter_expression / --filter_name ---
log "Starting TEST 1: --filter_expression / --filter_name (site-level JEXL filters)"

input_vcf="$test_data_dir/input.vcf"
cat <<EOF > "$input_vcf"
##fileformat=VCFv4.2
##contig=<ID=seq1,length=1000>
##INFO=<ID=QD,Number=1,Type=Float,Description="Variant Confidence/Quality by Depth">
##INFO=<ID=FS,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher's exact test to detect strand bias">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1
seq1	100	.	G	A	100	.	QD=20.0;FS=1.0	GT	0/1
seq1	200	.	G	A	100	.	QD=1.0;FS=1.0	GT	0/1
seq1	300	.	G	A	100	.	QD=1.0;FS=70.0	GT	0/1
EOF
check_file_exists "$input_vcf" "input VCF"

output_vcf="$meta_temp_dir/output1/filtered.vcf"
mkdir -p "$(dirname "$output_vcf")"

log "Executing $meta_name with --filter_expression/--filter_name pairs..."
"$meta_executable" \
  --variant "$input_vcf" \
  --reference "$ref_fasta" \
  --reference_fai "$ref_fai" \
  --reference_dict "$ref_dict" \
  --output "$output_vcf" \
  --filter_expression "QD < 2.0" \
  --filter_name "lowQD" \
  --filter_expression "FS > 60.0" \
  --filter_name "highFS"

log "Validating TEST 1 outputs..."
check_file_exists "$output_vcf" "filtered VCF"
check_file_not_empty "$output_vcf" "filtered VCF"

# Extract the FILTER column for each record so assertions can't be confused
# by filter names appearing in other records.
record1_filter="$meta_temp_dir/record1_filter.txt"
record2_filter="$meta_temp_dir/record2_filter.txt"
record3_filter="$meta_temp_dir/record3_filter.txt"

awk -F'\t' '$2==100 {print $7}' "$output_vcf" > "$record1_filter"
awk -F'\t' '$2==200 {print $7}' "$output_vcf" > "$record2_filter"
awk -F'\t' '$2==300 {print $7}' "$output_vcf" > "$record3_filter"

check_file_not_empty "$record1_filter" "record 1 (pos 100) FILTER field"
check_file_not_empty "$record2_filter" "record 2 (pos 200) FILTER field"
check_file_not_empty "$record3_filter" "record 3 (pos 300) FILTER field"

# Record 1: passes both filters -> PASS
log "Checking record 1 (QD=20.0, FS=1.0, passes both filters) is marked PASS..."
if [[ "$(cat "$record1_filter")" == "PASS" ]]; then
  log "✓ record 1 (pos 100) FILTER is PASS"
else
  log_error "✗ record 1 (pos 100) FILTER is '$(cat "$record1_filter")', expected PASS"
  exit 1
fi

# Record 2: fails only the QD filter -> lowQD, and never highFS
log "Checking record 2 (QD=1.0 only) is marked lowQD..."
check_file_contains "$record2_filter" "lowQD" "record 2 (pos 200) FILTER field"
check_file_not_contains "$record2_filter" "highFS" "record 2 (pos 200) FILTER field"

# Record 3: fails BOTH filters -> both names present (order not guaranteed)
log "Checking record 3 (QD=1.0, FS=70.0) fails both filters..."
check_file_contains "$record3_filter" "lowQD" "record 3 (pos 300) FILTER field"
check_file_contains "$record3_filter" "highFS" "record 3 (pos 300) FILTER field"

log "✅ TEST 1 completed successfully"

# --- Test Case 2: Mask and genotype-level filtering ---
log "Starting TEST 2: Mask filtering and genotype-level (FORMAT) filtering"

# Reuse create_test_known_sites_vcf as a stand-in for a "mask" 
mask_vcf="$test_data_dir/mask.vcf"
create_test_known_sites_vcf "$mask_vcf" seq1 1000 700
check_file_exists "$mask_vcf" "mask VCF"

input_vcf3="$test_data_dir/input3.vcf"
cat <<EOF > "$input_vcf3"
##fileformat=VCFv4.2
##contig=<ID=seq1,length=1000>
##INFO=<ID=QD,Number=1,Type=Float,Description="Variant Confidence/Quality by Depth">
##INFO=<ID=FS,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher's exact test to detect strand bias">
##INFO=<ID=SOR,Number=1,Type=Float,Description="Symmetric Odds Ratio of 2x2 contingency table to detect strand bias">
##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">
##INFO=<ID=MQRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities">
##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1
seq1	500	.	G	A	100	.	QD=20.0;FS=1.0;SOR=1.0;MQ=60.0;MQRankSum=0.0;ReadPosRankSum=0.0	GT:DP	0/1:5
seq1	600	.	G	A	100	.	QD=20.0;FS=1.0;SOR=1.0;MQ=60.0;MQRankSum=0.0;ReadPosRankSum=0.0	GT:DP	0/1:50
seq1	700	.	G	A	100	.	QD=20.0;FS=1.0;SOR=1.0;MQ=60.0;MQRankSum=0.0;ReadPosRankSum=0.0	GT:DP	0/1:50
EOF
check_file_exists "$input_vcf3" "input VCF for mask/genotype-filter test"

output_vcf3="$meta_temp_dir/output3/filtered.vcf"
mkdir -p "$(dirname "$output_vcf3")"

log "Executing $meta_name with --mask, --genotype_filter_expression and --set_filtered_genotype_to_no_call..."
"$meta_executable" \
  --variant "$input_vcf3" \
  --reference "$ref_fasta" \
  --reference_fai "$ref_fai" \
  --reference_dict "$ref_dict" \
  --output "$output_vcf3" \
  --mask "$mask_vcf" \
  --mask_name "MyMask" \
  --mask_extension 10 \
  --genotype_filter_expression "DP < 10" \
  --genotype_filter_name "lowDP" \
  --set_filtered_genotype_to_no_call

log "Validating TEST 2 outputs..."
check_file_exists "$output_vcf3" "filtered VCF (mask/genotype-filter test)"
check_file_not_empty "$output_vcf3" "filtered VCF (mask/genotype-filter test)"

record500="$meta_temp_dir/record500.txt"
record600="$meta_temp_dir/record600.txt"
record700_filter="$meta_temp_dir/record700_filter.txt"

awk -F'\t' '$2==500 {print $10}' "$output_vcf3" > "$record500"
awk -F'\t' '$2==600 {print $10}' "$output_vcf3" > "$record600"
awk -F'\t' '$2==700 {print $7}' "$output_vcf3" > "$record700_filter"

# Record at pos 500 (DP=5) fails the genotype filter: its genotype is set to
# no-call and the sample column carries the "lowDP" FT annotation.
log "Checking record at pos 500 (DP=5) genotype is set to no-call with FT=lowDP..."
check_file_contains "$record500" "\./\.:5:lowDP" "sample column at pos 500"

# Record at pos 600 (DP=50) passes the genotype filter: genotype is
# untouched and no FT annotation is added.
log "Checking record at pos 600 (DP=50) genotype is left untouched..."
check_file_contains "$record600" "0/1:50" "sample column at pos 600"

# Record at pos 700 coincides with the mask VCF's single record: its
# site-level FILTER is set to the configured mask name.
log "Checking record at pos 700 (overlaps mask) FILTER is MyMask..."
check_file_contains "$record700_filter" "MyMask" "record at pos 700 FILTER field"

log "✅ TEST 2 completed successfully"

# --- Test Case 3: Shared GATK engine options ---
log "Starting TEST 3: Shared GATK engine options"

output_vcf4="$meta_temp_dir/output4/filtered.vcf"
mkdir -p "$(dirname "$output_vcf4")"

log "Executing $meta_name with --interval_padding and --create_output_variant_index false..."
"$meta_executable" \
  --variant "$input_vcf" \
  --reference "$ref_fasta" \
  --reference_fai "$ref_fai" \
  --reference_dict "$ref_dict" \
  --output "$output_vcf4" \
  --interval_padding 10 \
  --create_output_variant_index false

log "Validating TEST 3 outputs..."
check_file_exists "$output_vcf4" "filtered VCF (engine options)"
check_file_not_empty "$output_vcf4" "filtered VCF (engine options)"
check_file_not_exists "$output_vcf4.idx" "output VCF index (should not exist with --create_output_variant_index false)"

log "✅ TEST 3 completed successfully"

# --- Test Case 4: VariantFiltration without a reference (optional per tool help) ---
log "Starting TEST 4: VariantFiltration without --reference"

output_vcf5="$meta_temp_dir/output5/filtered.vcf"
mkdir -p "$(dirname "$output_vcf5")"

log "Executing $meta_name without a reference..."
"$meta_executable" \
  --variant "$input_vcf" \
  --output "$output_vcf5" \
  --filter_expression "QD < 2.0" \
  --filter_name "lowQD"

log "Validating TEST 4 outputs..."
check_file_exists "$output_vcf5" "filtered VCF without reference"
check_file_not_empty "$output_vcf5" "filtered VCF without reference"

log "✅ TEST 4 completed successfully"

print_test_summary "All tests"
