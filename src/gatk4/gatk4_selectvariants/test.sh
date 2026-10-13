#!/bin/bash

## VIASH START
## VIASH END

# Source the centralized test helpers and the GATK4-specific helpers
source "$meta_resources_dir/test_helpers.sh"
source "$meta_resources_dir/gatk4/test_helpers.sh"

gatk4_test_init

# --- Shared fixtures ---

# A small VCF with 2 SNPs and 2 indels. indel1 (A -> ATG) is a 2bp insertion,
# indel2 (ATTG -> A) is a 3bp deletion.
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

# A reference that matches the seq1 contig of the test VCF
reference_fasta="$test_data_dir/reference.fasta"
create_test_reference "$reference_fasta" 1 1000

# Bgzipped --variant with its index in a different directory, and without an
# index
create_test_vcf_gz_with_separate_index "$test_data_dir/variants.vcf" \
  "$test_data_dir/variants.vcf.gz" "$test_data_dir/indexes/variants.vcf.gz.tbi"
bgzip -c "$test_data_dir/variants.vcf" > "$test_data_dir/variants_no_index.vcf.gz"

# A bgzipped comparison VCF with only the snp1 and indel2 records, with its
# index in a different directory
grep -v -e "snp2" -e "indel1" "$test_data_dir/variants.vcf" > "$test_data_dir/comparison.vcf"
create_test_vcf_gz_with_separate_index "$test_data_dir/comparison.vcf" \
  "$test_data_dir/comparison.vcf.gz" "$test_data_dir/indexes/comparison.vcf.gz.tbi"

# Interval files. The intersection of seq1:50-350 (BED is 0-based) and the
# sites in variants.vcf is snp1, snp2 and indel1. The exclusion of
# seq1:250-350 removes indel1.
mkdir -p "$test_data_dir/intervals"
printf 'seq1\t49\t350\n' > "$test_data_dir/intervals/include.bed"
printf 'seq1\t249\t350\n' > "$test_data_dir/intervals/exclude.bed"
create_test_vcf_gz_with_separate_index "$test_data_dir/variants.vcf" \
  "$test_data_dir/intervals/sites.vcf.gz" "$test_data_dir/intervals_index/sites.vcf.gz.tbi"
echo "seq1:50-250" > "$test_data_dir/intervals/region.intervals"

mkdir -p "$meta_temp_dir/indexes"

# --- Test Case 1: Basic run ---
log "Starting TEST 1: Plain --variant without options"

"$meta_executable" \
  --variant "$test_data_dir/variants.vcf" \
  --output "$meta_temp_dir/all.vcf" \
  --output_index "$meta_temp_dir/indexes/all.idx" 2> "$meta_temp_dir/test1.log"

log "Validating TEST 1 outputs..."
check_file_not_empty "$meta_temp_dir/all.vcf" "output VCF (all records)"
check_file_contains "$meta_temp_dir/all.vcf" "^##fileformat=VCF" "output VCF (all records) header"
check_vcf_column "$meta_temp_dir/all.vcf" 3 "snp1 snp2 indel1 indel2" "output VCF (all records)"
check_file_not_contains "$meta_temp_dir/test1.log" "no index was provided" "run log (plain VCF used without an index)"
check_file_not_exists "$test_data_dir/variants.vcf.idx" "index next to the plain input"
check_index_moved "$meta_temp_dir/indexes/all.idx" "$meta_temp_dir/all.vcf.idx" "output index file"
check_file_not_empty "$meta_temp_dir/indexes/all.idx" "output index file"

log "✅ TEST 1 completed successfully"

# --- Test Case 2: Inputs and indexes ---
log "Starting TEST 2: Input index handling and output options"

# Run 1: --concordance removes snp2 and indel1, SNP selection removes indel2
log "Executing $meta_name with indexed .vcf.gz --variant and --concordance, --reference, SNP selection and a .vcf.gz output..."
"$meta_executable" \
  --variant "$test_data_dir/variants.vcf.gz" \
  --variant_index "$test_data_dir/indexes/variants.vcf.gz.tbi" \
  --concordance "$test_data_dir/comparison.vcf.gz" \
  --concordance_index "$test_data_dir/indexes/comparison.vcf.gz.tbi" \
  --reference "$reference_fasta" \
  --reference_fai "${reference_fasta}.fai" \
  --reference_dict "${reference_fasta%.fasta}.dict" \
  --select_type_to_include SNP \
  --output "$meta_temp_dir/indexed.vcf.gz" \
  --output_index "$meta_temp_dir/indexes/indexed.tbi" 2> "$meta_temp_dir/test2_indexed.log"

log "Validating TEST 2 outputs (indexed inputs)..."
check_file_exists "$meta_temp_dir/indexed.vcf.gz" "bgzipped output VCF"
check_vcf_column "$meta_temp_dir/indexed.vcf.gz" 3 "snp1" "bgzipped output VCF (SNPs in the comparison VCF)"
check_file_not_contains "$meta_temp_dir/test2_indexed.log" "no index was provided" "run log (indexes were provided)"
check_file_not_exists "$test_data_dir/variants.vcf.gz.tbi" "index next to the input VCF"
check_file_not_exists "$test_data_dir/comparison.vcf.gz.tbi" "index next to the --concordance VCF"
check_index_moved "$meta_temp_dir/indexes/indexed.tbi" "$meta_temp_dir/indexed.vcf.gz.tbi" "output index file"
check_file_not_empty "$meta_temp_dir/indexes/indexed.tbi" "output index file"
# --reference has no effect on the selected records
zcat "$meta_temp_dir/indexed.vcf.gz" > "$meta_temp_dir/indexed.vcf"
check_file_contains "$meta_temp_dir/indexed.vcf" "^##GATKCommandLine=.* --reference [^ ]*/reference.fasta " "GATK command line (--reference)"

# Run 2: SNP selection removes the indels, --select_expressions "DP > 24"
# removes snp1
log "Executing $meta_name with unindexed .vcf.gz --variant, SNP selection, --select_expressions and engine options..."
"$meta_executable" \
  --variant "$test_data_dir/variants_no_index.vcf.gz" \
  --select_type_to_include SNP \
  --select_expressions "DP > 24" \
  --output "$meta_temp_dir/no_index.vcf" \
  --interval_padding 10 \
  --sites_only_vcf_output \
  --create_output_variant_index false 2> "$meta_temp_dir/test2_no_index.log"

log "Validating TEST 2 outputs (unindexed input, engine options)..."
check_vcf_column "$meta_temp_dir/no_index.vcf" 3 "snp2" "output VCF (SNPs with DP > 24)"
check_file_contains "$meta_temp_dir/test2_no_index.log" "Warning: no index was provided for '.*variants_no_index.vcf.gz'" "run log (copy of --variant indexed)"
check_file_not_exists "$test_data_dir/variants_no_index.vcf.gz.tbi" "index next to the input VCF (a copy is indexed instead)"
check_file_not_exists "$meta_temp_dir/no_index.vcf.idx" "output VCF index (should not exist with --create_output_variant_index false)"
check_file_contains "$meta_temp_dir/no_index.vcf" "^##GATKCommandLine=.* --interval-padding 10 " "GATK command line (--interval_padding)"

log "Checking that --sites_only_vcf_output dropped the FORMAT/genotype columns..."
chrom_line=$(grep "^#CHROM" "$meta_temp_dir/no_index.vcf")
if [[ "$chrom_line" == *"FORMAT"* || "$chrom_line" == *"sample1"* ]]; then
  log_error "✗ #CHROM header still contains FORMAT/sample columns despite --sites_only_vcf_output: $chrom_line"
  exit 1
else
  log "✓ #CHROM header has no FORMAT/sample columns (as expected): $chrom_line"
fi

log "✅ TEST 2 completed successfully"

# --- Test Case 3: Intervals ---
log "Starting TEST 3: Interval files, interval strings and --variant_output_filtering"

# Run 1: two --intervals files, INTERSECTION and --exclude_intervals
log "Executing $meta_name with two --intervals files, INTERSECTION and --exclude_intervals..."
"$meta_executable" \
  --variant "$test_data_dir/variants.vcf" \
  --intervals "$test_data_dir/intervals/include.bed" \
  --intervals "$test_data_dir/intervals/sites.vcf.gz" \
  --intervals_index "$test_data_dir/intervals_index/sites.vcf.gz.tbi" \
  --interval_set_rule INTERSECTION \
  --exclude_intervals "$test_data_dir/intervals/exclude.bed" \
  --output "$meta_temp_dir/intervals.vcf" \
  --output_index "$meta_temp_dir/intervals.vcf.idx" 2> "$meta_temp_dir/test3_intervals.log"

log "Validating TEST 3 outputs (interval files)..."
check_file_not_contains "$meta_temp_dir/test3_intervals.log" "no index was provided for '$test_data_dir/intervals/sites.vcf.gz'" "run log (no indexing of the interval VCF when --intervals_index is given)"
check_file_not_exists "$test_data_dir/intervals/sites.vcf.gz.tbi" "index next to the interval VCF"
check_vcf_column "$meta_temp_dir/intervals.vcf" 3 "snp1 snp2" "output VCF (records in the selected intervals)"

# Run 2: an --intervals file (seq1:50-250) and an interval string
# (seq1:380-420), minus an excluded string (seq1:190-210)
log "Executing $meta_name with an --intervals file, --interval_strings and --exclude_interval_strings..."
"$meta_executable" \
  --variant "$test_data_dir/variants.vcf" \
  --output "$meta_temp_dir/strings_mixed.vcf" \
  --output_index "$meta_temp_dir/strings_mixed.vcf.idx" \
  --intervals "$test_data_dir/intervals/region.intervals" \
  --interval_strings "seq1:380-420" \
  --exclude_interval_strings "seq1:190-210" 2> "$meta_temp_dir/test3_mixed.log"

log "Validating TEST 3 outputs (mixed intervals)..."
check_vcf_column "$meta_temp_dir/strings_mixed.vcf" 3 "snp1 indel2" "output VCF (interval file and string, minus the excluded string)"
check_file_contains "$meta_temp_dir/test3_mixed.log" "no index was provided" "run log (indexing warning with --intervals)"
check_file_not_exists "$test_data_dir/variants.vcf.idx" "index next to the plain input"

# Runs 3 and 4: indel2 (seq1:400-403) overlaps seq1:402-410 but does not start
# in it
log "Executing $meta_name with only --interval_strings and no --variant_output_filtering..."
"$meta_executable" \
  --variant "$test_data_dir/variants.vcf" \
  --output "$meta_temp_dir/filtering_default.vcf" \
  --create_output_variant_index false \
  --interval_strings "seq1:190-210" \
  --interval_strings "seq1:402-410" 2> "$meta_temp_dir/test3_strings.log"

log "Executing $meta_name with only --interval_strings and --variant_output_filtering STARTS_IN..."
"$meta_executable" \
  --variant "$test_data_dir/variants.vcf" \
  --output "$meta_temp_dir/filtering_starts_in.vcf" \
  --create_output_variant_index false \
  --interval_strings "seq1:190-210" \
  --interval_strings "seq1:402-410" \
  --variant_output_filtering STARTS_IN

log "Validating TEST 3 outputs (--variant_output_filtering)..."
check_vcf_column "$meta_temp_dir/filtering_default.vcf" 3 "snp2 indel2" "output VCF (no filtering mode, records that overlap the intervals)"
check_file_contains "$meta_temp_dir/test3_strings.log" "no index was provided" "run log (plain --variant indexed when interval strings are given)"
check_file_not_exists "$test_data_dir/variants.vcf.idx" "index next to the plain input"
check_vcf_column "$meta_temp_dir/filtering_starts_in.vcf" 3 "snp2" "output VCF (STARTS_IN, records that start in the intervals)"

log "✅ TEST 3 completed successfully"

# --- Test Case 4: Tool options ---
log "Starting TEST 4: Record selection options"

# Run 1: --keep_ids removes snp2 and indel1, --select "DP > 24" removes snp1
log "Executing $meta_name with --keep_ids 'snp1;indel2' and --select 'DP > 24'..."
"$meta_executable" \
  --variant "$test_data_dir/variants.vcf" \
  --output "$meta_temp_dir/keep_ids_select.vcf" \
  --create_output_variant_index false \
  --keep_ids "snp1;indel2" \
  --select "DP > 24"

log "Validating TEST 4 outputs (--keep_ids, --select)..."
check_vcf_column "$meta_temp_dir/keep_ids_select.vcf" 3 "indel2" "output VCF (kept IDs with DP > 24)"

# Run 2: INDEL selection removes the SNPs, --max_indel_size 2 removes indel2
log "Executing $meta_name with --select_type_to_include INDEL --max_indel_size 2..."
"$meta_executable" \
  --variant "$test_data_dir/variants.vcf" \
  --output "$meta_temp_dir/small_indels.vcf" \
  --create_output_variant_index false \
  --select_type_to_include INDEL \
  --max_indel_size 2

log "Validating TEST 4 outputs (INDEL selection, --max_indel_size)..."
check_vcf_column "$meta_temp_dir/small_indels.vcf" 3 "indel1" "output VCF (indels of 2bp or less)"

log "✅ TEST 4 completed successfully"

# --- Test Case 5: Argument errors ---
log "Starting TEST 5: Argument errors"

log "Executing $meta_name with --concordance_index and no --concordance..."
check_fails_with "Error: --concordance_index requires --concordance" "$meta_temp_dir/test5_concordance.log" \
  "$meta_executable" \
  --variant "$test_data_dir/variants.vcf" \
  --concordance_index "$test_data_dir/indexes/comparison.vcf.gz.tbi" \
  --output "$meta_temp_dir/index_only.vcf" \
  --output_index "$meta_temp_dir/index_only.vcf.idx"
check_file_not_exists "$meta_temp_dir/index_only.vcf" "output VCF (GATK did not run)"

log "Executing $meta_name without --output_index..."
check_fails_with "Error: --output_index is required unless --create_output_variant_index is false." "$meta_temp_dir/test5_output_index.log" \
  "$meta_executable" \
  --variant "$test_data_dir/variants.vcf" \
  --output "$meta_temp_dir/no_output_index.vcf"
check_file_not_exists "$meta_temp_dir/no_output_index.vcf" "output VCF (GATK did not run)"

log "✅ TEST 5 completed successfully"

print_test_summary "All tests"
