#!/bin/bash

## VIASH START
## VIASH END

# Source the centralized test helpers
source "$meta_resources_dir/test_helpers.sh"
source "$meta_resources_dir/gatk4/test_helpers.sh"

gatk4_test_init

# --- Shared fixtures ---
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

# Input for the site-level filters: record 1 passes both filters, record 2
# fails the QD filter, record 3 fails both
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
create_test_vcf_gz_with_separate_index "$input_vcf" \
  "$test_data_dir/input.vcf.gz" "$test_data_dir/indexes/input.vcf.gz.tbi"

# Input for the mask and genotype filters: the genotype at 500 has DP=5, the
# record at 700 is in the mask
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
create_test_vcf_gz_with_separate_index "$input_vcf3" \
  "$test_data_dir/input3.vcf.gz" "$test_data_dir/indexes/input3.vcf.gz.tbi"
mkdir -p "$test_data_dir/no_index"
bgzip -c "$input_vcf3" > "$test_data_dir/no_index/input3.vcf.gz"

# Masks at seq1:700: a bgzipped VCF with a separate index (the known-sites
# helper writes a one-record VCF), a plain BED without an index and a copy of
# the BED with a Tribble index in a different directory. GATK needs an index
# for a --mask file, also for a plain BED, so without --mask_index a copy of
# the BED is indexed.
create_test_known_sites_vcf --no-index "$test_data_dir/mask.vcf" seq1 1000 700
create_test_vcf_gz_with_separate_index "$test_data_dir/mask.vcf" \
  "$test_data_dir/mask.vcf.gz" "$test_data_dir/indexes/mask.vcf.gz.tbi"
printf 'seq1\t699\t700\n' > "$test_data_dir/mask.bed"
mkdir -p "$test_data_dir/bed_index"
cp "$test_data_dir/mask.bed" "$test_data_dir/bed_index/mask.bed"
gatk IndexFeatureFile --input "$test_data_dir/bed_index/mask.bed" --verbosity ERROR
check_file_exists "$test_data_dir/bed_index/mask.bed.idx" "separate index for the --mask BED"

# Input for the interval files: records at 100, 200, ..., 900. The area
# is the union of seq1:1-450 (BED) and seq1:800 (VCF), minus seq1:151-250
# (BED).
intervals_input_vcf="$test_data_dir/intervals_input.vcf"
{
  printf '##fileformat=VCFv4.2\n'
  printf '##contig=<ID=seq1,length=1000>\n'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
  for pos in 100 200 300 400 500 600 700 800 900; do
    printf 'seq1\t%d\t.\tG\tA\t100\t.\t.\n' "$pos"
  done
} > "$intervals_input_vcf"
printf 'seq1\t0\t450\n' > "$test_data_dir/intervals_a.bed"
printf 'seq1\t150\t250\n' > "$test_data_dir/exclude.bed"
create_test_known_sites_vcf --no-index "$test_data_dir/intervals_site.vcf" seq1 1000 800
create_test_vcf_gz_with_separate_index "$test_data_dir/intervals_site.vcf" \
  "$test_data_dir/intervals_site.vcf.gz" "$test_data_dir/interval_indexes/intervals_site.vcf.gz.tbi"

# Input for --variant_output_filtering
filtering_input_vcf="$test_data_dir/filtering_input.vcf"
{
  printf '##fileformat=VCFv4.2\n'
  printf '##contig=<ID=seq1,length=1000>\n'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
  printf 'seq1\t200\tsnp\tG\tA\t100\t.\t.\n'
  printf 'seq1\t400\tdeletion\tGATC\tG\t100\t.\t.\n'
} > "$filtering_input_vcf"

mkdir -p "$meta_temp_dir/indexes"

# --- Test Case 1: Basic run ---
log "Starting TEST 1: --filter_expression / --filter_name (site-level JEXL filters)"

output_vcf="$meta_temp_dir/output1/filtered.vcf"
mkdir -p "$(dirname "$output_vcf")"

"$meta_executable" \
  --variant "$input_vcf" \
  --reference "$ref_fasta" \
  --reference_fai "$ref_fai" \
  --reference_dict "$ref_dict" \
  --output "$output_vcf" \
  --output_index "$meta_temp_dir/indexes/filtered1.idx" \
  --filter_expression "QD < 2.0" \
  --filter_name "lowQD" \
  --filter_expression "FS > 60.0" \
  --filter_name "highFS" 2> "$meta_temp_dir/test1.log"

log "Validating TEST 1 outputs..."
check_file_not_empty "$output_vcf" "filtered VCF"
check_index_moved "$meta_temp_dir/indexes/filtered1.idx" "$output_vcf.idx" "output index file"
check_file_not_empty "$meta_temp_dir/indexes/filtered1.idx" "output index file"
check_file_not_contains "$meta_temp_dir/test1.log" "no index was provided" "run log (plain --variant used without an index)"
check_file_not_exists "$input_vcf.idx" "index next to the input VCF"

# Extract the FILTER column for each record so assertions can't be confused
# by filter names appearing in other records.
record2_filter="$meta_temp_dir/record2_filter.txt"
record3_filter="$meta_temp_dir/record3_filter.txt"
awk -F'\t' '$2==200 {print $7}' "$output_vcf" > "$record2_filter"
awk -F'\t' '$2==300 {print $7}' "$output_vcf" > "$record3_filter"
check_file_not_empty "$record2_filter" "record 2 (pos 200) FILTER field"
check_file_not_empty "$record3_filter" "record 3 (pos 300) FILTER field"

# Record 1: passes both filters -> PASS
check_vcf_has_record "$output_vcf" $'^seq1\t100\t([^\t]*\t){4}PASS\t' "filtered VCF (record 1 FILTER is PASS)"

# Record 2: fails only the QD filter -> lowQD, and never highFS
check_file_contains "$record2_filter" "lowQD" "record 2 (pos 200) FILTER field"
check_file_not_contains "$record2_filter" "highFS" "record 2 (pos 200) FILTER field"

# Record 3: fails BOTH filters -> both names present (order not guaranteed)
check_file_contains "$record3_filter" "lowQD" "record 3 (pos 300) FILTER field"
check_file_contains "$record3_filter" "highFS" "record 3 (pos 300) FILTER field"

log "✅ TEST 1 completed successfully"

# --- Test Case 2: Inputs and indexes ---
log "Starting TEST 2: Input and mask index handling, output options, mask and genotype filters"

# Run 1: bgzipped --variant with a separate index, no reference, and a
# .vcf.gz output with its index in a different directory
output_vcf_gz="$meta_temp_dir/output_gz/filtered.vcf.gz"
mkdir -p "$(dirname "$output_vcf_gz")"

log "Executing $meta_name with --variant .vcf.gz, --variant_index, no reference and a .vcf.gz output..."
"$meta_executable" \
  --variant "$test_data_dir/input.vcf.gz" \
  --variant_index "$test_data_dir/indexes/input.vcf.gz.tbi" \
  --output "$output_vcf_gz" \
  --output_index "$meta_temp_dir/indexes/filtered_gz.tbi" \
  --filter_expression "QD < 2.0" \
  --filter_name "lowQD" 2> "$meta_temp_dir/test2_gz.log"

log "Validating TEST 2 outputs (indexed .vcf.gz input and output)..."
check_file_exists "$output_vcf_gz" "bgzipped filtered VCF"
check_vcf_has_record "$output_vcf_gz" $'^seq1\t200\t([^\t]*\t){4}lowQD\t' "bgzipped filtered VCF (record 2 FILTER is lowQD)"
check_file_not_contains "$meta_temp_dir/test2_gz.log" "no index was provided" "run log (index was provided)"
check_file_not_exists "$test_data_dir/input.vcf.gz.tbi" "index next to the input VCF"
check_index_moved "$meta_temp_dir/indexes/filtered_gz.tbi" "$output_vcf_gz.tbi" "output index file"
check_file_not_empty "$meta_temp_dir/indexes/filtered_gz.tbi" "output index file"

# Run 2: unindexed bgzipped --variant and BED --mask, no output index and the
# --interval_padding engine option
output_vcf_no_index="$meta_temp_dir/output_no_index/filtered.vcf"
mkdir -p "$(dirname "$output_vcf_no_index")"

log "Executing $meta_name with unindexed .vcf.gz --variant, BED --mask without an index, --create_output_variant_index false and --interval_padding..."
"$meta_executable" \
  --variant "$test_data_dir/no_index/input3.vcf.gz" \
  --output "$output_vcf_no_index" \
  --mask "$test_data_dir/mask.bed" \
  --mask_name "MyMask" \
  --interval_padding 10 \
  --create_output_variant_index false 2> "$meta_temp_dir/test2_no_index.log"

log "Validating TEST 2 outputs (unindexed inputs)..."
check_file_contains "$output_vcf_no_index" $'^seq1\t700\t.*\tMyMask\t' "record at pos 700 (overlaps mask) of filtered VCF (BED mask)"
check_file_contains "$output_vcf_no_index" $'^seq1\t600\t.*\tPASS\t' "record at pos 600 (outside mask) of filtered VCF (BED mask)"
check_file_contains "$meta_temp_dir/test2_no_index.log" "Warning: no index was provided for '.*no_index/input3.vcf.gz'" "run log (copy of --variant indexed)"
check_file_contains "$meta_temp_dir/test2_no_index.log" "Warning: no index was provided for '.*mask.bed'" "run log (copy of --mask indexed)"
check_file_not_exists "$test_data_dir/no_index/input3.vcf.gz.tbi" "index next to the input VCF (a copy is indexed instead)"
check_file_not_exists "$test_data_dir/mask.bed.idx" "index next to the --mask BED (a copy is indexed instead)"
check_file_not_exists "$output_vcf_no_index.idx" "output VCF index (should not exist with --create_output_variant_index false)"
check_file_contains "$output_vcf_no_index" "^##GATKCommandLine=.* --interval-padding 10 " "GATK command line (--interval_padding)"

# Run 3: bgzipped --variant and VCF --mask with separate indexes, mask and
# genotype filters
output_vcf_mask="$meta_temp_dir/output_mask/filtered.vcf"
mkdir -p "$(dirname "$output_vcf_mask")"

log "Executing $meta_name with --mask .vcf.gz, --mask_index, --mask_extension, a genotype filter and --set_filtered_genotype_to_no_call..."
"$meta_executable" \
  --variant "$test_data_dir/input3.vcf.gz" \
  --variant_index "$test_data_dir/indexes/input3.vcf.gz.tbi" \
  --output "$output_vcf_mask" \
  --output_index "$output_vcf_mask.idx" \
  --mask "$test_data_dir/mask.vcf.gz" \
  --mask_index "$test_data_dir/indexes/mask.vcf.gz.tbi" \
  --mask_name "MyMask" \
  --mask_extension 10 \
  --genotype_filter_expression "DP < 10" \
  --genotype_filter_name "lowDP" \
  --set_filtered_genotype_to_no_call 2> "$meta_temp_dir/test2_mask.log"

log "Validating TEST 2 outputs (VCF mask, genotype filter)..."
check_file_not_empty "$output_vcf_mask" "filtered VCF (mask/genotype-filter test)"
check_file_contains "$output_vcf_mask" $'^seq1\t700\t.*\tMyMask\t' "record at pos 700 (overlaps mask) of filtered VCF (VCF mask)"
check_file_contains "$output_vcf_mask" $'^seq1\t600\t.*\tPASS\t' "record at pos 600 (outside mask) of filtered VCF (VCF mask)"
check_file_not_contains "$meta_temp_dir/test2_mask.log" "no index was provided" "run log (indexes were provided)"
check_file_not_exists "$test_data_dir/mask.vcf.gz.tbi" "index next to the --mask VCF"

record500="$meta_temp_dir/record500.txt"
record600="$meta_temp_dir/record600.txt"
awk -F'\t' '$2==500 {print $10}' "$output_vcf_mask" > "$record500"
awk -F'\t' '$2==600 {print $10}' "$output_vcf_mask" > "$record600"

# Record at pos 500 (DP=5) fails the genotype filter: its genotype is set to
# no-call and the sample column carries the "lowDP" FT annotation.
check_file_contains "$record500" "\./\.:5:lowDP" "sample column at pos 500"

# Record at pos 600 (DP=50) passes the genotype filter: genotype is
# untouched and no FT annotation is added.
check_file_contains "$record600" "0/1:50" "sample column at pos 600"

# --mask_extension has no visible effect on this data
check_file_contains "$output_vcf_mask" "^##GATKCommandLine=.* --mask-extension 10 " "GATK command line (--mask_extension)"

# Run 4: BED --mask with a separate Tribble index
output_vcf_bed_index="$meta_temp_dir/output_bed_index/filtered.vcf"
mkdir -p "$(dirname "$output_vcf_bed_index")"

log "Executing $meta_name with --mask .bed and --mask_index..."
"$meta_executable" \
  --variant "$test_data_dir/input3.vcf.gz" \
  --variant_index "$test_data_dir/indexes/input3.vcf.gz.tbi" \
  --output "$output_vcf_bed_index" \
  --output_index "$output_vcf_bed_index.idx" \
  --mask "$test_data_dir/mask.bed" \
  --mask_index "$test_data_dir/bed_index/mask.bed.idx" \
  --mask_name "MyMask" 2> "$meta_temp_dir/test2_bed_index.log"

log "Validating TEST 2 outputs (BED mask with index)..."
check_file_contains "$output_vcf_bed_index" $'^seq1\t700\t.*\tMyMask\t' "record at pos 700 (overlaps mask) of filtered VCF (indexed BED mask)"
check_file_not_contains "$meta_temp_dir/test2_bed_index.log" "no index was provided" "run log (indexes were provided)"

log "✅ TEST 2 completed successfully"

# --- Test Case 3: Intervals ---
log "Starting TEST 3: Interval files, interval strings and --variant_output_filtering"

# Run 1: a plain unindexed --variant, two --intervals files and
# --exclude_intervals
output_vcf_intervals="$meta_temp_dir/output_intervals/filtered.vcf"
mkdir -p "$(dirname "$output_vcf_intervals")"

log "Executing $meta_name with a plain unindexed --variant, two --intervals files, --intervals_index and --exclude_intervals..."
"$meta_executable" \
  --variant "$intervals_input_vcf" \
  --output "$output_vcf_intervals" \
  --output_index "$output_vcf_intervals.idx" \
  --intervals "$test_data_dir/intervals_a.bed" \
  --intervals "$test_data_dir/intervals_site.vcf.gz" \
  --intervals_index "$test_data_dir/interval_indexes/intervals_site.vcf.gz.tbi" \
  --interval_set_rule UNION \
  --exclude_intervals "$test_data_dir/exclude.bed" 2> "$meta_temp_dir/test3_intervals.log"

log "Validating TEST 3 outputs (interval files)..."
check_file_exists "$output_vcf_intervals.idx" "output index (multiple intervals)"
check_file_contains "$meta_temp_dir/test3_intervals.log" "Warning: no index was provided for '.*intervals_input.vcf'" "run log (copy of --variant indexed for traversal by intervals)"
check_file_not_exists "$intervals_input_vcf.idx" "index next to the input VCF (a copy is indexed instead)"
check_vcf_column "$output_vcf_intervals" 2 "100 300 400 800" "filtered VCF (records in the selected territory)"

# Runs 2 and 3: the deletion (seq1:400-403) overlaps seq1:402-410 but does not
# start in it
log "Executing $meta_name with --interval_strings and no --variant_output_filtering..."
"$meta_executable" \
  --variant "$filtering_input_vcf" \
  --output "$meta_temp_dir/filtering_default.vcf" \
  --create_output_variant_index false \
  --interval_strings "seq1:190-210" \
  --interval_strings "seq1:402-410"

log "Executing $meta_name with --interval_strings and --variant_output_filtering STARTS_IN..."
"$meta_executable" \
  --variant "$filtering_input_vcf" \
  --output "$meta_temp_dir/filtering_starts_in.vcf" \
  --create_output_variant_index false \
  --interval_strings "seq1:190-210" \
  --interval_strings "seq1:402-410" \
  --variant_output_filtering STARTS_IN

log "Validating TEST 3 outputs (--variant_output_filtering)..."
check_vcf_column "$meta_temp_dir/filtering_default.vcf" 3 "snp deletion" "output VCF (no filtering mode, records that overlap the intervals)"
check_vcf_column "$meta_temp_dir/filtering_starts_in.vcf" 3 "snp" "output VCF (STARTS_IN, records that start in the intervals)"

log "✅ TEST 3 completed successfully"

# --- Test Case 4: Argument errors ---
log "Starting TEST 4: Argument errors"

log "Executing $meta_name with --mask_index and no --mask..."
check_fails_with "Error: --mask_index requires --mask" "$meta_temp_dir/test4_mask_index.log" \
  "$meta_executable" \
  --variant "$input_vcf" \
  --mask_index "$test_data_dir/indexes/mask.vcf.gz.tbi" \
  --output "$meta_temp_dir/index_only.vcf" \
  --output_index "$meta_temp_dir/index_only.vcf.idx"
check_file_not_exists "$meta_temp_dir/index_only.vcf" "output VCF (GATK did not run)"

log "Executing $meta_name without --output_index..."
check_fails_with "Error: --output_index is required unless --create_output_variant_index is false." "$meta_temp_dir/test4_output_index.log" \
  "$meta_executable" \
  --variant "$input_vcf" \
  --output "$meta_temp_dir/no_output_index.vcf"
check_file_not_exists "$meta_temp_dir/no_output_index.vcf" "output VCF (GATK did not run)"

log "✅ TEST 4 completed successfully"

print_test_summary "All tests"
