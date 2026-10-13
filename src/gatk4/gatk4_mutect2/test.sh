#!/bin/bash

## VIASH START
## VIASH END

# Source the centralized test helpers
source "$meta_resources_dir/test_helpers.sh"
source "$meta_resources_dir/gatk4/test_helpers.sh"

# Set up the test environment and the test data directory
gatk4_test_init

# --- Build shared reference and BAM fixtures ---
# A 500bp reference with a known SNV at chr1:250 (T>A). The tumor BAM has the
# SNV at ~50% allele fraction, the normal BAM has only reference reads.
create_test_somatic_reference "$test_data_dir/reference.fasta"
create_test_somatic_bam "$test_data_dir/tumor.bam" tumor rg_tumor mixed
create_test_somatic_bam "$test_data_dir/normal.bam" normal rg_normal ref
check_file_exists "$test_data_dir/tumor.bai" "tumor BAM index"
check_file_exists "$test_data_dir/normal.bai" "normal BAM index"

# The same tumor reads with mapping quality 15, below the tool-default
# minimum of 20
awk 'BEGIN { OFS = "\t" } /^@/ { print; next } { $5 = 15; print }' "$test_data_dir/tumor.sam" > "$test_data_dir/tumor_mq15.sam"
sort_and_index_bam "$test_data_dir/tumor_mq15.sam" "$test_data_dir/tumor_mq15.bam"

reference_args=(
  --reference "$test_data_dir/reference.fasta"
  --reference_fai "$test_data_dir/reference.fasta.fai"
  --reference_dict "$test_data_dir/reference.dict"
)

# Regex for a VCF record of the somatic test SNV at chr1:250 (T>A)
snv_250=$'^chr1\t250\t.\tT\tA\t'

# --- Test Case 1: Tumor-normal mode ---
log "Starting TEST 1: Tumor-normal mode"

log "Executing $meta_name with a tumor and a normal BAM..."
"$meta_executable" \
  --input "$test_data_dir/tumor.bam" \
  --read_index "$test_data_dir/tumor.bai" \
  --input "$test_data_dir/normal.bam" \
  --read_index "$test_data_dir/normal.bai" \
  --normal_sample normal \
  "${reference_args[@]}" \
  --output "$meta_temp_dir/tumor_normal.vcf" \
  --output_stats "$meta_temp_dir/tumor_normal.stats" \
  --output_index "$meta_temp_dir/tumor_normal.idx" \
  --create_f1r2_tar_gz \
  --f1r2_tar_gz "$meta_temp_dir/tumor_normal_f1r2.tar.gz"

log "Validating TEST 1 outputs..."
check_file_exists "$meta_temp_dir/tumor_normal.vcf" "output VCF file"
check_file_contains "$meta_temp_dir/tumor_normal.vcf" "^##fileformat=VCF" "output VCF file header"
check_file_contains "$meta_temp_dir/tumor_normal.vcf" "^##normal_sample=normal" "output VCF normal sample header"
check_file_contains "$meta_temp_dir/tumor_normal.vcf" "^##tumor_sample=tumor" "output VCF tumor sample header"
check_vcf_has_record "$meta_temp_dir/tumor_normal.vcf" "$snv_250" "output VCF file (chr1:250 T>A call)"
check_file_exists "$meta_temp_dir/tumor_normal.stats" "output stats file"
check_file_contains "$meta_temp_dir/tumor_normal.stats" "callable" "output stats file"
check_file_not_exists "$meta_temp_dir/tumor_normal.vcf.stats" "stats file next to the output VCF (moved to --output_stats)"
check_index_moved "$meta_temp_dir/tumor_normal.idx" "$meta_temp_dir/tumor_normal.vcf.idx" "output index file (Tribble index for a .vcf)"
check_file_exists "$meta_temp_dir/tumor_normal_f1r2.tar.gz" "output F1R2 tarball"
check_file_not_empty "$meta_temp_dir/tumor_normal_f1r2.tar.gz" "output F1R2 tarball"

log "✅ TEST 1 completed successfully"

# --- Test Case 2: Tumor-only mode ---
log "Starting TEST 2: Tumor-only mode"

log "Executing $meta_name with only a tumor BAM, --output_stats next to --output and optional output paths without --create_*..."
mkdir -p "$meta_temp_dir/tumor_only"
"$meta_executable" \
  --input "$test_data_dir/tumor.bam" \
  --read_index "$test_data_dir/tumor.bai" \
  "${reference_args[@]}" \
  --output "$meta_temp_dir/tumor_only/tumor_only.vcf.gz" \
  --output_stats "$meta_temp_dir/tumor_only/tumor_only.vcf.gz.stats" \
  --output_index "$meta_temp_dir/tumor_only/tumor_only.tbi" \
  --f1r2_tar_gz "$meta_temp_dir/tumor_only/f1r2.tar.gz" \
  --bam_output "$meta_temp_dir/tumor_only/bamout.bam"

log "Validating TEST 2 outputs..."
check_file_exists "$meta_temp_dir/tumor_only/tumor_only.vcf.gz" "output VCF file"
zcat "$meta_temp_dir/tumor_only/tumor_only.vcf.gz" > "$meta_temp_dir/tumor_only/tumor_only.vcf"
check_file_not_contains "$meta_temp_dir/tumor_only/tumor_only.vcf" "^##normal_sample=" "tumor-only output VCF normal sample header"
check_file_contains "$meta_temp_dir/tumor_only/tumor_only.vcf" "^##tumor_sample=tumor" "output VCF tumor sample header"
check_vcf_has_record "$meta_temp_dir/tumor_only/tumor_only.vcf" "$snv_250" "output VCF file (chr1:250 T>A call)"
check_file_exists "$meta_temp_dir/tumor_only/tumor_only.vcf.gz.stats" "output stats file"
check_index_moved "$meta_temp_dir/tumor_only/tumor_only.tbi" "$meta_temp_dir/tumor_only/tumor_only.vcf.gz.tbi" "output index file (tabix index for a .vcf.gz)"

log "Checking that only the expected files were written..."
expected_files=$'tumor_only.tbi\ntumor_only.vcf\ntumor_only.vcf.gz\ntumor_only.vcf.gz.stats'
found_files=$(find "$meta_temp_dir/tumor_only" -type f -exec basename {} \; | sort)
if [[ "$found_files" == "$expected_files" ]]; then
  log "✓ Output directory contains only the expected files"
else
  log_error "✗ Output directory contains unexpected files: $(echo "$found_files" | tr '\n' ' ')"
  exit 1
fi

log "✅ TEST 2 completed successfully"

# --- Test Case 3: Germline resource, panel of normals and intervals ---
log "Starting TEST 3: Germline resource, panel of normals and intervals"

log "Writing germline resource and panel of normals VCFs..."
cat > "$test_data_dir/germline_resource.vcf" <<'VCF'
##fileformat=VCFv4.2
##contig=<ID=chr1,length=500>
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	250	.	T	A	.	.	AF=0.001
VCF
mkdir -p "$test_data_dir/indexed"
cat > "$test_data_dir/indexed/pon.vcf" <<'VCF'
##fileformat=VCFv4.2
##contig=<ID=chr1,length=500>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	400	.	A	G	.	.	.
VCF
gatk IndexFeatureFile --input "$test_data_dir/indexed/pon.vcf" --verbosity ERROR
printf 'chr1\t199\t300\n' > "$test_data_dir/intervals.bed"

log "Executing $meta_name with an unindexed --germline_resource, an indexed --panel_of_normals and --intervals..."
"$meta_executable" \
  --input "$test_data_dir/tumor.bam" \
  --read_index "$test_data_dir/tumor.bai" \
  --input "$test_data_dir/normal.bam" \
  --read_index "$test_data_dir/normal.bai" \
  --normal_sample normal \
  "${reference_args[@]}" \
  --germline_resource "$test_data_dir/germline_resource.vcf" \
  --panel_of_normals "$test_data_dir/indexed/pon.vcf" \
  --panel_of_normals_index "$test_data_dir/indexed/pon.vcf.idx" \
  --intervals "$test_data_dir/intervals.bed" \
  --output "$meta_temp_dir/resources.vcf" \
  --output_stats "$meta_temp_dir/resources.stats" \
  --output_index "$meta_temp_dir/resources.idx" 2>&1 | tee "$meta_temp_dir/test3.log"

log "Validating TEST 3 outputs..."
check_vcf_has_record "$meta_temp_dir/resources.vcf" "$snv_250" "output VCF file (chr1:250 T>A call)"
check_file_contains "$meta_temp_dir/resources.vcf" "germline-resource" "output VCF command line (germline resource used)"
check_file_contains "$meta_temp_dir/resources.vcf" "panel-of-normals" "output VCF command line (panel of normals used)"
check_file_contains "$meta_temp_dir/test3.log" "Warning: no index was provided for '.*germline_resource.vcf'" "run log (copy warning for the germline resource)"
check_file_not_contains "$meta_temp_dir/test3.log" "Warning: no index was provided for '.*pon.vcf'" "run log (no copy warning for the indexed panel of normals)"

log "✅ TEST 3 completed successfully"

# --- Test Case 4: Bamout and calling options ---
log "Starting TEST 4: Bamout and calling options"

log "Executing $meta_name with --bam_output and calling, annotation and assembly options..."
"$meta_executable" \
  --input "$test_data_dir/tumor.bam" \
  --read_index "$test_data_dir/tumor.bai" \
  --input "$test_data_dir/normal.bam" \
  --read_index "$test_data_dir/normal.bai" \
  --normal_sample normal \
  "${reference_args[@]}" \
  --output "$meta_temp_dir/options.vcf" \
  --output_stats "$meta_temp_dir/options.stats" \
  --output_index "$meta_temp_dir/options.idx" \
  --create_bam_output \
  --bam_output "$meta_temp_dir/bamout.bam" \
  --bam_output_index "$meta_temp_dir/bamout_index.bai" \
  --bam_writer_type ALL_POSSIBLE_HAPLOTYPES \
  --max_mnp_distance 0 \
  --tumor_lod_to_emit 3.0 \
  --normal_lod 2.2 \
  --min_base_quality_score 10 \
  --pcr_indel_model NONE \
  --kmer_size 10 \
  --kmer_size 25 \
  --annotation OrientationBiasReadCounts \
  --callable_depth 5 \
  --dont_use_soft_clipped_bases

log "Validating TEST 4 outputs..."
check_vcf_has_record "$meta_temp_dir/options.vcf" "$snv_250" "output VCF file (chr1:250 T>A call)"
check_file_exists "$meta_temp_dir/bamout.bam" "output bamout file"
check_file_not_empty "$meta_temp_dir/bamout.bam" "output bamout file"
check_index_moved "$meta_temp_dir/bamout_index.bai" "$meta_temp_dir/bamout.bai" "output bamout index file"
check_file_contains "$meta_temp_dir/options.vcf" "max-mnp-distance 0" "output VCF command line (--max_mnp_distance)"
check_file_contains "$meta_temp_dir/options.vcf" "pcr-indel-model NONE" "output VCF command line (--pcr_indel_model)"

log "✅ TEST 4 completed successfully"

# --- Test Case 5: Minimum mapping quality ---
log "Starting TEST 5: Minimum mapping quality"

# The read filter counts in the run log show the effect of the setting. (The
# MQ 15 reads do not give a call even when they are kept, because Mutect2 caps
# base qualities to the mapping quality.)
log "Executing $meta_name on MQ 15 reads with the default --minimum_mapping_quality (20)..."
"$meta_executable" \
  --input "$test_data_dir/tumor_mq15.bam" \
  --read_index "$test_data_dir/tumor_mq15.bai" \
  "${reference_args[@]}" \
  --output "$meta_temp_dir/mq_default.vcf" \
  --output_stats "$meta_temp_dir/mq_default.stats" \
  --output_index "$meta_temp_dir/mq_default.idx" > "$meta_temp_dir/test5_default.log" 2>&1

check_file_contains "$meta_temp_dir/test5_default.log" "16 read(s) filtered by: MappingQualityReadFilter" "run log (all 16 MQ 15 reads filtered)"

log "Executing $meta_name on MQ 15 reads with --minimum_mapping_quality 10..."
"$meta_executable" \
  --input "$test_data_dir/tumor_mq15.bam" \
  --read_index "$test_data_dir/tumor_mq15.bai" \
  "${reference_args[@]}" \
  --minimum_mapping_quality 10 \
  --output "$meta_temp_dir/mq10.vcf" \
  --output_stats "$meta_temp_dir/mq10.stats" \
  --output_index "$meta_temp_dir/mq10.idx" > "$meta_temp_dir/test5_mq10.log" 2>&1

check_file_contains "$meta_temp_dir/test5_mq10.log" "0 read(s) filtered by: MappingQualityReadFilter" "run log (no MQ 15 reads filtered)"

log "✅ TEST 5 completed successfully"

# --- Test Case 6: Invalid argument combinations ---
log "Starting TEST 6: Invalid argument combinations fail"

check_fails_with "Error: --read_index must be given once for each file in --input" "$meta_temp_dir/test6.log" \
  "$meta_executable" \
  --input "$test_data_dir/tumor.bam" \
  --input "$test_data_dir/normal.bam" \
  --read_index "$test_data_dir/tumor.bai" \
  "${reference_args[@]}" \
  --output "$meta_temp_dir/mismatch.vcf" \
  --output_stats "$meta_temp_dir/mismatch.stats" \
  --output_index "$meta_temp_dir/mismatch.idx"

check_fails_with "Error: --create_f1r2_tar_gz requires --f1r2_tar_gz" "$meta_temp_dir/test6_f1r2.log" \
  "$meta_executable" \
  --input "$test_data_dir/tumor.bam" \
  --read_index "$test_data_dir/tumor.bai" \
  "${reference_args[@]}" \
  --output "$meta_temp_dir/no_f1r2.vcf" \
  --output_stats "$meta_temp_dir/no_f1r2.stats" \
  --output_index "$meta_temp_dir/no_f1r2.idx" \
  --create_f1r2_tar_gz

# The arguments of a tumor-only run that the cases below add to
fail_args=(
  --input "$test_data_dir/tumor.bam"
  --read_index "$test_data_dir/tumor.bai"
  "${reference_args[@]}"
  --output "$meta_temp_dir/fail.vcf"
  --output_stats "$meta_temp_dir/fail.stats"
)

check_fails_with "Error: --output_index is required unless" "$meta_temp_dir/test6_output_index.log" \
  "$meta_executable" "${fail_args[@]}"
check_fails_with "Error: --create_bam_output requires --bam_output." "$meta_temp_dir/test6_bam_output.log" \
  "$meta_executable" "${fail_args[@]}" \
  --output_index "$meta_temp_dir/fail.idx" --create_bam_output
check_fails_with "Error: --bam_output_index is required unless --create_output_bam_index is false." "$meta_temp_dir/test6_bam_output_index.log" \
  "$meta_executable" "${fail_args[@]}" \
  --output_index "$meta_temp_dir/fail.idx" --create_bam_output --bam_output "$meta_temp_dir/fail.bam"
check_fails_with "Error: --bam_output must be a BAM file, GATK cannot write CRAM here." "$meta_temp_dir/test6_bam_output_cram.log" \
  "$meta_executable" "${fail_args[@]}" \
  --output_index "$meta_temp_dir/fail.idx" --create_bam_output --bam_output "$meta_temp_dir/fail.cram" \
  --bam_output_index "$meta_temp_dir/fail.crai"
check_file_not_contains "$meta_temp_dir/test6_bam_output_cram.log" "Mutect2 done" "run log (GATK did not run)"
check_fails_with "Error: --alleles_index requires --alleles." "$meta_temp_dir/test6_alleles_index.log" \
  "$meta_executable" "${fail_args[@]}" \
  --output_index "$meta_temp_dir/fail.idx" --alleles_index "$test_data_dir/indexed/pon.vcf.idx"

log "✅ TEST 6 completed successfully"

# --- Test Case 7: Output index without index creation ---
log "Starting TEST 7: --output_index is ignored when --create_output_variant_index is false"

"$meta_executable" \
  --input "$test_data_dir/tumor.bam" \
  --read_index "$test_data_dir/tumor.bai" \
  "${reference_args[@]}" \
  --output "$meta_temp_dir/no_index.vcf.bgz" \
  --output_stats "$meta_temp_dir/no_index.stats" \
  --output_index "$meta_temp_dir/no_index.tbi" \
  --create_output_variant_index false

check_file_exists "$meta_temp_dir/no_index.vcf.bgz" "output VCF file"
check_file_not_exists "$meta_temp_dir/no_index.tbi" "output index file"
check_file_not_exists "$meta_temp_dir/no_index.vcf.bgz.tbi" "index next to the output VCF"

log "Executing $meta_name with a .vcf.bgz output and --output_index..."
"$meta_executable" \
  --input "$test_data_dir/tumor.bam" \
  --read_index "$test_data_dir/tumor.bai" \
  "${reference_args[@]}" \
  --output "$meta_temp_dir/bgz.vcf.bgz" \
  --output_stats "$meta_temp_dir/bgz.stats" \
  --output_index "$meta_temp_dir/bgz.tbi"

check_file_exists "$meta_temp_dir/bgz.tbi" "output index file (tabix index for a .vcf.bgz)"

log "✅ TEST 7 completed successfully"

# --- Test Case 8: Force-calling alleles ---
log "Starting TEST 8: Force-calling alleles with --alleles"

# A site with no evidence in the reads, which is only in the output because it
# is force-called
cat > "$test_data_dir/alleles.vcf" <<'VCF'
##fileformat=VCFv4.2
##contig=<ID=chr1,length=500>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	300	.	A	C	.	.	.
VCF
cp "$test_data_dir/alleles.vcf" "$test_data_dir/indexed/alleles.vcf"
gatk IndexFeatureFile --input "$test_data_dir/indexed/alleles.vcf" --verbosity ERROR

log "Executing $meta_name with --alleles and --alleles_index..."
"$meta_executable" \
  --input "$test_data_dir/tumor.bam" \
  --read_index "$test_data_dir/tumor.bai" \
  "${reference_args[@]}" \
  --alleles "$test_data_dir/indexed/alleles.vcf" \
  --alleles_index "$test_data_dir/indexed/alleles.vcf.idx" \
  --output "$meta_temp_dir/alleles.vcf" \
  --output_stats "$meta_temp_dir/alleles.stats" \
  --output_index "$meta_temp_dir/alleles.idx"

log "Validating TEST 8 outputs..."
check_vcf_has_record "$meta_temp_dir/alleles.vcf" "$snv_250" "output VCF file (chr1:250 T>A call)"
check_vcf_has_record "$meta_temp_dir/alleles.vcf" $'^chr1\t300\t' "output VCF file (force-called chr1:300 allele)"

log "✅ TEST 8 completed successfully"

# --- Test Case 9: Interval padding and merging rule ---
log "Starting TEST 9: --interval_padding and --interval_merging_rule"

# Positions 241-245, which end 5bp before the SNV at position 250
printf 'chr1\t240\t245\n' > "$test_data_dir/near_250.bed"

log "Executing $meta_name with an --intervals file next to the SNV, without padding..."
"$meta_executable" \
  --input "$test_data_dir/tumor.bam" \
  --read_index "$test_data_dir/tumor.bai" \
  "${reference_args[@]}" \
  --intervals "$test_data_dir/near_250.bed" \
  --output "$meta_temp_dir/no_padding.vcf" \
  --output_stats "$meta_temp_dir/no_padding.stats" \
  --output_index "$meta_temp_dir/no_padding.idx"

check_vcf_lacks_record "$meta_temp_dir/no_padding.vcf" $'^chr1\t250\t' "output VCF file (chr1:250 SNV outside the interval)"

log "Executing $meta_name with the same --intervals file, --interval_padding 5 and --interval_merging_rule OVERLAPPING_ONLY..."
"$meta_executable" \
  --input "$test_data_dir/tumor.bam" \
  --read_index "$test_data_dir/tumor.bai" \
  "${reference_args[@]}" \
  --intervals "$test_data_dir/near_250.bed" \
  --interval_padding 5 \
  --interval_merging_rule OVERLAPPING_ONLY \
  --output "$meta_temp_dir/padding.vcf" \
  --output_stats "$meta_temp_dir/padding.stats" \
  --output_index "$meta_temp_dir/padding.idx"

check_vcf_has_record "$meta_temp_dir/padding.vcf" "$snv_250" "output VCF file (padded interval) (chr1:250 T>A call)"
check_file_contains "$meta_temp_dir/padding.vcf" "interval-padding 5" "output VCF command line (--interval_padding)"
check_file_contains "$meta_temp_dir/padding.vcf" "interval-merging-rule OVERLAPPING_ONLY" "output VCF command line (--interval_merging_rule)"

log "✅ TEST 9 completed successfully"

# --- Test Case 10: Tumor CRAM and normal BAM ---
log "Starting TEST 10: A tumor CRAM with a --read_index in another directory and a normal BAM"

create_test_cram "$test_data_dir/tumor.bam" "$test_data_dir/reference.fasta" \
  "$test_data_dir/cram/tumor.cram" "$test_data_dir/cram_index/tumor_index.crai"

"$meta_executable" \
  --input "$test_data_dir/cram/tumor.cram" \
  --read_index "$test_data_dir/cram_index/tumor_index.crai" \
  --input "$test_data_dir/normal.bam" \
  --read_index "$test_data_dir/normal.bai" \
  --normal_sample normal \
  "${reference_args[@]}" \
  --output "$meta_temp_dir/cram_normal.vcf" \
  --output_stats "$meta_temp_dir/cram_normal.stats" \
  --output_index "$meta_temp_dir/cram_normal.idx"

log "Validating TEST 10 outputs..."
check_file_contains "$meta_temp_dir/cram_normal.vcf" "^##normal_sample=normal" "output VCF normal sample header (CRAM input)"
check_file_contains "$meta_temp_dir/cram_normal.vcf" "^##tumor_sample=tumor" "output VCF tumor sample header (CRAM input)"
check_vcf_has_record "$meta_temp_dir/cram_normal.vcf" "$snv_250" "output VCF file (CRAM input) (chr1:250 T>A call)"
check_file_exists "$meta_temp_dir/cram_normal.idx" "output index file (CRAM input)"

log "✅ TEST 10 completed successfully"

# --- Test Case 11: Read length read filter ---
log "Starting TEST 11: --min_read_length"

# All tumor reads are 100bp, so a minimum length of 101 removes all of them
log "Executing $meta_name with --min_read_length 101..."
"$meta_executable" \
  --input "$test_data_dir/tumor.bam" \
  --read_index "$test_data_dir/tumor.bai" \
  "${reference_args[@]}" \
  --min_read_length 101 \
  --output "$meta_temp_dir/min_read_length.vcf" \
  --output_stats "$meta_temp_dir/min_read_length.stats" \
  --output_index "$meta_temp_dir/min_read_length.idx" > "$meta_temp_dir/test11.log" 2>&1

log "Validating TEST 11 outputs..."
check_file_contains "$meta_temp_dir/test11.log" "16 read(s) filtered by: ReadLengthReadFilter" "run log (all 16 reads filtered)"
check_vcf_record_count "$meta_temp_dir/min_read_length.vcf" 0 "output VCF file (all reads filtered)"

log "✅ TEST 11 completed successfully"

# --- Test Case 12: Other read filter, genotyping and annotation options ---
log "Starting TEST 12: --max_read_length, --maximum_mapping_quality, --genotype_germline_sites_fraction and annotation options"

log "Executing $meta_name with the options..."
"$meta_executable" \
  --input "$test_data_dir/tumor.bam" \
  --read_index "$test_data_dir/tumor.bai" \
  "${reference_args[@]}" \
  --min_read_length 50 \
  --max_read_length 100 \
  --maximum_mapping_quality 60 \
  --genotype_germline_sites \
  --genotype_germline_sites_fraction 0.5 \
  --enable_dynamic_read_disqualification_for_genotyping \
  --disable_tool_default_annotations \
  --output "$meta_temp_dir/options2.vcf" \
  --output_stats "$meta_temp_dir/options2.stats" \
  --output_index "$meta_temp_dir/options2.idx" > "$meta_temp_dir/test12.log" 2>&1

log "Validating TEST 12 outputs..."
check_vcf_has_record "$meta_temp_dir/options2.vcf" "$snv_250" "output VCF file (options) (chr1:250 T>A call)"
check_file_contains "$meta_temp_dir/test12.log" "0 read(s) filtered by: ReadLengthReadFilter" "run log (no reads filtered by length)"
check_file_contains "$meta_temp_dir/test12.log" "0 read(s) filtered by: MappingQualityReadFilter" "run log (no reads filtered by mapping quality)"
check_file_contains "$meta_temp_dir/options2.vcf" "max-read-length 100" "output VCF command line (--max_read_length)"
check_file_contains "$meta_temp_dir/options2.vcf" "maximum-mapping-quality 60" "output VCF command line (--maximum_mapping_quality)"
check_file_contains "$meta_temp_dir/options2.vcf" "genotype-germline-sites-fraction 0.5" "output VCF command line (--genotype_germline_sites_fraction)"
check_file_contains "$meta_temp_dir/options2.vcf" "enable-dynamic-read-disqualification-for-genotyping true" "output VCF command line (--enable_dynamic_read_disqualification_for_genotyping)"
check_file_contains "$meta_temp_dir/options2.vcf" "disable-tool-default-annotations true" "output VCF command line (--disable_tool_default_annotations)"

log "✅ TEST 12 completed successfully"

print_test_summary "All tests"
