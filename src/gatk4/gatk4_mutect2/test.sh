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

# Check that a VCF has the somatic test SNV at chr1:250
check_snv_called() {
  local vcf_path="$1"
  local description="$2"
  if grep -v '^#' "$vcf_path" | grep -q $'^chr1\t250\t.\tT\tA\t'; then
    log "✓ $description has the chr1:250 T>A call"
  else
    log_error "✗ $description does not have the chr1:250 T>A call"
    exit 1
  fi
}

# --- Test Case 1: Tumor-normal mode ---
log "Starting TEST 1: Tumor-normal mode"

log "Executing $meta_name with a tumor and a normal BAM..."
"$meta_executable" \
  --input "$test_data_dir/tumor.bam" \
  --bai "$test_data_dir/tumor.bai" \
  --input "$test_data_dir/normal.bam" \
  --bai "$test_data_dir/normal.bai" \
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
check_snv_called "$meta_temp_dir/tumor_normal.vcf" "output VCF file"
check_file_exists "$meta_temp_dir/tumor_normal.stats" "output stats file"
check_file_contains "$meta_temp_dir/tumor_normal.stats" "callable" "output stats file"
check_file_not_exists "$meta_temp_dir/tumor_normal.vcf.stats" "stats file next to the output VCF (moved to --output_stats)"
check_file_exists "$meta_temp_dir/tumor_normal.idx" "output index file (Tribble index for a .vcf)"
check_file_not_exists "$meta_temp_dir/tumor_normal.vcf.idx" "index next to the output VCF (moved to --output_index)"
check_file_exists "$meta_temp_dir/tumor_normal_f1r2.tar.gz" "output F1R2 tarball"
check_file_not_empty "$meta_temp_dir/tumor_normal_f1r2.tar.gz" "output F1R2 tarball"

log "✅ TEST 1 completed successfully"

# --- Test Case 2: Tumor-only mode ---
log "Starting TEST 2: Tumor-only mode"

log "Executing $meta_name with only a tumor BAM, --output_stats next to --output and optional output paths without --create_*..."
mkdir -p "$meta_temp_dir/tumor_only"
"$meta_executable" \
  --input "$test_data_dir/tumor.bam" \
  --bai "$test_data_dir/tumor.bai" \
  "${reference_args[@]}" \
  --output "$meta_temp_dir/tumor_only/tumor_only.vcf.gz" \
  --output_stats "$meta_temp_dir/tumor_only/tumor_only.vcf.gz.stats" \
  --output_index "$meta_temp_dir/tumor_only/tumor_only.tbi" \
  --f1r2_tar_gz "$meta_temp_dir/tumor_only/f1r2.tar.gz" \
  --bam_output "$meta_temp_dir/tumor_only/bamout.bam"

log "Validating TEST 2 outputs..."
check_file_exists "$meta_temp_dir/tumor_only/tumor_only.vcf.gz" "output VCF file"
if zcat "$meta_temp_dir/tumor_only/tumor_only.vcf.gz" | grep -q '^##normal_sample='; then
  log_error "✗ Tumor-only output VCF has a normal sample header"
  exit 1
fi
zcat "$meta_temp_dir/tumor_only/tumor_only.vcf.gz" > "$meta_temp_dir/tumor_only/tumor_only.vcf"
check_file_contains "$meta_temp_dir/tumor_only/tumor_only.vcf" "^##tumor_sample=tumor" "output VCF tumor sample header"
check_snv_called "$meta_temp_dir/tumor_only/tumor_only.vcf" "output VCF file"
check_file_exists "$meta_temp_dir/tumor_only/tumor_only.vcf.gz.stats" "output stats file"
check_file_exists "$meta_temp_dir/tumor_only/tumor_only.tbi" "output index file (tabix index for a .vcf.gz)"
check_file_not_exists "$meta_temp_dir/tumor_only/tumor_only.vcf.gz.tbi" "index next to the output VCF (moved to --output_index)"

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
  --bai "$test_data_dir/tumor.bai" \
  --input "$test_data_dir/normal.bam" \
  --bai "$test_data_dir/normal.bai" \
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
check_snv_called "$meta_temp_dir/resources.vcf" "output VCF file"
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
  --bai "$test_data_dir/tumor.bai" \
  --input "$test_data_dir/normal.bam" \
  --bai "$test_data_dir/normal.bai" \
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
check_snv_called "$meta_temp_dir/options.vcf" "output VCF file"
check_file_exists "$meta_temp_dir/bamout.bam" "output bamout file"
check_file_not_empty "$meta_temp_dir/bamout.bam" "output bamout file"
check_file_exists "$meta_temp_dir/bamout_index.bai" "output bamout index file"
check_file_not_exists "$meta_temp_dir/bamout.bai" "index next to the bamout file (moved to --bam_output_index)"
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
  --bai "$test_data_dir/tumor_mq15.bai" \
  "${reference_args[@]}" \
  --output "$meta_temp_dir/mq_default.vcf" \
  --output_stats "$meta_temp_dir/mq_default.stats" \
  --output_index "$meta_temp_dir/mq_default.idx" > "$meta_temp_dir/test5_default.log" 2>&1

check_file_contains "$meta_temp_dir/test5_default.log" "16 read(s) filtered by: MappingQualityReadFilter" "run log (all 16 MQ 15 reads filtered)"

log "Executing $meta_name on MQ 15 reads with --minimum_mapping_quality 10..."
"$meta_executable" \
  --input "$test_data_dir/tumor_mq15.bam" \
  --bai "$test_data_dir/tumor_mq15.bai" \
  "${reference_args[@]}" \
  --minimum_mapping_quality 10 \
  --output "$meta_temp_dir/mq10.vcf" \
  --output_stats "$meta_temp_dir/mq10.stats" \
  --output_index "$meta_temp_dir/mq10.idx" > "$meta_temp_dir/test5_mq10.log" 2>&1

check_file_contains "$meta_temp_dir/test5_mq10.log" "0 read(s) filtered by: MappingQualityReadFilter" "run log (no MQ 15 reads filtered)"

log "✅ TEST 5 completed successfully"

# --- Test Case 6: --input and --bai count mismatch ---
log "Starting TEST 6: Invalid argument combinations fail"

if "$meta_executable" \
  --input "$test_data_dir/tumor.bam" \
  --input "$test_data_dir/normal.bam" \
  --bai "$test_data_dir/tumor.bai" \
  "${reference_args[@]}" \
  --output "$meta_temp_dir/mismatch.vcf" \
  --output_stats "$meta_temp_dir/mismatch.stats" \
  --output_index "$meta_temp_dir/mismatch.idx" > "$meta_temp_dir/test6.log" 2>&1; then
  log_error "✗ $meta_name did not fail with 2 --input files and 1 --bai file"
  exit 1
fi
check_file_contains "$meta_temp_dir/test6.log" "Error: --input and --bai must be given the same number of times" "error message"

if "$meta_executable" \
  --input "$test_data_dir/tumor.bam" \
  --bai "$test_data_dir/tumor.bai" \
  "${reference_args[@]}" \
  --output "$meta_temp_dir/no_f1r2.vcf" \
  --output_stats "$meta_temp_dir/no_f1r2.stats" \
  --output_index "$meta_temp_dir/no_f1r2.idx" \
  --create_f1r2_tar_gz > "$meta_temp_dir/test6_f1r2.log" 2>&1; then
  log_error "✗ $meta_name did not fail with --create_f1r2_tar_gz and no --f1r2_tar_gz"
  exit 1
fi
check_file_contains "$meta_temp_dir/test6_f1r2.log" "Error: --create_f1r2_tar_gz requires --f1r2_tar_gz" "error message"

# Run $meta_executable with extra arguments and check that it fails with a
# message
check_fails_with() {
  local message="$1"
  local log_file="$2"
  shift 2
  if "$meta_executable" \
    --input "$test_data_dir/tumor.bam" \
    --bai "$test_data_dir/tumor.bai" \
    "${reference_args[@]}" \
    --output "$meta_temp_dir/fail.vcf" \
    --output_stats "$meta_temp_dir/fail.stats" \
    "$@" > "$log_file" 2>&1; then
    log_error "✗ $meta_name did not fail with: $*"
    exit 1
  fi
  check_file_contains "$log_file" "$message" "error message"
}

check_fails_with "Error: --output_index is required unless" "$meta_temp_dir/test6_output_index.log"
check_fails_with "Error: --create_bam_output requires --bam_output\." "$meta_temp_dir/test6_bam_output.log" \
  --output_index "$meta_temp_dir/fail.idx" --create_bam_output
check_fails_with "Error: --create_bam_output requires --bam_output_index" "$meta_temp_dir/test6_bam_output_index.log" \
  --output_index "$meta_temp_dir/fail.idx" --create_bam_output --bam_output "$meta_temp_dir/fail.bam"

log "✅ TEST 6 completed successfully"

# --- Test Case 7: Output index without index creation ---
log "Starting TEST 7: --output_index is ignored when --create_output_variant_index is false"

"$meta_executable" \
  --input "$test_data_dir/tumor.bam" \
  --bai "$test_data_dir/tumor.bai" \
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
  --bai "$test_data_dir/tumor.bai" \
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
  --bai "$test_data_dir/tumor.bai" \
  "${reference_args[@]}" \
  --alleles "$test_data_dir/indexed/alleles.vcf" \
  --alleles_index "$test_data_dir/indexed/alleles.vcf.idx" \
  --output "$meta_temp_dir/alleles.vcf" \
  --output_stats "$meta_temp_dir/alleles.stats" \
  --output_index "$meta_temp_dir/alleles.idx"

log "Validating TEST 8 outputs..."
check_snv_called "$meta_temp_dir/alleles.vcf" "output VCF file"
if grep -v '^#' "$meta_temp_dir/alleles.vcf" | grep -q $'^chr1\t300\t'; then
  log "✓ The force-called chr1:300 allele is in the output"
else
  log_error "✗ The force-called chr1:300 allele is not in the output"
  exit 1
fi

log "✅ TEST 8 completed successfully"

print_test_summary "All tests"
