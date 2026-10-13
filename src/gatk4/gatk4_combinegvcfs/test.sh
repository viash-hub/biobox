#!/bin/bash

## VIASH START
## VIASH END

# Source the centralized test helpers
source "$meta_resources_dir/test_helpers.sh"
source "$meta_resources_dir/gatk4/test_helpers.sh"

gatk4_test_init

# --- Shared fixtures ---
create_test_reference "$test_data_dir/reference.fasta" 1 2000
check_file_exists "$test_data_dir/reference.fasta" "test reference genome"
check_file_exists "$test_data_dir/reference.dict" "reference sequence dictionary"
check_file_exists "$test_data_dir/reference.fasta.fai" "reference FASTA index"

reference_args=(
  --reference "$test_data_dir/reference.fasta"
  --reference_fai "$test_data_dir/reference.fasta.fai"
  --reference_dict "$test_data_dir/reference.dict"
)

# Two per-sample GVCFs (plain VCF) from a direct HaplotypeCaller run
for sample in sample1 sample2; do
  create_test_sam_reads "$test_data_dir/$sample.sam" seq1 2000 "rg_$sample" "$sample"
  sort_and_index_bam "$test_data_dir/$sample.sam" "$test_data_dir/$sample.sorted.bam"
  check_file_exists "$test_data_dir/$sample.sorted.bai" "BAM index file for $sample"

  log "Running HaplotypeCaller on $sample to produce a GVCF..."
  gatk HaplotypeCaller \
    --input "$test_data_dir/$sample.sorted.bam" \
    --output "$test_data_dir/$sample.g.vcf" \
    --reference "$test_data_dir/reference.fasta" \
    --emit-ref-confidence GVCF \
    --verbosity ERROR
  check_file_not_empty "$test_data_dir/$sample.g.vcf" "$sample GVCF"
done

# A dbSNP VCF (only the bgzipped copies are used, so no index)
create_test_known_sites_vcf --no-index "$test_data_dir/dbsnp.vcf" seq1 2000 100

# Bgzipped inputs with their indexes in a different directory
for name in sample1.g sample2.g dbsnp; do
  create_test_vcf_gz_with_separate_index "$test_data_dir/$name.vcf" \
    "$test_data_dir/separate/$name.vcf.gz" "$test_data_dir/separate_index/$name.vcf.gz.tbi"
done

# Bgzipped inputs without indexes
mkdir -p "$test_data_dir/no_index"
for name in sample1.g sample2.g dbsnp; do
  bgzip -c "$test_data_dir/$name.vcf" > "$test_data_dir/no_index/$name.vcf.gz"
done

# Interval files: include seq1:100-200 (BED is 0-based) and seq1:1001-1020
# (the 20bp REF allele of the VCF record), exclude seq1:150-160
mkdir -p "$test_data_dir/intervals"
printf 'seq1\t99\t200\n' > "$test_data_dir/intervals/include.bed"
printf 'seq1\t149\t160\n' > "$test_data_dir/intervals/exclude.bed"
{
  echo "##fileformat=VCFv4.2"
  echo "##contig=<ID=seq1,length=2000>"
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
  printf 'seq1\t1001\t.\tATCGATCGATCGATCGATCG\tA\t50\tPASS\t.\n'
} > "$test_data_dir/intervals/include.vcf"
create_test_vcf_gz_with_separate_index "$test_data_dir/intervals/include.vcf" \
  "$test_data_dir/intervals/include.vcf.gz" "$test_data_dir/intervals_index/include.vcf.gz.tbi"

# A trio pedigree: sample1 is the child of sample2 and a parent that is not in
# the GVCFs
printf 'fam1\tsample1\tsample2\tmother\t1\t0\nfam1\tsample2\t0\t0\t1\t0\nfam1\tmother\t0\t0\t2\t0\n' > "$test_data_dir/trio.ped"

# --- Test Case 1: Basic run ---
log "Starting TEST 1: Combine two plain GVCFs with default options"

"$meta_executable" \
  --variant "$test_data_dir/sample1.g.vcf" \
  --variant "$test_data_dir/sample2.g.vcf" \
  "${reference_args[@]}" \
  --output "$meta_temp_dir/combined.g.vcf" \
  --output_index "$meta_temp_dir/combined.idx"

log "Validating TEST 1 outputs..."
check_file_not_empty "$meta_temp_dir/combined.g.vcf" "combined GVCF file"
check_file_contains "$meta_temp_dir/combined.g.vcf" "^##fileformat=VCF" "combined GVCF header"
check_file_contains "$meta_temp_dir/combined.g.vcf" "sample1" "combined GVCF sample1 column"
check_file_contains "$meta_temp_dir/combined.g.vcf" "sample2" "combined GVCF sample2 column"
check_index_moved "$meta_temp_dir/combined.idx" "$meta_temp_dir/combined.g.vcf.idx" "output index file (Tribble index for a .vcf)"
check_index_magic "$meta_temp_dir/combined.idx" "TIDX" "output index file"

log "✅ TEST 1 completed successfully"

# --- Test Case 2: Inputs and indexes ---
log "Starting TEST 2: Input index handling and output index options"

# Run 1: bgzipped inputs with separate indexes and a .g.vcf.gz output with
# its index in a different directory
mkdir -p "$meta_temp_dir/gz_output" "$meta_temp_dir/gz_output_index"

log "Executing $meta_name with separate input indexes and a .g.vcf.gz output..."
"$meta_executable" \
  --variant "$test_data_dir/separate/sample1.g.vcf.gz" \
  --variant "$test_data_dir/separate/sample2.g.vcf.gz" \
  --variant_index "$test_data_dir/separate_index/sample1.g.vcf.gz.tbi" \
  --variant_index "$test_data_dir/separate_index/sample2.g.vcf.gz.tbi" \
  --dbsnp "$test_data_dir/separate/dbsnp.vcf.gz" \
  --dbsnp_index "$test_data_dir/separate_index/dbsnp.vcf.gz.tbi" \
  "${reference_args[@]}" \
  --output "$meta_temp_dir/gz_output/combined.g.vcf.gz" \
  --output_index "$meta_temp_dir/gz_output_index/combined.g.vcf.gz.tbi" 2> "$meta_temp_dir/test2_separate.log"

log "Validating TEST 2 outputs (separate indexes, .g.vcf.gz output)..."
check_file_exists "$meta_temp_dir/gz_output/combined.g.vcf.gz" "combined GVCF file (bgzipped)"
zcat "$meta_temp_dir/gz_output/combined.g.vcf.gz" > "$meta_temp_dir/gz_output_index/combined.g.vcf"
check_file_contains "$meta_temp_dir/gz_output_index/combined.g.vcf" "^##fileformat=VCF" "combined GVCF (bgzipped) header"
check_file_contains "$meta_temp_dir/gz_output_index/combined.g.vcf" "sample1" "combined GVCF (bgzipped) sample1 column"
check_file_contains "$meta_temp_dir/gz_output_index/combined.g.vcf" "sample2" "combined GVCF (bgzipped) sample2 column"
check_file_not_contains "$meta_temp_dir/test2_separate.log" "no index was provided" "run log (no indexing when indexes are given)"
check_file_not_exists "$test_data_dir/separate/sample1.g.vcf.gz.tbi" "index next to the sample1 input"
check_file_not_exists "$test_data_dir/separate/sample2.g.vcf.gz.tbi" "index next to the sample2 input"
check_file_not_exists "$test_data_dir/separate/dbsnp.vcf.gz.tbi" "index next to the dbsnp input"
check_file_exists "$meta_temp_dir/gz_output_index/combined.g.vcf.gz.tbi" "output index file (tabix index for a .vcf.gz)"
check_index_magic "$meta_temp_dir/gz_output_index/combined.g.vcf.gz.tbi" "TBI" "output index file"

log "Checking that no index was left next to the output GVCF..."
found_files=$(find "$meta_temp_dir/gz_output" -type f -exec basename {} \; | sort)
if [[ "$found_files" == "combined.g.vcf.gz" ]]; then
  log "✓ Output directory contains only the output GVCF"
else
  log_error "✗ Output directory contains unexpected files: $(echo "$found_files" | tr '\n' ' ')"
  exit 1
fi

# Run 2: unindexed bgzipped inputs, no output index and the --interval_padding
# engine option
log "Executing $meta_name with unindexed .vcf.gz inputs, --create_output_variant_index false and --interval_padding..."
"$meta_executable" \
  --variant "$test_data_dir/no_index/sample1.g.vcf.gz" \
  --variant "$test_data_dir/no_index/sample2.g.vcf.gz" \
  --dbsnp "$test_data_dir/no_index/dbsnp.vcf.gz" \
  "${reference_args[@]}" \
  --output "$meta_temp_dir/combined_no_index.g.vcf" \
  --interval_padding 10 \
  --create_output_variant_index false 2> "$meta_temp_dir/test2_no_index.log"

log "Validating TEST 2 outputs (unindexed inputs, no output index)..."
check_file_contains "$meta_temp_dir/combined_no_index.g.vcf" "sample1" "combined GVCF (no indexes) sample1 column"
check_file_contains "$meta_temp_dir/combined_no_index.g.vcf" "sample2" "combined GVCF (no indexes) sample2 column"
check_file_contains "$meta_temp_dir/test2_no_index.log" "no index was provided" "run log (indexing warning)"
check_file_not_exists "$test_data_dir/no_index/sample1.g.vcf.gz.tbi" "index next to the sample1 input"
check_file_not_exists "$test_data_dir/no_index/dbsnp.vcf.gz.tbi" "index next to the dbsnp input"
check_file_not_exists "$meta_temp_dir/combined_no_index.g.vcf.idx" "combined GVCF index (should not exist with --create_output_variant_index false)"
check_file_contains "$meta_temp_dir/combined_no_index.g.vcf" "^##GATKCommandLine=.* --interval-padding 10 " "GATK command line (--interval_padding)"

log "✅ TEST 2 completed successfully"

# --- Test Case 3: Intervals ---
log "Starting TEST 3: Two --intervals files (BED and VCF), --exclude_intervals and --intervals_index"

"$meta_executable" \
  --variant "$test_data_dir/sample1.g.vcf" \
  --variant "$test_data_dir/sample2.g.vcf" \
  "${reference_args[@]}" \
  --output "$meta_temp_dir/combined_intervals.g.vcf" \
  --output_index "$meta_temp_dir/combined_intervals.g.vcf.idx" \
  --intervals "$test_data_dir/intervals/include.bed" \
  --intervals "$test_data_dir/intervals/include.vcf.gz" \
  --intervals_index "$test_data_dir/intervals_index/include.vcf.gz.tbi" \
  --exclude_intervals "$test_data_dir/intervals/exclude.bed" \
  --interval_set_rule UNION \
  --convert_to_base_pair_resolution 2> "$meta_temp_dir/test3.log"

log "Validating TEST 3 outputs..."
check_file_exists "$meta_temp_dir/combined_intervals.g.vcf" "combined GVCF file (intervals)"
check_file_not_contains "$meta_temp_dir/test3.log" "no index was provided for '$test_data_dir/intervals/include.vcf.gz'" "run log (no indexing of the interval VCF when --intervals_index is given)"
check_file_not_exists "$test_data_dir/intervals/include.vcf.gz.tbi" "index next to the interval VCF"

# CombineGVCFs only writes a reference block where it ends inside an
# interval, and warns when an interval cuts a block. Here each sample is one
# block (seq1:1-2000), so the records come only from
# --convert_to_base_pair_resolution. The first record of each interval has the
# position of the part of the block before it, not the interval start: POS 1
# (for 1-100), 150 (for 150-161) and 201 (for 201-1001).
check_file_contains "$meta_temp_dir/test3.log" "cuts in the middle of one or more gVCF blocks" "run log (GATK warning about cut blocks)"
expected_positions=$( { echo 1; seq 101 150; seq 162 201; seq 1002 1020; } | paste -sd ' ' -)
check_vcf_column "$meta_temp_dir/combined_intervals.g.vcf" 2 "$expected_positions" "combined GVCF (intervals) record positions"

log "✅ TEST 3 completed successfully"

# --- Test Case 4: Tool options ---
log "Starting TEST 4: Annotation, output-shaping, band, distance, de novo and RMS options"

# The --variant "a;b" list syntax gives the two files in one argument
"$meta_executable" \
  --variant "$test_data_dir/sample1.g.vcf;$test_data_dir/sample2.g.vcf" \
  "${reference_args[@]}" \
  --output "$meta_temp_dir/combined_options.g.vcf" \
  --output_index "$meta_temp_dir/combined_options.g.vcf.idx" \
  --annotation ChromosomeCounts \
  --annotation RMSMappingQuality \
  --annotation PossibleDeNovo \
  --annotation_group StandardAnnotation \
  --annotations_to_exclude InbreedingCoeff \
  --pedigree "$test_data_dir/trio.ped" \
  --sites_only_vcf_output \
  --variant_output_filtering ANYWHERE \
  --ignore_variants_starting_outside_interval \
  --break_bands_at_multiples_of 10 \
  --combine_variants_distance 10 \
  --max_distance 100 \
  --denovo_depth_threshold 5 \
  --denovo_parent_gq_threshold 10 \
  --allow_old_rms_mapping_quality_annotation_data

log "Validating TEST 4 outputs..."
output="$meta_temp_dir/combined_options.g.vcf"
check_file_contains "$output" "^##fileformat=VCF" "combined GVCF (options) header"
check_file_contains "$output" "^##GATKCommandLine=.* --variant [^ ]*sample1.g.vcf --variant [^ ]*sample2.g.vcf " "GATK command line (both files of --variant \"a;b\")"

# ChromosomeCounts adds AC, PossibleDeNovo with --pedigree adds the de novo INFO
# fields, --break_bands_at_multiples_of 10 cuts the single reference block of
# each sample at seq1:10, 20, ...
check_file_contains "$output" "##INFO=<ID=AC," "combined GVCF (options) AC INFO definition (ChromosomeCounts)"
check_file_contains "$output" "##INFO=<ID=hiConfDeNovo," "combined GVCF (options) hiConfDeNovo INFO definition (PossibleDeNovo)"
check_vcf_has_record "$output" $'^seq1\t10\t.*\tEND=19$' "combined GVCF (options) band seq1:10-19 (--break_bands_at_multiples_of)"

log "Checking that --sites_only_vcf_output dropped the FORMAT/genotype columns..."
chrom_line=$(grep "^#CHROM" "$output")
if [[ "$chrom_line" == *"FORMAT"* || "$chrom_line" == *"sample1"* || "$chrom_line" == *"sample2"* ]]; then
  log_error "✗ #CHROM header still contains FORMAT/sample columns despite --sites_only_vcf_output: $chrom_line"
  exit 1
else
  log "✓ #CHROM header has no FORMAT/sample columns (as expected): $chrom_line"
fi

# Pass-through options without a visible effect on this data: check the GATK
# command line in the output header
for expected in \
  "--annotation RMSMappingQuality" \
  "--annotation-group StandardAnnotation" \
  "--annotations-to-exclude InbreedingCoeff" \
  "--variant-output-filtering ANYWHERE" \
  "--ignore-variants-starting-outside-interval true" \
  "--combine-variants-distance 10" \
  "--max-distance 100" \
  "--denovo-depth-threshold 5" \
  "--denovo-parent-gq-threshold 10" \
  "--allow-old-rms-mapping-quality-annotation-data true"; do
  check_file_contains "$output" "^##GATKCommandLine=.* $expected " "GATK command line ($expected)"
done

log "✅ TEST 4 completed successfully"

# --- Test Case 5: Argument errors ---
log "Starting TEST 5: Argument errors"

log "Executing $meta_name with 2 --variant files and 1 --variant_index file..."
check_fails_with "Error: --variant_index must be given once for each file in --variant" "$meta_temp_dir/test5_mismatch.log" \
  "$meta_executable" \
  --variant "$test_data_dir/separate/sample1.g.vcf.gz" \
  --variant "$test_data_dir/separate/sample2.g.vcf.gz" \
  --variant_index "$test_data_dir/separate_index/sample1.g.vcf.gz.tbi" \
  "${reference_args[@]}" \
  --output "$meta_temp_dir/combined_mismatch.g.vcf" \
  --output_index "$meta_temp_dir/combined_mismatch.g.vcf.idx"

log "Executing $meta_name with --dbsnp_index and no --dbsnp..."
check_fails_with "Error: --dbsnp_index requires --dbsnp" "$meta_temp_dir/test5_dbsnp.log" \
  "$meta_executable" \
  --variant "$test_data_dir/sample1.g.vcf" \
  --dbsnp_index "$test_data_dir/separate_index/dbsnp.vcf.gz.tbi" \
  "${reference_args[@]}" \
  --output "$meta_temp_dir/combined_dbsnp_index.g.vcf" \
  --output_index "$meta_temp_dir/combined_dbsnp_index.g.vcf.idx"

log "Executing $meta_name without --output_index..."
check_fails_with "Error: --output_index is required unless --create_output_variant_index is false." "$meta_temp_dir/test5_output_index.log" \
  "$meta_executable" \
  --variant "$test_data_dir/sample1.g.vcf" \
  --variant "$test_data_dir/sample2.g.vcf" \
  "${reference_args[@]}" \
  --output "$meta_temp_dir/no_output_index.g.vcf"
check_file_not_exists "$meta_temp_dir/no_output_index.g.vcf" "output GVCF (no run without --output_index)"

log "✅ TEST 5 completed successfully"

print_test_summary "All tests"
