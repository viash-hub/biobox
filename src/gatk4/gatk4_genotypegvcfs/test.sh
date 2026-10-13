#!/bin/bash

## VIASH START
## VIASH END

# Source the centralized test helpers
source "$meta_resources_dir/test_helpers.sh"
source "$meta_resources_dir/gatk4/test_helpers.sh"

gatk4_test_init

# --- Build shared reference fixtures ---
create_test_reference "$test_data_dir/reference.fasta" 1 2000
check_file_exists "$test_data_dir/reference.fasta" "test reference genome"
check_file_exists "$test_data_dir/reference.dict" "reference sequence dictionary"
check_file_exists "$test_data_dir/reference.fasta.fai" "reference FASTA index"

reference_args=(
  --reference "$test_data_dir/reference.fasta"
  --reference_fai "$test_data_dir/reference.fasta.fai"
  --reference_dict "$test_data_dir/reference.dict"
)

# --- Build a GVCF fixture by running HaplotypeCaller directly ---
create_test_sam_reads "$test_data_dir/reads.sam" seq1 2000 rg1 sample1
check_file_exists "$test_data_dir/reads.sam" "hand-crafted SAM reads"

sort_and_index_bam "$test_data_dir/reads.sam" "$test_data_dir/reads.sorted.bam"
check_file_exists "$test_data_dir/reads.sorted.bam" "sorted BAM file"
check_file_exists "$test_data_dir/reads.sorted.bai" "BAM index file"

log "Running HaplotypeCaller in GVCF mode to build the input fixture..."
gatk HaplotypeCaller \
  --input "$test_data_dir/reads.sorted.bam" \
  --output "$test_data_dir/sample.g.vcf" \
  --reference "$test_data_dir/reference.fasta" \
  --emit-ref-confidence GVCF \
  --verbosity ERROR
check_file_exists "$test_data_dir/sample.g.vcf" "input GVCF fixture"
check_file_contains "$test_data_dir/sample.g.vcf" "<NON_REF>" "input GVCF fixture NON_REF placeholder allele"

# --- Test Case 1: Default parameters ---
log "Starting TEST 1: Default parameters"

log "Executing $meta_name with basic parameters..."
"$meta_executable" \
  --variant "$test_data_dir/sample.g.vcf" \
  "${reference_args[@]}" \
  --output "$meta_temp_dir/output.vcf" \
  --output_index "$meta_temp_dir/output.idx"

log "Validating TEST 1 outputs..."
check_file_exists "$meta_temp_dir/output.vcf" "output VCF file"
check_file_not_empty "$meta_temp_dir/output.vcf" "output VCF file"
check_file_contains "$meta_temp_dir/output.vcf" "^##fileformat=VCF" "output VCF file header"
check_file_contains "$meta_temp_dir/output.vcf" "sample1" "output VCF sample column"
check_file_not_contains "$meta_temp_dir/output.vcf" "<NON_REF>" "output VCF should not contain GVCF placeholder allele"
check_index_moved "$meta_temp_dir/output.idx" "$meta_temp_dir/output.vcf.idx" "output index file (Tribble index for a .vcf)"
check_index_magic "$meta_temp_dir/output.idx" "TIDX" "output index file"

log "✅ TEST 1 completed successfully"

# --- Test Case 2: Options ---
log "Starting TEST 2: Additional options"

log "Executing $meta_name with additional options (--stand_call_conf is an alternative name)..."
"$meta_executable" \
  --variant "$test_data_dir/sample.g.vcf" \
  "${reference_args[@]}" \
  --output "$meta_temp_dir/output_options.vcf" \
  --output_index "$meta_temp_dir/output_options.vcf.idx" \
  --include_non_variant_sites \
  --stand_call_conf 10.0 \
  --sample_ploidy 2

log "Validating TEST 2 outputs..."
check_file_exists "$meta_temp_dir/output_options.vcf" "output VCF file with options"
check_file_not_empty "$meta_temp_dir/output_options.vcf" "output VCF file with options"
check_file_contains "$meta_temp_dir/output_options.vcf" "^##fileformat=VCF" "output VCF file header"
check_file_not_contains "$meta_temp_dir/output_options.vcf" "<NON_REF>" "output VCF should not contain GVCF placeholder allele"

log "✅ TEST 2 completed successfully"

# --- Test Case 3: Annotation and output-shaping options ---
log "Starting TEST 3: Additional annotation and output-shaping options"

log "Executing $meta_name with --annotation, --annotations_to_exclude, --call_genotypes, --sites_only_vcf_output, and --variant_output_filtering..."
"$meta_executable" \
  --variant "$test_data_dir/sample.g.vcf" \
  "${reference_args[@]}" \
  --output "$meta_temp_dir/output_annotations.vcf" \
  --output_index "$meta_temp_dir/output_annotations.vcf.idx" \
  --annotation ChromosomeCounts \
  --annotations_to_exclude InbreedingCoeff \
  --call_genotypes \
  --sites_only_vcf_output \
  --variant_output_filtering ANYWHERE

log "Validating TEST 3 outputs..."
check_file_exists "$meta_temp_dir/output_annotations.vcf" "output VCF file with annotation options"
check_file_not_empty "$meta_temp_dir/output_annotations.vcf" "output VCF file with annotation options"
check_file_contains "$meta_temp_dir/output_annotations.vcf" "^##fileformat=VCF" "output VCF (annotations) header"
check_file_contains "$meta_temp_dir/output_annotations.vcf" "##INFO=<ID=AC," "output VCF (annotations) AC INFO definition"

log "Checking that --sites_only_vcf_output dropped the FORMAT/genotype columns..."
chrom_line=$(grep "^#CHROM" "$meta_temp_dir/output_annotations.vcf")
if [[ "$chrom_line" == *"FORMAT"* || "$chrom_line" == *"sample1"* ]]; then
  log_error "✗ #CHROM header still contains FORMAT/sample columns despite --sites_only_vcf_output: $chrom_line"
  exit 1
else
  log "✓ #CHROM header has no FORMAT/sample columns (as expected): $chrom_line"
fi

log "✅ TEST 3 completed successfully"

# --- Test Case 4: Shared GATK engine options ---
log "Starting TEST 4: Shared GATK engine options"

log "Executing $meta_name with --interval_padding and --create_output_variant_index false..."
"$meta_executable" \
  --variant "$test_data_dir/sample.g.vcf" \
  "${reference_args[@]}" \
  --output "$meta_temp_dir/output_engine.vcf" \
  --interval_padding 10 \
  --create_output_variant_index false

log "Validating TEST 4 outputs..."
check_file_exists "$meta_temp_dir/output_engine.vcf" "output VCF file (engine options)"
check_file_not_empty "$meta_temp_dir/output_engine.vcf" "output VCF file (engine options)"
check_file_not_exists "$meta_temp_dir/output_engine.vcf.idx" "output VCF index (should not exist with --create_output_variant_index false)"

log "✅ TEST 4 completed successfully"

# --- Test Case 5: Test GenomicsDB workspace input ---
log "Starting TEST 5: GenomicsDB workspace input"

log "Creating a GenomicsDB workspace from the input GVCF..."
printf 'sample1\t%s\n' "$test_data_dir/sample.g.vcf" > "$test_data_dir/sample_map.tsv"
gatk GenomicsDBImport \
  --genomicsdb-workspace-path "$meta_temp_dir/genomicsdb_workspace" \
  --batch-size 50 \
  --sample-name-map "$test_data_dir/sample_map.tsv" \
  --reader-threads 1 \
  --tmp-dir "$meta_temp_dir" \
  --intervals "seq1:1-2000" \
  --reference "$test_data_dir/reference.fasta"

log "Executing $meta_name with GenomicsDB workspace input..."
"$meta_executable" \
  --variant "$meta_temp_dir/genomicsdb_workspace" \
  "${reference_args[@]}" \
  --output "$meta_temp_dir/output_genomicsdb.vcf" \
  --output_index "$meta_temp_dir/output_genomicsdb.vcf.idx"

log "Validating TEST 5 outputs..."
check_file_exists "$meta_temp_dir/output_genomicsdb.vcf" "output VCF file (GenomicsDB input)"
check_file_not_empty "$meta_temp_dir/output_genomicsdb.vcf" "output VCF file (GenomicsDB input)"
check_file_contains "$meta_temp_dir/output_genomicsdb.vcf" "^##fileformat=VCF" "output VCF file header (GenomicsDB input)"
check_file_contains "$meta_temp_dir/output_genomicsdb.vcf" "sample1" "output VCF sample column (GenomicsDB input)"
check_file_not_contains "$meta_temp_dir/output_genomicsdb.vcf" "<NON_REF>" "output VCF should not contain GVCF placeholder allele (GenomicsDB input)"

log "✅ TEST 5 completed successfully"

# --- Test Case 6: Index files in a different directory ---
log "Starting TEST 6: --variant_index, --dbsnp_index and --population_callset_index in a different directory"

mkdir -p "$test_data_dir/separate"
create_test_known_sites_vcf "$test_data_dir/known_sites.vcf" seq1 2000 100
create_test_vcf_gz_with_separate_index "$test_data_dir/sample.g.vcf" \
  "$test_data_dir/separate/sample.g.vcf.gz" "$test_data_dir/separate_index/sample.g.vcf.gz.tbi"
create_test_vcf_gz_with_separate_index "$test_data_dir/known_sites.vcf" \
  "$test_data_dir/separate/dbsnp.vcf.gz" "$test_data_dir/separate_index/dbsnp.vcf.gz.tbi"
create_test_vcf_gz_with_separate_index "$test_data_dir/known_sites.vcf" \
  "$test_data_dir/separate/population.vcf.gz" "$test_data_dir/separate_index/population.vcf.gz.tbi"

log "Executing $meta_name with .vcf.gz inputs and separate .tbi index files..."
"$meta_executable" \
  --variant "$test_data_dir/separate/sample.g.vcf.gz" \
  --variant_index "$test_data_dir/separate_index/sample.g.vcf.gz.tbi" \
  --dbsnp "$test_data_dir/separate/dbsnp.vcf.gz" \
  --dbsnp_index "$test_data_dir/separate_index/dbsnp.vcf.gz.tbi" \
  --population_callset "$test_data_dir/separate/population.vcf.gz" \
  --population_callset_index "$test_data_dir/separate_index/population.vcf.gz.tbi" \
  "${reference_args[@]}" \
  --output "$meta_temp_dir/output_separate.vcf" \
  --output_index "$meta_temp_dir/output_separate.vcf.idx" 2> "$meta_temp_dir/test6.log"

log "Validating TEST 6 outputs..."
check_file_exists "$meta_temp_dir/output_separate.vcf" "output VCF file (separate indexes)"
check_file_contains "$meta_temp_dir/output_separate.vcf" "^##fileformat=VCF" "output VCF file header (separate indexes)"
check_file_contains "$meta_temp_dir/output_separate.vcf" "sample1" "output VCF sample column (separate indexes)"
check_file_not_contains "$meta_temp_dir/output_separate.vcf" "<NON_REF>" "output VCF should not contain GVCF placeholder allele (separate indexes)"
check_file_not_contains "$meta_temp_dir/test6.log" "no index was provided" "run log (no indexing when indexes are given)"
check_file_not_exists "$test_data_dir/separate/sample.g.vcf.gz.tbi" "index next to the --variant input"
check_file_not_exists "$test_data_dir/separate/dbsnp.vcf.gz.tbi" "index next to the --dbsnp input"
check_file_not_exists "$test_data_dir/separate/population.vcf.gz.tbi" "index next to the --population_callset input"

log "✅ TEST 6 completed successfully"

# --- Test Case 7: Unindexed .vcf.gz inputs ---
log "Starting TEST 7: .vcf.gz inputs without index files"

mkdir -p "$test_data_dir/no_index"
bgzip -c "$test_data_dir/sample.g.vcf" > "$test_data_dir/no_index/sample.g.vcf.gz"
bgzip -c "$test_data_dir/known_sites.vcf" > "$test_data_dir/no_index/dbsnp.vcf.gz"
bgzip -c "$test_data_dir/known_sites.vcf" > "$test_data_dir/no_index/population.vcf.gz"

log "Executing $meta_name with unindexed .vcf.gz inputs..."
"$meta_executable" \
  --variant "$test_data_dir/no_index/sample.g.vcf.gz" \
  --dbsnp "$test_data_dir/no_index/dbsnp.vcf.gz" \
  --population_callset "$test_data_dir/no_index/population.vcf.gz" \
  "${reference_args[@]}" \
  --output "$meta_temp_dir/output_no_index.vcf" \
  --output_index "$meta_temp_dir/output_no_index.vcf.idx" 2> "$meta_temp_dir/test7.log"

log "Validating TEST 7 outputs..."
check_file_exists "$meta_temp_dir/output_no_index.vcf" "output VCF file (no indexes)"
check_file_contains "$meta_temp_dir/output_no_index.vcf" "sample1" "output VCF sample column (no indexes)"
check_file_contains "$meta_temp_dir/test7.log" "no index was provided" "run log (indexing warning)"
check_file_not_exists "$test_data_dir/no_index/sample.g.vcf.gz.tbi" "index next to the --variant input"
check_file_not_exists "$test_data_dir/no_index/dbsnp.vcf.gz.tbi" "index next to the --dbsnp input"
check_file_not_exists "$test_data_dir/no_index/population.vcf.gz.tbi" "index next to the --population_callset input"

log "✅ TEST 7 completed successfully"

# --- Test Case 8: Index argument errors ---
log "Starting TEST 8: Index argument errors"

fail_args=(
  "${reference_args[@]}"
  --output "$meta_temp_dir/output_fail.vcf"
  --output_index "$meta_temp_dir/output_fail.vcf.idx"
)

log "Executing $meta_name with --variant_index and a GenomicsDB workspace --variant..."
check_fails_with "Error: --variant_index cannot be used when --variant is a GenomicsDB workspace directory" \
  "$meta_temp_dir/test8_workspace.log" \
  "$meta_executable" "${fail_args[@]}" \
  --variant "$meta_temp_dir/genomicsdb_workspace" \
  --variant_index "$test_data_dir/separate_index/sample.g.vcf.gz.tbi"

log "Executing $meta_name with --dbsnp_index and no --dbsnp..."
check_fails_with "Error: --dbsnp_index requires --dbsnp" \
  "$meta_temp_dir/test8_dbsnp.log" \
  "$meta_executable" "${fail_args[@]}" \
  --variant "$test_data_dir/sample.g.vcf" \
  --dbsnp_index "$test_data_dir/separate_index/dbsnp.vcf.gz.tbi"

log "Executing $meta_name with --population_callset_index and no --population_callset..."
check_fails_with "Error: --population_callset_index requires --population_callset" \
  "$meta_temp_dir/test8_population.log" \
  "$meta_executable" "${fail_args[@]}" \
  --variant "$test_data_dir/sample.g.vcf" \
  --population_callset_index "$test_data_dir/separate_index/population.vcf.gz.tbi"

log "✅ TEST 8 completed successfully"

# --- Test Case 9: Bgzipped output with an index in a different directory ---
log "Starting TEST 9: .vcf.gz output with --output_index in a different directory"

mkdir -p "$meta_temp_dir/gz_output" "$meta_temp_dir/gz_output_index"

log "Executing $meta_name with a .vcf.gz --output and --output_index in a different directory..."
"$meta_executable" \
  --variant "$test_data_dir/sample.g.vcf" \
  "${reference_args[@]}" \
  --output "$meta_temp_dir/gz_output/output.vcf.gz" \
  --output_index "$meta_temp_dir/gz_output_index/output.vcf.gz.tbi"

log "Validating TEST 9 outputs..."
check_file_exists "$meta_temp_dir/gz_output/output.vcf.gz" "output VCF file (bgzipped)"
zcat "$meta_temp_dir/gz_output/output.vcf.gz" > "$meta_temp_dir/gz_output_index/output.vcf"
check_file_contains "$meta_temp_dir/gz_output_index/output.vcf" "sample1" "output VCF (bgzipped) sample column"
check_file_exists "$meta_temp_dir/gz_output_index/output.vcf.gz.tbi" "output index file (tabix index for a .vcf.gz)"
check_index_magic "$meta_temp_dir/gz_output_index/output.vcf.gz.tbi" "TBI" "output index file"

log "Checking that no index was left next to the output VCF..."
found_files=$(find "$meta_temp_dir/gz_output" -type f -exec basename {} \; | sort)
if [[ "$found_files" == "output.vcf.gz" ]]; then
  log "✓ Output directory contains only the output VCF"
else
  log_error "✗ Output directory contains unexpected files: $(echo "$found_files" | tr '\n' ' ')"
  exit 1
fi

log "✅ TEST 9 completed successfully"

# --- Test Case 10: Missing --output_index ---
log "Starting TEST 10: A missing --output_index fails"

check_fails_with "Error: --output_index is required unless --create_output_variant_index is false." "$meta_temp_dir/test10.log" \
  "$meta_executable" \
  --variant "$test_data_dir/sample.g.vcf" \
  "${reference_args[@]}" \
  --output "$meta_temp_dir/no_output_index.vcf"
check_file_not_exists "$meta_temp_dir/no_output_index.vcf" "output VCF (no run without --output_index)"

log "✅ TEST 10 completed successfully"

# --- Test Case 11: Multiple intervals, exclusion and interval VCF index ---
log "Starting TEST 11: Two --intervals files (BED and VCF), --exclude_intervals and --intervals_index"

mkdir -p "$test_data_dir/intervals"
# Include seq1:100-200 (BED is 0-based) and seq1:1001-1020 (the 20bp REF
# allele of the VCF record), exclude seq1:150-160
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

log "Executing $meta_name with two --intervals files and --exclude_intervals..."
"$meta_executable" \
  --variant "$test_data_dir/sample.g.vcf" \
  "${reference_args[@]}" \
  --output "$meta_temp_dir/output_intervals.vcf" \
  --output_index "$meta_temp_dir/output_intervals.vcf.idx" \
  --intervals "$test_data_dir/intervals/include.bed" \
  --intervals "$test_data_dir/intervals/include.vcf.gz" \
  --intervals_index "$test_data_dir/intervals_index/include.vcf.gz.tbi" \
  --exclude_intervals "$test_data_dir/intervals/exclude.bed" \
  --interval_set_rule UNION \
  --include_non_variant_sites 2> "$meta_temp_dir/test11.log"

log "Validating TEST 11 outputs..."
check_file_exists "$meta_temp_dir/output_intervals.vcf" "output VCF file (intervals)"
check_file_not_contains "$meta_temp_dir/test11.log" "no index was provided for '$test_data_dir/intervals/include.vcf.gz'" "run log (no indexing of the interval VCF when --intervals_index is given)"
check_file_not_exists "$test_data_dir/intervals/include.vcf.gz.tbi" "index next to the interval VCF"

# With --include_non_variant_sites there is one record for each position
expected_positions=$({ seq 100 149; seq 161 200; seq 1001 1020; } | paste -sd ' ' -)
check_vcf_column "$meta_temp_dir/output_intervals.vcf" 2 "$expected_positions" "output VCF (intervals) record positions"

log "✅ TEST 11 completed successfully"

# --- Test Case 12: Forced output sites from a file and from strings ---
log "Starting TEST 12: --force_output_intervals and --force_output_interval_strings"

# Force seq1:600 with a BED file (0-based) and seq1:500 with a string
printf 'seq1\t599\t600\n' > "$test_data_dir/force_sites.bed"

log "Executing $meta_name with --force_output_intervals and --force_output_interval_strings..."
"$meta_executable" \
  --variant "$test_data_dir/sample.g.vcf" \
  "${reference_args[@]}" \
  --output "$meta_temp_dir/output_forced.vcf" \
  --create_output_variant_index false \
  --force_output_intervals "$test_data_dir/force_sites.bed" \
  --force_output_interval_strings "seq1:500"

log "Validating TEST 12 outputs..."
check_file_exists "$meta_temp_dir/output_forced.vcf" "output VCF file (forced sites)"
for pos in 500 600; do
  # The default output must not have the site, else the check below proves nothing
  check_vcf_lacks_record "$meta_temp_dir/output.vcf" $'^seq1\t'"$pos"$'\t' "default output VCF (TEST 1)"
  check_vcf_has_record "$meta_temp_dir/output_forced.vcf" $'^seq1\t'"$pos"$'\t' "output VCF (forced sites)"
done

log "✅ TEST 12 completed successfully"

# --- Test Case 13: De novo and RMS mapping quality options ---
log "Starting TEST 13: --denovo_depth_threshold, --denovo_parent_gq_threshold and --allow_old_rms_mapping_quality_annotation_data"

# A trio pedigree: sample1 is the child of two parents that are not in the GVCF
printf 'fam1\tsample1\tfather\tmother\t1\t0\nfam1\tfather\t0\t0\t1\t0\nfam1\tmother\t0\t0\t2\t0\n' > "$test_data_dir/trio.ped"

log "Executing $meta_name with the de novo and RMS mapping quality options..."
"$meta_executable" \
  --variant "$test_data_dir/sample.g.vcf" \
  "${reference_args[@]}" \
  --output "$meta_temp_dir/output_new_options.vcf" \
  --create_output_variant_index false \
  --annotation RMSMappingQuality \
  --annotation PossibleDeNovo \
  --pedigree "$test_data_dir/trio.ped" \
  --denovo_depth_threshold 5 \
  --denovo_parent_gq_threshold 10 \
  --allow_old_rms_mapping_quality_annotation_data \
  --standard_min_confidence_threshold_for_calling 10.0

log "Validating TEST 13 outputs..."
check_file_exists "$meta_temp_dir/output_new_options.vcf" "output VCF file (new options)"
check_file_contains "$meta_temp_dir/output_new_options.vcf" "^##fileformat=VCF" "output VCF (new options) header"
check_file_contains "$meta_temp_dir/output_new_options.vcf" "sample1" "output VCF (new options) sample1 column"

log "✅ TEST 13 completed successfully"

print_test_summary "All tests"
