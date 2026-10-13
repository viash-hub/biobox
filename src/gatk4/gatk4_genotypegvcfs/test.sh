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

# A GVCF (plain VCF) from a direct HaplotypeCaller run
create_test_sam_reads "$test_data_dir/reads.sam" seq1 2000 rg1 sample1
sort_and_index_bam "$test_data_dir/reads.sam" "$test_data_dir/reads.sorted.bam"
check_file_exists "$test_data_dir/reads.sorted.bai" "BAM index file"

log "Running HaplotypeCaller in GVCF mode to build the input fixture..."
gatk HaplotypeCaller \
  --input "$test_data_dir/reads.sorted.bam" \
  --output "$test_data_dir/sample.g.vcf" \
  --reference "$test_data_dir/reference.fasta" \
  --emit-ref-confidence GVCF \
  --verbosity ERROR
check_file_contains "$test_data_dir/sample.g.vcf" "<NON_REF>" "input GVCF fixture NON_REF placeholder allele"

# A GenomicsDB workspace made from the GVCF
log "Creating a GenomicsDB workspace from the input GVCF..."
printf 'sample1\t%s\n' "$test_data_dir/sample.g.vcf" > "$test_data_dir/sample_map.tsv"
gatk GenomicsDBImport \
  --genomicsdb-workspace-path "$test_data_dir/genomicsdb_workspace" \
  --batch-size 50 \
  --sample-name-map "$test_data_dir/sample_map.tsv" \
  --reader-threads 1 \
  --tmp-dir "$meta_temp_dir" \
  --intervals "seq1:1-2000" \
  --reference "$test_data_dir/reference.fasta"

# A known-sites VCF for --dbsnp and --population_callset (only the bgzipped
# copies are used, so no index)
create_test_known_sites_vcf --no-index "$test_data_dir/known_sites.vcf" seq1 2000 100

# Bgzipped inputs with their indexes in a different directory
create_test_vcf_gz_with_separate_index "$test_data_dir/sample.g.vcf" \
  "$test_data_dir/separate/sample.g.vcf.gz" "$test_data_dir/separate_index/sample.g.vcf.gz.tbi"
create_test_vcf_gz_with_separate_index "$test_data_dir/known_sites.vcf" \
  "$test_data_dir/separate/dbsnp.vcf.gz" "$test_data_dir/separate_index/dbsnp.vcf.gz.tbi"
create_test_vcf_gz_with_separate_index "$test_data_dir/known_sites.vcf" \
  "$test_data_dir/separate/population.vcf.gz" "$test_data_dir/separate_index/population.vcf.gz.tbi"

# Bgzipped inputs without indexes
mkdir -p "$test_data_dir/no_index"
bgzip -c "$test_data_dir/sample.g.vcf" > "$test_data_dir/no_index/sample.g.vcf.gz"
bgzip -c "$test_data_dir/known_sites.vcf" > "$test_data_dir/no_index/dbsnp.vcf.gz"
bgzip -c "$test_data_dir/known_sites.vcf" > "$test_data_dir/no_index/population.vcf.gz"

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

# A forced output site at seq1:600 (BED is 0-based)
printf 'seq1\t599\t600\n' > "$test_data_dir/force_sites.bed"

# A trio pedigree: sample1 is the child of two parents that are not in the GVCF
printf 'fam1\tsample1\tfather\tmother\t1\t0\nfam1\tfather\t0\t0\t1\t0\nfam1\tmother\t0\t0\t2\t0\n' > "$test_data_dir/trio.ped"

# --- Test Case 1: Basic run ---
log "Starting TEST 1: Default parameters"

"$meta_executable" \
  --variant "$test_data_dir/sample.g.vcf" \
  "${reference_args[@]}" \
  --output "$meta_temp_dir/output.vcf" \
  --output_index "$meta_temp_dir/output.idx"

log "Validating TEST 1 outputs..."
check_file_not_empty "$meta_temp_dir/output.vcf" "output VCF file"
check_file_contains "$meta_temp_dir/output.vcf" "^##fileformat=VCF" "output VCF file header"
check_file_contains "$meta_temp_dir/output.vcf" "sample1" "output VCF sample column"
check_file_not_contains "$meta_temp_dir/output.vcf" "<NON_REF>" "output VCF should not contain GVCF placeholder allele"
check_index_moved "$meta_temp_dir/output.idx" "$meta_temp_dir/output.vcf.idx" "output index file (Tribble index for a .vcf)"
check_index_magic "$meta_temp_dir/output.idx" "TIDX" "output index file"

log "✅ TEST 1 completed successfully"

# --- Test Case 2: Inputs and indexes ---
log "Starting TEST 2: Input index handling, output index options and GenomicsDB input"

# Run 1: bgzipped inputs with separate indexes and a .vcf.gz output with its
# index in a different directory
mkdir -p "$meta_temp_dir/gz_output" "$meta_temp_dir/gz_output_index"

log "Executing $meta_name with separate input indexes and a .vcf.gz output..."
"$meta_executable" \
  --variant "$test_data_dir/separate/sample.g.vcf.gz" \
  --variant_index "$test_data_dir/separate_index/sample.g.vcf.gz.tbi" \
  --dbsnp "$test_data_dir/separate/dbsnp.vcf.gz" \
  --dbsnp_index "$test_data_dir/separate_index/dbsnp.vcf.gz.tbi" \
  --population_callset "$test_data_dir/separate/population.vcf.gz" \
  --population_callset_index "$test_data_dir/separate_index/population.vcf.gz.tbi" \
  "${reference_args[@]}" \
  --output "$meta_temp_dir/gz_output/output.vcf.gz" \
  --output_index "$meta_temp_dir/gz_output_index/output.vcf.gz.tbi" 2> "$meta_temp_dir/test2_separate.log"

log "Validating TEST 2 outputs (separate indexes, .vcf.gz output)..."
check_file_exists "$meta_temp_dir/gz_output/output.vcf.gz" "output VCF file (bgzipped)"
zcat "$meta_temp_dir/gz_output/output.vcf.gz" > "$meta_temp_dir/gz_output_index/output.vcf"
check_file_contains "$meta_temp_dir/gz_output_index/output.vcf" "^##fileformat=VCF" "output VCF (bgzipped) header"
check_file_contains "$meta_temp_dir/gz_output_index/output.vcf" "sample1" "output VCF (bgzipped) sample column"
check_file_not_contains "$meta_temp_dir/gz_output_index/output.vcf" "<NON_REF>" "output VCF should not contain GVCF placeholder allele (separate indexes)"
check_file_not_contains "$meta_temp_dir/test2_separate.log" "no index was provided" "run log (no indexing when indexes are given)"
check_file_not_exists "$test_data_dir/separate/sample.g.vcf.gz.tbi" "index next to the --variant input"
check_file_not_exists "$test_data_dir/separate/dbsnp.vcf.gz.tbi" "index next to the --dbsnp input"
check_file_not_exists "$test_data_dir/separate/population.vcf.gz.tbi" "index next to the --population_callset input"
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

# Run 2: unindexed bgzipped inputs, no output index and the --interval_padding
# engine option
log "Executing $meta_name with unindexed .vcf.gz inputs, --create_output_variant_index false and --interval_padding..."
"$meta_executable" \
  --variant "$test_data_dir/no_index/sample.g.vcf.gz" \
  --dbsnp "$test_data_dir/no_index/dbsnp.vcf.gz" \
  --population_callset "$test_data_dir/no_index/population.vcf.gz" \
  "${reference_args[@]}" \
  --output "$meta_temp_dir/output_no_index.vcf" \
  --interval_padding 10 \
  --create_output_variant_index false 2> "$meta_temp_dir/test2_no_index.log"

log "Validating TEST 2 outputs (unindexed inputs, no output index)..."
check_file_contains "$meta_temp_dir/output_no_index.vcf" "sample1" "output VCF sample column (no indexes)"
check_file_contains "$meta_temp_dir/test2_no_index.log" "no index was provided" "run log (indexing warning)"
check_file_not_exists "$test_data_dir/no_index/sample.g.vcf.gz.tbi" "index next to the --variant input"
check_file_not_exists "$test_data_dir/no_index/dbsnp.vcf.gz.tbi" "index next to the --dbsnp input"
check_file_not_exists "$test_data_dir/no_index/population.vcf.gz.tbi" "index next to the --population_callset input"
check_file_not_exists "$meta_temp_dir/output_no_index.vcf.idx" "output VCF index (should not exist with --create_output_variant_index false)"
check_file_contains "$meta_temp_dir/output_no_index.vcf" "^##GATKCommandLine=.* --interval-padding 10 " "GATK command line (--interval_padding)"

# Run 3: a GenomicsDB workspace as --variant
log "Executing $meta_name with GenomicsDB workspace input..."
"$meta_executable" \
  --variant "$test_data_dir/genomicsdb_workspace" \
  "${reference_args[@]}" \
  --output "$meta_temp_dir/output_genomicsdb.vcf" \
  --output_index "$meta_temp_dir/output_genomicsdb.vcf.idx"

log "Validating TEST 2 outputs (GenomicsDB input)..."
check_file_not_empty "$meta_temp_dir/output_genomicsdb.vcf" "output VCF file (GenomicsDB input)"
check_file_contains "$meta_temp_dir/output_genomicsdb.vcf" "^##fileformat=VCF" "output VCF file header (GenomicsDB input)"
check_file_contains "$meta_temp_dir/output_genomicsdb.vcf" "sample1" "output VCF sample column (GenomicsDB input)"
check_file_not_contains "$meta_temp_dir/output_genomicsdb.vcf" "<NON_REF>" "output VCF should not contain GVCF placeholder allele (GenomicsDB input)"

log "✅ TEST 2 completed successfully"

# --- Test Case 3: Intervals ---
log "Starting TEST 3: Interval files and forced output sites"

# Run 1: two --intervals files and --exclude_intervals. With
# --include_non_variant_sites there is one record for each position.
log "Executing $meta_name with two --intervals files, --exclude_intervals and --include_non_variant_sites..."
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
  --include_non_variant_sites 2> "$meta_temp_dir/test3_intervals.log"

log "Validating TEST 3 outputs (interval files)..."
check_file_exists "$meta_temp_dir/output_intervals.vcf" "output VCF file (intervals)"
check_file_not_contains "$meta_temp_dir/test3_intervals.log" "no index was provided for '$test_data_dir/intervals/include.vcf.gz'" "run log (no indexing of the interval VCF when --intervals_index is given)"
check_file_not_exists "$test_data_dir/intervals/include.vcf.gz.tbi" "index next to the interval VCF"
check_file_not_contains "$meta_temp_dir/output_intervals.vcf" "<NON_REF>" "output VCF should not contain GVCF placeholder allele (intervals)"
expected_positions=$({ seq 100 149; seq 161 200; seq 1001 1020; } | paste -sd ' ' -)
check_vcf_column "$meta_temp_dir/output_intervals.vcf" 2 "$expected_positions" "output VCF (intervals) record positions"

# Run 2: force seq1:600 with a BED file and seq1:500 with a string
log "Executing $meta_name with --force_output_intervals and --force_output_interval_strings..."
"$meta_executable" \
  --variant "$test_data_dir/sample.g.vcf" \
  "${reference_args[@]}" \
  --output "$meta_temp_dir/output_forced.vcf" \
  --create_output_variant_index false \
  --force_output_intervals "$test_data_dir/force_sites.bed" \
  --force_output_interval_strings "seq1:500"

log "Validating TEST 3 outputs (forced sites)..."
check_file_exists "$meta_temp_dir/output_forced.vcf" "output VCF file (forced sites)"
for pos in 500 600; do
  # The default output must not have the site, else the check below proves nothing
  check_vcf_lacks_record "$meta_temp_dir/output.vcf" $'^seq1\t'"$pos"$'\t' "default output VCF (TEST 1)"
  check_vcf_has_record "$meta_temp_dir/output_forced.vcf" $'^seq1\t'"$pos"$'\t' "output VCF (forced sites)"
done

log "✅ TEST 3 completed successfully"

# --- Test Case 4: Tool options ---
log "Starting TEST 4: Annotation, output-shaping, calling, de novo and RMS options"

# --stand_call_conf is the alternative name of
# --standard_min_confidence_threshold_for_calling
"$meta_executable" \
  --variant "$test_data_dir/sample.g.vcf" \
  "${reference_args[@]}" \
  --output "$meta_temp_dir/output_options.vcf" \
  --output_index "$meta_temp_dir/output_options.vcf.idx" \
  --annotation ChromosomeCounts \
  --annotation RMSMappingQuality \
  --annotation PossibleDeNovo \
  --annotations_to_exclude InbreedingCoeff \
  --pedigree "$test_data_dir/trio.ped" \
  --call_genotypes \
  --sites_only_vcf_output \
  --variant_output_filtering ANYWHERE \
  --stand_call_conf 10.0 \
  --sample_ploidy 2 \
  --denovo_depth_threshold 5 \
  --denovo_parent_gq_threshold 10 \
  --allow_old_rms_mapping_quality_annotation_data

log "Validating TEST 4 outputs..."
output="$meta_temp_dir/output_options.vcf"
check_file_contains "$output" "^##fileformat=VCF" "output VCF (options) header"
check_file_not_contains "$output" "<NON_REF>" "output VCF (options) should not contain GVCF placeholder allele"

# Behaviour: ChromosomeCounts adds AC, PossibleDeNovo with --pedigree adds the
# de novo INFO fields
check_file_contains "$output" "##INFO=<ID=AC," "output VCF (options) AC INFO definition (ChromosomeCounts)"
check_file_contains "$output" "##INFO=<ID=hiConfDeNovo," "output VCF (options) hiConfDeNovo INFO definition (PossibleDeNovo)"

log "Checking that --sites_only_vcf_output dropped the FORMAT/genotype columns..."
chrom_line=$(grep "^#CHROM" "$output")
if [[ "$chrom_line" == *"FORMAT"* || "$chrom_line" == *"sample1"* ]]; then
  log_error "✗ #CHROM header still contains FORMAT/sample columns despite --sites_only_vcf_output: $chrom_line"
  exit 1
else
  log "✓ #CHROM header has no FORMAT/sample columns (as expected): $chrom_line"
fi

for expected in \
  "--annotation RMSMappingQuality" \
  "--annotations-to-exclude InbreedingCoeff" \
  "--call-genotypes true" \
  "--variant-output-filtering ANYWHERE" \
  "--standard-min-confidence-threshold-for-calling 10.0" \
  "--sample-ploidy 2" \
  "--denovo-depth-threshold 5" \
  "--denovo-parent-gq-threshold 10" \
  "--allow-old-rms-mapping-quality-annotation-data true"; do
  check_file_contains "$output" "^##GATKCommandLine=.* $expected " "GATK command line ($expected)"
done

log "✅ TEST 4 completed successfully"

# --- Test Case 5: Argument errors ---
log "Starting TEST 5: Argument errors"

fail_args=(
  "${reference_args[@]}"
  --output "$meta_temp_dir/output_fail.vcf"
  --output_index "$meta_temp_dir/output_fail.vcf.idx"
)

log "Executing $meta_name with --variant_index and a GenomicsDB workspace --variant..."
check_fails_with "Error: --variant_index cannot be used when --variant is a GenomicsDB workspace directory" \
  "$meta_temp_dir/test5_workspace.log" \
  "$meta_executable" "${fail_args[@]}" \
  --variant "$test_data_dir/genomicsdb_workspace" \
  --variant_index "$test_data_dir/separate_index/sample.g.vcf.gz.tbi"

log "Executing $meta_name with --dbsnp_index and no --dbsnp..."
check_fails_with "Error: --dbsnp_index requires --dbsnp" \
  "$meta_temp_dir/test5_dbsnp.log" \
  "$meta_executable" "${fail_args[@]}" \
  --variant "$test_data_dir/sample.g.vcf" \
  --dbsnp_index "$test_data_dir/separate_index/dbsnp.vcf.gz.tbi"

log "Executing $meta_name with --population_callset_index and no --population_callset..."
check_fails_with "Error: --population_callset_index requires --population_callset" \
  "$meta_temp_dir/test5_population.log" \
  "$meta_executable" "${fail_args[@]}" \
  --variant "$test_data_dir/sample.g.vcf" \
  --population_callset_index "$test_data_dir/separate_index/population.vcf.gz.tbi"

log "Executing $meta_name without --output_index..."
check_fails_with "Error: --output_index is required unless --create_output_variant_index is false." "$meta_temp_dir/test5_output_index.log" \
  "$meta_executable" \
  --variant "$test_data_dir/sample.g.vcf" \
  "${reference_args[@]}" \
  --output "$meta_temp_dir/no_output_index.vcf"
check_file_not_exists "$meta_temp_dir/no_output_index.vcf" "output VCF (no run without --output_index)"

log "✅ TEST 5 completed successfully"

print_test_summary "All tests"
