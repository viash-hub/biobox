#!/bin/bash

## VIASH START
## VIASH END

# Source the centralized test helpers and the gatk4-specific helpers
source "$meta_resources_dir/test_helpers.sh"
source "$meta_resources_dir/gatk4/test_helpers.sh"

gatk4_test_init

#############################################
# Shared fixtures: reference genome + BAM + BQSR recal table
#############################################

create_test_reference "$test_data_dir/reference.fasta" 1 2000
check_file_exists "$test_data_dir/reference.fasta" "test reference genome"
check_file_exists "$test_data_dir/reference.dict" "reference sequence dictionary"
check_file_exists "$test_data_dir/reference.fasta.fai" "reference FASTA index"

sam_file="$test_data_dir/reads.sam"
create_test_sam_reads "$sam_file" seq1 2000 rg1 sample1 20 100
check_file_exists "$sam_file" "test SAM reads"

sort_and_index_bam "$sam_file" "$test_data_dir/sample.sorted.bam"
check_file_exists "$test_data_dir/sample.sorted.bam" "sorted test BAM"
check_file_exists "$test_data_dir/sample.sorted.bai" "sorted test BAM index"

known_sites_vcf="$test_data_dir/known_sites.vcf"
create_test_known_sites_vcf "$known_sites_vcf" seq1 2000 500
check_file_exists "$known_sites_vcf" "known-sites VCF"

log "Running BaseRecalibrator to produce a recal table..."
recal_table="$test_data_dir/recal.table"
gatk BaseRecalibrator \
  --input "$test_data_dir/sample.sorted.bam" \
  --reference "$test_data_dir/reference.fasta" \
  --known-sites "$known_sites_vcf" \
  --output "$recal_table" \
  --verbosity ERROR
check_file_exists "$recal_table" "BQSR recalibration table"
check_file_not_empty "$recal_table" "BQSR recalibration table"

reference_args=(
  --reference "$test_data_dir/reference.fasta"
  --reference_fai "$test_data_dir/reference.fasta.fai"
  --reference_dict "$test_data_dir/reference.dict"
)

# --- Test Case 1: Basic ApplyBQSR ---
log "Starting TEST 1: Basic ApplyBQSR"

mkdir -p "$meta_temp_dir/test1" "$meta_temp_dir/test1_index"

log "Executing $meta_name with basic parameters..."
"$meta_executable" \
  --input "$test_data_dir/sample.sorted.bam" \
  --read_index "$test_data_dir/sample.sorted.bai" \
  "${reference_args[@]}" \
  --bqsr_recal_file "$recal_table" \
  --output "$meta_temp_dir/test1/output.bam" \
  --output_index "$meta_temp_dir/test1_index/recalibrated.bai"

log "Validating TEST 1 outputs..."
check_file_exists "$meta_temp_dir/test1/output.bam" "recalibrated BAM"
check_file_not_empty "$meta_temp_dir/test1/output.bam" "recalibrated BAM"
check_index_moved "$meta_temp_dir/test1_index/recalibrated.bai" "$meta_temp_dir/test1/output.bai" "output BAM index"
check_file_not_empty "$meta_temp_dir/test1_index/recalibrated.bai" "output BAM index"

log "Validating recalibrated BAM with ValidateSamFile..."
gatk ValidateSamFile \
  --INPUT "$meta_temp_dir/test1/output.bam" \
  --MODE SUMMARY \
  --VERBOSITY ERROR \
  > "$meta_temp_dir/test1/validate.log" 2>&1 || true
check_file_not_contains "$meta_temp_dir/test1/validate.log" "ERROR:" "ValidateSamFile summary (no ERROR-level issues)"

log "✅ TEST 1 completed successfully"

# --- Test Case 2: ApplyBQSR with options ---
log "Starting TEST 2: ApplyBQSR with intervals, static quantized quals and preserve-qscores-less-than"

mkdir -p "$meta_temp_dir/test2"

intervals_bed="$test_data_dir/intervals.bed"
echo -e "seq1\t0\t1000" > "$intervals_bed"
check_file_exists "$intervals_bed" "intervals BED file"

log "Executing $meta_name with advanced parameters..."
"$meta_executable" \
  --input "$test_data_dir/sample.sorted.bam" \
  --read_index "$test_data_dir/sample.sorted.bai" \
  "${reference_args[@]}" \
  --bqsr_recal_file "$recal_table" \
  --output "$meta_temp_dir/test2/output.bam" \
  --output_index "$meta_temp_dir/test2/output.bai" \
  --intervals "$intervals_bed" \
  --static_quantized_quals 10 \
  --static_quantized_quals 20 \
  --static_quantized_quals 30 \
  --preserve_qscores_less_than 20

log "Validating TEST 2 outputs..."
check_file_exists "$meta_temp_dir/test2/output.bam" "recalibrated BAM with options"
check_file_not_empty "$meta_temp_dir/test2/output.bam" "recalibrated BAM with options"
check_file_exists "$meta_temp_dir/test2/output.bai" "output BAM index (same path GATK writes)"

log "Validating recalibrated BAM with ValidateSamFile..."
gatk ValidateSamFile \
  --INPUT "$meta_temp_dir/test2/output.bam" \
  --MODE SUMMARY \
  --VERBOSITY ERROR \
  > "$meta_temp_dir/test2/validate.log" 2>&1 || true
check_file_not_contains "$meta_temp_dir/test2/validate.log" "ERROR:" "ValidateSamFile summary (no ERROR-level issues)"

log "✅ TEST 2 completed successfully"

# --- Test Case 3: ApplyBQSR without a reference ---
log "Starting TEST 3: ApplyBQSR without --reference"

mkdir -p "$meta_temp_dir/test3"

log "Executing $meta_name without a reference..."
"$meta_executable" \
  --input "$test_data_dir/sample.sorted.bam" \
  --read_index "$test_data_dir/sample.sorted.bai" \
  --bqsr_recal_file "$recal_table" \
  --output "$meta_temp_dir/test3/output.bam" \
  --output_index "$meta_temp_dir/test3/output.bai"

log "Validating TEST 3 outputs..."
check_file_exists "$meta_temp_dir/test3/output.bam" "recalibrated BAM without reference"
check_file_not_empty "$meta_temp_dir/test3/output.bam" "recalibrated BAM without reference"
check_file_exists "$meta_temp_dir/test3/output.bai" "output BAM index without reference"

log "✅ TEST 3 completed successfully"

# --- Test Case 4: Shared GATK engine options ---
log "Starting TEST 4: Shared GATK engine options"

mkdir -p "$meta_temp_dir/test4"

log "Executing $meta_name with --interval_padding, --read_validation_stringency, --lenient and --create_output_bam_index false..."
"$meta_executable" \
  --input "$test_data_dir/sample.sorted.bam" \
  --read_index "$test_data_dir/sample.sorted.bai" \
  "${reference_args[@]}" \
  --bqsr_recal_file "$recal_table" \
  --output "$meta_temp_dir/test4/output.bam" \
  --interval_padding 10 \
  --read_validation_stringency LENIENT \
  --lenient \
  --create_output_bam_index false

log "Validating TEST 4 outputs..."
check_file_exists "$meta_temp_dir/test4/output.bam" "recalibrated BAM (engine options)"
check_file_not_empty "$meta_temp_dir/test4/output.bam" "recalibrated BAM (engine options)"
check_file_not_exists "$meta_temp_dir/test4/output.bai" "output BAM index (should not exist with --create_output_bam_index false)"

log "✅ TEST 4 completed successfully"

# --- Test Case 5: Missing --output_index ---
log "Starting TEST 5: --output_index is required unless --create_output_bam_index is false"

mkdir -p "$meta_temp_dir/test5"

check_fails_with "Error: --output_index is required unless --create_output_bam_index is false." "$meta_temp_dir/test5/run.log" \
  "$meta_executable" \
  --input "$test_data_dir/sample.sorted.bam" \
  --read_index "$test_data_dir/sample.sorted.bai" \
  --bqsr_recal_file "$recal_table" \
  --output "$meta_temp_dir/test5/output.bam"
check_file_not_exists "$meta_temp_dir/test5/output.bam" "output BAM (should not be written)"

log "✅ TEST 5 completed successfully"

# --- Test Case 6: Multiple interval files with exclusions ---
log "Starting TEST 6: Multiple --intervals files (BED and bgzipped VCF) with --exclude_intervals"

mkdir -p "$meta_temp_dir/test6"

# Reads start at 1, 101, ..., 1901 (100bp each). The area is the union
# of seq1:1-600 (BED) and seq1:1550 (VCF), minus seq1:201-300 (BED).
printf 'seq1\t0\t600\n' > "$test_data_dir/intervals_a.bed"
printf 'seq1\t200\t300\n' > "$test_data_dir/exclude.bed"
create_test_known_sites_vcf "$test_data_dir/intervals_site.vcf" seq1 2000 1550
create_test_vcf_gz_with_separate_index "$test_data_dir/intervals_site.vcf" \
  "$test_data_dir/intervals_site.vcf.gz" "$test_data_dir/interval_indexes/intervals_site.vcf.gz.tbi"

log "Executing $meta_name with two --intervals files, --intervals_index and --exclude_intervals..."
"$meta_executable" \
  --input "$test_data_dir/sample.sorted.bam" \
  --read_index "$test_data_dir/sample.sorted.bai" \
  "${reference_args[@]}" \
  --bqsr_recal_file "$recal_table" \
  --output "$meta_temp_dir/test6/output.bam" \
  --output_index "$meta_temp_dir/test6/output.bai" \
  --intervals "$test_data_dir/intervals_a.bed" \
  --intervals "$test_data_dir/intervals_site.vcf.gz" \
  --intervals_index "$test_data_dir/interval_indexes/intervals_site.vcf.gz.tbi" \
  --interval_set_rule UNION \
  --exclude_intervals "$test_data_dir/exclude.bed"

log "Validating TEST 6 outputs..."
check_file_exists "$meta_temp_dir/test6/output.bam" "recalibrated BAM (multiple intervals)"
check_file_exists "$meta_temp_dir/test6/output.bai" "output BAM index (multiple intervals)"

# ApplyBQSR writes the reads that overlap the area
read_starts=$(samtools view "$meta_temp_dir/test6/output.bam" | cut -f4 | tr '\n' ' ')
expected_starts="1 101 301 401 501 1501 "
if [[ "$read_starts" == "$expected_starts" ]]; then
  log "✓ output BAM only has the reads in the selected territory: $read_starts"
else
  log_error "✗ output BAM read starts are '$read_starts', expected '$expected_starts'"
  exit 1
fi

log "✅ TEST 6 completed successfully"

# --- Test Case 7: CRAM input with a separate index ---
log "Starting TEST 7: CRAM --input with its .crai in another directory"

mkdir -p "$meta_temp_dir/test7"

create_test_cram "$test_data_dir/sample.sorted.bam" "$test_data_dir/reference.fasta" \
  "$test_data_dir/cram/sample.cram" "$test_data_dir/cram_index/sample.cram.crai"
check_file_not_exists "$test_data_dir/cram/sample.cram.crai" "index next to the CRAM input"

log "Executing $meta_name with a CRAM --input and --read_index..."
"$meta_executable" \
  --input "$test_data_dir/cram/sample.cram" \
  --read_index "$test_data_dir/cram_index/sample.cram.crai" \
  "${reference_args[@]}" \
  --bqsr_recal_file "$recal_table" \
  --output "$meta_temp_dir/test7/output.bam" \
  --output_index "$meta_temp_dir/test7/output.bai"

log "Validating TEST 7 outputs..."
check_file_exists "$meta_temp_dir/test7/output.bam" "recalibrated BAM from CRAM input"
check_file_exists "$meta_temp_dir/test7/output.bai" "output BAM index from CRAM input"

gatk ValidateSamFile \
  --INPUT "$meta_temp_dir/test7/output.bam" \
  --MODE SUMMARY \
  --VERBOSITY ERROR \
  > "$meta_temp_dir/test7/validate.log" 2>&1 || true
check_file_not_contains "$meta_temp_dir/test7/validate.log" "ERROR:" "ValidateSamFile summary (no ERROR-level issues)"

read_count=$(samtools view -c "$meta_temp_dir/test7/output.bam")
if [[ "$read_count" == "20" ]]; then
  log "✓ output BAM has all $read_count reads of the CRAM input"
else
  log_error "✗ output BAM has $read_count reads, expected 20"
  exit 1
fi

log "✅ TEST 7 completed successfully"

# --- Test Case 8: CRAM output ---
log "Starting TEST 8: CRAM --output with --output_index in another directory"

mkdir -p "$meta_temp_dir/test8" "$meta_temp_dir/test8_index"

log "Executing $meta_name with a CRAM --output..."
"$meta_executable" \
  --input "$test_data_dir/sample.sorted.bam" \
  --read_index "$test_data_dir/sample.sorted.bai" \
  "${reference_args[@]}" \
  --bqsr_recal_file "$recal_table" \
  --output "$meta_temp_dir/test8/output.cram" \
  --output_index "$meta_temp_dir/test8_index/output.cram.bai"

log "Validating TEST 8 outputs..."
check_file_exists "$meta_temp_dir/test8/output.cram" "recalibrated CRAM"
check_index_moved "$meta_temp_dir/test8_index/output.cram.bai" "$meta_temp_dir/test8/output.cram.bai" "output CRAM index"
check_file_not_exists "$meta_temp_dir/test8/output.cram.crai" "CRAI index next to the output CRAM"
check_file_not_exists "$meta_temp_dir/test8/output.bai" "BAI index next to the output CRAM"

if [[ "$(head -c 4 "$meta_temp_dir/test8/output.cram")" == "CRAM" ]]; then
  log "✓ output has the CRAM magic bytes"
else
  log_error "✗ output does not have the CRAM magic bytes"
  exit 1
fi
if [[ "$(head -c 4 "$meta_temp_dir/test8_index/output.cram.bai")" == $'BAI\x01' ]]; then
  log "✓ output index has the BAI magic bytes"
else
  log_error "✗ output index does not have the BAI magic bytes"
  exit 1
fi

gatk ValidateSamFile \
  --INPUT "$meta_temp_dir/test8/output.cram" \
  --REFERENCE_SEQUENCE "$test_data_dir/reference.fasta" \
  --MODE SUMMARY \
  --VERBOSITY ERROR \
  > "$meta_temp_dir/test8/validate.log" 2>&1 || true
check_file_not_contains "$meta_temp_dir/test8/validate.log" "ERROR:" "ValidateSamFile summary (no ERROR-level issues)"

log "✅ TEST 8 completed successfully"

# --- Test Case 9: CRAM input without a reference ---
log "Starting TEST 9: CRAM --input without --reference"

mkdir -p "$meta_temp_dir/test9"

check_fails_with "Error: --reference is required for CRAM input or output." "$meta_temp_dir/test9/run.log" \
  "$meta_executable" \
  --input "$test_data_dir/cram/sample.cram" \
  --read_index "$test_data_dir/cram_index/sample.cram.crai" \
  --bqsr_recal_file "$recal_table" \
  --output "$meta_temp_dir/test9/output.bam" \
  --output_index "$meta_temp_dir/test9/output.bai"
check_file_not_exists "$meta_temp_dir/test9/output.bam" "output BAM (should not be written)"

log "✅ TEST 9 completed successfully"

# --- Test Case 10: Inverted read filter ---
log "Starting TEST 10: --inverted_read_filter"

mkdir -p "$meta_temp_dir/test10"

# All test reads are mapped, so the inverted MappedReadFilter removes all of them
log "Executing $meta_name with --inverted_read_filter MappedReadFilter..."
"$meta_executable" \
  --input "$test_data_dir/sample.sorted.bam" \
  --read_index "$test_data_dir/sample.sorted.bai" \
  "${reference_args[@]}" \
  --bqsr_recal_file "$recal_table" \
  --output "$meta_temp_dir/test10/output.bam" \
  --output_index "$meta_temp_dir/test10/output.bai" \
  --inverted_read_filter MappedReadFilter

log "Validating TEST 10 outputs..."
read_count=$(samtools view -c "$meta_temp_dir/test10/output.bam")
if [[ "$read_count" -eq 0 ]]; then
  log "✓ Output BAM has no reads with the inverted MappedReadFilter"
else
  log_error "✗ Output BAM has $read_count reads, expected 0 with the inverted MappedReadFilter"
  exit 1
fi

log "✅ TEST 10 completed successfully"

# --- Test Case 11: Two input files ---
log "Starting TEST 11: Two --input files give one output BAM"

mkdir -p "$meta_temp_dir/test11"

# A second BAM with another read group and sample (10 reads at 51, 151, ...,
# 951), and a recalibration table made from both BAMs
positions_rg2=()
for ((i = 0; i < 10; i++)); do
  positions_rg2+=("$((i * 100 + 51))")
done
create_test_sam "$test_data_dir/reads_rg2.sam" seq1 2000 rg2 sample2 "${positions_rg2[@]}"
sed -i 's/^read/rg2_read/' "$test_data_dir/reads_rg2.sam"
sort_and_index_bam "$test_data_dir/reads_rg2.sam" "$test_data_dir/sample_rg2.sorted.bam"
check_file_exists "$test_data_dir/sample_rg2.sorted.bai" "second BAM index"

recal_table_two="$test_data_dir/recal_two.table"
gatk BaseRecalibrator \
  --input "$test_data_dir/sample.sorted.bam" \
  --input "$test_data_dir/sample_rg2.sorted.bam" \
  --reference "$test_data_dir/reference.fasta" \
  --known-sites "$known_sites_vcf" \
  --output "$recal_table_two" \
  --verbosity ERROR
check_file_exists "$recal_table_two" "BQSR recalibration table for two BAMs"

log "Executing $meta_name with two --input and --read_index files..."
"$meta_executable" \
  --input "$test_data_dir/sample.sorted.bam" \
  --read_index "$test_data_dir/sample.sorted.bai" \
  --input "$test_data_dir/sample_rg2.sorted.bam" \
  --read_index "$test_data_dir/sample_rg2.sorted.bai" \
  "${reference_args[@]}" \
  --bqsr_recal_file "$recal_table_two" \
  --output "$meta_temp_dir/test11/output.bam" \
  --output_index "$meta_temp_dir/test11/output.bai"

log "Validating TEST 11 outputs..."
check_file_exists "$meta_temp_dir/test11/output.bam" "recalibrated BAM (two inputs)"
check_file_exists "$meta_temp_dir/test11/output.bai" "output BAM index (two inputs)"

read_count=$(samtools view -c "$meta_temp_dir/test11/output.bam")
if [[ "$read_count" == "30" ]]; then
  log "✓ output BAM has the reads of both inputs ($read_count reads)"
else
  log_error "✗ output BAM has $read_count reads, expected 30"
  exit 1
fi
read_groups=$(samtools view -H "$meta_temp_dir/test11/output.bam" | grep '^@RG' | cut -f2 | sort | tr '\n' ' ')
if [[ "$read_groups" == "ID:rg1 ID:rg2 " ]]; then
  log "✓ output BAM header has the read groups of both inputs: $read_groups"
else
  log_error "✗ output BAM header read groups are '$read_groups', expected 'ID:rg1 ID:rg2 '"
  exit 1
fi

log "✅ TEST 11 completed successfully"

# --- Test Case 12: Read group and sample read filters ---
log "Starting TEST 12: --keep_read_group and --sample read filter values"

mkdir -p "$meta_temp_dir/test12"

# Check that a BAM has the expected number of reads and that all of them
# have the expected read group
check_read_group_reads() {
  local bam="$1" expected_count="$2" expected_rg="$3"
  local count other_count
  count=$(samtools view -c "$bam")
  other_count=$(samtools view "$bam" | grep -vc "RG:Z:${expected_rg}$" || true)
  if [[ "$count" == "$expected_count" && "$other_count" == "0" ]]; then
    log "✓ $(basename "$bam") has only the $count reads of read group $expected_rg"
  else
    log_error "✗ $(basename "$bam") has $count reads ($other_count not in $expected_rg), expected $expected_count reads of $expected_rg"
    exit 1
  fi
}

log "Executing $meta_name with --read_filter ReadGroupReadFilter --keep_read_group rg2..."
"$meta_executable" \
  --input "$test_data_dir/sample.sorted.bam" \
  --read_index "$test_data_dir/sample.sorted.bai" \
  --input "$test_data_dir/sample_rg2.sorted.bam" \
  --read_index "$test_data_dir/sample_rg2.sorted.bai" \
  --bqsr_recal_file "$recal_table_two" \
  --output "$meta_temp_dir/test12/keep_rg2.bam" \
  --output_index "$meta_temp_dir/test12/keep_rg2.bai" \
  --read_filter ReadGroupReadFilter \
  --keep_read_group rg2

check_read_group_reads "$meta_temp_dir/test12/keep_rg2.bam" 10 rg2

log "Executing $meta_name with --read_filter SampleReadFilter --sample sample1..."
"$meta_executable" \
  --input "$test_data_dir/sample.sorted.bam" \
  --read_index "$test_data_dir/sample.sorted.bai" \
  --input "$test_data_dir/sample_rg2.sorted.bam" \
  --read_index "$test_data_dir/sample_rg2.sorted.bai" \
  --bqsr_recal_file "$recal_table_two" \
  --output "$meta_temp_dir/test12/sample1.bam" \
  --output_index "$meta_temp_dir/test12/sample1.bai" \
  --read_filter SampleReadFilter \
  --sample sample1

check_read_group_reads "$meta_temp_dir/test12/sample1.bam" 20 rg1

log "✅ TEST 12 completed successfully"

# --- Test Case 13: --read_index count mismatch ---
log "Starting TEST 13: Fewer --read_index files than --input files fails"

mkdir -p "$meta_temp_dir/test13"

check_fails_with "Error: --read_index must be given once for each file in --input" "$meta_temp_dir/test13/run.log" \
  "$meta_executable" \
  --input "$test_data_dir/sample.sorted.bam" \
  --input "$test_data_dir/sample_rg2.sorted.bam" \
  --read_index "$test_data_dir/sample.sorted.bai" \
  --bqsr_recal_file "$recal_table_two" \
  --output "$meta_temp_dir/test13/output.bam" \
  --output_index "$meta_temp_dir/test13/output.bai"
check_file_not_exists "$meta_temp_dir/test13/output.bam" "output BAM (should not be written)"

log "✅ TEST 13 completed successfully"

print_test_summary "All tests"
