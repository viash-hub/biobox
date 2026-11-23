#!/bin/bash

## VIASH START
## VIASH END

# Source the centralized test helpers and the gatk4-specific helpers
source "$meta_resources_dir/test_helpers.sh"
source "$meta_resources_dir/gatk4/test_helpers.sh"

gatk4_test_init

#############################################
# Shared fixtures
#############################################

create_test_reference "$test_data_dir/reference.fasta" 1 2000
check_file_exists "$test_data_dir/reference.dict" "reference sequence dictionary"
check_file_exists "$test_data_dir/reference.fasta.fai" "reference FASTA index"
reference_args=(
  --reference "$test_data_dir/reference.fasta"
  --reference_fai "$test_data_dir/reference.fasta.fai"
  --reference_dict "$test_data_dir/reference.dict"
)

# First BAM: 20 reads (rg1, sample1) at 1, 101, ..., 1901
create_test_sam_reads "$test_data_dir/reads.sam" seq1 2000 rg1 sample1 20 100
sort_and_index_bam "$test_data_dir/reads.sam" "$test_data_dir/sample.sorted.bam"
check_file_exists "$test_data_dir/sample.sorted.bai" "sorted test BAM index"

# The same reads as a CRAM, with the .crai in another directory
create_test_cram "$test_data_dir/sample.sorted.bam" "$test_data_dir/reference.fasta" \
  "$test_data_dir/cram/sample.cram" "$test_data_dir/cram_index/sample.cram.crai"
check_file_not_exists "$test_data_dir/cram/sample.cram.crai" "index next to the CRAM input"

# Second BAM: another read group and sample (10 reads at 51, 151, ..., 951)
positions_rg2=()
for ((i = 0; i < 10; i++)); do
  positions_rg2+=("$((i * 100 + 51))")
done
create_test_sam "$test_data_dir/reads_rg2.sam" seq1 2000 rg2 sample2 "${positions_rg2[@]}"
sed -i 's/^read/rg2_read/' "$test_data_dir/reads_rg2.sam"
sort_and_index_bam "$test_data_dir/reads_rg2.sam" "$test_data_dir/sample_rg2.sorted.bam"
check_file_exists "$test_data_dir/sample_rg2.sorted.bai" "second BAM index"

# Recalibration tables for the first BAM, and for both BAMs
known_sites_vcf="$test_data_dir/known_sites.vcf"
create_test_known_sites_vcf "$known_sites_vcf" seq1 2000 500

log "Running BaseRecalibrator to produce the recal tables..."
recal_table="$test_data_dir/recal.table"
gatk BaseRecalibrator \
  --input "$test_data_dir/sample.sorted.bam" \
  --reference "$test_data_dir/reference.fasta" \
  --known-sites "$known_sites_vcf" \
  --output "$recal_table" \
  --verbosity ERROR
check_file_not_empty "$recal_table" "BQSR recalibration table"

recal_table_two="$test_data_dir/recal_two.table"
gatk BaseRecalibrator \
  --input "$test_data_dir/sample.sorted.bam" \
  --input "$test_data_dir/sample_rg2.sorted.bam" \
  --reference "$test_data_dir/reference.fasta" \
  --known-sites "$known_sites_vcf" \
  --output "$recal_table_two" \
  --verbosity ERROR
check_file_not_empty "$recal_table_two" "BQSR recalibration table for two BAMs"

# Check that a BAM or CRAM has no ERROR-level issues with ValidateSamFile
#
# Usage: check_validate_sam "/path/to/output.bam" [ValidateSamFile arguments...]
check_validate_sam() {
  local file="$1"
  shift
  gatk ValidateSamFile \
    --INPUT "$file" \
    "$@" \
    --MODE SUMMARY \
    --VERBOSITY ERROR \
    > "${file}.validate.log" 2>&1 || true
  check_file_not_contains "${file}.validate.log" "ERROR:" "ValidateSamFile summary of $(basename "$file") (no ERROR-level issues)"
}

# Check that a BAM has the expected number of reads and that all of them
# have the expected read group
#
# Usage: check_read_group_reads "/path/to/output.bam" count read_group
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

# --- Test Case 1: Basic run (defaults) ---
log "Starting TEST 1: Basic run (defaults)"

mkdir -p "$meta_temp_dir/test1" "$meta_temp_dir/test1_index"

"$meta_executable" \
  --input "$test_data_dir/sample.sorted.bam" \
  --read_index "$test_data_dir/sample.sorted.bai" \
  "${reference_args[@]}" \
  --bqsr_recal_file "$recal_table" \
  --output "$meta_temp_dir/test1/output.bam" \
  --output_index "$meta_temp_dir/test1_index/recalibrated.bai"

log "Validating TEST 1 outputs..."
check_file_not_empty "$meta_temp_dir/test1/output.bam" "recalibrated BAM"
check_index_moved "$meta_temp_dir/test1_index/recalibrated.bai" "$meta_temp_dir/test1/output.bai" "output BAM index"
check_file_not_empty "$meta_temp_dir/test1_index/recalibrated.bai" "output BAM index"
check_validate_sam "$meta_temp_dir/test1/output.bam"

log "✅ TEST 1 completed successfully"

# --- Test Case 2: Inputs and indexes ---
log "Starting TEST 2: Inputs and indexes"

test2_dir="$meta_temp_dir/test2"
mkdir -p "$test2_dir" "$meta_temp_dir/test2_index"

# Run A: CRAM input with a separate .crai, CRAM output with the index in
# another directory
log "Executing $meta_name with a CRAM --input and a CRAM --output..."
"$meta_executable" \
  --input "$test_data_dir/cram/sample.cram" \
  --read_index "$test_data_dir/cram_index/sample.cram.crai" \
  "${reference_args[@]}" \
  --bqsr_recal_file "$recal_table" \
  --output "$test2_dir/output.cram" \
  --output_index "$meta_temp_dir/test2_index/output.cram.bai"

check_file_exists "$test2_dir/output.cram" "recalibrated CRAM"
check_index_moved "$meta_temp_dir/test2_index/output.cram.bai" "$test2_dir/output.cram.bai" "output CRAM index"
check_file_not_exists "$test2_dir/output.cram.crai" "CRAI index next to the output CRAM"
check_file_not_exists "$test2_dir/output.bai" "BAI index next to the output CRAM"
if [[ "$(head -c 4 "$test2_dir/output.cram")" == "CRAM" ]]; then
  log "✓ output has the CRAM magic bytes"
else
  log_error "✗ output does not have the CRAM magic bytes"
  exit 1
fi
if [[ "$(head -c 4 "$meta_temp_dir/test2_index/output.cram.bai")" == $'BAI\x01' ]]; then
  log "✓ output index has the BAI magic bytes"
else
  log_error "✗ output index does not have the BAI magic bytes"
  exit 1
fi
read_count=$(samtools view -c -T "$test_data_dir/reference.fasta" "$test2_dir/output.cram")
if [[ "$read_count" == "20" ]]; then
  log "✓ output CRAM has all $read_count reads of the CRAM input"
else
  log_error "✗ output CRAM has $read_count reads, expected 20"
  exit 1
fi
check_validate_sam "$test2_dir/output.cram" --REFERENCE_SEQUENCE "$test_data_dir/reference.fasta"

# Run B: two BAM inputs, static quantized qualities and
# --preserve_qscores_less_than, index at the path GATK writes
log "Executing $meta_name with two --input files, --static_quantized_quals and --preserve_qscores_less_than..."
"$meta_executable" \
  --input "$test_data_dir/sample.sorted.bam" \
  --read_index "$test_data_dir/sample.sorted.bai" \
  --input "$test_data_dir/sample_rg2.sorted.bam" \
  --read_index "$test_data_dir/sample_rg2.sorted.bai" \
  "${reference_args[@]}" \
  --bqsr_recal_file "$recal_table_two" \
  --output "$test2_dir/two_inputs.bam" \
  --output_index "$test2_dir/two_inputs.bai" \
  --static_quantized_quals 10 \
  --static_quantized_quals 20 \
  --static_quantized_quals 30 \
  --preserve_qscores_less_than 20

check_file_exists "$test2_dir/two_inputs.bam" "recalibrated BAM (two inputs)"
check_file_exists "$test2_dir/two_inputs.bai" "output BAM index (same path GATK writes)"
read_count=$(samtools view -c "$test2_dir/two_inputs.bam")
if [[ "$read_count" == "30" ]]; then
  log "✓ output BAM has the reads of both inputs ($read_count reads)"
else
  log_error "✗ output BAM has $read_count reads, expected 30"
  exit 1
fi
read_groups=$(samtools view -H "$test2_dir/two_inputs.bam" | grep '^@RG' | cut -f2 | sort | tr '\n' ' ')
if [[ "$read_groups" == "ID:rg1 ID:rg2 " ]]; then
  log "✓ output BAM header has the read groups of both inputs: $read_groups"
else
  log_error "✗ output BAM header read groups are '$read_groups', expected 'ID:rg1 ID:rg2 '"
  exit 1
fi
# The rg1 reads match the reference, so their recalibrated base qualities
# can only have the static values 10, 20 and 30 (Phred+33 '+', '5' and '?').
# The rg2 reads do not match the reference (they start out of phase with the
# ATCG repeat), and GATK gives them a quality below the static values.
other_quals=$(samtools view -r rg1 "$test2_dir/two_inputs.bam" | cut -f11 | tr -d '+5?\n' | wc -c)
if [[ "$other_quals" -eq 0 ]]; then
  log "✓ rg1 reads have only the static quantized base qualities (--static_quantized_quals)"
else
  log_error "✗ rg1 reads have $other_quals base qualities that are not 10, 20 or 30"
  exit 1
fi
# --preserve_qscores_less_than has no visible effect on these reads (all
# input qualities are 40), so check the GATK command line in the header
samtools view -H "$test2_dir/two_inputs.bam" > "$test2_dir/two_inputs_header.sam"
check_file_contains "$test2_dir/two_inputs_header.sam" "preserve-qscores-less-than 20" "@PG command line (--preserve_qscores_less_than)"

log "✅ TEST 2 completed successfully"

# --- Test Case 3: Intervals ---
log "Starting TEST 3: Multiple --intervals files (BED and bgzipped VCF) with --exclude_intervals"

mkdir -p "$meta_temp_dir/test3"

# Reads start at 1, 101, ..., 1901 (100bp each). The area is the union
# of seq1:1-600 (BED) and seq1:1550 (VCF), minus seq1:201-300 (BED).
printf 'seq1\t0\t600\n' > "$test_data_dir/intervals_a.bed"
printf 'seq1\t200\t300\n' > "$test_data_dir/exclude.bed"
create_test_known_sites_vcf --no-index "$test_data_dir/intervals_site.vcf" seq1 2000 1550
create_test_vcf_gz_with_separate_index "$test_data_dir/intervals_site.vcf" \
  "$test_data_dir/intervals_site.vcf.gz" "$test_data_dir/interval_indexes/intervals_site.vcf.gz.tbi"

"$meta_executable" \
  --input "$test_data_dir/sample.sorted.bam" \
  --read_index "$test_data_dir/sample.sorted.bai" \
  "${reference_args[@]}" \
  --bqsr_recal_file "$recal_table" \
  --output "$meta_temp_dir/test3/output.bam" \
  --output_index "$meta_temp_dir/test3/output.bai" \
  --intervals "$test_data_dir/intervals_a.bed" \
  --intervals "$test_data_dir/intervals_site.vcf.gz" \
  --intervals_index "$test_data_dir/interval_indexes/intervals_site.vcf.gz.tbi" \
  --interval_set_rule UNION \
  --exclude_intervals "$test_data_dir/exclude.bed"

log "Validating TEST 3 outputs..."
check_file_exists "$meta_temp_dir/test3/output.bam" "recalibrated BAM (multiple intervals)"
check_file_exists "$meta_temp_dir/test3/output.bai" "output BAM index (same path GATK writes)"

# ApplyBQSR writes the reads that overlap the area
read_starts=$(samtools view "$meta_temp_dir/test3/output.bam" | cut -f4 | tr '\n' ' ')
expected_starts="1 101 301 401 501 1501 "
if [[ "$read_starts" == "$expected_starts" ]]; then
  log "✓ output BAM only has the reads in the selected territory: $read_starts"
else
  log_error "✗ output BAM read starts are '$read_starts', expected '$expected_starts'"
  exit 1
fi

log "✅ TEST 3 completed successfully"

# --- Test Case 4: Read filters and index options ---
log "Starting TEST 4: Read filters and index options"

test4_dir="$meta_temp_dir/test4"
mkdir -p "$test4_dir"

# Run A: all test reads are mapped, so the inverted MappedReadFilter removes
# all of them. No index and no --output_index with
# --create_output_bam_index false. The engine options are pass-through.
log "Executing $meta_name with --inverted_read_filter, --create_output_bam_index false and engine options..."
"$meta_executable" \
  --input "$test_data_dir/sample.sorted.bam" \
  --read_index "$test_data_dir/sample.sorted.bai" \
  "${reference_args[@]}" \
  --bqsr_recal_file "$recal_table" \
  --output "$test4_dir/inverted.bam" \
  --inverted_read_filter MappedReadFilter \
  --interval_padding 10 \
  --read_validation_stringency LENIENT \
  --lenient \
  --create_output_bam_index false

check_file_not_empty "$test4_dir/inverted.bam" "recalibrated BAM (inverted read filter)"
check_file_not_exists "$test4_dir/inverted.bai" "output BAM index (should not exist with --create_output_bam_index false)"
read_count=$(samtools view -c "$test4_dir/inverted.bam")
if [[ "$read_count" -eq 0 ]]; then
  log "✓ Output BAM has no reads with the inverted MappedReadFilter"
else
  log_error "✗ Output BAM has $read_count reads, expected 0 with the inverted MappedReadFilter"
  exit 1
fi

# Run B: BAM inputs without --reference, keep one read group
log "Executing $meta_name without --reference, with --read_filter ReadGroupReadFilter --keep_read_group rg2..."
"$meta_executable" \
  --input "$test_data_dir/sample.sorted.bam" \
  --read_index "$test_data_dir/sample.sorted.bai" \
  --input "$test_data_dir/sample_rg2.sorted.bam" \
  --read_index "$test_data_dir/sample_rg2.sorted.bai" \
  --bqsr_recal_file "$recal_table_two" \
  --output "$test4_dir/keep_rg2.bam" \
  --output_index "$test4_dir/keep_rg2.bai" \
  --read_filter ReadGroupReadFilter \
  --keep_read_group rg2

check_file_exists "$test4_dir/keep_rg2.bai" "output BAM index without reference"
check_read_group_reads "$test4_dir/keep_rg2.bam" 10 rg2

# Run C: BAM inputs without --reference, keep one sample
log "Executing $meta_name without --reference, with --read_filter SampleReadFilter --sample sample1..."
"$meta_executable" \
  --input "$test_data_dir/sample.sorted.bam" \
  --read_index "$test_data_dir/sample.sorted.bai" \
  --input "$test_data_dir/sample_rg2.sorted.bam" \
  --read_index "$test_data_dir/sample_rg2.sorted.bai" \
  --bqsr_recal_file "$recal_table_two" \
  --output "$test4_dir/sample1.bam" \
  --output_index "$test4_dir/sample1.bai" \
  --read_filter SampleReadFilter \
  --sample sample1

check_read_group_reads "$test4_dir/sample1.bam" 20 rg1

log "✅ TEST 4 completed successfully"

# --- Test Case 5: Argument errors ---
log "Starting TEST 5: Argument errors"

test5_dir="$meta_temp_dir/test5"
mkdir -p "$test5_dir"

log "Executing $meta_name without --output_index..."
check_fails_with "Error: --output_index is required unless --create_output_bam_index is false." "$test5_dir/no_output_index.log" \
  "$meta_executable" \
  --input "$test_data_dir/sample.sorted.bam" \
  --read_index "$test_data_dir/sample.sorted.bai" \
  --bqsr_recal_file "$recal_table" \
  --output "$test5_dir/no_output_index.bam"
check_file_not_exists "$test5_dir/no_output_index.bam" "output BAM (should not be written)"

log "Executing $meta_name with a CRAM --input and no --reference..."
check_fails_with "Error: --reference is required for CRAM input or output." "$test5_dir/cram_no_reference.log" \
  "$meta_executable" \
  --input "$test_data_dir/cram/sample.cram" \
  --read_index "$test_data_dir/cram_index/sample.cram.crai" \
  --bqsr_recal_file "$recal_table" \
  --output "$test5_dir/cram_no_reference.bam" \
  --output_index "$test5_dir/cram_no_reference.bai"
check_file_not_exists "$test5_dir/cram_no_reference.bam" "output BAM (should not be written)"

log "Executing $meta_name with fewer --read_index than --input files..."
check_fails_with "Error: --read_index must be given once for each file in --input" "$test5_dir/index_count.log" \
  "$meta_executable" \
  --input "$test_data_dir/sample.sorted.bam" \
  --input "$test_data_dir/sample_rg2.sorted.bam" \
  --read_index "$test_data_dir/sample.sorted.bai" \
  --bqsr_recal_file "$recal_table_two" \
  --output "$test5_dir/index_count.bam" \
  --output_index "$test5_dir/index_count.bai"
check_file_not_exists "$test5_dir/index_count.bam" "output BAM (should not be written)"

log "✅ TEST 5 completed successfully"

print_test_summary "All tests"
