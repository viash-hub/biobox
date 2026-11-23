#!/bin/bash

## VIASH START
## VIASH END

# Source the centralized test helpers
source "$meta_resources_dir/test_helpers.sh"

# Initialize test environment with strict error handling
setup_test_env

#############################################
# Helper functions
#############################################

# Reverse complement a DNA sequence
revcomp() {
  local seq="$1" out="" i
  for (( i=${#seq}-1; i>=0; i-- )); do
    out+="${seq:i:1}"
  done
  echo "$out" | tr 'ACGT' 'TGCA'
}

# Append a FASTQ record with constant high base qualities
write_fastq_record() {
  local file="$1" name="$2" seq="$3"
  printf "@%s\n%s\n+\n%s\n" "$name" "$seq" "$(printf "%*s" "${#seq}" "" | tr ' ' 'I')" >> "$file"
}

# Count SAM records that are mapped (FLAG bit 0x4 unset)
count_mapped() {
  local sam="$1" qname flag rest n=0
  while IFS=$'\t' read -r qname flag rest; do
    if [[ "$qname" != @* ]] && (( (flag & 4) == 0 )); then
      n=$((n + 1))
    fi
  done < "$sam"
  echo "$n"
}

check_count() {
  local actual="$1" expected="$2" description="$3"
  if [[ "$actual" -ne "$expected" ]]; then
    log_error "Expected $expected $description, found $actual"
    exit 1
  fi
  log "✓ $actual $description"
}

#############################################
# Test data
#############################################

log "Starting tests for $meta_name"

test_data_dir="$meta_temp_dir/test_data"
mkdir -p "$test_data_dir"

# Generate test reference genome
log "Generating test reference genome..."
create_test_fasta "$test_data_dir/reference.fasta" 1 500
check_file_exists "$test_data_dir/reference.fasta" "test reference genome"

# Build a BWA index directory, laid out like the bwa_mem3_index output
log "Building BWA index for alignment tests..."
mkdir -p "$test_data_dir/index"
bwa-mem3 index -p "$test_data_dir/index/reference" "$test_data_dir/reference.fasta" >/dev/null 2>&1

# Verify index was created (bwa-mem3 does not write .0123 by default)
for suffix in amb ann bwt.2bit.64 pac; do
  check_file_exists "$test_data_dir/index/reference.$suffix" "BWA index file reference.$suffix"
done

# Take the reads from the reference, so they align to it
log "Generating test FASTQ files from the reference..."
ref_seq=$(grep -v '^>' "$test_data_dir/reference.fasta" | tr -d '\n')
read_length=60
insert_size=200
num_reads=10

: > "$test_data_dir/reads_single.fastq"
: > "$test_data_dir/reads_R1.fastq"
: > "$test_data_dir/reads_R2.fastq"
: > "$test_data_dir/reads_meth.fastq"

for (( i=0; i<num_reads; i++ )); do
  pos=$(( i * 30 ))

  # Single-end: forward-strand fragment
  write_fastq_record "$test_data_dir/reads_single.fastq" "read$i" "${ref_seq:pos:read_length}"

  # Paired-end (FR): R1 forward from the fragment start, R2 reverse complement of the fragment end
  r2_seq=$(revcomp "${ref_seq:pos+insert_size-read_length:read_length}")
  write_fastq_record "$test_data_dir/reads_R1.fastq" "pair$i" "${ref_seq:pos:read_length}"
  write_fastq_record "$test_data_dir/reads_R2.fastq" "pair$i" "$r2_seq"

  # EM-seq: fully unmethylated, so every C on the forward strand reads as T
  meth_seq="${ref_seq:pos:read_length}"
  write_fastq_record "$test_data_dir/reads_meth.fastq" "meth$i" "${meth_seq//C/T}"
done

check_file_exists "$test_data_dir/reads_single.fastq" "single-end reads"
check_file_exists "$test_data_dir/reads_R1.fastq" "paired-end R1 reads"
check_file_exists "$test_data_dir/reads_R2.fastq" "paired-end R2 reads"
check_file_exists "$test_data_dir/reads_meth.fastq" "methylation reads"

#############################################
# Tests
#############################################

# --- Test Case 1: Single-end alignment ---
log "Starting TEST 1: Single-end BWA MEM alignment"

"$meta_executable" \
  --index "$test_data_dir/index" \
  --reads1 "$test_data_dir/reads_single.fastq" \
  --output "$meta_temp_dir/single_end.sam"

log "Validating TEST 1 outputs..."
check_file_exists "$meta_temp_dir/single_end.sam" "single-end SAM output"
check_file_contains "$meta_temp_dir/single_end.sam" "@SQ"
check_count "$(count_mapped "$meta_temp_dir/single_end.sam")" "$num_reads" "mapped single-end reads"

log "✅ TEST 1 completed successfully"

# --- Test Case 2: Paired-end alignment ---
log "Starting TEST 2: Paired-end BWA MEM alignment"

"$meta_executable" \
  --index "$test_data_dir/index" \
  --reads1 "$test_data_dir/reads_R1.fastq" \
  --reads2 "$test_data_dir/reads_R2.fastq" \
  --output "$meta_temp_dir/paired_end.sam"

log "Validating TEST 2 outputs..."
check_file_exists "$meta_temp_dir/paired_end.sam" "paired-end SAM output"
check_file_contains "$meta_temp_dir/paired_end.sam" "@SQ"
check_count "$(count_mapped "$meta_temp_dir/paired_end.sam")" "$(( num_reads * 2 ))" "mapped paired-end reads"

log "✅ TEST 2 completed successfully"

# --- Test Case 3: BAM output ---
log "Starting TEST 3: BWA MEM with BAM output"

"$meta_executable" \
  --index "$test_data_dir/index" \
  --reads1 "$test_data_dir/reads_R1.fastq" \
  --reads2 "$test_data_dir/reads_R2.fastq" \
  --output "$meta_temp_dir/paired_end.bam" \
  --bam_compression 6

log "Validating TEST 3 outputs..."
check_file_exists "$meta_temp_dir/paired_end.bam" "paired-end BAM output"
check_file_not_empty "$meta_temp_dir/paired_end.bam" "paired-end BAM output"

# BAM is BGZF-compressed, so it starts with the gzip magic bytes 1f 8b
if [[ "$(head -c 2 "$meta_temp_dir/paired_end.bam" | od -An -tx1 | tr -d ' \n')" == "1f8b" ]]; then
  log "✓ Output is BGZF-compressed BAM"
else
  log_error "Output is not BGZF-compressed BAM"
  exit 1
fi

log "✅ TEST 3 completed successfully"

# --- Test Case 4: Methylation (--meth) alignment ---
log "Starting TEST 4: BWA MEM methylation alignment"

# --meth needs the original index plus the seed index <fa>.meth.* in the same directory
log "Building BWA methylation index..."
mkdir -p "$test_data_dir/meth_index"
cp "$test_data_dir/reference.fasta" "$test_data_dir/meth_index/"
bwa-mem3 index --meth "$test_data_dir/meth_index/reference.fasta" >/dev/null 2>&1

for suffix in amb ann bwt.2bit.64 pac meth.amb meth.ann meth.bwt.2bit.64 meth.pac; do
  check_file_exists "$test_data_dir/meth_index/reference.fasta.$suffix" "BWA meth index file reference.fasta.$suffix"
done

"$meta_executable" \
  --index "$test_data_dir/meth_index" \
  --reads1 "$test_data_dir/reads_meth.fastq" \
  --output "$meta_temp_dir/meth.sam" \
  --meth emseq

log "Validating TEST 4 outputs..."
check_file_exists "$meta_temp_dir/meth.sam" "methylation SAM output"
check_count "$(count_mapped "$meta_temp_dir/meth.sam")" "$num_reads" "mapped methylation reads"

# Mapped --meth records carry the Bismark-style methylation call, genome and read conversion tags
for tag in XM XG XR; do
  check_file_contains "$meta_temp_dir/meth.sam" "$(printf '\t%s:Z:' "$tag")"
done

log "✅ TEST 4 completed successfully"

print_test_summary "All tests completed successfully"
