#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

source "$meta_resources_dir/test_helpers.sh"
setup_test_env

log "Starting tests for $meta_name"

test_dir="$meta_temp_dir/test_data"
mkdir -p "$test_dir"

# --- Test data -------------------------------------------------------------
# A 20 kb reference of deterministic pseudo-random sequence (MINSTD LCG in awk;
# the engine image has no python3). The reference has to be this long and this
# non-repetitive because the long-read presets index with a large k and w -
# map-hifi uses k=19, w=19 - and find no usable minimizers in a short,
# low-complexity sequence.
awk 'BEGIN {
  bases = "ACGT"
  x = 42
  line = ""
  print ">seq1"
  for (i = 0; i < 20000; i++) {
    x = (x * 16807) % 2147483647
    line = line substr(bases, int(x / 2147483647 * 4) + 1, 1)
    if (length(line) == 70) { print line; line = "" }
  }
  if (length(line) > 0) print line
}' > "$test_dir/ref.fasta"

# Reads are 1 kb slices of the reference with a single substitution each, so
# the alignments are realistic rather than trivial exact matches.
ref_seq=$(grep -v '^>' "$test_dir/ref.fasta" | tr -d '\n')
: > "$test_dir/reads.fastq"
read_nr=0
for offset in 500 5000 12000; do
  read_nr=$((read_nr + 1))
  read_seq="${ref_seq:$offset:1000}"
  read_seq="${read_seq:0:500}A${read_seq:501}"
  quality=$(printf 'I%.0s' $(seq 1 ${#read_seq}))
  printf '@read%d\n%s\n+\n%s\n' "$read_nr" "$read_seq" "$quality" >> "$test_dir/reads.fastq"
done

check_file_exists "$test_dir/ref.fasta" "reference FASTA"
check_file_exists "$test_dir/reads.fastq" "query FASTQ"

# --- TEST 1: PAF output ---
log "TEST 1: alignment to PAF"
"$meta_executable" \
  --reference "$test_dir/ref.fasta" \
  --query "$test_dir/reads.fastq" \
  --output "$meta_temp_dir/out.paf"

check_file_exists "$meta_temp_dir/out.paf" "PAF output"
check_file_not_empty "$meta_temp_dir/out.paf" "PAF output"
check_file_contains "$meta_temp_dir/out.paf" "seq1" "PAF alignment to seq1"
log "OK: TEST 1 passed"

# --- TEST 2: sorted + indexed BAM output ---
log "TEST 2: alignment to sorted BAM"
"$meta_executable" \
  --reference "$test_dir/ref.fasta" \
  --query "$test_dir/reads.fastq" \
  --bam \
  --output "$meta_temp_dir/out.bam"

check_file_exists "$meta_temp_dir/out.bam" "BAM output"
check_file_not_empty "$meta_temp_dir/out.bam" "BAM output"
check_file_exists "$meta_temp_dir/out.bam.bai" "BAM index"

if ! samtools view "$meta_temp_dir/out.bam" | grep -q "seq1"; then
  log_error "BAM alignment to seq1 missing"
  exit 1
fi
log "- BAM contains alignment to seq1"

if ! samtools view -H "$meta_temp_dir/out.bam" | grep -q "SO:coordinate"; then
  log_error "BAM is not marked as coordinate-sorted"
  exit 1
fi
log "- BAM is coordinate-sorted"
log "OK: TEST 2 passed"

# --- TEST 3: BAM with an explicitly declared index path ---
# Without a declared output the .bai is invisible to viash: it is not published
# by the Nextflow runner and is left owned by the container user.
log "TEST 3: BAM with explicit --output_index"
"$meta_executable" \
  --reference "$test_dir/ref.fasta" \
  --query "$test_dir/reads.fastq" \
  --bam \
  --output "$meta_temp_dir/named.bam" \
  --output_index "$meta_temp_dir/named.bam.bai"

check_file_not_empty "$meta_temp_dir/named.bam" "BAM output"
check_file_not_empty "$meta_temp_dir/named.bam.bai" "explicitly named BAM index"
log "OK: TEST 3 passed"

# --- TEST 4: preset the old engine image could not run ---
# map-hifi was added in minimap2 2.19; the component used to ship 2.17, where
# this preset failed with "unknown preset".
log "TEST 4: map-hifi preset"
"$meta_executable" \
  --reference "$test_dir/ref.fasta" \
  --query "$test_dir/reads.fastq" \
  --preset map-hifi \
  --cigar_paf \
  --output "$meta_temp_dir/hifi.paf"

check_file_not_empty "$meta_temp_dir/hifi.paf" "map-hifi PAF output"
check_file_contains "$meta_temp_dir/hifi.paf" "cg:Z:" "CIGAR string in PAF (--cigar_paf)"
log "OK: TEST 4 passed"

# --- TEST 5: contradictory flags are rejected ---
log "TEST 5: mutually exclusive flags are rejected"
if "$meta_executable" \
  --reference "$test_dir/ref.fasta" \
  --query "$test_dir/reads.fastq" \
  --cigar_bam \
  --output "$meta_temp_dir/bad.paf" > /dev/null 2>&1; then
  log_error "--cigar_bam without --bam should have failed"
  exit 1
fi
log "- --cigar_bam without --bam is rejected"

if "$meta_executable" \
  --reference "$test_dir/ref.fasta" \
  --query "$test_dir/reads.fastq" \
  --output "$meta_temp_dir/bad2.paf" \
  --output_index "$meta_temp_dir/bad2.paf.bai" > /dev/null 2>&1; then
  log_error "--output_index without --bam should have failed"
  exit 1
fi
log "- --output_index without --bam is rejected"
log "OK: TEST 5 passed"

print_test_summary "$meta_name tests passed"
