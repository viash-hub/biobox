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
# the alignments are realistic rather than trivial exact matches. The read
# name carries the 0-based offset it was cut from, which is what the placement
# checks below compare against.
ref_seq=$(grep -v '^>' "$test_dir/ref.fasta" | tr -d '\n')
: > "$test_dir/reads.fastq"
offsets=(500 5000 12000)
for offset in "${offsets[@]}"; do
  read_seq="${ref_seq:$offset:1000}"
  read_seq="${read_seq:0:500}A${read_seq:501}"
  quality=$(printf 'I%.0s' $(seq 1 ${#read_seq}))
  printf '@read_%d\n%s\n+\n%s\n' "$offset" "$read_seq" "$quality" >> "$test_dir/reads.fastq"
done
n_reads=${#offsets[@]}

check_file_exists "$test_dir/ref.fasta" "reference FASTA"
check_file_exists "$test_dir/reads.fastq" "query FASTQ"

# --- Helpers ---------------------------------------------------------------
# SAM/BAM is read with samtools, never with grep over the raw file, so no check
# can be satisfied by a header line.
# Primary mapped records: drop unmapped (0x4), secondary (0x100), supplementary (0x800)
count_primary() { samtools view -c -F 0x904 "$1"; }
# QNAME/FLAG/RNAME/POS/MAPQ/CIGAR of every record, sorted: the comparison key
# used to assert two runs produced the same alignments.
alignment_key() { samtools view "$1" | cut -f1-6 | sort; }

# Every read must map to seq1 within 50 bp of the offset it was cut from.
# Input is "<qname> <target> <0-based start>" per line, one line per read.
check_placement() {
  local placements="$1" label="$2" n_ok
  n_ok=$(awk -v n="$n_reads" '
    { split($1, a, "_"); d = $3 - a[2]; if (d < 0) d = -d
      if ($2 == "seq1" && d <= 50) ok[$1] = 1 }
    END { c = 0; for (k in ok) c++; print c }' "$placements")
  if [[ "$n_ok" -ne "$n_reads" ]]; then
    log_error "$label: only $n_ok of $n_reads reads placed at their source offset"
    cat "$placements" >&2
    exit 1
  fi
  log "- $label: all $n_reads reads placed at their source offset"
}
# PAF columns 1/6/8 = query name, target name, 0-based target start; primary
# alignments only (tp:A:P)
paf_placement() { awk '/\ttp:A:P/ { print $1, $6, $8 }' "$1" > "$2"; }
# SAM POS is 1-based
bam_placement() { samtools view -F 0x904 "$1" | awk '{ print $1, $3, $4 - 1 }' > "$2"; }

# Assert the file is really BAM. samtools quickcheck accepts plain SAM too, so
# ask htslib what the file is: "BAM version 1 ..." vs "SAM version 1.6 ...".
check_is_bam() {
  local file="$1" label="$2" format
  if ! samtools quickcheck "$file"; then
    log_error "$label failed samtools quickcheck"
    exit 1
  fi
  format=$(htsfile "$file")
  if [[ "$format" != *"BAM version"* ]]; then
    log_error "$label is not BAM: $format"
    exit 1
  fi
  log "- $label is a valid BAM file"
}

# A negative test that only checks the exit code passes on any failure at all,
# so each one is pinned to the component's own error message.
check_rejected() {
  local log_file="$1" message="$2" label="$3"
  shift 3
  if "$meta_executable" "$@" > "$log_file" 2>&1; then
    log_error "$label should have failed"
    exit 1
  fi
  check_file_contains "$log_file" "$message" "error message ($label)"
}

# --- TEST 1: PAF output ---
log "TEST 1: alignment to PAF"
"$meta_executable" \
  --reference "$test_dir/ref.fasta" \
  --query "$test_dir/reads.fastq" \
  --output "$meta_temp_dir/out.paf"

check_file_not_empty "$meta_temp_dir/out.paf" "PAF output"
paf_placement "$meta_temp_dir/out.paf" "$meta_temp_dir/out.paf.pos"
check_placement "$meta_temp_dir/out.paf.pos" "PAF"
log "OK: TEST 1 passed"

# --- TEST 2: sorted + indexed BAM output ---
log "TEST 2: alignment to sorted BAM"
"$meta_executable" \
  --reference "$test_dir/ref.fasta" \
  --query "$test_dir/reads.fastq" \
  --bam \
  --output "$meta_temp_dir/out.bam"

check_is_bam "$meta_temp_dir/out.bam" "BAM output"
check_file_not_empty "$meta_temp_dir/out.bam.bai" "BAM index (default location)"
check_file_contains <(samtools view -H "$meta_temp_dir/out.bam") "SO:coordinate" "BAM header (sort order)"
t2_primary=$(count_primary "$meta_temp_dir/out.bam")
if [[ "$t2_primary" -ne "$n_reads" ]]; then
  log_error "BAM has $t2_primary primary alignments, expected $n_reads"
  exit 1
fi
bam_placement "$meta_temp_dir/out.bam" "$meta_temp_dir/out.bam.pos"
check_placement "$meta_temp_dir/out.bam.pos" "BAM"
# The index must actually serve region queries, not just exist
t2_region=$(samtools view -c "$meta_temp_dir/out.bam" seq1:4000-6000)
if [[ "$t2_region" -ne 1 ]]; then
  log_error "Region query seq1:4000-6000 returned $t2_region records, expected 1"
  exit 1
fi
log "- BAM index serves region queries"
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
  --output_index "$meta_temp_dir/custom_index.bai"

check_is_bam "$meta_temp_dir/named.bam" "BAM output"
check_file_not_empty "$meta_temp_dir/custom_index.bai" "explicitly named BAM index"
check_file_not_exists "$meta_temp_dir/named.bam.bai" "default-location index (should not be written)"
log "OK: TEST 3 passed"

# --- TEST 4: --bam writes BAM even when --output ends in .sam ---
# Regression guard: samtools sort infers its output format from the file
# extension, so without an explicit -O bam this run wrote plain SAM and
# samtools index then failed.
log "TEST 4: --bam with a .sam output path"
"$meta_executable" \
  --reference "$test_dir/ref.fasta" \
  --query "$test_dir/reads.fastq" \
  --bam \
  --output "$meta_temp_dir/out.sam"

check_is_bam "$meta_temp_dir/out.sam" "BAM output (.sam extension)"
check_file_not_empty "$meta_temp_dir/out.sam.bai" "BAM index (.sam extension)"
log "OK: TEST 4 passed"

# --- TEST 5: --bam under a cpu/memory allocation ---
# ---cpus / ---memory are what populate meta_cpus / meta_memory_mb, which is
# the only way the -t / -@ / -m branches of the script run at all.
log "TEST 5: --bam under ---cpus 2 ---memory 2gb"
"$meta_executable" \
  --reference "$test_dir/ref.fasta" \
  --query "$test_dir/reads.fastq" \
  --bam \
  --output "$meta_temp_dir/threads.bam" \
  ---cpus 2 \
  ---memory 2gb

check_is_bam "$meta_temp_dir/threads.bam" "BAM output (with cpu/memory)"
check_file_not_empty "$meta_temp_dir/threads.bam.bai" "BAM index (with cpu/memory)"
if ! diff -q <(alignment_key "$meta_temp_dir/out.bam") \
             <(alignment_key "$meta_temp_dir/threads.bam") > /dev/null; then
  log_error "Allocating cpus changed the alignments"
  exit 1
fi
log "- same alignments as the unconstrained run"
# Same alignments alone would also hold if the threads were never passed on,
# so check the @PG command lines that minimap2 and samtools sort record.
# 2gb = 2000 MB; sort gets half of it over 1 + 2 threads -> -m 333M.
samtools view -H "$meta_temp_dir/threads.bam" | grep "^@PG" > "$meta_temp_dir/threads.pg"
grep "ID:samtools" "$meta_temp_dir/threads.pg" > "$meta_temp_dir/threads_sort.pg"
check_file_contains "$meta_temp_dir/threads.pg" "CL:minimap2 -t 2" "@PG header (minimap2 -t)"
check_file_contains "$meta_temp_dir/threads_sort.pg" " -@ 2 " "@PG header (samtools sort -@)"
check_file_contains "$meta_temp_dir/threads_sort.pg" " -m 333M " "@PG header (samtools sort -m)"

# A small allocation must not drive -m below the 128M floor.
"$meta_executable" \
  --reference "$test_dir/ref.fasta" \
  --query "$test_dir/reads.fastq" \
  --bam \
  --output "$meta_temp_dir/lowmem.bam" \
  ---cpus 2 \
  ---memory 512mb
check_is_bam "$meta_temp_dir/lowmem.bam" "BAM output (low memory)"
samtools view -H "$meta_temp_dir/lowmem.bam" | grep "ID:samtools" > "$meta_temp_dir/lowmem_sort.pg"
check_file_contains "$meta_temp_dir/lowmem_sort.pg" " -m 128M " "@PG header (samtools sort -m floor)"
log "OK: TEST 5 passed"

# --- TEST 6: presets the old engine image could not run ---
# map-hifi (added in minimap2 2.19) and lr:hq were both rejected with "unknown
# preset" by the minimap2 2.17 the component used to ship.
for preset in map-hifi lr:hq; do
  log "TEST 6: $preset preset"
  "$meta_executable" \
    --reference "$test_dir/ref.fasta" \
    --query "$test_dir/reads.fastq" \
    --preset "$preset" \
    --cigar_paf \
    --output "$meta_temp_dir/preset_${preset/:/_}.paf"
  paf="$meta_temp_dir/preset_${preset/:/_}.paf"
  check_file_not_empty "$paf" "$preset PAF output"
  # every alignment line carries a CIGAR, not just one of them
  n_cigar=$(grep -c "cg:Z:" "$paf")
  if [[ "$n_cigar" -ne "$(wc -l < "$paf")" ]]; then
    log_error "$preset: $n_cigar of $(wc -l < "$paf") PAF lines carry cg:Z: (--cigar_paf)"
    exit 1
  fi
  paf_placement "$paf" "$paf.pos"
  check_placement "$paf.pos" "$preset"
done
log "OK: TEST 6 passed"

# --- TEST 7: pre-built index built for the same preset ---
# The preset's indexing parameters (-k/-w/-H) are baked into a .mmi. An index
# built for the same preset must reproduce the FASTA run exactly and without
# minimap2's override warning; minimap2_index --preset exists to produce it.
# The index is built with minimap2 directly because this test cannot call the
# sibling component.
log "TEST 7: map-hifi on a map-hifi .mmi"
minimap2 -x map-hifi -d "$test_dir/ref_hifi.mmi" "$test_dir/ref.fasta" 2> /dev/null
"$meta_executable" \
  --reference "$test_dir/ref_hifi.mmi" \
  --query "$test_dir/reads.fastq" \
  --preset map-hifi \
  --cigar_paf \
  --output "$meta_temp_dir/mmi_hifi.paf" \
  > "$meta_temp_dir/mmi_hifi.log" 2>&1

check_file_not_contains "$meta_temp_dir/mmi_hifi.log" "overridden" "matching .mmi log (override warning)"
if ! diff -q <(sort "$meta_temp_dir/preset_map-hifi.paf") \
             <(sort "$meta_temp_dir/mmi_hifi.paf") > /dev/null; then
  log_error "map-hifi .mmi produced different alignments than the FASTA reference"
  exit 1
fi
log "- identical to the FASTA reference run"

# And the trap the matching index avoids: a default (k=15, w=10) .mmi silently
# keeps its own parameters when aligned with map-hifi.
minimap2 -d "$test_dir/ref_default.mmi" "$test_dir/ref.fasta" 2> /dev/null
"$meta_executable" \
  --reference "$test_dir/ref_default.mmi" \
  --query "$test_dir/reads.fastq" \
  --preset map-hifi \
  --output "$meta_temp_dir/mmi_default.paf" \
  > "$meta_temp_dir/mmi_default.log" 2>&1
check_file_contains "$meta_temp_dir/mmi_default.log" "overridden by parameters used in the prebuilt index" "mismatched .mmi log (override warning)"
log "OK: TEST 7 passed"

# --- TEST 8: contradictory flags are rejected ---
log "TEST 8: mutually exclusive flags are rejected"
common=(--reference "$test_dir/ref.fasta" --query "$test_dir/reads.fastq")
check_rejected "$meta_temp_dir/bad1.log" "cigar_bam is only valid together with" \
  "cigar_bam without bam" "${common[@]}" --cigar_bam --output "$meta_temp_dir/bad1.paf"
check_rejected "$meta_temp_dir/bad2.log" "cigar_paf is only valid without" \
  "cigar_paf with bam" "${common[@]}" --bam --cigar_paf --output "$meta_temp_dir/bad2.bam"
log "OK: TEST 8 passed"

# --- TEST 9: --output_index without --bam is ignored, not rejected ---
# The Nextflow runner fills in a default path for every output file argument,
# so a PAF run there always arrives with --output_index set; rejecting it made
# every PAF run fail under Nextflow. Under the executable runner viash itself
# still reports the index as a missing output afterwards, so only the script's
# own behaviour is checked here: it must align, not reject.
log "TEST 9: --output_index without --bam still produces PAF"
"$meta_executable" \
  "${common[@]}" \
  --output "$meta_temp_dir/paf_with_index.paf" \
  --output_index "$meta_temp_dir/paf_with_index.bai" \
  > "$meta_temp_dir/paf_with_index.log" 2>&1 || true

check_file_not_contains "$meta_temp_dir/paf_with_index.log" "only valid together with" "PAF with --output_index log (rejection)"
check_file_not_empty "$meta_temp_dir/paf_with_index.paf" "PAF output"
paf_placement "$meta_temp_dir/paf_with_index.paf" "$meta_temp_dir/paf_with_index.pos"
check_placement "$meta_temp_dir/paf_with_index.pos" "PAF with --output_index"
if [[ -e "$meta_temp_dir/paf_with_index.bai" ]]; then
  log_error "--output_index was written although no BAM was produced"
  exit 1
fi
log "OK: TEST 9 passed"

print_test_summary "$meta_name tests passed"
