#!/bin/bash

## VIASH START
## VIASH END

# Source the centralized test helpers
source "$meta_resources_dir/test_helpers.sh"

# Initialize test environment with strict error handling (set -euo pipefail)
setup_test_env

log "Starting tests for $meta_name"

# Test fixtures are pre-generated and shipped in test_data/ (see
# test_data/script.sh for how they are built): one unique contig plus one
# contig carrying two tandem-repeat arrays of differing copy number, and reads
# sampled FROM that reference so they actually align. The shared
# create_test_fasta/create_test_fastq helpers cannot be used here - their pure
# ATCGATCG... repeat collapses to a single distinct minimizer, so meryl reports
# "no. of kmers read=0" and the weighted-minimizer path Winnowmap exists for is
# never exercised.
test_data_dir="$meta_resources_dir/test_data"
read_length=2000

# SAM/BAM is read with samtools throughout, never with grep/awk over the raw
# file: samtools skips the header itself, so no check can be fooled by an @SQ
# line, and the same helper works on the SAM and the BAM runs alike.
# Alignment records (headers excluded by samtools view)
count_records() { samtools view -c "$1"; }
# Mapped records: -F 4 drops the unmapped flag
count_mapped() { samtools view -c -F 4 "$1"; }
# Distinct reads that mapped. Repeat-array reads pick up supplementary
# alignments, so record count exceeds read count - count QNAMEs instead.
count_mapped_reads() { samtools view -F 4 "$1" | cut -f1 | sort -u | wc -l; }
# Records whose RNAME is a given contig, via an htslib filter expression
count_on_contig() { samtools view -c -e "rname==\"$2\"" "$1"; }
# QNAME/FLAG/RNAME/POS/MAPQ/CIGAR of every record, sorted: the comparison key
# used to assert two runs produced the same alignments.
alignment_key() { samtools view "$1" | cut -f1-6 | sort; }

# Assert at least one alignment landed on a contig
check_contig_aligned() {
  local sam="$1" contig="$2" label="$3" n
  n=$(count_on_contig "$sam" "$contig")
  if [[ "$n" -lt 1 ]]; then
    log_error "$label: no alignment records with RNAME $contig"
    exit 1
  fi
  log "$label: $n records on $contig"
}

check_file_exists "$test_data_dir/reference.fasta" "test reference genome"
check_file_not_empty "$test_data_dir/reference.fasta" "test reference genome"
check_file_contains "$test_data_dir/reference.fasta" "chr_repeat" "reference (repeat contig)"
for f in reads.fastq reads.fastq.gz reads.fasta; do
  check_file_not_empty "$test_data_dir/$f" "test query reads ($f)"
done

n_reads=$(( $(wc -l < "$test_data_dir/reads.fastq") / 4 ))
log "Using $n_reads pre-generated test reads of ${read_length}bp"

##############################################################
log "Starting TEST 1: basic SAM output with automatic meryl k-mer computation"
##############################################################
"$meta_executable" \
  --reference "$test_data_dir/reference.fasta" \
  --query "$test_data_dir/reads.fastq" \
  --preset map-ont \
  --kmer_size 15 \
  --output "$meta_temp_dir/test1.sam"

check_file_exists "$meta_temp_dir/test1.sam" "SAM output"
check_file_not_empty "$meta_temp_dir/test1.sam" "SAM output"
check_file_contains <(samtools view -H "$meta_temp_dir/test1.sam") "@SQ" "SAM header (sequence dictionary)"

# A SAM containing only headers would pass a naive "has @ lines" check, so
# assert real alignment records and that the reads actually mapped.
t1_records=$(count_records "$meta_temp_dir/test1.sam")
t1_mapped_reads=$(count_mapped_reads "$meta_temp_dir/test1.sam")
log "TEST 1: $t1_records alignment records, $t1_mapped_reads of $n_reads reads mapped"
if [[ "$t1_records" -lt 1 ]]; then
  log_error "SAM output contains headers but no alignment records"
  exit 1
fi
if [[ "$t1_mapped_reads" -lt "$n_reads" ]]; then
  log_error "Only $t1_mapped_reads of $n_reads reads mapped; all reads were sampled from the reference so all should map"
  exit 1
fi
# Both contigs must carry alignments, asserted on RNAME rather than on the
# contig name appearing anywhere in the file (which the @SQ lines guarantee).
check_contig_aligned "$meta_temp_dir/test1.sam" "chr_unique" "TEST 1"
check_contig_aligned "$meta_temp_dir/test1.sam" "chr_repeat" "TEST 1"
log "TEST 1 completed successfully"

##############################################################
log "Starting TEST 2: pre-computed repetitive k-mers give identical alignments"
##############################################################
meryl_db="$meta_temp_dir/merylDB"
repetitive_kmers="$meta_temp_dir/repetitive_k15.txt"
meryl count k=15 output "$meryl_db" "$test_data_dir/reference.fasta"
meryl print greater-than distinct=0.9998 "$meryl_db" > "$repetitive_kmers"

check_file_exists "$repetitive_kmers" "meryl repetitive k-mers"
# The tandem-repeat arrays must make this non-empty, otherwise the weighted
# minimizer path is not being tested at all.
check_file_not_empty "$repetitive_kmers" "meryl repetitive k-mers"
log "TEST 2: meryl reported $(wc -l < "$repetitive_kmers") repetitive k-mers"

"$meta_executable" \
  --reference "$test_data_dir/reference.fasta" \
  --query "$test_data_dir/reads.fastq" \
  --repetitive_kmers "$repetitive_kmers" \
  --preset map-ont \
  --output "$meta_temp_dir/test2.sam"

check_file_exists "$meta_temp_dir/test2.sam" "SAM output (pre-computed k-mers)"
check_file_not_empty "$meta_temp_dir/test2.sam" "SAM output (pre-computed k-mers)"

# Supplying the same k-mer set the component would have computed itself must
# yield the same alignments; this is what proves -W is actually wired up.
alignment_key "$meta_temp_dir/test1.sam" > "$meta_temp_dir/t1.key"
alignment_key "$meta_temp_dir/test2.sam" > "$meta_temp_dir/t2.key"
if ! diff -q "$meta_temp_dir/t1.key" "$meta_temp_dir/t2.key" > /dev/null; then
  log_error "Pre-computed k-mers produced different alignments than auto-computed k-mers"
  diff "$meta_temp_dir/t1.key" "$meta_temp_dir/t2.key" | head -10
  exit 1
fi
log "TEST 2 completed successfully"

##############################################################
log "Starting TEST 3: sorted BAM output with index"
##############################################################
# Run unconditionally: samtools is bundled in the component's engine image, so
# skipping this test on a missing samtools would hide a broken --bam feature.
"$meta_executable" \
  --reference "$test_data_dir/reference.fasta" \
  --query "$test_data_dir/reads.fastq" \
  --preset map-ont \
  --bam \
  --kmer_size 15 \
  --output "$meta_temp_dir/test3.bam"

check_file_exists "$meta_temp_dir/test3.bam" "BAM output"
check_file_not_empty "$meta_temp_dir/test3.bam" "BAM output"
check_file_exists "$meta_temp_dir/test3.bam.bai" "BAM index (default location)"

if ! samtools quickcheck "$meta_temp_dir/test3.bam"; then
  log_error "BAM output failed samtools quickcheck"
  exit 1
fi

# The BAM must hold the same alignments as the SAM run, just sorted.
t3_records=$(samtools view -c "$meta_temp_dir/test3.bam")
if [[ "$t3_records" -ne "$t1_records" ]]; then
  log_error "BAM has $t3_records records but SAM had $t1_records"
  exit 1
fi
# Coordinate-sorted output must declare so in the header
check_file_contains <(samtools view -H "$meta_temp_dir/test3.bam") "SO:coordinate" "BAM header (sort order)"
log "TEST 3 completed successfully"

##############################################################
log "Starting TEST 3b: --bam writes BAM even when --output ends in .sam"
##############################################################
# Regression guard: samtools sort infers its output format from the file
# extension, so without an explicit -O bam this run wrote plain SAM and
# samtools index then failed with "not a BGZF file". The Nextflow runner
# derives the output name from the config `example:` (alignment.sam), so this
# is the default path there, not an edge case.
"$meta_executable" \
  --reference "$test_data_dir/reference.fasta" \
  --query "$test_data_dir/reads.fastq" \
  --preset map-ont \
  --bam \
  --kmer_size 15 \
  --output "$meta_temp_dir/test3b.sam" \
  --output_index "$meta_temp_dir/test3b.sam.bai"

check_file_not_empty "$meta_temp_dir/test3b.sam" "BAM output (.sam extension)"
check_file_not_empty "$meta_temp_dir/test3b.sam.bai" "BAM index (.sam extension)"
if ! samtools quickcheck "$meta_temp_dir/test3b.sam"; then
  log_error "--bam with a .sam output path did not produce a valid BAM"
  exit 1
fi
# quickcheck accepts plain SAM too, so ask htslib what the file actually is:
# "BAM version 1 compressed sequence data" vs "SAM version 1.6 sequence text".
# (`file` is not installed in the engine image; `htsfile` ships with samtools.)
t3b_format=$(htsfile "$meta_temp_dir/test3b.sam")
if [[ "$t3b_format" != *"BAM version"* ]]; then
  log_error "Output is not BAM; samtools sort fell back to SAM: $t3b_format"
  exit 1
fi
log "TEST 3b completed successfully"

##############################################################
log "Starting TEST 3c: --bam under a cpu/memory allocation"
##############################################################
# ---cpus / ---memory are what populate meta_cpus / meta_memory_mb, which is
# the only way the `samtools sort -m <per-thread>` branch in the script runs at
# all. Without this the whole memory-splitting block is dead code in the tests.
# The assertion is that the run works, not what -m was set to.
"$meta_executable" \
  --reference "$test_data_dir/reference.fasta" \
  --query "$test_data_dir/reads.fastq" \
  --preset map-ont \
  --bam \
  --kmer_size 15 \
  --output "$meta_temp_dir/test3c.bam" \
  ---cpus 2 \
  ---memory 4gb

check_file_not_empty "$meta_temp_dir/test3c.bam" "BAM output (with cpu/memory)"
check_file_not_empty "$meta_temp_dir/test3c.bam.bai" "BAM index (with cpu/memory)"
if ! samtools quickcheck "$meta_temp_dir/test3c.bam"; then
  log_error "BAM produced under a cpu/memory allocation failed samtools quickcheck"
  exit 1
fi
# Same alignments as the unconstrained BAM run: -t/-@/-m must not change output
if ! diff -q <(alignment_key "$meta_temp_dir/test3.bam") \
             <(alignment_key "$meta_temp_dir/test3c.bam") > /dev/null; then
  log_error "Allocating cpus/memory changed the alignments"
  exit 1
fi
log "TEST 3c completed successfully"

##############################################################
log "Starting TEST 4: BAM index written to an explicit --output_index path"
##############################################################
"$meta_executable" \
  --reference "$test_data_dir/reference.fasta" \
  --query "$test_data_dir/reads.fastq" \
  --preset map-ont \
  --bam \
  --kmer_size 15 \
  --output "$meta_temp_dir/test4.bam" \
  --output_index "$meta_temp_dir/test4_custom.bam.bai"

check_file_exists "$meta_temp_dir/test4.bam" "BAM output (explicit index path)"
check_file_exists "$meta_temp_dir/test4_custom.bam.bai" "BAM index (explicit path)"
check_file_not_empty "$meta_temp_dir/test4_custom.bam.bai" "BAM index (explicit path)"
log "TEST 4 completed successfully"

##############################################################
log "Starting TEST 4b: --output_index without --bam is ignored, not an error"
##############################################################
# The Nextflow runner auto-fills every declared output argument, so
# --output_index is always set there even in the default SAM mode. Rejecting
# that combination made every default-mode Nextflow run fail. `must_exist: false`
# on the argument is what keeps viash's own post-run output check quiet when no
# index is produced.
"$meta_executable" \
  --reference "$test_data_dir/reference.fasta" \
  --query "$test_data_dir/reads.fastq" \
  --preset map-ont \
  --kmer_size 15 \
  --output "$meta_temp_dir/test4b.sam" \
  --output_index "$meta_temp_dir/test4b.sam.bai"

check_file_not_empty "$meta_temp_dir/test4b.sam" "SAM output (--output_index given without --bam)"
if [[ "$(count_mapped_reads "$meta_temp_dir/test4b.sam")" -lt 1 ]]; then
  log_error "SAM run with --output_index mapped no reads"
  exit 1
fi
log "TEST 4b completed successfully"

##############################################################
log "Starting TEST 5: custom meryl k-mer size"
##############################################################
# Regression guard: winnowmap aborts with "input list of k-mers and winnowmap
# parameter k are inconsistent" unless --kmer_size drives winnowmap's -k as
# well as the meryl count step.
"$meta_executable" \
  --reference "$test_data_dir/reference.fasta" \
  --query "$test_data_dir/reads.fastq" \
  --preset map-ont \
  --kmer_size 19 \
  --output "$meta_temp_dir/test5.sam"

check_file_exists "$meta_temp_dir/test5.sam" "SAM output (kmer_size 19)"
check_file_not_empty "$meta_temp_dir/test5.sam" "SAM output (kmer_size 19)"
if [[ "$(count_records "$meta_temp_dir/test5.sam")" -lt 1 ]]; then
  log_error "kmer_size 19 run produced no alignment records"
  exit 1
fi
if [[ "$(count_mapped_reads "$meta_temp_dir/test5.sam")" -lt 1 ]]; then
  log_error "kmer_size 19 run mapped no reads"
  exit 1
fi
log "TEST 5 completed successfully"

##############################################################
log "Starting TEST 5b: pre-computed k=19 k-mer list with matching --kmer_size"
##############################################################
meryl_db19="$meta_temp_dir/merylDB19"
repetitive_k19="$meta_temp_dir/repetitive_k19.txt"
meryl count k=19 output "$meryl_db19" "$test_data_dir/reference.fasta"
meryl print greater-than distinct=0.9998 "$meryl_db19" > "$repetitive_k19"
check_file_not_empty "$repetitive_k19" "meryl repetitive k-mers (k=19)"

"$meta_executable" \
  --reference "$test_data_dir/reference.fasta" \
  --query "$test_data_dir/reads.fastq" \
  --repetitive_kmers "$repetitive_k19" \
  --kmer_size 19 \
  --preset map-ont \
  --output "$meta_temp_dir/test5b.sam"

check_file_not_empty "$meta_temp_dir/test5b.sam" "SAM output (pre-computed k=19)"
if [[ "$(count_mapped_reads "$meta_temp_dir/test5b.sam")" -lt 1 ]]; then
  log_error "pre-computed k=19 list mapped no reads"
  exit 1
fi
log "TEST 5b completed successfully"

##############################################################
log "Starting TEST 5c: mismatched k must fail, not mis-align silently"
##############################################################
# A k=19 list with the default k=15 is exactly the inconsistency winnowmap
# guards against; the component must surface that rather than swallow it.
if "$meta_executable" \
     --reference "$test_data_dir/reference.fasta" \
     --query "$test_data_dir/reads.fastq" \
     --repetitive_kmers "$repetitive_k19" \
     --kmer_size 15 \
     --preset map-ont \
     --output "$meta_temp_dir/test5c.sam" > "$meta_temp_dir/test5c.log" 2>&1; then
  log_error "A k=19 k-mer list was accepted with --kmer_size 15; expected failure"
  exit 1
fi
# A negative test that does not look at the message passes on any failure at
# all, so pin it to winnowmap's own inconsistency check.
check_file_contains "$meta_temp_dir/test5c.log" "inconsistent" "error message (k mismatch)"
log "TEST 5c completed successfully (k mismatch rejected)"

##############################################################
log "Starting TEST 5d: --kmer_size above winnowmap's ceiling is rejected by the config"
##############################################################
# max: 28 in the config turns this into an argument-parsing failure, instead of
# a full meryl count over the reference followed by a winnowmap abort.
if "$meta_executable" \
     --reference "$test_data_dir/reference.fasta" \
     --query "$test_data_dir/reads.fastq" \
     --kmer_size 30 \
     --output "$meta_temp_dir/test5d.sam" > "$meta_temp_dir/test5d.log" 2>&1; then
  log_error "--kmer_size 30 was accepted; winnowmap allows at most 28"
  exit 1
fi
# Must be the config's `max:` talking, before meryl ever starts
check_file_contains "$meta_temp_dir/test5d.log" "less than or equal to 28" "error message (kmer_size ceiling)"
log "TEST 5d completed successfully (out-of-range kmer_size rejected)"

##############################################################
log "Starting TEST 5e: meryl path without --kmer_size must fail"
##############################################################
# No k is defaulted in the script: meryl has to be told which k to count, and
# hardcoding winnowmap's current default here would drift the day upstream
# changes it. So the component refuses instead of guessing.
if "$meta_executable" \
     --reference "$test_data_dir/reference.fasta" \
     --query "$test_data_dir/reads.fastq" \
     --preset map-ont \
     --output "$meta_temp_dir/test5e.sam" > "$meta_temp_dir/test5e.log" 2>&1; then
  log_error "--kmer_size was omitted while meryl had to build the list; expected failure"
  exit 1
fi
# Leading "--" is dropped from the pattern on purpose: check_file_contains
# passes it straight to grep, which would read it as an option.
check_file_contains "$meta_temp_dir/test5e.log" "kmer_size is required" "error message (missing kmer_size)"
log "TEST 5e completed successfully (missing kmer_size rejected)"

##############################################################
log "Starting TEST 6: gzip-compressed FASTQ query"
##############################################################
"$meta_executable" \
  --reference "$test_data_dir/reference.fasta" \
  --query "$test_data_dir/reads.fastq.gz" \
  --preset map-ont \
  --kmer_size 15 \
  --output "$meta_temp_dir/test6.sam"

check_file_exists "$meta_temp_dir/test6.sam" "SAM output (gzipped query)"
# Compressed input must give the same alignments as the plain-text input
alignment_key "$meta_temp_dir/test6.sam" > "$meta_temp_dir/t6.key"
if ! diff -q "$meta_temp_dir/t1.key" "$meta_temp_dir/t6.key" > /dev/null; then
  log_error "Gzipped query produced different alignments than uncompressed query"
  exit 1
fi
log "TEST 6 completed successfully"

##############################################################
log "Starting TEST 7: FASTA query input"
##############################################################
"$meta_executable" \
  --reference "$test_data_dir/reference.fasta" \
  --query "$test_data_dir/reads.fasta" \
  --preset map-ont \
  --kmer_size 15 \
  --output "$meta_temp_dir/test7.sam"

check_file_exists "$meta_temp_dir/test7.sam" "SAM output (FASTA query)"
if [[ "$(count_mapped "$meta_temp_dir/test7.sam")" -lt 1 ]]; then
  log_error "FASTA query produced no mapped reads"
  exit 1
fi
log "TEST 7 completed successfully"

##############################################################
log "Starting TEST 8: map-pb preset (PacBio HiFi)"
##############################################################
# Winnowmap's own help documents map-pb as "hifi-to-ref"; there is no map-hifi.
"$meta_executable" \
  --reference "$test_data_dir/reference.fasta" \
  --query "$test_data_dir/reads.fastq" \
  --preset map-pb \
  --kmer_size 15 \
  --output "$meta_temp_dir/test8.sam"

check_file_exists "$meta_temp_dir/test8.sam" "SAM output (map-pb)"
if [[ "$(count_mapped "$meta_temp_dir/test8.sam")" -lt 1 ]]; then
  log_error "map-pb preset produced no mapped reads"
  exit 1
fi
log "TEST 8 completed successfully"

##############################################################
log "Starting TEST 8b: asm5 preset, whose own k must not override --kmer_size"
##############################################################
# asm5/asm10/asm20 set k=19 themselves, so they are the presets that can
# desynchronise winnowmap's k from the k the -W list was built with
# ("input list of k-mers and winnowmap parameter k are inconsistent"). No asm
# preset was covered before. Verified with winnowmap 2.03: an explicit -k wins
# over the preset regardless of where -x sits on the command line.
"$meta_executable" \
  --reference "$test_data_dir/reference.fasta" \
  --query "$test_data_dir/reads.fasta" \
  --preset asm5 \
  --kmer_size 15 \
  --output "$meta_temp_dir/test8b.sam"

check_file_not_empty "$meta_temp_dir/test8b.sam" "SAM output (asm5)"
if [[ "$(count_mapped "$meta_temp_dir/test8b.sam")" -lt 1 ]]; then
  log_error "asm5 preset produced no mapped reads"
  exit 1
fi

# And the same preset with an explicit k: the user's value must win, both for
# meryl and for winnowmap.
"$meta_executable" \
  --reference "$test_data_dir/reference.fasta" \
  --query "$test_data_dir/reads.fasta" \
  --preset asm5 \
  --kmer_size 19 \
  --output "$meta_temp_dir/test8c.sam"

check_file_not_empty "$meta_temp_dir/test8c.sam" "SAM output (asm5, kmer_size 19)"
if [[ "$(count_mapped "$meta_temp_dir/test8c.sam")" -lt 1 ]]; then
  log_error "asm5 preset with --kmer_size 19 produced no mapped reads"
  exit 1
fi
log "TEST 8b completed successfully"

##############################################################
log "Starting TEST 9: no preset at all (winnowmap defaults)"
##############################################################
# --preset has no default, matching winnowmap's own `-x STR ... []`, so this
# exercises the no -x code path rather than repeating TEST 1.
"$meta_executable" \
  --reference "$test_data_dir/reference.fasta" \
  --query "$test_data_dir/reads.fastq" \
  --kmer_size 15 \
  --output "$meta_temp_dir/test9.sam"

check_file_exists "$meta_temp_dir/test9.sam" "SAM output (no preset)"
if [[ "$(count_mapped "$meta_temp_dir/test9.sam")" -lt 1 ]]; then
  log_error "run without a preset produced no mapped reads"
  exit 1
fi
log "TEST 9 completed successfully"

##############################################################
log "Starting TEST 10: an unsupported preset must fail loudly"
##############################################################
# `choices:` in the config rejects this during argument parsing, before the
# meryl count step; guard with "if !" because setup_test_env enabled `set -e`.
if "$meta_executable" \
     --reference "$test_data_dir/reference.fasta" \
     --query "$test_data_dir/reads.fastq" \
     --preset not-a-real-preset \
     --kmer_size 15 \
     --output "$meta_temp_dir/test10.sam" > "$meta_temp_dir/test10.log" 2>&1; then
  log_error "An invalid preset was accepted; the component should have failed"
  exit 1
fi
# Must be the config's `choices:` talking, not a downstream winnowmap error
check_file_contains "$meta_temp_dir/test10.log" "not in the list of allowed values" "error message (invalid preset)"
log "TEST 10 completed successfully (invalid preset rejected)"

##############################################################
log "Starting TEST 11: pass-through of the wider winnowmap option surface"
##############################################################
"$meta_executable" \
  --reference "$test_data_dir/reference.fasta" \
  --query "$test_data_dir/reads.fastq" \
  --preset map-ont \
  --read_group '@RG\tID:rg1\tSM:sample1' \
  --md_tag \
  --eqx \
  --cigar_bam \
  --soft_clipping \
  --window_size 30 \
  --secondary_ratio 0.5 \
  --min_chaining_score 30 \
  --kmer_size 15 \
  --output "$meta_temp_dir/test11.sam"

check_file_not_empty "$meta_temp_dir/test11.sam" "SAM output (extra options)"
# -R lands in the header, and every record must carry the RG tag
check_file_contains <(samtools view -H "$meta_temp_dir/test11.sam") "@RG" "SAM header (read group)"
check_file_contains <(samtools view -H "$meta_temp_dir/test11.sam") "SM:sample1" "SAM header (read group sample)"
if [[ "$(samtools view -c -e '[RG]=="rg1"' "$meta_temp_dir/test11.sam")" -lt 1 ]]; then
  log_error "--read_group did not tag any alignment records"
  exit 1
fi
# --MD adds MD:Z:, --eqx replaces M operators with =/X in the CIGAR
if [[ "$(samtools view -c -e '[MD]' "$meta_temp_dir/test11.sam")" -lt 1 ]]; then
  log_error "--md_tag produced no MD tags"
  exit 1
fi
if [[ "$(samtools view -c -e 'cigar =~ "="' "$meta_temp_dir/test11.sam")" -lt 1 ]]; then
  log_error "--eqx produced no =/X CIGAR operators"
  exit 1
fi
log "TEST 11 completed successfully"

##############################################################
log "Starting TEST 11b: --cs_tag"
##############################################################
"$meta_executable" \
  --reference "$test_data_dir/reference.fasta" \
  --query "$test_data_dir/reads.fastq" \
  --preset map-ont \
  --cs_tag short \
  --kmer_size 15 \
  --output "$meta_temp_dir/test11b.sam"

if [[ "$(samtools view -c -e '[cs]' "$meta_temp_dir/test11b.sam")" -lt 1 ]]; then
  log_error "--cs_tag short produced no cs tags"
  exit 1
fi
log "TEST 11b completed successfully"

##############################################################
log "Starting TEST 11c: the meryl scratch directory is cleaned up"
##############################################################
# The scratch directory holds the meryl database, about 3x the size of the
# reference, and it is created under the component's meta_temp_dir: a shared
# root on the host that nothing else clears, written as root under the docker
# engine. A run that leaves it behind fills the host with directories the
# caller cannot delete.
#
# meta_temp_dir of the component is VIASH_TEMP, not this test's meta_temp_dir
# (the -W paths in the logs above are /tmp/<name>_XXXXXX inside the component's
# own container). So point VIASH_TEMP at a directory this test can inspect.
scratch_root="$meta_temp_dir/scratch_root"
mkdir -p "$scratch_root"
count_scratch() { find "$scratch_root" -maxdepth 1 -type d -name "${meta_name}_*" | wc -l; }

VIASH_TEMP="$scratch_root" "$meta_executable" \
  --reference "$test_data_dir/reference.fasta" \
  --query "$test_data_dir/reads.fastq" \
  --preset map-ont \
  --kmer_size 15 \
  --output "$meta_temp_dir/test11c.sam"

check_file_not_empty "$meta_temp_dir/test11c.sam" "SAM output (scratch-root run)"
if [[ "$(count_scratch)" -ne 0 ]]; then
  log_error "scratch directory left behind after a successful run:"
  find "$scratch_root" -maxdepth 1 -type d -name "${meta_name}_*"
  exit 1
fi

# Same on the failure path: the scratch directory is created before winnowmap
# runs, so the k-mismatch abort below has to leave nothing behind either. (A
# meryl-step failure would be the other half of this, but meryl exits 0 on
# every malformed reference tried here, so there is nothing to trigger it with.)
if VIASH_TEMP="$scratch_root" "$meta_executable" \
     --reference "$test_data_dir/reference.fasta" \
     --query "$test_data_dir/reads.fastq" \
     --repetitive_kmers "$repetitive_k19" \
     --kmer_size 15 \
     --preset map-ont \
     --output "$meta_temp_dir/test11c_fail.sam" > /dev/null 2>&1; then
  log_error "expected the k mismatch to fail this run"
  exit 1
fi
if [[ "$(count_scratch)" -ne 0 ]]; then
  log_error "scratch directory left behind after a failed run:"
  find "$scratch_root" -maxdepth 1 -type d -name "${meta_name}_*"
  exit 1
fi
log "TEST 11c completed successfully"

##############################################################
log "Starting TEST 12: two concurrent runs sharing one temp root"
##############################################################
# The meryl database and k-mer list are built under meta_temp_dir, which is a
# shared root (VIASH_TEMP, /tmp inside the container) rather than a
# per-invocation directory - the -W paths in the logs above show it. With fixed
# names there, parallel runs on the same host write the same merylDB and the
# same k-mer list: one run reads what the other is still writing, giving
# silently wrong weighting or a mid-run abort.
#
# The two runs below use different k so their meryl databases genuinely differ;
# running them with the same reference and k would collide on byte-identical
# content and prove nothing.
"$meta_executable" \
  --reference "$test_data_dir/reference.fasta" \
  --query "$test_data_dir/reads.fastq" \
  --preset map-ont \
  --kmer_size 15 \
  --output "$meta_temp_dir/test12a.sam" > "$meta_temp_dir/test12a.log" 2>&1 &
pid_a=$!
"$meta_executable" \
  --reference "$test_data_dir/reference.fasta" \
  --query "$test_data_dir/reads.fastq" \
  --preset map-ont \
  --kmer_size 19 \
  --output "$meta_temp_dir/test12b.sam" > "$meta_temp_dir/test12b.log" 2>&1 &
pid_b=$!

status=0
wait "$pid_a" || status=1
wait "$pid_b" || status=1
if [[ "$status" -ne 0 ]]; then
  log_error "Concurrent runs failed; the meryl scratch files are colliding"
  cat "$meta_temp_dir/test12a.log" "$meta_temp_dir/test12b.log"
  exit 1
fi

# Each must reproduce exactly what the same arguments produced sequentially
for pair in "test12a:test1" "test12b:test5"; do
  conc="${pair%%:*}"
  seq="${pair##*:}"
  alignment_key "$meta_temp_dir/$conc.sam" > "$meta_temp_dir/$conc.key"
  alignment_key "$meta_temp_dir/$seq.sam" > "$meta_temp_dir/$seq.key"
  if ! diff -q "$meta_temp_dir/$conc.key" "$meta_temp_dir/$seq.key" > /dev/null; then
    log_error "Concurrent run $conc differs from the sequential run $seq"
    diff "$meta_temp_dir/$seq.key" "$meta_temp_dir/$conc.key" | head -10
    exit 1
  fi
done
log "TEST 12 completed successfully"

print_test_summary "All tests completed successfully"
