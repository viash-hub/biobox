#!/bin/bash

## VIASH START
## VIASH END

# Source the centralized test helpers
source "$meta_resources_dir/test_helpers.sh"

# Initialize test environment with strict error handling (set -euo pipefail)
setup_test_env

log "Starting tests for $meta_name"

test_data_dir="$meta_temp_dir/test_data"
mkdir -p "$test_data_dir"

##############################################################
# Test data generation
#
# The shared create_test_fasta helper emits a pure ATCGATCG... tandem repeat,
# which collapses to a single distinct minimizer: meryl then reports
# "no. of kmers read=0" and the weighted-minimizer code path that Winnowmap
# exists for is never exercised. create_test_fastq likewise emits reads that
# are unrelated to the reference, so nothing aligns.
#
# So build fixtures locally instead: one unique contig plus one contig
# carrying real tandem-repeat arrays, and reads sampled FROM that reference.
# awk with a fixed seed keeps this deterministic (no python3 in this container).
##############################################################

read_length=2000

# `meryl print greater-than distinct=0.9998` keeps only the top 0.02% of
# distinct k-mers by count, so the repeat structure has to be GRADED to produce
# a non-empty set. A single tandem array is not enough: all of its k-mers share
# one count, so the quantile threshold lands exactly on them and
# "greater-than" excludes every one. Two arrays with different copy numbers put
# the threshold at the lower count and let the higher-count k-mers through.
# Sizes are chosen so the top 0.02% (~38 k-mers of ~190k distinct) is wider
# than motif A's ~20 distinct k-mers.
generate_reference() {
  local out="$1"
  awk '
    function emit(c) { buf = buf c; if (length(buf) == 60) { print buf; buf = "" } }
    function flush_buf() { if (length(buf) > 0) { print buf; buf = "" } }
    function random_seq(n,   i) { for (i = 0; i < n; i++) emit(substr(b, int(rand() * 4) + 1, 1)) }
    function make_motif(n,   i, m) { m = ""; for (i = 0; i < n; i++) m = m substr(b, int(rand() * 4) + 1, 1); return m }
    function tandem(motif, copies,   r, i) {
      for (r = 0; r < copies; r++) for (i = 1; i <= length(motif); i++) emit(substr(motif, i, 1))
    }
    BEGIN {
      srand(42); b = "ACGT"; buf = ""
      # contig 1: unique background
      print ">chr_unique"
      random_seq(120000)
      flush_buf()

      # contig 2: unique flanks around two tandem arrays of differing copy number
      motif_a = make_motif(20)   # high copy  -> survives the 0.9998 cutoff
      motif_b = make_motif(31)   # lower copy -> sets the cutoff
      print ">chr_repeat"
      random_seq(20000)
      tandem(motif_a, 600)
      random_seq(10000)
      tandem(motif_b, 250)
      random_seq(20000)
      flush_buf()
    }' > "$out"
}

# Sample reads from the reference so they actually align. fmt = fastq | fasta
generate_reads() {
  local ref="$1" out="$2" fmt="$3" rl="$4"
  awk -v fmt="$fmt" -v rl="$rl" '
    BEGIN { srand(7); b = "ACGT" }
    /^>/ { if (name != "") seqs[name] = s; name = substr($0, 2); s = ""; next }
    { s = s $0 }
    END {
      seqs[name] = s
      q = ""
      for (i = 0; i < rl; i++) q = q "I"
      n = 0
      # deterministic contig order
      split("chr_unique chr_repeat", order, " ")
      for (o = 1; o <= 2; o++) {
        nm = order[o]
        L = length(seqs[nm])
        if (L <= rl) continue
        for (k = 0; k < 8; k++) {
          st = int(rand() * (L - rl)) + 1
          r = substr(seqs[nm], st, rl)
          # ~1% substitutions so alignment is not a trivial exact match
          out = ""
          for (i = 1; i <= rl; i++) {
            c = substr(r, i, 1)
            if (rand() < 0.01) c = substr(b, int(rand() * 4) + 1, 1)
            out = out c
          }
          n++
          if (fmt == "fastq") printf("@read%d_%s_%d\n%s\n+\n%s\n", n, nm, st, out, q)
          else printf(">read%d_%s_%d\n%s\n", n, nm, st, out)
        }
      }
    }' "$ref" > "$out"
}

# Alignment records only (drop @ header lines)
count_records() { grep -vc "^@" "$1" || true; }
# Mapped records: RNAME is not "*"
count_mapped() { awk '!/^@/ && $3 != "*"' "$1" | wc -l; }
# Distinct reads that mapped. Repeat-array reads pick up supplementary
# alignments, so record count exceeds read count - count QNAMEs instead.
count_mapped_reads() { awk '!/^@/ && $3 != "*" { seen[$1] = 1 } END { print length(seen) }' "$1"; }
# Records whose RNAME is a given contig. Grepping the file for the contig name
# instead would match the @SQ header line and pass even with nothing aligned.
count_on_contig() { awk -v c="$2" '!/^@/ && $3 == c' "$1" | wc -l; }

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

log "Generating test reference (unique contig + 20bp x600 and 31bp x250 tandem arrays)..."
generate_reference "$test_data_dir/reference.fasta"
check_file_exists "$test_data_dir/reference.fasta" "test reference genome"
check_file_not_empty "$test_data_dir/reference.fasta" "test reference genome"
check_file_contains "$test_data_dir/reference.fasta" "chr_repeat" "reference (repeat contig)"

log "Generating test reads sampled from the reference..."
generate_reads "$test_data_dir/reference.fasta" "$test_data_dir/reads.fastq" fastq "$read_length"
check_file_exists "$test_data_dir/reads.fastq" "test query reads"
check_file_not_empty "$test_data_dir/reads.fastq" "test query reads"

n_reads=$(( $(wc -l < "$test_data_dir/reads.fastq") / 4 ))
log "Generated $n_reads test reads of ${read_length}bp"

##############################################################
log "Starting TEST 1: basic SAM output with automatic meryl k-mer computation"
##############################################################
"$meta_executable" \
  --reference "$test_data_dir/reference.fasta" \
  --query "$test_data_dir/reads.fastq" \
  --preset map-ont \
  --output "$meta_temp_dir/test1.sam"

check_file_exists "$meta_temp_dir/test1.sam" "SAM output"
check_file_not_empty "$meta_temp_dir/test1.sam" "SAM output"
check_file_contains "$meta_temp_dir/test1.sam" "@SQ" "SAM output (sequence dictionary header)"

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
awk '!/^@/ {print $1, $2, $3, $4, $5, $6}' "$meta_temp_dir/test1.sam" | sort > "$meta_temp_dir/t1.key"
awk '!/^@/ {print $1, $2, $3, $4, $5, $6}' "$meta_temp_dir/test2.sam" | sort > "$meta_temp_dir/t2.key"
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
  --output "$meta_temp_dir/test3b.sam" \
  --output_index "$meta_temp_dir/test3b.sam.bai"

check_file_not_empty "$meta_temp_dir/test3b.sam" "BAM output (.sam extension)"
check_file_not_empty "$meta_temp_dir/test3b.sam.bai" "BAM index (.sam extension)"
if ! samtools quickcheck "$meta_temp_dir/test3b.sam"; then
  log_error "--bam with a .sam output path did not produce a valid BAM"
  exit 1
fi
# BGZF magic: a plain-text SAM starts with "@HD"/"@SQ" instead
if [[ "$(head -c 2 "$meta_temp_dir/test3b.sam" | od -An -tx1 | tr -d ' ')" != "1f8b" ]]; then
  log_error "Output is not BGZF-compressed; samtools sort fell back to SAM"
  exit 1
fi
log "TEST 3b completed successfully"

##############################################################
log "Starting TEST 4: BAM index written to an explicit --output_index path"
##############################################################
"$meta_executable" \
  --reference "$test_data_dir/reference.fasta" \
  --query "$test_data_dir/reads.fastq" \
  --preset map-ont \
  --bam \
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
     --output "$meta_temp_dir/test5c.sam" > /dev/null 2>&1; then
  log_error "A k=19 k-mer list was accepted with --kmer_size 15; expected failure"
  exit 1
fi
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
     --output "$meta_temp_dir/test5d.sam" > /dev/null 2>&1; then
  log_error "--kmer_size 30 was accepted; winnowmap allows at most 28"
  exit 1
fi
log "TEST 5d completed successfully (out-of-range kmer_size rejected)"

##############################################################
log "Starting TEST 6: gzip-compressed FASTQ query"
##############################################################
gzip -c "$test_data_dir/reads.fastq" > "$test_data_dir/reads.fastq.gz"
"$meta_executable" \
  --reference "$test_data_dir/reference.fasta" \
  --query "$test_data_dir/reads.fastq.gz" \
  --preset map-ont \
  --output "$meta_temp_dir/test6.sam"

check_file_exists "$meta_temp_dir/test6.sam" "SAM output (gzipped query)"
# Compressed input must give the same alignments as the plain-text input
awk '!/^@/ {print $1, $2, $3, $4, $5, $6}' "$meta_temp_dir/test6.sam" | sort > "$meta_temp_dir/t6.key"
if ! diff -q "$meta_temp_dir/t1.key" "$meta_temp_dir/t6.key" > /dev/null; then
  log_error "Gzipped query produced different alignments than uncompressed query"
  exit 1
fi
log "TEST 6 completed successfully"

##############################################################
log "Starting TEST 7: FASTA query input"
##############################################################
generate_reads "$test_data_dir/reference.fasta" "$test_data_dir/reads.fasta" fasta "$read_length"
"$meta_executable" \
  --reference "$test_data_dir/reference.fasta" \
  --query "$test_data_dir/reads.fasta" \
  --preset map-ont \
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
     --output "$meta_temp_dir/test10.sam" > /dev/null 2>&1; then
  log_error "An invalid preset was accepted; the component should have failed"
  exit 1
fi
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
  --output "$meta_temp_dir/test11.sam"

check_file_not_empty "$meta_temp_dir/test11.sam" "SAM output (extra options)"
# -R lands in the header, and every record must carry the RG tag
check_file_contains "$meta_temp_dir/test11.sam" "@RG" "SAM header (read group)"
check_file_contains "$meta_temp_dir/test11.sam" "SM:sample1" "SAM header (read group sample)"
if [[ "$(awk '!/^@/ && /RG:Z:rg1/' "$meta_temp_dir/test11.sam" | wc -l)" -lt 1 ]]; then
  log_error "--read_group did not tag any alignment records"
  exit 1
fi
# --MD adds MD:Z:, --eqx replaces M operators with =/X in the CIGAR
if [[ "$(awk '!/^@/ && /MD:Z:/' "$meta_temp_dir/test11.sam" | wc -l)" -lt 1 ]]; then
  log_error "--md_tag produced no MD tags"
  exit 1
fi
if [[ "$(awk '!/^@/ && $6 ~ /=/' "$meta_temp_dir/test11.sam" | wc -l)" -lt 1 ]]; then
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
  --output "$meta_temp_dir/test11b.sam"

if [[ "$(awk '!/^@/ && /cs:Z:/' "$meta_temp_dir/test11b.sam" | wc -l)" -lt 1 ]]; then
  log_error "--cs_tag short produced no cs tags"
  exit 1
fi
log "TEST 11b completed successfully"

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
  awk '!/^@/ {print $1, $2, $3, $4, $5, $6}' "$meta_temp_dir/$conc.sam" | sort > "$meta_temp_dir/$conc.key"
  awk '!/^@/ {print $1, $2, $3, $4, $5, $6}' "$meta_temp_dir/$seq.sam" | sort > "$meta_temp_dir/$seq.key"
  if ! diff -q "$meta_temp_dir/$conc.key" "$meta_temp_dir/$seq.key" > /dev/null; then
    log_error "Concurrent run $conc differs from the sequential run $seq"
    diff "$meta_temp_dir/$seq.key" "$meta_temp_dir/$conc.key" | head -10
    exit 1
  fi
done
log "TEST 12 completed successfully"

print_test_summary "All tests completed successfully"
