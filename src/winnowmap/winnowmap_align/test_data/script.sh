#!/bin/bash

# Generates the winnowmap_align test fixtures. Run once from the repository
# root; the generated files are committed alongside this script.
#
#   bash src/winnowmap/winnowmap_align/test_data/script.sh
#
# The shared create_test_fasta helper emits a pure ATCGATCG... tandem repeat,
# which collapses to a single distinct minimizer: meryl then reports
# "no. of kmers read=0" and the weighted-minimizer code path that Winnowmap
# exists for is never exercised. create_test_fastq likewise emits reads that
# are unrelated to the reference, so nothing aligns. So build fixtures here
# instead: one unique contig plus one contig carrying real tandem-repeat
# arrays, and reads sampled FROM that reference.

set -eo pipefail

OUT_DIR=src/winnowmap/winnowmap_align/test_data
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

echo "Generating reference (unique contig + 20bp x600 and 31bp x250 tandem arrays)..."
generate_reference "$OUT_DIR/reference.fasta"

echo "Generating reads sampled from the reference..."
generate_reads "$OUT_DIR/reference.fasta" "$OUT_DIR/reads.fastq" fastq "$read_length"
generate_reads "$OUT_DIR/reference.fasta" "$OUT_DIR/reads.fasta" fasta "$read_length"
gzip -c "$OUT_DIR/reads.fastq" > "$OUT_DIR/reads.fastq.gz"

echo "Wrote $(( $(wc -l < "$OUT_DIR/reads.fastq") / 4 )) reads of ${read_length}bp to $OUT_DIR"
