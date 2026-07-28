#!/bin/bash

# Regenerates the test files in this directory. Not used by the component or
# its tests at runtime. Run manually with:
#   bash src/mosdepth/mosdepth/test_data/script.sh
# Requires Docker.

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

SAMTOOLS_IMAGE="quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1"
samtools() {
  docker run --rm -v "$PWD":/data -w /data "$SAMTOOLS_IMAGE" samtools "$@"
}

# --- Reference FASTA -------------------------------------------------------
# A single 1000bp synthetic "chr1" contig with a repeating ATCG pattern
{
  echo ">chr1"
  head -c 1000 /dev/zero | tr '\0' 'A' | sed 's/A/ATCG/g' | head -c 1000
  echo
} > reference.fasta
samtools faidx reference.fasta

# --- Alignment (SAM) -------------------------------------------------------
# Six synthetic, properly-paired read pairs (MAPQ 60, insert size 150bp) at
# staggered start positions, so per-base depth builds up from 0 to 5 and back
# down across chr1:100-300
SEQ50=$(head -c 50 /dev/zero | tr '\0' 'A' | sed 's/A/ACGT/g' | head -c 50)
QUAL50=$(printf '%*s' 50 '' | tr ' ' 'I')

{
  echo -e "@HD\tVN:1.6\tSO:coordinate"
  echo -e "@SQ\tSN:chr1\tLN:1000"
  echo -e "@PG\tID:mosdepth_test_data\tPN:mosdepth_test_data\tVN:1.0"
  for i in 1 2 3 4 5 6; do
    r1_start=$((100 + (i - 1) * 10))
    r2_start=$((r1_start + 100))
    tlen=150
    # flag 99  = paired + proper_pair + mate_reverse + first_in_pair
    # flag 147 = paired + proper_pair + read_reverse  + second_in_pair
    echo -e "pair${i}\t99\tchr1\t${r1_start}\t60\t50M\t=\t${r2_start}\t${tlen}\t${SEQ50}\t${QUAL50}"
    echo -e "pair${i}\t147\tchr1\t${r2_start}\t60\t50M\t=\t${r1_start}\t-${tlen}\t${SEQ50}\t${QUAL50}"
  done
} > test.sam

samtools sort -O bam -o test.paired_end.sorted.bam test.sam
samtools index test.paired_end.sorted.bam
rm test.sam

# --- CRAM (same alignment, for --fasta/CRAM-input tests) -------------------
samtools view -O cram -T reference.fasta -o test.cram test.paired_end.sorted.bam
samtools index test.cram
