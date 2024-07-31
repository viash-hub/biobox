#!/bin/bash

TMP_DIR=/tmp/bd_rhapsody_make_reference
OUT_DIR=src/bd_rhapsody/bd_rhapsody_make_reference/test_data

# check if seqkit is installed
if ! command -v seqkit &> /dev/null; then
  echo "seqkit could not be found"
  exit 1
fi

# create temporary directory and clean up on exit
mkdir -p $TMP_DIR
function clean_up {
    rm -rf "$TMP_DIR"
}
trap clean_up EXIT

# fetch reference
ORIG_FA=$TMP_DIR/reference.fa.gz
if [ ! -f $ORIG_FA ]; then
  wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/GRCh38.primary_assembly.genome.fa.gz \
    -O $ORIG_FA
fi

ORIG_GTF=$TMP_DIR/reference.gtf.gz
if [ ! -f $ORIG_GTF ]; then
  wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.annotation.gtf.gz \
    -O $ORIG_GTF
fi

# create small reference
START=30000
END=31500
CHR=chr1

# subset to small region
seqkit grep -r -p "^$CHR\$" "$ORIG_FA" | \
  seqkit subseq -r "$START:$END" > $OUT_DIR/reference_small.fa

zcat "$ORIG_GTF" | \
  awk -v FS='\t' -v OFS='\t' "
    \$1 == \"$CHR\" && \$4 >= $START && \$5 <= $END {
      \$4 = \$4 - $START + 1;
      \$5 = \$5 - $START + 1;
      print;
    }" > $OUT_DIR/reference_small.gtf
