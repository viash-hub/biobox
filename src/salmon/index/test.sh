#!/bin/bash

set -e

echo "==============================================================================="
echo "> Prepare test data"

dir_in="test_data"
mkdir -p "$dir_in"

cat > "$dir_in/transcriptome.fasta" <<'EOF'
>contig1
AGCTCCAGATTCGCTCAGGCCCTTGATCATCAGTCGTCGTCGTCTTCGATTTGCCAGAGG
AGTTTAGATGAAGAATGTCAAGGATGTTCCTCCCTGCCCTCCCATCTAGCCAAGAACATT
TCCAAGAAGATAAAACTGTCACTGAGACAGGTCTGGATGCGCCCTAGGGGCAAATAGAGA
>contig2
AGGCCTTTACCACATTGCTGCTGGCTATAGGAAGTCCCAGGTACTAGCCTGAAACAGCTG
ATATTTGGGGCTGTCACAGACAATATGGCCACCCCTTGGTCTTTATGCATGAAGATTATG
TAAAGGTTTTTATTAAAAAATATATATATATATATAAATGATCTAGATTATTTTCCTCTT
TCTGAAGTACTTTCTTAAAAAAATAAAATTAAATGTTTATAGTATTCCCGGT
EOF

echo "==============================================================================="
echo "> Run salmon_index"
"$meta_executable" \
  --transcripts "$dir_in/transcriptome.fasta" \
  --index index

echo ">>> Checking whether output exists"
[ ! -d "index" ] && echo "'index' does not exist!" && exit 1
[ -z "$(ls -A 'index')" ] && echo "'index' is empty!" && exit 1
[ ! -f "index/info.json" ] && echo "Salmon index does not contain 'info.json'! Not all files were generated correctly!" && exit 1

# TODO: check contents of 'index/info.json'

echo "==============================================================================="
echo "> Test successful"
