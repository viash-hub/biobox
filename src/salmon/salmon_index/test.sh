#!/bin/bash

set -e

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

printf ">>> Run salmon_index"
"$meta_executable" \
  --transcripts $dir_in/transcriptome.fasta \
  --index index \
  --kmer_len 31

printf ">>> Checking whether output exists"
[ ! -d "index" ] && echo "'index' does not exist!" && exit 1
[ -z "$(ls -A 'index')" ] && echo "'index' is empty!" && exit 1
[ ! -f "index/info.json" ] && echo "Salmon index does not contain 'info.json'! Not all files were generated correctly!" && exit 1
[ $(grep '"k": [0-9]*' index/info.json | cut -d':' -f 2) != '31,' ] && printf "The generated Salmon index seems to be incorrect!" && exit 1

echo "All tests succeeded!"
exit 0