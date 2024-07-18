#!/bin/bash

echo ">>> Testing $meta_functionality_name"

gunzip "$meta_resources_dir/genes.gtf"

echo ">>> Prepare RSEM reference"
"$meta_executable" \
  --fasta "$meta_resources_dir/genome.fasta" \
  --gtf "$meta_resources_dir/genes.gtf" \
  --star true \
  --rsem RSEM_index \
  --transcript_fasta transcripts.fasta 

echo ">>> Checking whether output exists"
[ ! -d "RSEM_index" ] && echo "RSEM index does not exist!" && exit 1
[ -z "$(ls -A 'RSEM_index')" ] && echo "RSEM index is empty!" && exit 1
[ ! -f "transcripts.fasta" ] && echo "Transcripts FASTA file does not exist!" && exit 1
[ ! -s "transcripts.fasta" ] && echo "Transcripts FASTA file is empty!" && exit 1

echo ">>> Make transcripts FASTA file"
"$meta_executable" \
  --fasta "$meta_resources_dir/genome.fasta" \
  --gtf "$meta_resources_dir/genes.gtf" \
  --star false \
  --transcript_fasta transcripts.fasta 

echo ">>> Checking whether output exists"
[ ! -f "transcripts.fasta" ] && echo "Transcripts FASTA file does not exist!" && exit 1
[ ! -s "transcripts.fasta" ] && echo "Transcripts FASTA file is empty!" && exit 1

echo "All tests succeeded!"
exit 0
