#!/bin/bash

set -e

dir_in="$meta_resources_dir/test_data"

echo "> Run salmon quant for single-end reads"
"$meta_executable" \
  --lib_type "A" \
  --index "$dir_in/salmon/transcriptome_index" \
  --unmated_reads "$dir_in/reads/a_se.fq.gz" \
  --output "quant_se_results"

echo ">> Checking output"
[ ! -d "quant_se_results" ] && echo "Output directory quant_se_results does not exist" && exit 
[ ! -s "quant_se_results/quant.sf" ] && echo "Output file quant_se_results/quant.sf is empty" && exit 1

echo "> Run salmon quant for paired-end reads"
"$meta_executable" \
  --lib_type "A" \
  --index "$dir_in/salmon/transcriptome_index" \
  --mates1 "$dir_in/reads/a_1.fq.gz" \
  --mates2 "$dir_in/reads/a_2.fq.gz" \
  --output "quant_pe_results"

echo ">> Checking output"
[ ! -d "quant_pe_results" ] && echo "Output directory quant_pe_results does not exist" && exit 
[ ! -s "quant_pe_results/quant.sf" ] && echo "Output file quant_pe_results/quant.sf is empty" && exit 1

echo "> Run salmon quant for paired-end reads with technical replicates"
"$meta_executable" \
  --lib_type "A" \
  --index "$dir_in/salmon/transcriptome_index" \
  --mates1 $dir_in/reads/a_1.fq.gz $dir_in/reads/b_1.fq.gz \
  --mates2 $dir_in/reads/a_2.fq.gz $dir_in/reads/b_2.fq.gz \
  --output "quant_pe_rep_results"

echo ">> Checking output"
[ ! -d "quant_pe_rep_results" ] && echo "Output directory quant_pe_rep_results does not exist" && exit 
[ ! -s "quant_pe_rep_results/quant.sf" ] && echo "Output file quant_pe_rep_results/quant.sf is empty" && exit 1

echo "> Test successful"