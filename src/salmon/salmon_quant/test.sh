#!/bin/bash

set -e

dir_in="$meta_resources_dir/test_data"

echo "> Run salmon quant for single-end reads"
"$meta_executable" \
  --lib_type "A" \
  --index "$dir_in/salmon/transcriptome_index" \
  --unmated_reads "$dir_in/reads/a_se.fq.gz" \
  --output "quant_se_results" \
  --quant_results "quant_se.sf"

echo ">> Checking output"
[ ! -d "quant_se_results" ] && echo "Output directory quant_se_results does not exist" && exit 1
[ ! -f "quant_se.sf" ] && echo "Output file quant_se.sf does not exist!" && exit 1
[ ! -s "quant_se.sf" ] && echo "Output file quant_se.sf is empty!" && exit 1
grep -q "Name	Length	EffectiveLength	TPM	NumReads" "quant_se.sf" || (echo "Output file quant_se.sf does not have the right format!" && exit 1)

echo "> Run salmon quant for paired-end reads"
"$meta_executable" \
  --lib_type "A" \
  --index "$dir_in/salmon/transcriptome_index" \
  --mates1 "$dir_in/reads/a_1.fq.gz" \
  --mates2 "$dir_in/reads/a_2.fq.gz" \
  --output "quant_pe_results" \
  --quant_results "quant_pe.sf"

echo ">> Checking output"
[ ! -d "quant_pe_results" ] && echo "Output directory quant_pe_results does not exist" && exit 
[ ! -f "quant_pe.sf" ] && echo "Output file quant_pe.sf does not exist!" && exit 1
[ ! -s "quant_pe.sf" ] && echo "Output file quant_pe.sf is empty!" && exit 1
grep -q "Name	Length	EffectiveLength	TPM	NumReads" "quant_pe.sf" || (echo "Output file quant_pe.sf does not have the right format!" && exit 1)

echo "> Run salmon quant for paired-end reads with technical replicates"
"$meta_executable" \
  --lib_type "A" \
  --index "$dir_in/salmon/transcriptome_index" \
  --mates1 "$dir_in/reads/a_1.fq.gz;$dir_in/reads/b_1.fq.gz" \
  --mates2 "$dir_in/reads/a_2.fq.gz;$dir_in/reads/b_2.fq.gz" \
  --output "quant_pe_rep_results" \
  --quant_results "quant_pe_rep.sf"

echo ">> Checking output"
[ ! -d "quant_pe_rep_results" ] && echo "Output directory quant_pe_rep_results does not exist" && exit 
[ ! -f "quant_pe_rep.sf" ] && echo "Output file quant_pe_rep.sf does not exist!" && exit 1
[ ! -s "quant_pe_rep.sf" ] && echo "Output file quant_pe_rep.sf is empty!" && exit 1
grep -q "Name	Length	EffectiveLength	TPM	NumReads" "quant_pe_rep.sf" || (echo "Output file quant_pe_rep.sf does not have the right format!" && exit 1)

echo "> Test successful"