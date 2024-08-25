#!/bin/bash

set -ep pipefail

mkdir -p tmp/genes
# cut -f 1,2 `ls $par_count_genes/*` | head -n 1` > gene_ids.txt
for file_id in ${par_count_genes[*]}; do
    samplename=`basename $file_id | sed s/\\.genes.results\$//g`
    echo $samplename > tmp/genes/${samplename}.counts.txt
    cut -f 5 ${file_id} | tail -n+2 >> tmp/genes/${samplename}.counts.txt
    echo $samplename > tmp/genes/${samplename}.tpm.txt
    cut -f 6 ${file_id} | tail -n+2 >> tmp/genes/${samplename}.tpm.txt
done

mkdir -p tmp/isoforms
# cut -f 1,2 `ls $par_counts_transcripts/*` | head -n 1` > transcript_ids.txt
for file_id in ${par_counts_transcripts[*]}; do
    samplename=`basename $file_id | sed s/\\.isoforms.results\$//g`
    echo $samplename > tmp/isoforms/${samplename}.counts.txt
    cut -f 5 ${file_id} | tail -n+2 >> tmp/isoforms/${samplename}.counts.txt
    echo $samplename > tmp/isoforms/${samplename}.tpm.txt
    cut -f 6 ${file_id} | tail -n+2 >> tmp/isoforms/${samplename}.tpm.txt
done

paste gene_ids.txt tmp/genes/*.counts.txt > $par_merged_gene_counts
paste gene_ids.txt tmp/genes/*.tpm.txt > $par_merged_gene_tpm
paste transcript_ids.txt tmp/isoforms/*.counts.txt > $par_merged_transcript_counts
paste transcript_ids.txt tmp/isoforms/*.tpm.txt > $par_merged_transcript_tpm
