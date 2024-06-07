#!/bin/bash

set -eo pipefail

tmp_dir=$(mktemp -d -p "$meta_temp_dir" qualimap_XXXXXXXXX)


if [ -n "$par_output_file" ]; then
    outfile=$(basename "$par_output_file")
fi

if [ -n "$par_output_counts" ]; then
    counts=$(basename "$par_output_counts")
fi

qualimap rnaseq \
    --java-mem-size=$par_java_memory_size \
    --algorithm $par_algorithm \
    --num-pr-bases $par_pr_bases \
    --num-tr-bias $par_tr_bias \
    --sequencing-protocol $par_sequencing_protocol \
    -bam $par_input \
    -gtf $par_gtf \
    -outdir "$tmp_dir" \
    -outformat pdf \
    ${par_output_file:+-outfile "$outfile"} \
    ${par_output_counts:+-oc "$counts"}

mv "$tmp_dir/rnaseq_qc_results.txt" "$par_output"

if [ -n "$par_output_file" ]; then
    mv "$tmp_dir/$outfile" "$par_output_file"
fi

if [ -n "$par_output_counts" ]; then
    mv "$tmp_dir/$counts" "$par_output_counts"
fi
