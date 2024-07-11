#!/bin/bash

set -eo pipefail

tmp_dir=$(mktemp -d -p "$meta_temp_dir" qualimap_XXXXXXXXX)

# Handle output parameters
if [ -n "$par_output_report" ]; then
    outfile=$(basename "$par_output_report")
    report_extension="${outfile##*.}"
fi

if [ -n "$par_output_counts" ]; then
    counts=$(basename "$par_output_counts")
fi

# disable flags
[[ "$par_paired" == "false" ]] && unset par_paired
[[ "$par_sorted" == "false" ]] && unset par_sorted

# Run qualimap
qualimap rnaseq \
    --java-mem-size=$par_java_memory_size \
    --algorithm $par_algorithm \
    --num-pr-bases $par_pr_bases \
    --num-tr-bias $par_tr_bias \
    --sequencing-protocol $par_sequencing_protocol \
    -bam $par_input \
    -gtf $par_gtf \
    -outdir "$tmp_dir" \
    --java-mem-size=$par_java_memory_size \
    ${par_output_report:+-outformat $report_extension} \
    ${par_paired:+--paired} \
    ${par_sorted:+--sorted} \
    ${par_output_report:+-outfile "$outfile"} \
    ${par_output_counts:+-oc "$counts"}

# Move output files
mv "$tmp_dir/rnaseq_qc_results.txt" "$par_output"

if [ -n "$par_output_report" ] && [ $report_extension = "html" ]; then
    mv "$tmp_dir/qualimapReport.html" "$par_output_report"
fi

if [ -n "$par_output_report" ] && [ $report_extension = "pdf" ]; then
    mv "$tmp_dir/$outfile" "$par_output_report"
fi

if [ -n "$par_output_counts" ]; then
    mv "$tmp_dir/$counts" "$par_output_counts"
fi
