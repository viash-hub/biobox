#!/bin/bash

set -eo pipefail

tmp_dir=$(mktemp -d -p "$meta_temp_dir" qualimap_XXXXXXXXX)

# Handle output parameters
if [ -n "$par_report" ]; then
    outfile=$(basename "$par_report")
    report_extension="${outfile##*.}"
fi

if [ -n "$par_counts" ]; then
    counts=$(basename "$par_counts")
fi

# disable flags
[[ "$par_paired" == "false" ]] && unset par_paired
[[ "$par_sorted" == "false" ]] && unset par_sorted

# Run qualimap
qualimap rnaseq \
    ${meta_memory_mb:+--java-mem-size=${meta_memory_mb}M} \
    ${par_algorithm:+--algorithm $par_algorithm} \
    ${par_sequencing_protocol:+--sequencing-protocol $par_sequencing_protocol} \
    -bam $par_bam \
    -gtf $par_gtf \
    -outdir "$tmp_dir" \
    ${par_num_pr_bases:+--num-pr-bases $par_num_pr_bases} \
    ${par_num_tr_bias:+--num-tr-bias $par_num_tr_bias} \
    ${par_report:+-outformat $report_extension} \
    ${par_paired:+--paired} \
    ${par_sorted:+--sorted} \
    ${par_report:+-outfile "$outfile"} \
    ${par_counts:+-oc "$counts"}

# Move output files
mv "$tmp_dir/rnaseq_qc_results.txt" "$par_qc_results"

if [ -n "$par_report" ] && [ $report_extension = "html" ]; then
    mv "$tmp_dir/qualimapReport.html" "$par_report"
fi

if [ -n "$par_report" ] && [ $report_extension = "pdf" ]; then
    mv "$tmp_dir/$outfile" "$par_report"
fi

if [ -n "$par_counts" ]; then
    mv "$tmp_dir/$counts" "$par_counts"
fi
