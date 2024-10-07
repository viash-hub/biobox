#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

# Unset flags
[[ "$par_verbose" == "false" ]] && unset par_verbose
[[ "$par_store" == "false" ]] && unset par_store
[[ "$par_raw" == "false" ]] && unset par_raw
[[ "$par_huge" == "false" ]] && unset par_huge
[[ "$par_no_static" == "true" ]] && unset par_no_static
[[ "$par_tsv_stats" == "false" ]] && unset par_tsv_stats
[[ "$par_only_report" == "false" ]] && unset par_only_report
[[ "$par_info_in_report" == "false" ]] && unset par_info_in_report
[[ "$par_drop_outliers" == "true" ]] && unset par_drop_outliers
[[ "$par_loglength" == "false" ]] && unset par_loglength
[[ "$par_percentqual" == "false" ]] && unset par_percentqual
[[ "$par_alength" == "false" ]] && unset par_alength
[[ "$par_barcoded" == "false" ]] && unset par_barcoded
[[ "$par_no_supplementary" == "true" ]] && unset par_no_supplementary
[[ "$par_listcolors" == "false" ]] && unset par_listcolors
[[ "$par_listcolormaps" == "false" ]] && unset par_listcolormaps
[[ "$par_no_N50" == "true" ]] && unset par_no_N50
[[ "$par_N50" == "false" ]] && unset par_N50
[[ "$par_hide_stats" == "true" ]] && unset par_hide_stats

# Multiple inputs - replace ';' with ' ' (space)
par_fastq=$(echo $par_fastq | tr ';' ' ')
par_fasta=$(echo $par_fasta | tr ';' ' ')
par_fastq_rich=$(echo $par_fastq_rich | tr ';' ' ')
par_fastq_minimal=$(echo $par_fastq_minimal | tr ';' ' ')
par_summary=$(echo $par_summary | tr ';' ' ')
par_bam=$(echo $par_bam | tr ';' ' ')
par_ubam=$(echo $par_ubam | tr ';' ' ')
par_cram=$(echo $par_cram | tr ';' ' ')
par_pickle=$(echo $par_pickle | tr ';' ' ')
par_feather=$(echo $par_feather | tr ';' ' ')

# Run NanoPlot
NanoPlot \
    ${par_fastq:+--fastq "$par_fastq"} \
    ${par_fasta:+--fasta "$par_fasta"} \
    ${par_fastq_rich:+--fastq_rich "$par_fastq_rich"} \
    ${par_fastq_minimal:+--fastq_minimal "$par_fastq_minimal"} \
    ${par_summary:+--summary "$par_summary"} \
    ${par_bam:+--bam "$par_bam"} \
    ${par_ubam:+--ubam "$par_ubam"} \
    ${par_cram:+--cram "$par_cram"} \
    ${par_pickle:+--pickle "$par_pickle"} \
    ${par_feather:+--feather "$par_feather"} \
    ${par_verbose:+--verbose} \
    ${par_store:+--store} \
    ${par_raw:+--raw} \
    ${par_huge:+--huge} \
    ${par_no_static:+--no_static} \
    ${par_prefix:+--prefix "$par_prefix"} \
    ${par_tsv_stats:+--tsv_stats} \
    ${par_only_report:+--only-report} \
    ${par_info_in_report:+--info_in_report} \
    ${par_maxlength:+--maxlength "$par_maxlength"} \
    ${par_minlength:+--minlength "$par_minlength"} \
    ${par_drop_outliers:+--drop_outliers} \
    ${par_downsample:+--downsample "$par_downsample"} \
    ${par_loglength:+--loglength} \
    ${par_percentqual:+--percentqual} \
    ${par_alength:+--alength} \
    ${par_minqual:+--minqual "$par_minqual"} \
    ${par_runtime_until:+--runtime_until "$par_runtime_until"} \
    ${par_readtype:+--readtype "$par_readtype"} \
    ${par_barcoded:+--barcoded} \
    ${par_no_supplementary:+--no_supplementary} \
    ${par_color:+--color "$par_color"} \
    ${par_colormap:+--colormap "$par_colormap"} \
    ${par_format:+--format "$par_format"} \
    ${par_plots:+--plots "$par_plots"} \
    ${par_legacy:+--legacy "$par_legacy"} \
    ${par_listcolors:+--listcolors} \
    ${par_listcolormaps:+--listcolormaps} \
    ${par_no_N50:+--no-N50} \
    ${par_N50:+--N50} \
    ${par_title:+--title "$par_title"} \
    ${par_font_scale:+--font_scale "$par_font_scale"} \
    ${par_dpi:+--dpi "$par_dpi"} \
    ${par_hide_stats:+--hide_stats} \
    ${meta_cpus:+--threads "$meta_cpus"} \
    --outdir "$par_outdir"