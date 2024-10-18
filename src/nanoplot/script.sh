#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

# Unset flags
unset_if_false=( 
    par_verbose
    par_store
    par_raw
    par_huge
    par_no_static
    par_tsv_stats
    par_only_report
    par_info_in_report
    par_drop_outliers
    par_loglength
    par_percentqual
    par_alength
    par_barcoded
    par_no_supplementary
    par_listcolors
    par_listcolormaps
    par_no_N50
    par_N50
    par_hide_stats
)

for var in "${unset_if_false[@]}"; do
    test_val="${!var}"
    [[ "$test_val" == "false" ]] && unset $var
done

par_fastq="${par_fastq//;/ }"
par_fasta="${par_fasta//;/ }"
par_fastq_rich="${par_fastq_rich//;/ }"
par_fastq_minimal="${par_fastq_minimal//;/ }"
par_summary="${par_summary//;/ }"
par_bam="${par_bam//;/ }"
par_ubam="${par_ubam//;/ }"
par_cram="${par_cram//;/ }"
par_pickle="${par_pickle//;/ }"
par_feather="${par_feather//;/ }"


inputs=( 
    "$par_fastq" 
    "$par_fasta"
    "$par_fastq_rich"
    "$par_fastq_minimal"
    "$par_summary"
    "$par_bam"
    "$par_ubam"
    "$par_cram"
    "$par_pickle"
    "$par_feather"
)

one_input=false
for var in "${inputs[@]}"; do
    if [ -n "$var" ]; then # if the parameter is not empty
        if [ "$one_input" = "false" ]; then
            one_input=true
        else # Multiple input file types specified
            echo "Error: Multiple input file types specified."
            exit 1
        fi
    fi
done

if [ ! "$one_input" ]; then
    echo "Error: No input file type specified."
    exit 1
fi



# Run NanoPlot
NanoPlot \
    ${par_fastq:+--fastq $par_fastq} \
    ${par_fasta:+--fasta $par_fasta} \
    ${par_fastq_rich:+--fastq_rich $par_fastq_rich} \
    ${par_fastq_minimal:+--fastq_minimal $par_fastq_minimal} \
    ${par_summary:+--summary $par_summary} \
    ${par_bam:+--bam $par_bam} \
    ${par_ubam:+--ubam $par_ubam} \
    ${par_cram:+--cram $par_cram} \
    ${par_pickle:+--pickle $par_pickle} \
    ${par_feather:+--feather $par_feather} \
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

exit 0
