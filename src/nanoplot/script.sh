#!/bin/bash

## VIASH START
par_output='/output'
## VIASH END

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

NanoPlot \
    "$par_input_fastq" \
    > "$par_output"