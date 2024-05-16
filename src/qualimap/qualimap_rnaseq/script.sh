#!/bin/bash

sleep 10000

set -eo pipefail

mkdir -p $par_output_dir

qualimap rnaseq \
    # --algorithm=$par_algorithm \
    --java-mem-size=$par_java_memory_size \
    --num-pr-bases $par_pr_bases \
    --num-tr-bias $par_tr_bias \
    --sequencing-protocol $par_sequencing_protocol \
    -bam $par_input \
    -gtf $par_gtf \
    ${par_paired:+-pe} \
    ${par_sorted:+-s} \
    -outdir $par_output_dir \
    # -oc $par_output_counts \
    # -outfile $par_output_file \
    -outformat $par_output_format 
