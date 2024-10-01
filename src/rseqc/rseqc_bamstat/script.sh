#!/bin/bash


set -eo pipefail 

bam_stat.py \
    --input-file "${par_input_file}" \
    ${par_mapq:+--mapq "${par_mapq}"} \
> $par_output
