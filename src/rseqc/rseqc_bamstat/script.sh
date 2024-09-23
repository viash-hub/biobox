#!/bin/bash

set -eo pipefail 

bam_stat.py \
    --input $par_input \
    --mapq $par_map_qual \
> $par_output
