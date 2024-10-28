#!/bin/bash

set -eo pipefail 

infer_experiment.py \
    -i $par_input_file \
    -r $par_refgene \
    ${par_sample_size:+-s "${par_sample_size}"} \
    ${par_mapq:+-q "${par_mapq}"} \
> $par_output
