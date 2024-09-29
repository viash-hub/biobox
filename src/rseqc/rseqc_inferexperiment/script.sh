#!/bin/bash

set -eo pipefail 

infer_experiment.py \
    -i $par_input \
    -r $par_refgene \
    ${par_sample_size:+-s "${par_sample_size}"} \
    ${par_map_qual:+-q "${par_map_qual}"} \
> $par_output
