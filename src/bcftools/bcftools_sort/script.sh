#!/bin/bash

## VIASH START
## VIASH END

# Exit on error
set -eo pipefail

# Execute bedtools bamtofastq with the provided arguments
bcftools sort \
    -o "$par_output" \
    ${par_output_type:+-O "$par_output_type"} \
    $par_input \

