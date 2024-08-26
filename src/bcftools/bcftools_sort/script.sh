#!/bin/bash

## VIASH START
## VIASH END

# Exit on error
set -eo pipefail

# Unset parameters
# [[ "$par_tags" == "false" ]] && unset par_tags

# Execute bedtools bamtofastq with the provided arguments
bcftools sort \

