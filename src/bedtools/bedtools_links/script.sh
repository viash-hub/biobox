#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# Execute bedtools links
bedtools links \
    ${par_base_url:+-base "$par_base_url"} \
    ${par_organism:+-org "$par_organism"} \
    ${par_database:+-db "$par_database"} \
    -i "$par_input" \
    > "$par_output"
