#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

source "$meta_resources_dir/gatk4/script_helpers.sh"

# Build command arguments array
cmd_args=(
  -R "$par_input"
  -O "$par_output"
  ${par_alt_names:+--ALT_NAMES "$par_alt_names"}
  ${par_genome_assembly:+--GENOME_ASSEMBLY "$par_genome_assembly"}
  ${par_num_sequences:+--NUM_SEQUENCES "$par_num_sequences"}
  ${par_species:+--SPECIES "$par_species"}
  ${par_truncate_names_at_whitespace:+--TRUNCATE_NAMES_AT_WHITESPACE "$par_truncate_names_at_whitespace"}
  ${par_uri:+--URI "$par_uri"}
)

# Run GATK CreateSequenceDictionary
run_gatk CreateSequenceDictionary "${cmd_args[@]}"
