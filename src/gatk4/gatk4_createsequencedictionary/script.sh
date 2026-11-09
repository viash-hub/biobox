#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# Compute available memory for the JVM (80% of allocated memory, fallback to 3072MB)
avail_mem_mb=$(( ${meta_memory_mb:-3072} * 8 / 10 ))

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

# Run gatk CreateSequenceDictionary
gatk --java-options "-Xmx${avail_mem_mb}M -XX:-UsePerfData" CreateSequenceDictionary "${cmd_args[@]}"
