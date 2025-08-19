#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# unset flags
[[ "$par_preload_index" == "false" ]] && unset par_preload_index
[[ "$par_disable_smith_waterman" == "false" ]] && unset par_disable_smith_waterman
[[ "$par_disable_insert_size_estimate" == "false" ]] && unset par_disable_insert_size_estimate

# Build the command
cmd_args=(
    # Pairing options
    ${par_max_insert_size:+-a "$par_max_insert_size"}
    ${par_max_occ_one_end:+-o "$par_max_occ_one_end"}
    ${par_max_hits_paired:+-n "$par_max_hits_paired"}
    ${par_max_hits_discordant:+-N "$par_max_hits_discordant"}
    ${par_chimeric_rate:+-c "$par_chimeric_rate"}
    
    # Output options
    ${par_output:+-f "$par_output"}
    ${par_read_group:+-r "$par_read_group"}
    
    # Algorithm options
    ${par_preload_index:+-P}
    ${par_disable_smith_waterman:+-s}
    ${par_disable_insert_size_estimate:+-A}
    
    # Required arguments: index, SAI files, FASTQ files
    "$par_index"
    "$par_sai1"
    "$par_sai2"
    "$par_reads1"
    "$par_reads2"
)

# Run bwa sampe
bwa sampe "${cmd_args[@]}"
