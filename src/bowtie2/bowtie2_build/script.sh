#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# unset flags
[[ "$par_large_index" == "false" ]] && unset par_large_index
[[ "$par_noauto" == "false" ]] && unset par_noauto
[[ "$par_packed" == "false" ]] && unset par_packed
[[ "$par_nodc" == "false" ]] && unset par_nodc
[[ "$par_noref" == "false" ]] && unset par_noref
[[ "$par_justref" == "false" ]] && unset par_justref
[[ "$par_quiet" == "false" ]] && unset par_quiet
[[ "$par_fasta" == "false" ]] && unset par_fasta
[[ "$par_cmdline" == "false" ]] && unset par_cmdline

# Create output directory
mkdir -p "$par_output"

# Determine the index base name for the output
if [ -n "$par_index_name" ]; then
    index_basename="$par_index_name"
else
    index_basename=$(basename "$par_input" .fasta)
    index_basename=$(basename "$index_basename" .fa)
    index_basename=$(basename "$index_basename" .fna)
fi

# Set the full path for the index
index_path="$par_output/$index_basename"

# Build the command arguments
cmd_args=(
    ${par_fasta:+-f}
    ${par_cmdline:+-c}
    ${par_large_index:+--large-index}
    ${par_noauto:+-a}
    ${par_packed:+-p}
    ${par_bmax:+--bmax "$par_bmax"}
    ${par_bmaxdivn:+--bmaxdivn "$par_bmaxdivn"}
    ${par_dcv:+--dcv "$par_dcv"}
    ${par_nodc:+--nodc}
    ${par_noref:+-r}
    ${par_justref:+-3}
    ${par_offrate:+-o "$par_offrate"}
    ${par_ftabchars:+-t "$par_ftabchars"}
    ${par_seed:+--seed "$par_seed"}
    ${par_quiet:+-q}
    ${meta_cpus:+--threads "$meta_cpus"}
    "$par_input"
    "$index_path"
)

# Run bowtie2-build
bowtie2-build "${cmd_args[@]}"
