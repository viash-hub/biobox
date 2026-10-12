#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# Create output directory
mkdir -p "$par_output"

# Determine the index base name for the output
index_basename=$(basename "$par_input" .gz)
index_basename=$(basename "$index_basename" .fasta)
index_basename=$(basename "$index_basename" .fa)
index_basename=$(basename "$index_basename" .fna)

# Unset boolean_true parameters
[[ "$par_emit_unpacked_ref" == "false" ]] && unset par_emit_unpacked_ref
[[ "$par_meth" == "false" ]] && unset par_meth

# Create temporary directory
tmp_dir=$(mktemp -d "$meta_temp_dir/$meta_name-XXXXXXXX")
function clean_up {
  rm -rf "$tmp_dir"
}
trap clean_up EXIT

# Set prefix to write directly to output directory
if [ -n "$par_prefix" ]; then
  # Use custom prefix in output directory
  output_prefix="$par_output/$par_prefix"
else
  # Use input filename (without extension) as prefix in output directory
  output_prefix="$par_output/$index_basename"
fi

# Build the command
cmd_args=(
  ${par_meth:+--meth}
  ${par_emit_unpacked_ref:+--emit-unpacked-ref}
  ${par_sa_sample_rate:+-u "$par_sa_sample_rate"}
  ${meta_cpus:+-t "$meta_cpus"}
  ${meta_memory_gb:+--max-memory "${meta_memory_gb}G"}
  --tmp-dir "$tmp_dir"
)

if [ -n "$par_meth" ]; then
  # --meth does not accept prefix parameter (i.e. '-p') and writes <in.fasta>.* and <in.fasta>.meth.*
  # next to the input, so index a copy named after the output prefix
  cp "$par_input" "$output_prefix"
  bwa-mem3 index "${cmd_args[@]}" "$output_prefix"
  rm "$output_prefix"
else
  bwa-mem3 index -p "$output_prefix" "${cmd_args[@]}" "$par_input"
fi
