#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# unset flags
[[ "$par_nosecondary" == "false" ]] && unset par_nosecondary
[[ "$par_dry" == "false" ]] && unset par_dry
[[ "$par_nopreflight" == "false" ]] && unset par_nopreflight

# Validate required arguments
if [[ -z "$par_id" ]]; then
  echo "Error: --id is required" >&2
  exit 1
fi

if [[ -z "$par_csv" ]]; then
  echo "Error: --csv is required" >&2
  exit 1
fi

if [[ -z "$par_output_dir" ]]; then
  echo "Error: --output_dir is required" >&2
  exit 1
fi

# Make absolute paths
par_csv=$(realpath "$par_csv")
par_output_dir=$(realpath "$par_output_dir")

# if memory is defined, subtract 2GB from memory
if [[ "$meta_memory_gb" != "" ]]; then
  # if memory is less than 2gb, unset it
  if [[ "$meta_memory_gb" -lt 2 ]]; then
    echo "WARNING: Memory is less than 2GB, unsetting memory requirements"
    unset meta_memory_gb
  else
    meta_memory_gb=$((meta_memory_gb-2))
  fi
fi

# Create temporary directory
tmp_dir=$(mktemp -d "$meta_temp_dir/$meta_name-XXXXXXXX")
function clean_up {
    rm -rf "$tmp_dir"
}
trap clean_up EXIT

# Change to temp directory
cd "$tmp_dir"

# Run cellranger aggr
echo "> Running cellranger aggr"
cellranger aggr \
  --id="$par_id" \
  --csv="$par_csv" \
  --disable-ui \
  ${meta_cpus:+--localcores="$meta_cpus"} \
  ${meta_memory_gb:+--localmem="$meta_memory_gb"} \
  ${par_description:+--description="$par_description"} \
  ${par_normalize:+--normalize="$par_normalize"} \
  ${par_nosecondary:+--nosecondary} \
  ${par_dry:+--dry} \
  ${par_min_crispr_umi:+--min-crispr-umi="$par_min_crispr_umi"} \
  ${par_enable_tsne:+--enable-tsne="$par_enable_tsne"} \
  ${par_nopreflight:+--nopreflight}

# Copy output
echo "> Copying output"
if [ -d "$par_id/outs/" ]; then
  if [ ! -d "$par_output_dir" ]; then
    mkdir -p "$par_output_dir"
  fi
  mv "$par_id/outs/"* "$par_output_dir"
fi
