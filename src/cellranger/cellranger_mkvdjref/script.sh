#!/bin/bash

set -eo pipefail

## VIASH START
par_genome="vdj_reference"
par_fasta="test_data/genome.fa"
par_genes="test_data/genes.gtf"
par_output="vdj_reference"
## VIASH END

# create temporary directory
tmp_dir=$(mktemp -d "$meta_temp_dir/$meta_name-XXXXXXXX")
function clean_up {
    rm -rf "$tmp_dir"
}
trap clean_up EXIT

## PROCESS INPUT FILES
# Cell Ranger reads the FASTA and GTF inputs as plain text, so gzipped files
# have to be decompressed first. Echoes the path to use to stdout.
function decompress_if_needed {
  local input
  input=$(realpath "$1")
  if file --dereference --brief "$input" | grep -q "^gzip compressed"; then
    # a directory per file, so that inputs sharing a basename don't collide
    local output_dir
    output_dir=$(mktemp -d "$tmp_dir/decompressed-XXXXXXXX")
    local output="$output_dir/$(basename "${input%.gz}")"
    echo "> Decompressing $(basename "$input")" >&2
    unpigz -c "$input" > "$output"
    input="$output"
  fi
  echo "$input"
}

# We change into the tempdir later, so we need absolute paths.
par_output=$(realpath "$par_output")

# `--genes` can be specified multiple times, so split the `;`-separated values
# into one `--genes` argument each
genes_args=()
if [[ -n "$par_genes" ]]; then
  IFS=';' read -ra genes <<< "$par_genes"
  for gtf in "${genes[@]}"; do
    genes_args+=("--genes=$(decompress_if_needed "$gtf")")
  done
fi

## PROCESS PARAMETERS
# if memory is defined, subtract 2GB from memory
if [[ -n "$meta_memory_gb" ]]; then
  # if memory is less than 2gb, unset it
  if [[ "$meta_memory_gb" -lt 2 ]]; then
    echo "WARNING: Memory is less than 2GB, unsetting memory requirements"
    unset meta_memory_gb
  else
    meta_memory_gb=$((meta_memory_gb-2))
  fi
fi

cmd_args=(
  --genome="$par_genome"
  --disable-ui
  ${par_fasta:+--fasta="$(decompress_if_needed "$par_fasta")"}
  ${par_seqs:+--seqs="$(decompress_if_needed "$par_seqs")"}
  ${par_rm_transcripts:+--rm-transcripts="$(realpath "$par_rm_transcripts")"}
  ${par_ref_version:+--ref-version="$par_ref_version"}
  ${meta_memory_gb:+--memgb="$meta_memory_gb"}
  "${genes_args[@]}"
)

## RUN CELLRANGER MKVDJREF
# cellranger creates the reference in a folder named after the genome, relative
# to the working directory, next to the pipestance folder
echo "> Running cellranger mkvdjref"
( cd "$tmp_dir" && cellranger mkvdjref "${cmd_args[@]}" )

echo "> Copying output"
mkdir -p "$par_output"
mv "$tmp_dir/$par_genome/"* "$par_output"