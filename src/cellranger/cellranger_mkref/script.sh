#!/bin/bash

set -eo pipefail

## VIASH START
par_genome_fasta="test_data/reference_small.fa.gz"
par_transcriptome_gtf="test_data/reference_small.gtf.gz"
par_output="output.tar.gz"
## VIASH END

# create temporary directory
tmp_dir=$(mktemp -d "$meta_temp_dir/$meta_name-XXXXXXXX")
function clean_up {
    rm -rf "$tmp_dir"
}
trap clean_up EXIT

# We change into the tempdir later, so we need absolute paths.
par_genome_fasta=$(realpath $par_genome_fasta)
par_transcriptome_gtf=$(realpath $par_transcriptome_gtf)
par_output=$(realpath $par_output)

# if memory is defined, subtract 2GB from memory
if [[ "$meta_memory_gb" != "" ]]; then
  # if memory is less than 2gb, unset it
  if [[ "$meta_memory_gb" -lt 2 ]]; then
    unset meta_memory_gb
  else
    meta_memory_gb=$((meta_memory_gb-2))
  fi
fi

echo "> Unzipping input files"
unpigz -c "$par_genome_fasta" > "$tmp_dir/genome.fa"

echo "> Building star index"
cd "$tmp_dir"
cellranger mkref \
  --fasta "$tmp_dir/genome.fa" \
  --genes "$par_transcriptome_gtf" \
  --genome output \
  ${par_reference_version:+--ref-version $par_reference_version} \
  ${meta_cpus:+--nthreads $meta_cpus} \
  ${meta_memory_gb:+--memgb ${meta_memory_gb}}

echo "> Creating archive"
tar --use-compress-program="pigz -k " -cf "$par_output" -C "$tmp_dir/output" .

exit 0
