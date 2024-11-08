#!/bin/bash

set -eo pipefail

## VIASH START
par_genome_fasta="test_data/reference_small.fa.gz"
par_transcriptome_gtf="test_data/reference_small.gtf.gz"
par_output="output.tar.gz"
## VIASH END

# create temporary directory
tmpdir=$(mktemp -d "$meta_temp_dir/$meta_name-XXXXXXXX")
function clean_up {
    rm -rf "$tmpdir"
}
trap clean_up EXIT

# We change into the tempdir later, so we need absolute paths.
par_genome_fasta=$(realpath $par_genome_fasta)
par_transcriptome_gtf=$(realpath $par_transcriptome_gtf)
par_output=$(realpath $par_output)


echo "> Unzipping input files"
unpigz -c "$par_genome_fasta" > "$tmpdir/genome.fa"

echo "> Building star index"
cd "$tmpdir"
cellranger mkref \
  --fasta "$tmpdir/genome.fa" \
  --genes "$par_transcriptome_gtf" \
  --genome output \
  ${par_reference_version:+--ref-version $par_reference_version} \
  ${meta_cpus:+--nthreads $meta_cpus} \
  ${meta_memory_gb:+--memgb $(($meta_memory_gb-2))} # always keep 2 gb for the OS itself

echo "> Creating archive"
tar --use-compress-program="pigz -k " -cf "$par_output" -C "$tmpdir/output" .
