#!/bin/bash

set -eo pipefail

## VIASH START
par_input='/opt/cellranger-8.0.0/lib/python/cellranger-tiny-fastq'
par_reference='/opt/cellranger-8.0.0/lib/python/cellranger-tiny-ref'
par_output='test_data/bam'
par_chemistry="auto"
par_expect_cells="3000"
par_secondary_analysis="false"
## VIASH END

## PROCESS INPUT FILES
# We change into the tempdir later, so we need absolute paths.
par_reference=$(realpath $par_reference)
par_output=$(realpath $par_output)

# create temporary directory
tmp_dir=$(mktemp -d "$meta_temp_dir/$meta_name-XXXXXXXX")
function clean_up {
    rm -rf "$tmp_dir"
}
trap clean_up EXIT

# process inputs
# for every fastq file found, make a symlink into the tempdir
fastq_dir="$tmp_dir/fastqs"
mkdir -p "$fastq_dir"
IFS=";"
for var in $par_input; do
  unset IFS
  abs_path=$(realpath $var)
  if [ -d "$abs_path" ]; then
    find "$abs_path" -name *.fastq.gz -exec ln -s {} "$fastq_dir" \;
  else
    ln -s "$abs_path" "$fastq_dir"
  fi
done

# process reference
if file "${par_reference}" | grep -q 'gzip compressed data'; then
  echo "> Untarring genome"
  ref_dir="${tmp_dir}/reference"
  mkdir -p "$ref_dir"
  tar -xvf "${par_reference}" -C "$ref_dir"
  par_reference="${ref_dir}"
fi

## PROCESS PARAMETERS
# unset flags
[[ "$par_secondary_analysis" == "false" ]] && unset par_secondary_analysis

# change ifs from ; to ,
par_lanes=${par_lanes//;/,}

# if memory is defined, subtract 2GB from memory
if [[ "$meta_memory_gb" != "" ]]; then
  # if memory is less than 2gb, unset it
  if [[ "$meta_memory_gb" -lt 2 ]]; then
    unset meta_memory_gb
  else
    meta_memory_gb=$((meta_memory_gb-2))
  fi
fi

## RUN CELLRANGER COUNT
echo "> Running cellranger count"
cd "$tmp_dir"
id=myoutput
cellranger count \
  --id="$id" \
  --fastqs="${fastq_dir}" \
  --transcriptome="${par_reference}" \
  --include-introns="${par_include_introns}" \
  ${meta_cpus:+"--localcores=${meta_cpus}"} \
  ${meta_memory_gb:+"--localmem=${meta_memory_gb}"} \
  ${par_expect_cells:+"--expect-cells=${par_expect_cells}"} \
  ${par_force_cells:+"--force-cells=${par_force_cells}"} \
  ${par_chemistry:+"--chemistry=${par_chemistry}"} \
  ${par_generate_bam:+"--create-bam=${par_generate_bam}"} \
  ${no_secondary_analysis:+--nosecondary} \
  ${par_r1_length:+"--r1-length=${par_r1_length}"} \
  ${par_r2_length:+"--r2-length=${par_r2_length}"} \
  ${par_lanes:+"--lanes=${par_lanes}"} \
  ${par_library_compatibility_check:+"--check-library-compatibility=${par_library_compatibility_check}"}\
  --disable-ui

echo "> Copying output"
if [ -d "$id/outs/" ]; then
  if [ ! -d "${par_output}" ]; then
    mkdir -p "${par_output}"
  fi
  mv "$id/outs/"* "${par_output}"
fi

exit 0
