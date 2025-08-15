#!/bin/bash

set -eo pipefail

## VIASH START
par_fastqs='/opt/cellranger-9.0.1/lib/python/cellranger-tiny-fastq'
par_transcriptome='/opt/cellranger-9.0.1/lib/python/cellranger-tiny-ref'
par_output='test_data/bam'
par_chemistry="auto"
par_expect_cells="3000"
par_secondary_analysis="false"
## VIASH END

## PROCESS INPUT FILES
# We change into the tempdir later, so we need absolute paths.
par_transcriptome=$(realpath $par_transcriptome)
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
for var in $par_fastqs; do
  unset IFS
  abs_path=$(realpath $var)
  if [ -d "$abs_path" ]; then
    find "$abs_path" -name *.fastq.gz -exec ln -s {} "$fastq_dir" \;
  else
    ln -s "$abs_path" "$fastq_dir"
  fi
done

# process reference
# Note: should we do this?
if file "${par_transcriptome}" | grep -q 'gzip compressed data'; then
  echo "> Untarring transcriptome"
  ref_dir="${tmp_dir}/reference"
  mkdir -p "${ref_dir}"
  tar -xvf "${par_transcriptome}" -C "${ref_dir}"
  par_transcriptome="${ref_dir}"
fi

## PROCESS PARAMETERS
# unset flags
[[ "$par_no_secondary" == "false" ]] && unset par_no_secondary
[[ "$par_no_libraries" == "false" ]] && unset par_no_libraries
[[ "$par_dry" == "false" ]] && unset par_dry


# change ifs from ; to ,
par_lanes=${par_lanes//;/,}

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

## RUN CELLRANGER COUNT
echo "> Running cellranger count"
cd "$tmp_dir"
id=run
cellranger count \
  --id="$id" \
  --fastqs="${fastq_dir}" \
  --transcriptome="${par_transcriptome}" \
  --disable-ui \
  ${meta_cpus:+"--localcores=${meta_cpus}"} \
  ${meta_memory_gb:+"--localmem=${meta_memory_gb}"} \
  ${par_description:+"--description=${par_description}"} \
  ${par_sample:+"--sample=${par_sample}"} \
  ${par_lanes:+"--lanes=${par_lanes}"} \
  ${par_libraries:+"--libraries=${par_libraries}"} \
  ${par_feature_ref:+"--feature-ref=${par_feature_ref}"} \
  ${par_expect_cells:+"--expect-cells=${par_expect_cells}"} \
  ${par_force_cells:+"--force-cells=${par_force_cells}"} \
  ${par_create_bam:+"--create-bam=${par_create_bam}"} \
  ${par_no_secondary:+--nosecondary} \
  ${par_r1_length:+"--r1-length=${par_r1_length}"} \
  ${par_r2_length:+"--r2-length=${par_r2_length}"} \
  ${par_include_introns:+--include-introns=${par_include_introns}} \
  ${par_chemistry:+"--chemistry=${par_chemistry}"} \
  ${par_no_libraries:+--no-libraries} \
  ${par_check_library_compatibility:+"--check-library-compatibility=${par_check_library_compatibility}"} \
  ${par_cell_annotation_model:+"--cell-annotation-model=${par_cell_annotation_model}"} \
  ${par_min_crispr_umi:+"--min-crispr-umi=${par_min_crispr_umi}"} \
  ${par_tenx_cloud_token:+"--tenx-cloud-token-path=${par_tenx_cloud_token}"} \
  ${par_dry:+--dry}

echo "> Copying output"
if [ -d "$id/outs/" ]; then
  if [ ! -d "${par_output}" ]; then
    mkdir -p "${par_output}"
  fi
  mv "$id/outs/"* "${par_output}"
fi

exit 0
