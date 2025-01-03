#!/bin/bash

set -eo pipefail

unset_if_false=(
    par_override_id
    par_nosecondary
)

for par in ${unset_if_false[@]}; do
    test_val="${!par}"
    [[ "$test_val" == "false" ]] && unset $par
done

spaceranger count \
  ${par_id:+--id="$par_id"} \
  ${par_transcriptome:+--transcriptome="$par_transcriptome"} \
  ${par_fastqs:+--fastqs="$par_fastqs"} \
  ${par_probe_set:+--probe-set="$par_probe_set"} \
  ${par_cytaimage:+--cytaimage="$par_cytaimage"} \
  ${par_image:+--image="$par_image"} \
  ${par_slide:+--slide="$par_slide"} \
  ${par_area:+--area="$par_area"} \
  ${par_unknown_slide:+--unknown-slide="$par_unknown_slide"} \
  ${par_slidefile:+--slidefile="$par_slidefile"} \
  ${par_override_id:+--override-id} \
  ${par_darkimage:+--darkimage="$par_darkimage"} \
  ${par_colorizedimage:+--colorizedimage="$par_colorizedimage"} \
  ${par_dapi_index:+--dapi-index="$par_dapi_index"} \
  ${par_image_scale:+--image-scale="$par_image_scale"} \
  ${par_reorient_images:+--reorient-images="$par_reorient_images"} \
  ${par_create_bam:+--create-bam="$par_create_bam"} \
  ${par_nosecondary:+--nosecondary} \
  ${par_r1_length:+--r1-length="$par_r1_length"} \
  ${par_r2_length:+--r2-length="$par_r2_length"} \
  ${par_filter_probes:+--filter-probes="$par_filter_probes"} \
  ${par_custom_bin_size:+--custom-bin-size="$par_custom_bin_size"} \
  ${par_project:+--project="$par_project"} \
  ${par_sample:+--sample="$par_sample"} \
  ${par_lanes:+--lanes="$par_lanes"} \
  ${meta_cpus:+--localcores="$meta_cpus"} \
  --localmem=16
  # --localcores=4 \
  # ${meta_memory_gb:+--localmem="$meta_memory_gb"}


# spaceranger count \
#   --id=test_count \
#   --transcriptome=src/spaceranger/spaceranger_count/test_data/GRCh38 \
#   --fastqs=src/spaceranger/spaceranger_count/test_data/subsampled_fastqs \
#   --probe-set=src/spaceranger/spaceranger_count/test_data/Visium_HD_Human_Pancreas_probe_set.csv \
#   --cytaimage=src/spaceranger/spaceranger_count/test_data/Visium_HD_Human_Pancreas_image.tif \
#   --image=src/spaceranger/spaceranger_count/test_data/Visium_HD_Human_Pancreas_tissue_image.btf \
#   --create-bam=false \
#   --lanes=1 \
#   --slide=H1-HBNMBMC \
#   --area=D1 \
#   --nosecondary \
#   --localcores=4 \
#   --localmem=16
  # --custom-bin-size 16 \
