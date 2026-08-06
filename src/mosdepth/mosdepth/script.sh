#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# Unset flags
[[ "$par_use_median" == "false" ]] && unset par_use_median
[[ "$par_no_per_base" == "false" ]] && unset par_no_per_base
[[ "$par_fast_mode" == "false" ]] && unset par_fast_mode
[[ "$par_fragment_mode" == "false" ]] && unset par_fragment_mode

tmp_dir=$(mktemp -d)
trap 'rm -rf "$tmp_dir"' EXIT

# If an input index is provided, stage it to a temporary directory with the main
# file so they are available alongside each other. If only the main input is
# provided, we don't stage it and assume the index is available in the existing
# location.
if [ -n "$par_input_index" ]; then
  staged_input="$tmp_dir/$(basename "$par_input")"
  ln -s "$(realpath "$par_input")" "$staged_input"
  ln -s "$(realpath "$par_input_index")" "${staged_input}.${par_input_index##*.}"
else
  staged_input="$par_input"
fi

if [ -n "$par_fasta" ]; then
  if [ -n "$par_fasta_index" ]; then
    staged_fasta="$tmp_dir/$(basename "$par_fasta")"
    ln -s "$(realpath "$par_fasta")" "$staged_fasta"
    ln -s "$(realpath "$par_fasta_index")" "${staged_fasta}.${par_fasta_index##*.}"
  else
    staged_fasta="$par_fasta"
  fi
fi

# Custom quantize bin labels are set via MOSDEPTH_Q0..MOSDEPTH_QN env vars
if [ -n "$par_quantize_labels" ]; then
  IFS=',' read -ra quantize_labels <<<"$par_quantize_labels"
  for i in "${!quantize_labels[@]}"; do
    export "MOSDEPTH_Q${i}=${quantize_labels[$i]}"
  done
fi

# Distribution decimal precision is set via the MOSDEPTH_PRECISION env var
[ -n "$par_dist_precision" ] && export MOSDEPTH_PRECISION="$par_dist_precision"

# mosdepth names output files using a fixed prefix. We stage output to a
# temporary directory and rename it at the end of the script.
output_prefix="$tmp_dir/mosdepth"

# Select either --by_bed or --by_window, but not both depending on which is set
if [ -n "$par_by_bed" ] && [ -n "$par_by_window" ]; then
  echo "Error: --by_bed and --by_window cannot be set at the same time." >&2
  exit 1
elif [ -n "$par_by_bed" ]; then
  par_by="$par_by_bed"
elif [ -n "$par_by_window" ]; then
  par_by="$par_by_window"
fi

cmd_args=(
  ${meta_cpus:+--threads "$meta_cpus"}
  ${par_chrom:+--chrom "$par_chrom"}
  ${par_by:+--by "$par_by"}
  ${par_no_per_base:+--no-per-base}
  ${staged_fasta:+--fasta "$staged_fasta"}
  ${par_flag:+--flag "$par_flag"}
  ${par_include_flag:+--include-flag "$par_include_flag"}
  ${par_fast_mode:+--fast-mode}
  ${par_fragment_mode:+--fragment-mode}
  ${par_quantize:+--quantize "$par_quantize"}
  ${par_mapq:+--mapq "$par_mapq"}
  ${par_min_frag_len:+--min-frag-len "$par_min_frag_len"}
  ${par_max_frag_len:+--max-frag-len "$par_max_frag_len"}
  ${par_thresholds:+--thresholds "$par_thresholds"}
  ${par_use_median:+--use-median}
  ${par_read_groups:+--read-groups "$par_read_groups"}
  "$output_prefix"
  "$staged_input"
)

mosdepth "${cmd_args[@]}"

# Move output files to their final destination. The optional third argument
# is the destination for the file's .csi index, if it produces one.
move_output() {
  local src="$1"
  local dest="$2"
  local dest_index="${3:-}"
  if [ -n "$dest" ] && [ -f "$src" ]; then
    mv "$src" "$dest"
    if [ -n "$dest_index" ] && [ -f "${src}.csi" ]; then
      mv "${src}.csi" "$dest_index"
    fi
  fi
}

move_output "${output_prefix}.mosdepth.summary.txt" "$par_output_summary"
move_output "${output_prefix}.mosdepth.global.dist.txt" "$par_output_global_dist"
move_output "${output_prefix}.mosdepth.region.dist.txt" "$par_output_region_dist"
move_output "${output_prefix}.per-base.bed.gz" "$par_output_per_base" "$par_output_per_base_index"
move_output "${output_prefix}.regions.bed.gz" "$par_output_regions" "$par_output_regions_index"
move_output "${output_prefix}.quantized.bed.gz" "$par_output_quantized" "$par_output_quantized_index"
move_output "${output_prefix}.thresholds.bed.gz" "$par_output_thresholds" "$par_output_thresholds_index"
