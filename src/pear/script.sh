#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

[[ "$par_emperical_freqs" == "false" ]] && unset par_emperical_freqs
[[ "$par_nbase" == "false" ]] && unset par_nbase

if [[ "${par_forward_fastq##*.}" == "gz" ]]; then
  gunzip $par_forward_fastq
  par_forward_fastq=${par_forward_fastq%.*}
fi
if [[ "${par_reverse_fastq##*.}" == "gz" ]]; then
  gunzip $par_reverse_fastq
  par_reverse_fastq=${par_reverse_fastq%.*}
fi

output_dir=$(mktemp -d -p "$meta_temp_dir" "pear.XXXXXX")

pear \
  -f "$par_forward_fastq" \
  -r "$par_reverse_fastq" \
  -o "$output_dir" \
  ${par_p_value:+-p "${par_p_value}"} \
  ${par_min_overlap:+-v "${par_min_overlap}"} \
  ${par_max_assembly_length:+-m "${par_max_assembly_length}"} \
  ${par_min_assembly_length:+-n "${par_min_assembly_length}"} \
  ${par_min_trim_length:+-t "${par_min_trim_length}"} \
  ${par_quality_threshold:+-q "${par_quality_threshold}"} \
  ${par_max_uncalled_base:+-u "${par_max_uncalled_base}"} \
  ${par_test_method:+-g "${par_test_method}"} \
  ${par_score_method:+-s "${par_score_method}"} \
  ${par_phred_base:+-b "${par_phred_base}"} \
  ${meta_memory_mb:+--memory "${meta_memory_mb}M"} \
  ${par_cap:+-c "${par_cap}"} \
  ${meta_cpus:+-j "${meta_cpus}"} \
  ${par_emperical_freqs:+-e} \
  ${par_nbase:+-z}


if [[ "${par_assembled##*.}" == "gz" ]]; then
  gzip -9 -c ${output_dir}.assembled.fastq > ${par_assembled}
else
  mv ${output_dir}.assembled.fastq ${par_assembled}
fi

if [[ "${par_unassembled_forward##*.}" == "gz" ]]; then
  gzip -9 -c ${output_dir}.unassembled.forward.fastq > ${par_unassembled_forward}
else
  mv ${output_dir}.unassembled.forward.fastq ${par_unassembled_forward}
fi

if [[ "${par_unassembled_reverse##*.}" == "gz" ]]; then
  gzip -9 -c ${output_dir}.unassembled.reverse.fastq > ${par_unassembled_reverse}
else
  mv ${output_dir}.unassembled.reverse.fastq ${par_unassembled_reverse}
fi

if [[ "${par_discarded##*.}" == "gz" ]]; then
  gzip -9 -c ${output_dir}.discarded.fastq > ${par_discarded}
else
  mv ${output_dir}.discarded.fastq ${par_discarded}
fi
