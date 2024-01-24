#!/bin/bash

## VIASH START
## VIASh END

[["$par_emperical_freqs" == "false"]] && unset par_emperical_freqs
[["$par_nbase" == "false"]] && unset par_nbase

pear \
  -f "$par_forward" \
  -r "$par_reverse" \
  -o "$par_output" \
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
  ${par_memory:+-y "${par_memory}"} \
  ${par_cap:+-c "${par_cap}"} \
  ${par_threads:+-j "${par_threads}"} \
  ${par_emperical_freqs:+-e} \
  ${par_nbase:+-z}
