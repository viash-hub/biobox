#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

[[ "$par_flat_output" == "false" ]] && unset par_flat_output
[[ "$par_paired" == "false" ]] && unset par_paired
[[ "$par_skip_dup_check" == "false" ]] && unset par_skip_dup_check
[[ "$par_quiet" == "false" ]] && unset par_quiet
[[ "$par_verbose" == "false" ]] && unset par_verbose
[[ "$par_skip_tin" == "false" ]] && unset par_skip_tin
[[ "$par_skip_read_duplication" == "false" ]] && unset par_skip_read_duplication
[[ "$par_skip_preseq" == "false" ]] && unset par_skip_preseq

IFS=';' read -ra input_files <<< "$par_input"

cmd_args=(
  # Input / Output:
  --gtf "$par_gtf"
  ${par_reference:+--reference "$par_reference"}
  --outdir "$par_outdir"
  ${par_sample_name:+--sample-name "$par_sample_name"}
  ${par_flat_output:+--flat-output}
  ${par_config:+--config "$par_config"}
  ${par_json_summary:+--json-summary "$par_json_summary"}
  # Library:
  ${par_stranded:+--stranded "$par_stranded"}
  ${par_paired:+--paired}
  # General:
  ${meta_cpus:+--threads "$meta_cpus"}
  ${par_mapq:+--mapq "$par_mapq"}
  ${par_biotype_attribute:+--biotype-attribute "$par_biotype_attribute"}
  ${par_skip_dup_check:+--skip-dup-check}
  ${par_quiet:+--quiet}
  ${par_verbose:+--verbose}
  # Tool parameters:
  ${par_infer_experiment_sample_size:+--infer-experiment-sample-size "$par_infer_experiment_sample_size"}
  ${par_min_intron:+--min-intron "$par_min_intron"}
  ${par_junction_saturation_seed:+--junction-saturation-seed "$par_junction_saturation_seed"}
  ${par_junction_saturation_min_coverage:+--junction-saturation-min-coverage "$par_junction_saturation_min_coverage"}
  ${par_junction_saturation_percentile_floor:+--junction-saturation-percentile-floor "$par_junction_saturation_percentile_floor"}
  ${par_junction_saturation_percentile_ceiling:+--junction-saturation-percentile-ceiling "$par_junction_saturation_percentile_ceiling"}
  ${par_junction_saturation_percentile_step:+--junction-saturation-percentile-step "$par_junction_saturation_percentile_step"}
  ${par_inner_distance_sample_size:+--inner-distance-sample-size "$par_inner_distance_sample_size"}
  ${par_inner_distance_lower_bound:+--inner-distance-lower-bound "$par_inner_distance_lower_bound"}
  ${par_inner_distance_upper_bound:+--inner-distance-upper-bound "$par_inner_distance_upper_bound"}
  ${par_inner_distance_step:+--inner-distance-step "$par_inner_distance_step"}
  ${par_tin_seed:+--tin-seed "$par_tin_seed"}
  ${par_skip_tin:+--skip-tin}
  ${par_skip_read_duplication:+--skip-read-duplication}
  ${par_skip_preseq:+--skip-preseq}
  ${par_preseq_seed:+--preseq-seed "$par_preseq_seed"}
  ${par_preseq_max_extrap:+--preseq-max-extrap "$par_preseq_max_extrap"}
  ${par_preseq_step_size:+--preseq-step-size "$par_preseq_step_size"}
  ${par_preseq_n_bootstraps:+--preseq-n-bootstraps "$par_preseq_n_bootstraps"}
  ${par_preseq_seg_len:+--preseq-seg-len "$par_preseq_seg_len"}
  "${input_files[@]}"
)

rustqc rna "${cmd_args[@]}"
