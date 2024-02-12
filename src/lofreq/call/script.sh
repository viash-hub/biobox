#!/bin/bash

## VIASH START
## VIASH END

# Unset all parameters that are set to "false"
[[ "$par_no_baq" == "false" ]] && unset par_no_baq
[[ "$par_no_idaq" == "false" ]] && unset par_no_idaq
[[ "$par_del_baq" == "false" ]] && unset par_del_baq
[[ "$par_no_ext_baq" == "false" ]] && unset par_no_ext_baq
[[ "$par_no_mq" == "false" ]] && unset par_no_mq
[[ "$par_call_indels" == "false" ]] && unset par_call_indels
[[ "$par_only_indels" == "false" ]] && unset par_only_indels
[[ "$par_src_qual" == "false" ]] && unset par_src_qual
[[ "$par_illumina_13" == "false" ]] && unset par_illumina_13
[[ "$par_use_orphan" == "false" ]] && unset par_use_orphan
[[ "$par_plp_summary_only" == "false" ]] && unset par_plp_summary_only
[[ "$par_no_default_filter" == "false" ]] && unset par_no_default_filter
[[ "$par_force_overwrite" == "false" ]] && unset par_force_overwrite
[[ "$par_verbose" == "false" ]] && unset par_verbose
[[ "$par_debug" == "false" ]] && unset par_debug

# Run lofreq call
lofreq call \
  -f "$par_ref" \
  -o "$par_out" \
  ${par_region:+-r "${par_region}"} \
  ${par_bed:+-l "${par_bed}"} \
  ${par_min_bq:+-q "${par_min_bq}"} \
  ${par_min_alt_bq:+-Q "${par_min_alt_bq}"} \
  ${par_def_alt_bq:+-R "${par_def_alt_bq}"} \
  ${par_min_jq:+-j "${par_min_jq}"} \
  ${par_alt_jq:+-K "${par_alt_jq}"} \
  ${par_no_baq:+-B} \
  ${par_no_idaq:+-A} \
  ${par_del_baq:+-D} \
  ${par_no_ext_baq:+-e} \
  ${par_min_mq:+-m "${par_min_mq}"} \
  ${par_max_mq:+-M "${par_max_mq}"} \
  ${par_no_mq:+-N} \
  ${par_call_indels:+--call-indels} \
  ${par_only_indels:+--only-indels} \
  ${par_src_qual:+-s} \
  ${par_ign_vcf:+-S "${par_ign_vcf}"} \
  ${par_def_nm_q:+-T "${par_def_nm_q}"} \
  ${par_sig:+-a "${par_sig}"} \
  ${par_bonf:+-b "${par_bonf}"} \
  ${par_min_cov:+-C "${par_min_cov}"} \
  ${par_max_depth:+-d "${par_max_depth}"} \
  ${par_illumina_13:+--illumina-1.3} \
  ${par_use_orphan:+--use-orphan} \
  ${par_plp_summary_only:+--plp-summary-only} \
  ${par_no_default_filter:+--no-default-filter} \
  ${par_force_overwrite:+--force-overwrite} \
  ${par_verbose:+--verbose} \
  ${par_debug:+--debug} \
  "$par_input"