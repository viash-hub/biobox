#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# Unset false boolean parameters (biobox standard)
[[ "$par_vcf" == "false" ]] && unset par_vcf
[[ "$par_cache" == "false" ]] && unset par_cache
[[ "$par_offline" == "false" ]] && unset par_offline
[[ "$par_everything" == "false" ]] && unset par_everything
[[ "$par_canonical" == "false" ]] && unset par_canonical
[[ "$par_ccds" == "false" ]] && unset par_ccds
[[ "$par_protein" == "false" ]] && unset par_protein
[[ "$par_symbol" == "false" ]] && unset par_symbol
[[ "$par_hgvs" == "false" ]] && unset par_hgvs
[[ "$par_af" == "false" ]] && unset par_af
[[ "$par_af_1kg" == "false" ]] && unset par_af_1kg
[[ "$par_af_gnomad" == "false" ]] && unset par_af_gnomad
[[ "$par_max_af" == "false" ]] && unset par_max_af
[[ "$par_pick" == "false" ]] && unset par_pick
[[ "$par_pick_allele" == "false" ]] && unset par_pick_allele
[[ "$par_flag_pick" == "false" ]] && unset par_flag_pick
[[ "$par_per_gene" == "false" ]] && unset par_per_gene
[[ "$par_most_severe" == "false" ]] && unset par_most_severe
[[ "$par_summary" == "false" ]] && unset par_summary
[[ "$par_filter_common" == "false" ]] && unset par_filter_common
[[ "$par_no_check_variants_order" == "false" ]] && unset par_no_check_variants_order
[[ "$par_allow_non_variant" == "false" ]] && unset par_allow_non_variant

# Build command array (preferred pattern)
cmd_args=(
  vep
  --input_file "$par_input_file"
  --output_file "$par_output_file"
  ${par_format:+--format "$par_format"}
  ${par_vcf:+--vcf}
  ${par_species:+--species "$par_species"}
  ${par_assembly:+--assembly "$par_assembly"}
  ${par_cache:+--cache}
  ${par_dir:+--dir "$par_dir"}
  ${par_cache_version:+--cache_version "$par_cache_version"}
  ${par_offline:+--offline}
  ${par_everything:+--everything}
  ${par_canonical:+--canonical}
  ${par_ccds:+--ccds}
  ${par_protein:+--protein}
  ${par_symbol:+--symbol}
  ${par_hgvs:+--hgvs}
  ${par_sift:+--sift "$par_sift"}
  ${par_polyphen:+--polyphen "$par_polyphen"}
  ${par_af:+--af}
  ${par_af_1kg:+--af_1kg}
  ${par_af_gnomad:+--af_gnomad}
  ${par_max_af:+--max_af}
  ${par_pick:+--pick}
  ${par_pick_allele:+--pick_allele}
  ${par_flag_pick:+--flag_pick}
  ${par_per_gene:+--per_gene}
  ${par_pick_order:+--pick_order "$par_pick_order"}
  ${par_most_severe:+--most_severe}
  ${par_summary:+--summary}
  ${par_filter_common:+--filter_common}
  ${par_buffer_size:+--buffer_size "$par_buffer_size"}
  ${par_no_check_variants_order:+--no_check_variants_order}
  ${par_allow_non_variant:+--allow_non_variant}
  ${meta_cpus:+--fork "$meta_cpus"}
)

# Execute command
"${cmd_args[@]}"