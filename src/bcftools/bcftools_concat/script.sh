#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# Unset false boolean parameters
unset_if_false=(
  par_allow_overlaps
  par_compact_ps
  par_remove_duplicates
  par_drop_genotypes
  par_ligate
  par_ligate_force
  par_ligate_warn
  par_no_version
  par_naive
  par_naive_force
)

for par in ${unset_if_false[@]}; do
  test_val="${!par}"
  [[ "$test_val" == "false" ]] && unset $par
done

# Check that either input files or file_list is provided
if [[ -z "${par_input}" && -z "${par_file_list}" ]]; then
  echo "Error: One of the parameters '--input' or '--file_list' must be used."
  exit 1
fi

# Handle multiple input files (semicolon-separated from Viash)
if [[ -n "$par_input" ]]; then
  IFS=';' read -ra input_files <<< "$par_input"
fi

# Build command array
cmd_args=(
  bcftools concat
  ${par_allow_overlaps:+--allow-overlaps}
  ${par_compact_ps:+--compact-PS}
  ${par_rm_dups:+--rm-dups "$par_rm_dups"}
  ${par_remove_duplicates:+--remove-duplicates}
  ${par_drop_genotypes:+--drop-genotypes}
  ${par_ligate:+--ligate}
  ${par_ligate_force:+--ligate-force}
  ${par_ligate_warn:+--ligate-warn}
  ${par_no_version:+--no-version}
  ${par_naive:+--naive}
  ${par_naive_force:+--naive-force}
  ${par_output_type:+--output-type "$par_output_type"}
  ${par_min_pq:+--min-PQ "$par_min_pq"}
  ${par_regions:+--regions "$par_regions"}
  ${par_regions_file:+--regions-file "$par_regions_file"}
  ${par_regions_overlap:+--regions-overlap "$par_regions_overlap"}
  ${meta_cpus:+--threads "$meta_cpus"}
  ${par_verbosity:+--verbosity "$par_verbosity"}
  ${par_write_index:+--write-index="$par_write_index"}
  ${par_output:+--output "$par_output"}
  ${par_file_list:+--file-list "$par_file_list"}
)

# Add input files to command array
if [[ -n "$par_input" ]]; then
  cmd_args+=("${input_files[@]}")
fi

# Execute command
"${cmd_args[@]}"