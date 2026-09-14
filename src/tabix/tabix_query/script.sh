#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

[[ "$par_print_header" == "false" ]] && unset par_print_header
[[ "$par_only_header" == "false" ]] && unset par_only_header
[[ "$par_list_chroms" == "false" ]] && unset par_list_chroms
[[ "$par_unique" == "false" ]] && unset par_unique
[[ "$par_separate_regions" == "false" ]] && unset par_separate_regions
[[ "$par_skip_index_download" == "false" ]] && unset par_skip_index_download

# Stage the input and index together so tabix can find them
work_dir=$(mktemp -d "$meta_temp_dir/tabix_query.XXXXXX")
trap 'rm -rf "$work_dir"' EXIT

input_name=$(basename "$par_input")
staged_input="$work_dir/$input_name"
ln -s "$(realpath "$par_input")" "$staged_input"

index_ext="${par_index##*.}"
ln -s "$(realpath "$par_index")" "${staged_input}.${index_ext}"

staged_regions_file=""
if [[ -n "${par_regions_file:-}" ]]; then
  staged_regions_file="$work_dir/$(basename "$par_regions_file")"
  ln -s "$(realpath "$par_regions_file")" "$staged_regions_file"
fi

staged_targets_file=""
if [[ -n "${par_targets_file:-}" ]]; then
  staged_targets_file="$work_dir/$(basename "$par_targets_file")"
  ln -s "$(realpath "$par_targets_file")" "$staged_targets_file"
fi

cmd_args=(
  tabix
  ${par_print_header:+--print-header}
  ${par_only_header:+--only-header}
  ${par_list_chroms:+--list-chroms}
  ${par_unique:+--unique}
  ${par_separate_regions:+--separate-regions}
  ${par_cache:+--cache "$par_cache"}
  ${par_verbosity:+--verbosity "$par_verbosity"}
  ${par_skip_index_download:+-D}
  ${meta_cpus:+--threads "$meta_cpus"}
  ${staged_regions_file:+--regions "$staged_regions_file"}
  ${staged_targets_file:+--targets "$staged_targets_file"}
  "$staged_input"
)

if [[ -n "${par_region:-}" ]]; then
  IFS=';' read -ra region_array <<<"$par_region"
  cmd_args+=("${region_array[@]}")
fi

"${cmd_args[@]}" >"$par_output"
