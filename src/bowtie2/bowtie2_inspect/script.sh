#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# unset flags
[[ "$par_summary" == "false" ]] && unset par_summary
[[ "$par_names" == "false" ]] && unset par_names
[[ "$par_verbose_inspect" == "false" ]] && unset par_verbose_inspect
[[ "$par_debug" == "false" ]] && unset par_debug
[[ "$par_sanitized" == "false" ]] && unset par_sanitized
[[ "$par_verbose" == "false" ]] && unset par_verbose
[[ "$par_large_index" == "false" ]] && unset par_large_index

# Resolve the index path prefix from the index directory. When --large_index is
# set, only a large (.bt2l) index is considered.
index_dir="${par_index%/}"
if [[ ! -d "$index_dir" ]]; then
  echo "Error: --index must be a directory containing the bowtie2 index files" >&2
  exit 1
fi

# a large index is always a candidate, a small one only when not overridden
index_suffixes=(bt2l)
if [[ -z "$par_large_index" ]]; then
  index_suffixes+=(bt2)
fi

index_files=()
for index_suffix in "${index_suffixes[@]}"; do
  # -L follows symlinks, as a workflow engine may stage the index files as such.
  # The .rev.1 files are excluded, as they belong to the same index.
  mapfile -d '' -t index_files < <(
    find -L "$index_dir" -maxdepth 1 -type f \
      -name "*.1.$index_suffix" ! -name "*.rev.1.$index_suffix" -print0
  )
  [[ "${#index_files[@]}" -gt 0 ]] && break
done

if [[ "${#index_files[@]}" -eq 0 ]]; then
  if [[ -n "$par_large_index" ]]; then
    echo "Error: no large bowtie2 index files (.bt2l) found in '$par_index'." \
      "Omit --large_index to inspect a small (.bt2) index." >&2
  else
    echo "Error: no bowtie2 index files (.bt2 or .bt2l) found in '$par_index'" >&2
  fi
  exit 1
fi

if [[ "${#index_files[@]}" -gt 1 ]]; then
  echo "Error: multiple bowtie2 indices found in '$par_index': ${index_files[*]##*/}." \
    "The index directory must contain exactly one index." >&2
  exit 1
fi

index_path="${index_files[0]%.1.$index_suffix}"

# Build the command arguments
cmd_args=(
    ${par_summary:+-s}
    ${par_names:+-n}
    ${par_across:+-a "$par_across"}
    ${par_verbose_inspect:+-v}
    ${par_debug:+--debug}
    ${par_sanitized:+--sanitized}
    ${par_verbose:+--verbose}
    ${par_large_index:+--large-index}

    # The index must come last, as documented: the bowtie2-inspect wrapper reads
    # the basename from the final argument to decide whether to run the small or
    # the large binary. Passing it earlier silently breaks that detection.
    "$index_path"
)

# Run bowtie2-inspect
if [[ -n "$par_output" ]]; then
  bowtie2-inspect "${cmd_args[@]}" > "$par_output"
else
  bowtie2-inspect "${cmd_args[@]}"
fi
