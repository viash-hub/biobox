#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# Unset boolean flags that are false
[[ "$par_same_strand" == "false" ]] && unset par_same_strand
[[ "$par_opposite_strand" == "false" ]] && unset par_opposite_strand
[[ "$par_use_names" == "false" ]] && unset par_use_names
[[ "$par_use_scores" == "false" ]] && unset par_use_scores
[[ "$par_use_intervals" == "false" ]] && unset par_use_intervals

# Convert semicolon-separated files to array
IFS=';' read -ra files_array <<< "$par_files"

# Convert semicolon-separated labels to array if provided
if [ -n "${par_labels:-}" ]; then
  IFS=';' read -ra labels_array <<< "$par_labels"
fi

# Validate that if labels are provided, files and labels have the same length
# Labels are required unless using names/scores/intervals options
if [ -n "${par_labels:-}" ]; then
  if [ ${#files_array[@]} -ne ${#labels_array[@]} ]; then
    echo "Error: Number of files (${#files_array[@]}) must match number of labels (${#labels_array[@]})" >&2
    exit 1
  fi
elif [ -z "${par_use_names:-}" ] && [ -z "${par_use_scores:-}" ] && [ -z "${par_use_intervals:-}" ]; then
  echo "Error: Must provide --labels unless using --use_names, --use_scores, or --use_intervals" >&2
  exit 1
fi

# Build command arguments array
cmd_args=(
  -i "$par_input"
  -files "${files_array[@]}"
  ${par_labels:+-labels "${labels_array[@]}"}
  ${par_min_overlap:+-f "$par_min_overlap"}
  ${par_same_strand:+-s}
  ${par_opposite_strand:+-S}
  ${par_tag_name:+-tag "$par_tag_name"}
  ${par_use_names:+-names}
  ${par_use_scores:+-scores}
  ${par_use_intervals:+-intervals}
)

# Execute bedtools tag and redirect output to the specified output file
bedtools tag "${cmd_args[@]}" > "$par_output"
