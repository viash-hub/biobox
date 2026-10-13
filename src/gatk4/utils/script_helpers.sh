#!/bin/bash

# GATK4-specific runtime script helper functions for biobox components
#
# Source this file (alongside the component's own script.sh) via:
#   source "$meta_resources_dir/gatk4/script_helpers.sh"

# Split a semicolon-separated `multiple: true` argument value into a bash
# array of repeated "--flag value" pairs, e.g. "a;b;c" with flag
# "--read-filter" becomes ("--read-filter" "a" "--read-filter" "b"
# "--read-filter" "c"). Produces an empty array if the value is unset/empty.
#
# Usage: split_multiple_to_flags "$par_value" "--flag-name" result_array_name
split_multiple_to_flags() {
  local value="$1"
  local flag="$2"
  local -n result_ref="$3"

  result_ref=()
  if [[ -n "$value" ]]; then
    local items=()
    IFS=';' read -ra items <<< "$value"
    local item
    for item in "${items[@]}"; do
      result_ref+=("$flag" "$item")
    done
  fi
}

# Return success if a boolean argument value is set and not "false"
#
# Usage: if is_true "$par_flag"; then ...; fi
is_true() {
  [[ -n "${1:-}" && "$1" != "false" ]]
}

# Unset `boolean_true` argument variables that are "false", so that
# `${par_x:+--flag}` only adds the flag when the argument is set.
# Only use it for `boolean_true` arguments.
#
# Usage: unset_if_false par_a par_b ...
unset_if_false() {
  local par
  for par in "$@"; do
    if [[ "${!par:-}" == "false" ]]; then
      unset "$par"
    fi
  done
}

# Create a temp dir for the component in $meta_temp_dir and remove it when the
# script exits. Sets the `tmp_dir` variable. Call it directly, not in $(...),
# so that the trap is set in the main shell.
#
# Usage: setup_tmp_dir
setup_tmp_dir() {
  tmp_dir=$(mktemp -d "$meta_temp_dir/${meta_name}.XXXXXX") || exit 1
  trap 'rm -rf "$tmp_dir"' EXIT
}

# Run a GATK (or Picard) tool with the JVM heap set to 80% of the memory given
# to the task (fallback 3072 MB, as in the nf-core GATK4 modules). The other 20%
# is left for memory outside the heap.
#
# Usage: run_gatk Tool [arguments...]
run_gatk() {
  local tool="$1"
  shift
  local avail_mem_mb=$(( ${meta_memory_mb:-3072} * 8 / 10 ))
  gatk --java-options "-Xmx${avail_mem_mb}M -XX:-UsePerfData" "$tool" "$@"
}

# Stop with an error if an output that is only written on request is requested
# (`--create_*` is set) but its path is not given.
#
# Usage: require_output_if_created "$par_create_x" "$par_x" --create_x --x
require_output_if_created() {
  local create="$1"
  local output="$2"
  local create_arg="$3"
  local output_arg="$4"

  if is_true "$create" && [[ -z "$output" ]]; then
    echo "Error: $create_arg requires $output_arg." >&2
    exit 1
  fi
}

# Stop with an error if an index output is not given although GATK will write
# the index. For the GATK index options (`--create_output_variant_index`,
# `--create_output_bam_index`), the index is written unless the option is
# "false".
#
# Usage: require_output_index "$par_output_index" "$par_create_output_variant_index" \
#          --output_index --create_output_variant_index
require_output_index() {
  local output_index="$1"
  local create_index="$2"
  local index_arg="$3"
  local create_arg="$4"

  if [[ "$create_index" != "false" && -z "$output_index" ]]; then
    echo "Error: $index_arg is required unless $create_arg is false." >&2
    exit 1
  fi
}

# Move an index that GATK wrote next to an output to the path of the index
# output argument. Does nothing if no index path is given or if both paths are
# the same file.
#
# Usage: move_output_index "$gatk_index" "$par_output_index"
move_output_index() {
  local gatk_index="$1"
  local output_index="$2"

  if [[ -n "$gatk_index" && -n "$output_index" && ! "$gatk_index" -ef "$output_index" ]]; then
    mv "$gatk_index" "$output_index"
  fi
}

# Stop with an error if a CRAM file is used and no reference is given.
#
# Usage: require_reference_for_cram "$par_reference" --reference "$par_input" "$par_output" ...
require_reference_for_cram() {
  local reference="$1"
  local reference_arg="$2"
  shift 2

  local value file files
  for value in "$@"; do
    IFS=';' read -ra files <<< "$value"
    for file in "${files[@]}"; do
      if [[ "$file" == *.cram && -z "$reference" ]]; then
        echo "Error: $reference_arg is required for CRAM input or output." >&2
        exit 1
      fi
    done
  done
}

# Symlink a reference FASTA and its .fai/.dict into a temp dir under matching
# basenames, so GATK can find them
#
# Usage: staged=$(stage_reference_trio "$tmp_dir" "$par_reference" "$par_reference_fai" "$par_reference_dict")
stage_reference_trio() {
  local tmp_dir="$1"
  local reference="$2"
  local reference_fai="$3"
  local reference_dict="$4"

  if [[ -z "$reference_fai" || -z "$reference_dict" ]]; then
    echo "Error: --reference_fai and --reference_dict must both be provided when --reference is set." >&2
    exit 1
  fi

  # Create links, explicitly exit on failure
  ln -s "$(readlink -f "$reference")" "$tmp_dir/reference.fasta" || exit 1
  ln -s "$(readlink -f "$reference_fai")" "$tmp_dir/reference.fasta.fai" || exit 1
  ln -s "$(readlink -f "$reference_dict")" "$tmp_dir/reference.dict" || exit 1

  echo "$tmp_dir/reference.fasta"
}

# Stage an optional reference trio with stage_reference_trio and build the
# GATK flags for it. The result is empty if no reference is given. Stops with
# an error if the .fai or .dict is given without the reference.
#
# Usage: stage_optional_reference "$tmp_dir" "$par_reference" "$par_reference_fai" \
#          "$par_reference_dict" result_array_name
stage_optional_reference() {
  local tmp_dir="$1"
  local reference="$2"
  local reference_fai="$3"
  local reference_dict="$4"
  local -n reference_args_ref="$5"

  reference_args_ref=()
  if [[ -z "$reference" ]]; then
    if [[ -n "$reference_fai" || -n "$reference_dict" ]]; then
      echo "Error: --reference_fai and --reference_dict require --reference." >&2
      exit 1
    fi
    return 0
  fi
  local staged_reference
  staged_reference=$(stage_reference_trio "$tmp_dir" "$reference" "$reference_fai" "$reference_dict") || exit 1
  reference_args_ref=(--reference "$staged_reference")
}

# Build the GATK flags for one or more read files (BAM or CRAM) and their
# indexes. The indexes are passed with --read-index, so they can have any name
# and do not need to be staged next to the read files. GATK matches them to the
# read files by order, and a wrong order gives an unclear error, so this stops
# with an error if the number of read files and indexes is different.
#
# Usage: gatk_read_args "$par_input" "$par_read_index" result_array_name
gatk_read_args() {
  local input_files="$1"
  local index_files="$2"
  local -n result_ref="$3"

  local inputs=()
  local indexes=()
  IFS=';' read -ra inputs <<< "$input_files"
  IFS=';' read -ra indexes <<< "$index_files"

  if [[ ${#inputs[@]} -ne ${#indexes[@]} ]]; then
    echo "Error: --read_index must be given once for each file in --input (got ${#inputs[@]} read file(s) and ${#indexes[@]} index file(s))." >&2
    exit 1
  fi

  result_ref=()
  local file
  for file in "${inputs[@]}"; do
    result_ref+=(--input "$file")
  done
  for file in "${indexes[@]}"; do
    result_ref+=(--read-index "$file")
  done
}

# Symlink a reference FASTA and its .fai into a temp dir under matching
# basenames. For tools that do not need the .dict (for example Picard tools,
# which use the reference only for CRAM files).
#
# Usage: staged=$(stage_reference_fai "$tmp_dir" "$par_reference" "$par_reference_fai" arg_name)
stage_reference_fai() {
  local tmp_dir="$1"
  local reference="$2"
  local reference_fai="$3"
  local arg_name="$4"

  if [[ -z "$reference_fai" ]]; then
    echo "Error: ${arg_name}_fai must be provided when $arg_name is set." >&2
    exit 1
  fi

  # Create links, explicitly exit on failure
  ln -s "$(readlink -f "$reference")" "$tmp_dir/reference.fasta" || exit 1
  ln -s "$(readlink -f "$reference_fai")" "$tmp_dir/reference.fasta.fai" || exit 1

  echo "$tmp_dir/reference.fasta"
}

# Stage a VCF (or other feature file) into a temp dir so that GATK can find
# its index. If an index is given, the file and the index are symlinked under
# matching names (`.tbi` for tabix indexes of `.vcf.gz`/`.vcf.bgz` files,
# `.idx` for Tribble indexes of plain `.vcf` files). If no index is given, the
# file is copied and indexed with IndexFeatureFile.
#
# For some input, GATK can read an uncompressed VCF without an index. For these,
# set index_optional to "true" and the file will not be copied and indexed. Only
# set it when GATK does not need random access to the input (check this in the container).
#
# Usage: staged=$(stage_vcf_with_index "$tmp_dir" "$par_vcf" "$par_vcf_index" prefix [index_optional])
stage_vcf_with_index() {
  local tmp_dir="$1"
  local vcf="$2"
  local vcf_index="$3"
  local prefix="$4"
  local index_optional="${5:-false}"
  local staged_vcf="$tmp_dir/${prefix}_$(basename "$vcf")"

  if [[ -z "$vcf_index" && "$index_optional" == "true" && "$vcf" != *.gz && "$vcf" != *.bgz ]]; then
    echo "$vcf"
    return
  fi

  if [[ -n "$vcf_index" ]]; then
    local index_ext
    case "$vcf_index" in
      *.tbi) index_ext="tbi" ;;
      *.idx) index_ext="idx" ;;
      *)
        echo "Error: index file '$vcf_index' for '$vcf' must have a .tbi or .idx extension." >&2
        exit 1
        ;;
    esac
    # Create links, explicitly exit on failure
    ln -s "$(readlink -f "$vcf")" "$staged_vcf" || exit 1
    ln -s "$(readlink -f "$vcf_index")" "${staged_vcf}.${index_ext}" || exit 1
  else
    echo "Warning: no index was provided for '$vcf'. Copying and indexing it, which can be slow for a large file. Provide an index to skip this step." >&2
    cp "$(readlink -f "$vcf")" "$staged_vcf" || exit 1
    if ! run_gatk IndexFeatureFile --input "$staged_vcf" --verbosity ERROR >&2; then
      echo "Error: could not index '$vcf'. A .vcf.gz file must be compressed with bgzip." >&2
      exit 1
    fi
  fi

  echo "$staged_vcf"
}

# Stage an optional VCF (or other feature file) input with stage_vcf_with_index
# and add the GATK flag for it to an array. The GATK flag is the argument name
# with dashes (for example `--dbsnp` or `--population-callset`). Stops with an
# error if the index is given without the input.
#
# Usage: stage_optional_vcf "$tmp_dir" "$par_x" "$par_x_index" --x result_array_name [index_optional]
stage_optional_vcf() {
  local tmp_dir="$1"
  local vcf="$2"
  local vcf_index="$3"
  local arg_name="$4"
  local -n vcf_args_ref="$5"
  local index_optional="${6:-false}"

  if [[ -z "$vcf" ]]; then
    if [[ -n "$vcf_index" ]]; then
      echo "Error: ${arg_name}_index requires $arg_name." >&2
      exit 1
    fi
    return 0
  fi
  local staged_vcf
  staged_vcf=$(stage_vcf_with_index "$tmp_dir" "$vcf" "$vcf_index" "${arg_name#--}" "$index_optional") || exit 1
  vcf_args_ref+=("${arg_name//_/-}" "$staged_vcf")
}

# Stage one or more VCFs given as a semicolon-separated `multiple: true`
# argument value with stage_vcf_with_index, and build the repeated GATK flags
# for them. The index files are matched by order. If index files are given,
# there must be one for each file. If no index files are given, every file is
# copied and indexed, unless index_optional is "true"
# (see stage_vcf_with_index).
#
# Usage: stage_vcfs_with_indexes "$tmp_dir" "$par_vcfs" "$par_vcfs_index" \
#          arg_name gatk_flag result_array_name [index_optional]
stage_vcfs_with_indexes() {
  local tmp_dir="$1"
  local vcf_files="$2"
  local vcf_indexes="$3"
  local arg_name="$4"
  local gatk_flag="$5"
  local -n result_ref="$6"
  local index_optional="${7:-false}"

  local vcfs=()
  local indexes=()
  IFS=';' read -ra vcfs <<< "$vcf_files"
  IFS=';' read -ra indexes <<< "$vcf_indexes"

  if [[ ${#indexes[@]} -gt 0 && ${#indexes[@]} -ne ${#vcfs[@]} ]]; then
    echo "Error: ${arg_name}_index must be given once for each file in $arg_name (got ${#indexes[@]} index file(s) and ${#vcfs[@]} file(s))." >&2
    exit 1
  fi

  result_ref=()
  local prefix="${arg_name#--}"
  local i
  for i in "${!vcfs[@]}"; do
    local staged_vcf
    staged_vcf=$(stage_vcf_with_index "$tmp_dir" "${vcfs[$i]}" "${indexes[$i]:-}" "${prefix}_$((i + 1))" "$index_optional") || exit 1
    result_ref+=("$gatk_flag" "$staged_vcf")
  done
}

# Stage one or more interval files (BED, interval_list, VCF) and build the
# repeated GATK interval flags for them. VCF files are staged with
# stage_vcf_with_index because GATK needs an index for a block-compressed VCF
# (a plain VCF works without one, so it is not copied).
# The index files are matched by order to the VCF files only, other interval
# files are passed directly. If index files are given, there must be one for
# each VCF file.
#
# Usage: stage_interval_files "$tmp_dir" "$par_intervals" "$par_intervals_index" \
#          arg_name gatk_flag result_array_name
stage_interval_files() {
  local tmp_dir="$1"
  local interval_files="$2"
  local interval_indexes="$3"
  local arg_name="$4"
  local gatk_flag="$5"
  local -n result_ref="$6"

  local files=()
  local indexes=()
  IFS=';' read -ra files <<< "$interval_files"
  IFS=';' read -ra indexes <<< "$interval_indexes"

  local vcf_count=0
  local file
  for file in "${files[@]}"; do
    [[ "$file" == *.vcf || "$file" == *.vcf.gz || "$file" == *.vcf.bgz ]] && vcf_count=$((vcf_count + 1))
  done
  if [[ ${#indexes[@]} -gt 0 && ${#indexes[@]} -ne $vcf_count ]]; then
    echo "Error: ${arg_name}_index must be given once for each VCF file in $arg_name (got ${#indexes[@]} index file(s) and $vcf_count VCF file(s))." >&2
    exit 1
  fi

  result_ref=()
  local prefix="${arg_name#--}"
  local vcf_i=0
  local i
  for i in "${!files[@]}"; do
    file="${files[$i]}"
    if [[ "$file" == *.vcf || "$file" == *.vcf.gz || "$file" == *.vcf.bgz ]]; then
      local staged_file
      staged_file=$(stage_vcf_with_index "$tmp_dir" "$file" "${indexes[$vcf_i]:-}" "${prefix}_$((i + 1))" true) || exit 1
      result_ref+=("$gatk_flag" "$staged_file")
      vcf_i=$((vcf_i + 1))
    else
      result_ref+=("$gatk_flag" "$file")
    fi
  done
}

# Stage the interval files and build the GATK flags for all the arguments in
# the shared "GATK Interval Options" group (interval_argument_group.yaml). The
# argument values are read from the `par_*` variables of that group.
#
# Usage: stage_interval_arguments "$tmp_dir" result_array_name
stage_interval_arguments() {
  local tmp_dir="$1"
  local -n interval_args_ref="$2"

  local intervals_args=()
  local exclude_intervals_args=()
  local interval_strings_args=()
  local exclude_interval_strings_args=()
  stage_interval_files "$tmp_dir" "$par_intervals" "$par_intervals_index" \
    --intervals --intervals intervals_args
  stage_interval_files "$tmp_dir" "$par_exclude_intervals" "$par_exclude_intervals_index" \
    --exclude_intervals --exclude-intervals exclude_intervals_args
  split_multiple_to_flags "$par_interval_strings" "--intervals" interval_strings_args
  split_multiple_to_flags "$par_exclude_interval_strings" "--exclude-intervals" exclude_interval_strings_args

  # GATK accepts a mixture of interval files and strings, files are given first
  interval_args_ref=(
    "${intervals_args[@]}"
    "${interval_strings_args[@]}"
    "${exclude_intervals_args[@]}"
    "${exclude_interval_strings_args[@]}"
    ${par_interval_set_rule:+--interval-set-rule "$par_interval_set_rule"}
    ${par_interval_merging_rule:+--interval-merging-rule "$par_interval_merging_rule"}
    ${par_interval_padding:+--interval-padding "$par_interval_padding"}
    ${par_interval_exclusion_padding:+--interval-exclusion-padding "$par_interval_exclusion_padding"}
  )
}

# Build the GATK flags for the shared "GATK Engine Options" group
# (engine_argument_group.yaml) and the temp dir. The argument values are read
# from the `par_*` variables of that group.
#
# Usage: gatk_engine_args "$tmp_dir" result_array_name
gatk_engine_args() {
  local tmp_dir="$1"
  local -n engine_args_ref="$2"

  engine_args_ref=(--tmp-dir "$tmp_dir")
  if is_true "${par_disable_sequence_dictionary_validation:-}"; then
    engine_args_ref+=(--disable-sequence-dictionary-validation)
  fi
  if is_true "${par_lenient:-}"; then
    engine_args_ref+=(--lenient)
  fi
}

# Build the GATK flags for the shared "GATK VCF Output Options" group
# (vcf_output_argument_group.yaml).
#
# Usage: gatk_vcf_output_args result_array_name
gatk_vcf_output_args() {
  local -n vcf_output_args_ref="$1"

  vcf_output_args_ref=(
    ${par_create_output_variant_index:+--create-output-variant-index "$par_create_output_variant_index"}
  )
  if is_true "${par_sites_only_vcf_output:-}"; then
    vcf_output_args_ref+=(--sites-only-vcf-output)
  fi
}

# Build the GATK flags for the shared "GATK BAM Output Options" group
# (bam_output_argument_group.yaml).
#
# Usage: gatk_bam_output_args result_array_name
gatk_bam_output_args() {
  local -n bam_output_args_ref="$1"

  bam_output_args_ref=(
    ${par_create_output_bam_index:+--create-output-bam-index "$par_create_output_bam_index"}
    ${par_output_cram_version:+--output-cram-version "$par_output_cram_version"}
  )
}

# Build the GATK flags for all the arguments in the shared "GATK Read Filter
# Options" group (read_argument_group.yaml). The argument values are read from
# the `par_*` variables of that group.
#
# Usage: gatk_read_filter_args result_array_name
gatk_read_filter_args() {
  local -n read_filter_args_ref="$1"

  local read_filter_flags=()
  local disable_read_filter_flags=()
  local inverted_read_filter_flags=()
  local sample_flags=()
  local library_flags=()
  local read_group_black_list_flags=()
  split_multiple_to_flags "$par_read_filter" "--read-filter" read_filter_flags
  split_multiple_to_flags "$par_disable_read_filter" "--disable-read-filter" disable_read_filter_flags
  split_multiple_to_flags "$par_inverted_read_filter" "--inverted-read-filter" inverted_read_filter_flags
  split_multiple_to_flags "$par_sample" "--sample" sample_flags
  split_multiple_to_flags "$par_library" "--library" library_flags
  split_multiple_to_flags "$par_read_group_black_list" "--read-group-black-list" read_group_black_list_flags

  read_filter_args_ref=(
    "${read_filter_flags[@]}"
    "${disable_read_filter_flags[@]}"
    "${inverted_read_filter_flags[@]}"
    ${par_read_validation_stringency:+--read-validation-stringency "$par_read_validation_stringency"}
    ${par_keep_read_group:+--keep-read-group "$par_keep_read_group"}
    "${sample_flags[@]}"
    "${library_flags[@]}"
    "${read_group_black_list_flags[@]}"
  )
  if is_true "${par_disable_tool_default_read_filters:-}"; then
    read_filter_args_ref+=(--disable-tool-default-read-filters)
  fi
}

# Print "true" if any of the shared interval arguments select intervals, and
# "false" otherwise. GATK needs random access to the main input when it
# traverses by intervals, so an unindexed plain VCF must be indexed.
#
# Usage: uses_intervals=$(gatk_uses_intervals)
gatk_uses_intervals() {
  if [[ -n "$par_intervals" || -n "$par_interval_strings" || -n "$par_exclude_intervals" || -n "$par_exclude_interval_strings" ]]; then
    echo "true"
  else
    echo "false"
  fi
}

# Stop with an error if no interval is given, for tools that require one.
#
# Usage: require_intervals
require_intervals() {
  if [[ -z "${par_intervals:-}" && -z "${par_interval_strings:-}" ]]; then
    echo "Error: --intervals or --interval_strings is required." >&2
    exit 1
  fi
}

# Print the path of the index that GATK writes next to a VCF output: a tabix
# index for a bgzipped VCF (`.vcf.gz` or `.vcf.bgz`) and a Tribble index for
# any other extension.
#
# Usage: gatk_index=$(gatk_output_vcf_index_path "$par_output")
gatk_output_vcf_index_path() {
  local output="$1"
  case "$output" in
    *.vcf.gz|*.vcf.bgz) echo "${output}.tbi" ;;
    *) echo "${output}.idx" ;;
  esac
}

# Print the path of the index that GATK writes next to a BAM or CRAM output:
# `<output without .bam>.bai` for a BAM file and `<output>.bai` for a CRAM
# file (a BAI format index, not a .crai).
#
# Usage: gatk_index=$(gatk_output_bam_index_path "$par_bam_output")
gatk_output_bam_index_path() {
  local output="$1"
  case "$output" in
    *.cram) echo "${output}.bai" ;;
    *) echo "${output%.bam}.bai" ;;
  esac
}
