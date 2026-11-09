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

# Symlink a BAM and its .bai companion into a temp dir under matching
# basenames, so GATK can find them
#
# Usage: staged=$(stage_bam_bai "$tmp_dir" "$par_input" "$par_bai")
stage_bam_bai() {
  local tmp_dir="$1"
  local input_bam="$2"
  local input_bai="$3"

  # Create links, explicitly exit on failure
  ln -s "$(readlink -f "$input_bam")" "$tmp_dir/sample.bam" || exit 1
  ln -s "$(readlink -f "$input_bai")" "$tmp_dir/sample.bai" || exit 1

  echo "$tmp_dir/sample.bam"
}

# Symlink one or more BAMs and their .bai companions into a temp dir under
# matching basenames, so GATK can find them. The BAMs and BAIs are given as
# semicolon-separated `multiple: true` argument values and are matched by
# order. Stops with an error if the number of BAMs and BAIs is different.
#
# Usage: stage_bams_bais "$tmp_dir" "$par_input" "$par_bai" result_array_name
stage_bams_bais() {
  local tmp_dir="$1"
  local input_bams="$2"
  local input_bais="$3"
  local -n result_ref="$4"

  local bams=()
  local bais=()
  IFS=';' read -ra bams <<< "$input_bams"
  IFS=';' read -ra bais <<< "$input_bais"

  if [[ ${#bams[@]} -ne ${#bais[@]} ]]; then
    echo "Error: --input and --bai must be given the same number of times (got ${#bams[@]} BAM file(s) and ${#bais[@]} BAI file(s))." >&2
    exit 1
  fi

  # Stage to the original basename with a prefix to avoid collisions
  result_ref=()
  local i
  for i in "${!bams[@]}"; do
    local staged_bam
    staged_bam="$tmp_dir/input_$((i + 1))_$(basename "${bams[$i]}")"
    ln -s "$(readlink -f "${bams[$i]}")" "$staged_bam"
    ln -s "$(readlink -f "${bais[$i]}")" "${staged_bam%.bam}.bai"
    result_ref+=("$staged_bam")
  done
}

# Stage a VCF (or other feature file) into a temp dir so that GATK can find
# its index. If an index is given, the file and the index are symlinked under
# matching names (`.tbi` for tabix indexes of `.vcf.gz`/`.vcf.bgz` files,
# `.idx` for Tribble indexes of plain `.vcf` files). If no index is given, the
# file is copied and indexed with IndexFeatureFile.
#
# Usage: staged=$(stage_vcf_with_index "$tmp_dir" "$par_vcf" "$par_vcf_index" prefix)
stage_vcf_with_index() {
  local tmp_dir="$1"
  local vcf="$2"
  local vcf_index="$3"
  local prefix="$4"
  local staged_vcf="$tmp_dir/${prefix}_$(basename "$vcf")"

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
    if ! gatk IndexFeatureFile --input "$staged_vcf" --verbosity ERROR >&2; then
      echo "Error: could not index '$vcf'. A .vcf.gz file must be compressed with bgzip." >&2
      exit 1
    fi
  fi

  echo "$staged_vcf"
}

# Stage one or more interval files (BED, interval_list, VCF) and build the
# repeated GATK interval flags for them. VCF files are staged with
# stage_vcf_with_index because GATK needs an index for a block-compressed VCF.
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
      staged_file=$(stage_vcf_with_index "$tmp_dir" "$file" "${indexes[$vcf_i]:-}" "${prefix}_$((i + 1))")
      result_ref+=("$gatk_flag" "$staged_file")
      vcf_i=$((vcf_i + 1))
    else
      result_ref+=("$gatk_flag" "$file")
    fi
  done
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

# Print the path of the index that GATK writes next to a BAM output
# (`<output without .bam>.bai`).
#
# Usage: gatk_index=$(gatk_output_bam_index_path "$par_bam_output")
gatk_output_bam_index_path() {
  local output="$1"
  echo "${output%.bam}.bai"
}
