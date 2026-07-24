#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

source "$meta_resources_dir/gatk4/script_helpers.sh"

# unset flags
unset_if_false=(
  par_remove_duplicates
  par_remove_sequencing_duplicates
  par_tag_duplicate_set_members
  par_create_index
  par_duplex_umi
)

for par in "${unset_if_false[@]}"; do
  test_val="${!par}"
  [[ "$test_val" == "false" ]] && unset $par
done

# Create a temporary directory for GATK's own scratch space
tmp_dir=$(mktemp -d)
trap 'rm -rf "$tmp_dir"' EXIT

# Compute available memory for the JVM (80% of allocated memory, fallback to 3072MB)
avail_mem_mb=$(( ${meta_memory_mb:-3072} * 8 / 10 ))

# Convert semicolon-separated multi-value arguments into repeated flags
split_multiple_to_flags "$par_input" "--INPUT" input_args
split_multiple_to_flags "$par_comment" "--COMMENT" comment_args

# Build command arguments array
cmd_args=(
  "${input_args[@]}"
  --OUTPUT "$par_output"
  --METRICS_FILE "$par_metrics"
  --TMP_DIR "$tmp_dir"
  ${par_add_pg_tag_to_reads:+--ADD_PG_TAG_TO_READS "$par_add_pg_tag_to_reads"}
  ${par_assume_sort_order:+--ASSUME_SORT_ORDER "$par_assume_sort_order"}
  ${par_barcode_tag:+--BARCODE_TAG "$par_barcode_tag"}
  ${par_clear_dt:+--CLEAR_DT "$par_clear_dt"}
  "${comment_args[@]}"
  ${par_create_index:+--CREATE_INDEX true}
  ${par_duplex_umi:+--DUPLEX_UMI true}
  ${par_duplicate_scoring_strategy:+--DUPLICATE_SCORING_STRATEGY "$par_duplicate_scoring_strategy"}
  ${par_max_file_handles_for_read_ends_map:+--MAX_FILE_HANDLES_FOR_READ_ENDS_MAP "$par_max_file_handles_for_read_ends_map"}
  ${par_max_optical_duplicate_set_size:+--MAX_OPTICAL_DUPLICATE_SET_SIZE "$par_max_optical_duplicate_set_size"}
  ${par_max_records_in_ram:+--MAX_RECORDS_IN_RAM "$par_max_records_in_ram"}
  ${par_molecular_identifier_tag:+--MOLECULAR_IDENTIFIER_TAG "$par_molecular_identifier_tag"}
  ${par_optical_duplicate_pixel_distance:+--OPTICAL_DUPLICATE_PIXEL_DISTANCE "$par_optical_duplicate_pixel_distance"}
  ${par_program_group_command_line:+--PROGRAM_GROUP_COMMAND_LINE "$par_program_group_command_line"}
  ${par_program_group_name:+--PROGRAM_GROUP_NAME "$par_program_group_name"}
  ${par_program_group_version:+--PROGRAM_GROUP_VERSION "$par_program_group_version"}
  ${par_program_record_id:+--PROGRAM_RECORD_ID "$par_program_record_id"}
  ${par_read_name_regex:+--READ_NAME_REGEX "$par_read_name_regex"}
  ${par_read_one_barcode_tag:+--READ_ONE_BARCODE_TAG "$par_read_one_barcode_tag"}
  ${par_read_two_barcode_tag:+--READ_TWO_BARCODE_TAG "$par_read_two_barcode_tag"}
  ${par_reference_sequence:+--REFERENCE_SEQUENCE "$par_reference_sequence"}
  ${par_remove_duplicates:+--REMOVE_DUPLICATES true}
  ${par_remove_sequencing_duplicates:+--REMOVE_SEQUENCING_DUPLICATES true}
  ${par_sorting_collection_size_ratio:+--SORTING_COLLECTION_SIZE_RATIO "$par_sorting_collection_size_ratio"}
  ${par_tag_duplicate_set_members:+--TAG_DUPLICATE_SET_MEMBERS true}
  ${par_tagging_policy:+--TAGGING_POLICY "$par_tagging_policy"}
  ${par_validation_stringency:+--VALIDATION_STRINGENCY "$par_validation_stringency"}
)

# Run gatk MarkDuplicates
gatk --java-options "-Xmx${avail_mem_mb}M -XX:-UsePerfData" MarkDuplicates "${cmd_args[@]}"
