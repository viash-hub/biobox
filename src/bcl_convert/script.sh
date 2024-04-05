#!/bin/bash

set -eo pipefail

$(which bcl-convert) \
  ${par_bcl_input_directory:+ --bcl-input-directory "$par_bcl_input_directory"} \
  ${par_sample_sheet:+ --sample-sheet "$par_sample_sheet"} \
  ${par_run_info:+ --run-info "$par_run_info"} \
  ${par_bcl_only_lane:+ --bcl-only-lane "$par_bcl_only_lane"} \
  ${par_first_tile_only:+ --first-tile-only "$par_first_tile_only"} \
  ${par_tiles:+ --tiles "$par_tiles"} \
  ${par_exclude_tiles:+ --exclude-tiles "$par_exclude_tiles"} \
  ${par_shared_thread_odirect_output:+ --shared-thread-odirect-output "$par_shared_thread_odirect_output"} \
  ${par_bcl_num_parallel_tiles:+ --bcl-num-parallel-tiles "$par_bcl_num_parallel_tiles"} \
  ${par_bcl_num_conversion_threads:+ --bcl-num-conversion-threads "$par_bcl_num_conversion_threads"} \
  ${par_bcl_num_compression_threads:+ --bcl-num-compression-threads "$par_bcl_num_compression_threads"} \
  ${par_bcl_num_decompression_threads:+ --bcl-num-decompression-threads "$par_bcl_num_decompression_threads"} \
  ${par_bcl_only_matched_reads:+ --bcl-only-matched-reads "$par_bcl_only_matched_reads"} \
  ${par_no_lane_splitting:+ --no-lane-splitting "$par_no_lane_splitting"} \
  ${par_num_unknown_barcodes_reported:+ --num-unknown-barcodes-reported "$par_num_unknown_barcodes_reported"} \
  ${par_bcl_validate_sample_sheet_only:+ --bcl-validate-sample-sheet-only "$par_bcl_validate_sample_sheet_only"} \
  ${par_strict_mode:+ --strict-mode "$par_strict_mode"} \
  ${par_sample_name_column_enabled:+ --sample-name-column-enabled "$par_sample_name_column_enabled"} \
  ${par_output_directory:+ --output-directory "$par_output_directory"} \
  ${par_bcl_sampleproject_subdirectories:+ --bcl-sampleproject-subdirectories "$par_bcl_sampleproject_subdirectories"} \
  ${par_fastq_gzip_compression_level:+ --fastq-gzip-compression-level "$par_fastq_gzip_compression_level"}

if [ ! -z "$par_reports" ]; then
  echo "Moving reports to their own location"
  mv "${par_output_directory}/Reports" "$par_reports"
else
  echo "Leaving reports alone"
fi

if [ ! -z "$par_logs" ]; then
  echo "Moving logs to their own location"
  mv "${par_output_directory}/Logs" "$par_logs"
else
  echo "Leaving logs alone"
fi
