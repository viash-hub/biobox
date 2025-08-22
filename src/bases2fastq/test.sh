#!/bin/bash

set -eou pipefail

# Source centralized test helpers
source "$meta_resources_dir/test_helpers.sh"

# Example output
# Note that the format of the fastq file names and organization into subfolders
# can differ based on the arguments provided to bases2fastq

# |-- 20230404-Bases2Fastq-Sim_QC.html
# |-- IndexAssignment.csv
# |-- Metrics.csv
# |-- RunManifest.csv
# |-- RunManifest.json
# |-- RunParameters.json
# |-- RunStats.json
# |-- Samples
# |   |-- DefaultProject
# |   |   |-- DefaultProject_IndexAssignment.csv
# |   |   |-- DefaultProject_Metrics.csv
# |   |   |-- DefaultProject_QC.html
# |   |   |-- DefaultProject_RunStats.json
# |   |   |-- sample_0
# |   |   |   |-- sample_0_L1_R1.fastq.gz
# |   |   |   |-- sample_0_L1_R2.fastq.gz
# |   |   |   |-- sample_0_L2_R1.fastq.gz
# |   |   |   |-- sample_0_L2_R2.fastq.gz
# |   |   |   `-- sample_0_stats.json
# |   |   |-- sample_1
# |   |   |   |-- sample_1_L1_R1.fastq.gz
# |   |   |   |-- sample_1_L1_R2.fastq.gz
# |   |   |   |-- sample_1_L2_R1.fastq.gz
# |   |   |   |-- sample_1_L2_R2.fastq.gz
# |   |   |   `-- sample_1_stats.json
# |   |   |-- sample_2
# |   |   |   |-- sample_2_L1_R1.fastq.gz
# |   |   |   |-- sample_2_L1_R2.fastq.gz
# |   |   |   |-- sample_2_L2_R1.fastq.gz
# |   |   |   |-- sample_2_L2_R2.fastq.gz
# |   |   |   `-- sample_2_stats.json
# |   |   |-- sample_3
# |   |   |   |-- sample_3_L1_R1.fastq.gz
# |   |   |   |-- sample_3_L1_R2.fastq.gz
# |   |   |   |-- sample_3_L2_R1.fastq.gz
# |   |   |   |-- sample_3_L2_R2.fastq.gz
# |   |   |   `-- sample_3_stats.json
# |   |   `-- sample_4
# |   |       |-- sample_4_L1_R1.fastq.gz
# |   |       |-- sample_4_L1_R2.fastq.gz
# |   |       |-- sample_4_L2_R1.fastq.gz
# |   |       |-- sample_4_L2_R2.fastq.gz
# |   |       `-- sample_4_stats.json
# |   `-- Unassigned
# |       |-- Unassigned_L1_R1.fastq.gz
# |       |-- Unassigned_L1_R2.fastq.gz
# |       |-- Unassigned_L2_R1.fastq.gz
# |       `-- Unassigned_L2_R2.fastq.gz
# |-- UnassignedSequences.csv
# `-- info
#     |-- Bases2Fastq.log
#     `-- RunManifestErrors.json


# create temporary directory and clean up on exit
TMPDIR=$(mktemp -d "$meta_temp_dir/$meta_name-XXXXXX")
function clean_up {
 [[ -d "$TMPDIR" ]] && rm -rf "$TMPDIR"
}
trap clean_up EXIT

log_info "Downloading and extracting test data"

# Unpack test input files
log_info "Downloading test data from Element Biosciences"
TAR_DIR="$TMPDIR/tar"
mkdir -p "$TAR_DIR"
wget http://element-public-data.s3.amazonaws.com/bases2fastq-share/bases2fastq-v2/20230404-bases2fastq-sim-151-151-9-9.tar.gz \
-O "$TAR_DIR/20230404-bases2fastq-sim-151-151-9-9.tar.gz"

log_info "Extracting test data"
BCL_DIR="$TMPDIR/bcl"
mkdir "$BCL_DIR"
tar -xzf "$TAR_DIR/20230404-bases2fastq-sim-151-151-9-9.tar.gz" -C "$BCL_DIR"

log_info "Running test 1 with multiple options"
mkdir "$TMPDIR/test1" && pushd "$TMPDIR/test1" > /dev/null
expected_out_dir="$TMPDIR/test1/out"
expected_report="$TMPDIR/report.html"
expected_logs="$TMPDIR/logs"
"$meta_executable" \
    --analysis_directory "$BCL_DIR/20230404-bases2fastq-sim-151-151-9-9" \
    --output_directory "$expected_out_dir" \
    --logs "$expected_logs" \
    --report "$expected_report" \
    --include_tile "L1R02C01S1;L2R21C01S1;L1R02C01S2;L2R21C01S2;L1R03C01S2;L2R20C01S2" \
    --exclude_tile "L1R04C01S1" \
    --chemistry_version 2 \
    --i1_cycles 10 \
    --i2_cycles 10 \
    --r1_cycles 152 \
    --r2_cycles 152 \
    --kit_configuration "300Cycles" \
    --detect_adapters \
    --error_on_missing \
    --flowcell_id foo \
    --force_index_orientation \
    --group_fastq \
    --legacy_fastq \
    --log_level DEBUG \
    --no_projects \
    --num_unassigned 30 \
    --run_manifest "$BCL_DIR/20230404-bases2fastq-sim-151-151-9-9/RunManifest.csv"

log_info "Validating test 1 outputs"
check_dir_exists "$expected_out_dir" "Output directory"
check_dir_exists "$expected_logs" "Logs directory"
check_file_exists "$expected_report" "HTML report"
check_file_not_empty "$expected_report" "HTML report (should contain data)"

expected_samples=(
  Undetermined_S0
  sample_0_S1
  sample_1_S2
  sample_2_S3
  sample_3_S4
  sample_4_S5
)

log_info "Checking FASTQ files for all samples and lanes"
for sample in "${expected_samples[@]}"; do
  for lane in "L001" "L002"; do 
    for orientation in "R1" "R2"; do
      check_file_exists "$expected_out_dir/${sample}_${lane}_${orientation}_001.fastq.gz" "FASTQ file for ${sample}_${lane}_${orientation}"
    done
  done
done
popd > /dev/null

log_info "Running test 3 with basic options"
mkdir "$TMPDIR/test3" && pushd "$TMPDIR/test3" > /dev/null
expected_out_dir="$TMPDIR/test3/out"
"$meta_executable" \
    --analysis_directory "$BCL_DIR/20230404-bases2fastq-sim-151-151-9-9" \
    --output_directory "$expected_out_dir"

expected_samples=(
  sample_0
  sample_1
  sample_2
  sample_3
  sample_4
)
log_info "Inspecting output directory structure:"
find "$expected_out_dir" -name "*.fastq.gz" | head -10

log_info "Checking sample FASTQ files"
for sample in "${expected_samples[@]}"; do
    for orientation in "R1" "R2"; do
      check_file_exists "$expected_out_dir/DefaultProject/${sample}/${sample}_${orientation}.fastq.gz" "Sample ${sample} ${orientation} FASTQ file"
  done
done
check_file_exists "$expected_out_dir/Unassigned/Unassigned_R1.fastq.gz" "Unassigned R1 FASTQ file"
check_file_exists "$expected_out_dir/Unassigned/Unassigned_R2.fastq.gz" "Unassigned R2 FASTQ file"
popd > /dev/null

log_info "Running test 4 with split lanes option"
mkdir "$TMPDIR/test4" && pushd "$TMPDIR/test4" > /dev/null
expected_out_dir="$TMPDIR/test4/out"
"$meta_executable" \
  --analysis_directory "$BCL_DIR/20230404-bases2fastq-sim-151-151-9-9" \
  --output_directory "$expected_out_dir" \
  --split_lanes

expected_samples=(
  "Unassigned/Unassigned"
  DefaultProject/sample_0/sample_0
  DefaultProject/sample_1/sample_1
  DefaultProject/sample_2/sample_2
  DefaultProject/sample_3/sample_3
  DefaultProject/sample_4/sample_4
)
log_info "Inspecting split lanes output directory:"
find "$expected_out_dir" -name "*.fastq.gz" | head -10

log_info "Checking split lane FASTQ files"
for sample in "${expected_samples[@]}"; do
  for lane in "L1" "L2"; do 
    for orientation in "R1" "R2"; do
      check_file_exists "$expected_out_dir/${sample}_${lane}_${orientation}.fastq.gz" "Split lane FASTQ file ${sample}_${lane}_${orientation}"
    done
  done
done
popd > /dev/null

log_info "All tests completed successfully"
