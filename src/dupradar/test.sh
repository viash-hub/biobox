#!/bin/bash

source "$meta_resources_dir/test_helpers.sh"
setup_test_env

test_data_dir="$meta_resources_dir/test_data"
input_bam="$test_data_dir/sample.bam"
input_gtf="$test_data_dir/genes.gtf"

# A number as written by R, e.g. 12, 0.114, 5.8e-09 or 1e+05
num="[0-9.e+-]+"

# Check all seven outputs in a directory.
# Usage: check_outputs "/path/to/dir" dupmatrix intercept_mqc boxplot densplot denscurve_mqc histogram intercept_slope
check_outputs() {
  local dir="$1"
  local dupmatrix="$dir/$2"
  local intercept_mqc="$dir/$3"
  local boxplot="$dir/$4"
  local densplot="$dir/$5"
  local denscurve_mqc="$dir/$6"
  local histogram="$dir/$7"
  local intercept_slope="$dir/$8"

  log "Checking whether output exists"
  check_file_exists "$dupmatrix" "duplication matrix"
  check_file_not_empty "$dupmatrix" "duplication matrix"
  check_file_exists "$intercept_mqc" "MultiQC intercept table"
  check_file_not_empty "$intercept_mqc" "MultiQC intercept table"
  check_file_exists "$boxplot" "duplication rate boxplot"
  check_file_not_empty "$boxplot" "duplication rate boxplot"
  check_file_exists "$densplot" "density scatter plot"
  check_file_not_empty "$densplot" "density scatter plot"
  check_file_exists "$denscurve_mqc" "MultiQC density curve"
  check_file_not_empty "$denscurve_mqc" "MultiQC density curve"
  check_file_exists "$histogram" "expression histogram"
  check_file_not_empty "$histogram" "expression histogram"
  check_file_exists "$intercept_slope" "intercept and slope"
  check_file_not_empty "$intercept_slope" "intercept and slope"
  check_file_not_exists "$dir/Rplots.pdf" "stray R default plot device file"

  log "Checking MultiQC headers"
  check_file_contains "$intercept_mqc" "#id: DupInt" "MultiQC intercept table"
  check_file_contains "$denscurve_mqc" "#id: dupradar" "MultiQC density curve"

  log "Checking output structure"
  check_file_contains "$dupmatrix" "dupRate" "duplication matrix header"
  check_file_matches_regex "$intercept_mqc" "^test $num$" "MultiQC intercept row"
  check_file_matches_regex "$denscurve_mqc" "^$num $num$" "MultiQC density curve data rows"
  check_file_matches_regex "$intercept_slope" "^test - dupRadar Int .*: $num" "intercept line"
  check_file_matches_regex "$intercept_slope" "^test - dupRadar Sl .*: $num" "slope line"
}

log "Run dupRadar on single-end test data, with explicit output paths"
test_1_output="$meta_temp_dir/output1"
mkdir -p "$test_1_output"

"$meta_executable" \
  --input_bam "$input_bam" \
  --input_gtf "$input_gtf" \
  --id "test" \
  --strandedness 1 \
  --output_dupmatrix "$test_1_output/dup_matrix.txt" \
  --output_dup_intercept_mqc "$test_1_output/dup_intercept_mqc.txt" \
  --output_duprate_exp_boxplot "$test_1_output/duprate_exp_boxplot.pdf" \
  --output_duprate_exp_densplot "$test_1_output/duprate_exp_densityplot.pdf" \
  --output_duprate_exp_denscurve_mqc "$test_1_output/duprate_exp_density_curve_mqc.txt" \
  --output_expression_histogram "$test_1_output/expression_hist.pdf" \
  --output_intercept_slope "$test_1_output/intercept_slope.txt"

check_outputs "$test_1_output" \
  dup_matrix.txt \
  dup_intercept_mqc.txt \
  duprate_exp_boxplot.pdf \
  duprate_exp_densityplot.pdf \
  duprate_exp_density_curve_mqc.txt \
  expression_hist.pdf \
  intercept_slope.txt

log "Cleanup"
rm -rf "$test_1_output"

log "Run dupRadar on single-end test data, without output paths (default <id>_* names)"
test_2_output="$meta_temp_dir/output2"
mkdir -p "$test_2_output"

(
  cd "$test_2_output"
  "$meta_executable" \
    --input_bam "$input_bam" \
    --input_gtf "$input_gtf" \
    --id "test" \
    --strandedness 1
)

check_outputs "$test_2_output" \
  test_dupMatrix.txt \
  test_dup_intercept_mqc.txt \
  test_duprateExpBoxplot.pdf \
  test_duprate_exp_densplot.pdf \
  test_duprateExpDensCurve_mqc.txt \
  test_expressionHist.pdf \
  test_intercept_slope.txt

log "Cleanup"
rm -rf "$test_2_output"

log "Run dupRadar on single-end test data, with only some output paths (rest default <id>_* names)"
test_3_output="$meta_temp_dir/output3"
mkdir -p "$test_3_output"

(
  cd "$test_3_output"
  "$meta_executable" \
    --input_bam "$input_bam" \
    --input_gtf "$input_gtf" \
    --id "test" \
    --strandedness 1 \
    --output_dupmatrix "$test_3_output/custom_dup_matrix.txt" \
    --output_intercept_slope "$test_3_output/custom_intercept_slope.txt"
)

check_outputs "$test_3_output" \
  custom_dup_matrix.txt \
  test_dup_intercept_mqc.txt \
  test_duprateExpBoxplot.pdf \
  test_duprate_exp_densplot.pdf \
  test_duprateExpDensCurve_mqc.txt \
  test_expressionHist.pdf \
  custom_intercept_slope.txt

log "Checking that given output paths replace the default names"
check_file_not_exists "$test_3_output/test_dupMatrix.txt" "default duplication matrix"
check_file_not_exists "$test_3_output/test_intercept_slope.txt" "default intercept and slope"

log "Cleanup"
rm -rf "$test_3_output"

print_test_summary "dupRadar tests"
