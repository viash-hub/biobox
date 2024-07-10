#!/bin/bash

# define input and output for script
input_bam="${meta_resources_dir}/test_data/sample.bam"
input_gtf="${meta_resources_dir}/test_data/genes.gtf"

output_dupmatrix="dup_matrix.txt"
output_dup_intercept_mqc="dup_intercept_mqc.txt"
output_duprate_exp_boxplot="duprate_exp_boxplot.pdf"
output_duprate_exp_densplot="duprate_exp_densityplot.pdf"
output_duprate_exp_denscurve_mqc="duprate_exp_density_curve_mqc.pdf"
output_expression_histogram="expression_hist.pdf"
output_intercept_slope="intercept_slope.txt"

# Run executable
echo "> Running $meta_functionality_name for unpaired reads, writing to tmpdir $tmpdir."

"$meta_executable" \
    --input "$input_bam" \
    --id "test" \
    --gtf_annotation "$input_gtf" \
    --strandedness "forward" \
    --paired false \
    --output_dupmatrix $output_dupmatrix \
    --output_dup_intercept_mqc $output_dup_intercept_mqc \
    --output_duprate_exp_boxplot $output_duprate_exp_boxplot \
    --output_duprate_exp_densplot $output_duprate_exp_densplot \
    --output_duprate_exp_denscurve_mqc $output_duprate_exp_denscurve_mqc \
    --output_expression_histogram $output_expression_histogram \
    --output_intercept_slope $output_intercept_slope

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

echo ">> asserting output has been created for paired read input"
[ ! -f "$output_dupmatrix" ] && echo "$output_dupmatrix was not created" && exit 1
[ ! -s "$output_dupmatrix" ] && echo "$output_dupmatrix is empty" && exit 1
[ ! -f "$output_dup_intercept_mqc" ] && echo "$output_dup_intercept_mqc was not created" && exit 1
[ ! -s "$output_dup_intercept_mqc" ] && echo "$output_dup_intercept_mqc is empty" && exit 1
[ ! -f "$output_duprate_exp_boxplot" ] && echo "$output_duprate_exp_boxplot was not created" && exit 1
[ ! -s "$output_duprate_exp_boxplot" ] && echo "$output_duprate_exp_boxplot is empty" && exit 1
[ ! -f "$output_duprate_exp_densplot" ] && echo "$output_duprate_exp_densplot was not created" && exit 1
[ ! -s "$output_duprate_exp_densplot" ] && echo "$output_duprate_exp_densplot is empty" && exit 1
[ ! -f "$output_duprate_exp_denscurve_mqc" ] && echo "$output_duprate_exp_denscurve_mqc was not created" && exit 1
[ ! -s "$output_duprate_exp_denscurve_mqc" ] && echo "$output_duprate_exp_denscurve_mqc is empty" && exit 1
[ ! -f "$output_expression_histogram" ] && echo "$output_expression_histogram was not created" && exit 1
[ ! -s "$output_expression_histogram" ] && echo "$output_expression_histogram is empty" && exit 1
[ ! -f "$output_intercept_slope" ] && echo "$output_intercept_slope was not created" && exit 1
[ ! -s "$output_intercept_slope" ] && echo "$output_intercept_slope is empty" && exit 1


echo ">> Check if output is empty"
[ ! -s "$output_dupmatrix" ] \
    && echo "Output file $output_dupmatrix is empty" && exit 1
[ ! -s "$output_dup_intercept_mqc" ] \
    && echo "Output file $output_dup_intercept_mqc is empty" && exit 1
[ ! -s "$output_intercept_slope" ] \
    && echo "Output file $output_intercept_slope is empty" && exit 1

echo ">> Check if output is correct"
cat "$output_dupmatrix"
cat "$output_dup_intercept_mqc"
cat "$output_intercept_slope"
# diff ignoring white spaces
diff -B -b "test_dupMatrix.pdf" "${meta_resources_dir}/test_data/test_dupMatrix.txt" || 
    (echo "Output file $output_dupmatrix is not correct" && exit 1)
diff -B -b "$output_dup_intercept_mqc" "${meta_resources_dir}/test_data/test_duprateExpDensCurve_mqc.txt" || \
    (echo "Output file $output_dup_intercept_mqc is not correct" && exit 1)
diff -B -b "$output_intercept_slope" "${meta_resources_dir}/test_data/test_intercept_slope.txt" || \
    (echo "Output file $output_intercept_slope is not correct" && exit 1)
exit 0