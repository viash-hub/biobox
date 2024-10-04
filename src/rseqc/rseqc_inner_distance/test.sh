#!/bin/bash


# define input and output for script
input_bam="$meta_resources_dir/test_data/test.paired_end.sorted.bam"
input_bed="$meta_resources_dir/test_data/test.bed12"

output_stats="inner_distance_stats.txt"
output_dist="inner_distance.txt"
output_plot="inner_distance_plot.pdf"
output_plot_r="inner_distance_plot.r"
output_freq="inner_distance_freq.txt"

# Run executable
echo "> Running $meta_functionality_name"

"$meta_executable" \
    --input_file $input_bam \
    --refgene $input_bed \
    --output_prefix "test" \
    --output_stats $output_stats \
    --output_dist $output_dist \
    --output_plot $output_plot \
    --output_plot_r $output_plot_r \
    --output_freq $output_freq

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

echo ">> Check whether output is present and not empty"

[[ -f "$output_stats" ]] || { echo "$output_stats was not created"; exit 1; }
[[ -s "$output_stats" ]] || { echo "$output_stats is empty"; exit 1; }
[[ -f "$output_dist" ]] || { echo "$output_dist was not created"; exit 1; }
[[ -s "$output_dist" ]] || { echo "$output_dist is empty"; exit 1; }
[[ -f "$output_plot" ]] || { echo "$output_plot was not created"; exit 1; }
[[ -s "$output_plot" ]] || { echo "$output_plot is empty"; exit 1; }
[[ -f "$output_plot_r" ]] || { echo "$output_plot_r was not created"; exit 1; }
[[ -s "$output_plot_r" ]] || { echo "$output_plot_r is empty"; exit 1; }
[[ -f "$output_freq" ]] || { echo "$output_freq was created"; exit 1; }
[[ -s "$output_freq" ]] || { echo "$output_freq is empty"; exit 1; }

echo ">> Check whether output is correct"
diff "$output_freq" "$meta_resources_dir/test_data/test1.inner_distance_freq.txt" || { echo "Output is not correct"; exit 1; }
diff "$output_dist" "$meta_resources_dir/test_data/test1.inner_distance.txt" || { echo "Output is not correct"; exit 1; }

# clean up
rm "$output_stats" "$output_dist" "$output_plot" "$output_plot_r" "$output_freq"
################################################################################

echo "> Running $meta_functionality_name with non-default parameters and default output file names"
"$meta_executable" \
    --input_file $input_bam \
    --refgene $input_bed \
    --output_prefix "test" \
    --sample_size 4 \
    --mapq 10

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

echo ">> Check whether output is present and not empty"

[[ -f "test.inner_distance.txt" ]] || { echo "test.inner_distance.txt was not created"; exit 1; }
[[ -s "test.inner_distance.txt" ]] || { echo "test.inner_distance.txt is empty"; exit 1; }
[[ -f "test.inner_distance_plot.pdf" ]] || { echo "test.inner_distance_plot.pdf was not created"; exit 1; }
[[ -s "test.inner_distance_plot.pdf" ]] || { echo "test.inner_distance_plot.pdf is empty"; exit 1; }
[[ -f "test.inner_distance_plot.r" ]] || { echo "test.inner_distance_plot.r was not created"; exit 1; }
[[ -s "test.inner_distance_plot.r" ]] || { echo "test.inner_distance_plot.r is empty"; exit 1; }
[[ -f "test.inner_distance_freq.txt" ]] || { echo "test.inner_distance_freq.txt was created"; exit 1; }
[[ -s "test.inner_distance_freq.txt" ]] || { echo "test.inner_distance_freq.txt is empty"; exit 1; }

echo ">> Check whether output is correct"
diff "test.inner_distance_freq.txt" "$meta_resources_dir/test_data/test2.inner_distance_freq.txt" || { echo "Output is not correct"; exit 1; }
diff "test.inner_distance.txt" "$meta_resources_dir/test_data/test2.inner_distance.txt" || { echo "Output is not correct"; exit 1; }

exit 0