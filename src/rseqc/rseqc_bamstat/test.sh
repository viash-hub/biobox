#!/bin/bash

# define input and output for script

input_bam="sample.bam"
output_summary="mapping_quality.txt"

# run executable and tests
echo "> Running $meta_functionality_name."

"$meta_executable" \
    --input_file "$meta_resources_dir/test_data/$input_bam" \
    --output "$output_summary"

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

echo ">> Checking whether output is present"
[ ! -f "$output_summary" ] && echo "$output_summary file missing" && exit 1
[ ! -s "$output_summary" ] && echo "$output_summary file is empty" && exit 1

echo ">> Checking whether output is correct"
diff "$meta_resources_dir/test_data/ref_output.txt" "$meta_resources_dir/$output_summary" || { echo "Output is not correct"; exit 1; }

#############################################################################

echo ">>> Test 2: Test with non-default mapping quality threshold"

output_summary="mapping_quality_mapq_50.txt"

# run executable and tests
echo "> Running $meta_functionality_name."

"$meta_executable" \
    --input_file "$meta_resources_dir/test_data/$input_bam" \
    --output "$output_summary" \
    --mapq 50

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

echo ">> Checking whether output is present"
[ ! -f "$output_summary" ] && echo "$output_summary file missing" && exit 1
[ ! -s "$output_summary" ] && echo "$output_summary file is empty" && exit 1

echo ">> Checking whether output is correct"
diff "$meta_resources_dir/test_data/ref_output_mapq.txt" "$meta_resources_dir/$output_summary" || { echo "Output is not correct"; exit 1; }

exit 0