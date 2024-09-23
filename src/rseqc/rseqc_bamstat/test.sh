#!/bin/bash

# define input and output for script

input_bam="test.paired_end.sorted.bam"
output_summary="mapping_quality.txt"

# run executable and tests
echo "> Running $meta_functionality_name."

"$meta_executable" \
    --input "$meta_resources_dir/$input_bam" \
    --output "$output_summary"

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

echo ">> Checking whether output can be found and has content"

[ ! -f "$output_summary" ] && echo "$output_summary file missing" && exit 1
[ ! -s "$output_summary" ] && echo "$output_summary file is empty" && exit 1

exit 0