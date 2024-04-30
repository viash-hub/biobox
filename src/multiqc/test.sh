#!/bin/bash

echo ">>> Testing input/output handling"

"$meta_executable" \
  --input "$meta_resources_dir/test_data/" \
  --output_report test1.html \
  --output_data data1 \
  --output_plots plots1 \
  --quiet

[ ! -f test1.html ] && echo "MultiQC report does not exist!" && exit 1
[ ! -d data1 ] && echo "MultiQC data directory does not exist!" && exit 1
[ ! -d plots1 ] && echo "MultiQC plots directory does not exist!" && exit 1

echo ">>> Testing module exclusion"

"$meta_executable" \
  --input "$meta_resources_dir/test_data/" \
  --output_report test2.html \
  --output_data data2 \
  --output_plots plots2 \
  --exclude_modules samtools \
  --quiet

[ -f test2.html ] && echo "MultiQC report should not exist!" && exit 1
[ -d data2 ] && echo "MultiQC data directory should not exist!" && exit 1
[ -d plots2 ] && echo "MultiQC plots directory should not exist!" && exit 1

echo ">>> Testing sample exclusion"

"$meta_executable" \
  --input "$meta_resources_dir/test_data/" \
  --output_report test3.html \
  --output_data data3 \
  --ignore_samples a \
  --quiet

key_to_check=".report_general_stats_data[0].a"
json_file="data3/multiqc_data.json"
[[ $(jq -r "$key_to_check" "$json_file") != null ]] && echo "$key_to_check should not be present in $json_file" && exit 1

echo "All tests succeeded!"
exit 0