#!/bin/bash

test_dir="$meta_resources_dir/test_data"

"${meta_executable}" \
    --input "$test_dir/test_dedup.sam" \
    --output "$test_dir/test_output.sam" \
    --sam


echo ">>> Check if output is present"
[[ ! -f "$test_dir/test_output.sam" ]] && echo "Output file not found" && exit 1
[[ ! -s "$test_dir/test_output.sam" ]] && echo "Output file is empty" && exit 1

echo ">>> Check if output is correct"
# use diff but ignoring the header lines (which start with @) as they may differ slightly
diff <(grep -v "^@" "$test_dir/test_output.sam") <(grep -v "^@" "$test_dir/test.sam") && echo "Output is correct" || (echo "Output is incorrect" && exit 1)

echo ">>> All test succeeded"
exit 0