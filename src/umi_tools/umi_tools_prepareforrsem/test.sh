#!/bin/bash

test_dir="$meta_resources_dir/test_data"

################################################################################
echo ">>> Test 1: with --sam:"

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

################################################################################
echo ">>> Test 2: without --sam:"

"${meta_executable}" \
    --input "$test_dir/test_dedup.bam" \
    --output "$test_dir/test_output.bam"

echo ">>> Check if output is present"
[[ ! -f "$test_dir/test_output.bam" ]] && echo "Output file not found" && exit 1
[[ ! -s "$test_dir/test_output.bam" ]] && echo "Output file is empty" && exit 1

echo ">>> Check if output is correct"
apt-get update && apt-get install -y samtools
diff <(samtools view "$test_dir/test_output.bam") <(samtools view "$test_dir/test.bam") || (echo "Output is incorrect" && exit 1)

################################################################################
echo ">>> Test 3: with --log:"

"${meta_executable}" \
    --log "$test_dir/test_log.log" \
    --input "$test_dir/test_dedup.sam" \
    --output "$test_dir/test_output.sam" \
    --sam

echo ">>> Check if output is present"
[[ ! -f "$test_dir/test_output.sam" ]] && echo "Output file not found" && exit 1
[[ ! -s "$test_dir/test_output.sam" ]] && echo "Output file is empty" && exit 1
[[ ! -f "$test_dir/test_log.log" ]] && echo "Log file not found" && exit 1
[[ ! -s "$test_dir/test_log.log" ]] && echo "Log file is empty" && exit 1

echo ">>> Check if log file is correct"
diff <(grep -v '^#' "$test_dir/test_log.log" | sed 's/^[0-9-]* [0-9:]*,[0-9]\{3\} //') <(grep -v '^#' "$test_dir/log.log" | sed 's/^[0-9-]* [0-9:]*,[0-9]\{3\} //') || (echo "Log file is incorrect" && exit 1)

echo ">>> All test succeeded"
exit 0