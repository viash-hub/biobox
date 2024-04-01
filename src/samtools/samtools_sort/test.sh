#!/bin/bash

test_dir="${meta_resources_dir}test_data"

echo ">>> Test 1: ******"

"$meta_executable" \
  --bam "$test_dir/a.bam" \
  --output "$test_dir/a.sorted.bam"

############################################################################################

echo ">>> Test 2: ******"

"$meta_executable" \
  --bam "$test_dir/test.paired_end.bam" \
  --output "$test_dir/test.paired_end.sorted.bam"


############################################################################################

echo "All tests succeeded!"
exit 0