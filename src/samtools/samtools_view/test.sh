#!/bin/bash

test_dir="${meta_resources_dir}/test_data"
temp_dir="${meta_resources_dir}/out"

############################################################################################

echo ">>> Test 1: Import SAM to BAM when @SQ lines are present in the header"
"$meta_executable" \
   --bam \
   --output "$temp_dir/a.bam" \
   --input "$test_dir/a.sam"

echo ">>> Checking whether output exists"
[ ! -f "$temp_dir/a.bam" ] && echo "File 'a.bam' does not exist!" && exit 1

echo ">>> Checking whether output is non-empty"
[ ! -s "$temp_dir/a.bam" ] && echo "File 'a.bam' is empty!" && exit 1

echo ">>> Checking whether output is correct"
# compare output of "samtools view" for both files
diff <(samtools view "$temp_dir/a.bam") <(samtools view "$test_dir/a.bam") || \
    (echo "Output file a.bam does not match expected output" && exit 1)

############################################################################################

echo ">>> Test 2: ${meta_functionality_name} with CRAM format output"

"$meta_executable" \
   --cram \
   --output "$temp_dir/a.cram" \
   --input "$test_dir/a.sam"

echo ">>> Checking whether output exists"
[ ! -f "$temp_dir/a.cram" ] && echo "File 'a.cram' does not exist!" && exit 1

echo ">>> Checking whether output is non-empty"
[ ! -s "$temp_dir/a.cram" ] && echo "File 'a.cram' is empty!" && exit 1

echo ">>> Checking whether output is correct"
# compare output of "samtools view" for both files
diff <(samtools view "$temp_dir/a.cram") <(samtools view "$test_dir/a.cram") || \
    (echo "Output file a.cram does not match expected output" && exit 1)

############################################################################################

echo ">>> Test 3: ${meta_functionality_name} with --count option"

"$meta_executable" \
   --count \
   --output "$temp_dir/a.count" \
   --input "$test_dir/a.sam"

echo ">>> Checking whether output exists"
[ ! -f "$temp_dir/a.count" ] && echo "File 'a.count' does not exist!" && exit 1

echo ">>> Checking whether output is non-empty"
[ ! -s "$temp_dir/a.count" ] && echo "File 'a.count' is empty!" && exit 1

echo ">>> Checking whether output is correct"
diff "$temp_dir/a.count" "$test_dir/a.count" || \
    (echo "Output file a.count does not match expected output" && exit 1)

############################################################################################

echo ">>> Test 4: ${meta_functionality_name} including only the forward reads from read pairs"

"$meta_executable" \
   --output "$temp_dir/a.forward" \
   --excl_flags "0x80" \
   --input "$test_dir/a.sam"

echo ">>> Checking whether output exists"
[ ! -f "$temp_dir/a.forward" ] && echo "File 'a.forward' does not exist!" && exit 1

echo ">>> Checking whether output is non-empty"
[ ! -s "$temp_dir/a.forward" ] && echo "File 'a.forward' is empty!" && exit 1

echo ">>> Checking whether output is correct"
diff "$temp_dir/a.forward" "$test_dir/a.forward" || \
    (echo "Output file a.forward does not match expected output" && exit 1)

############################################################################################

echo ">>> All test passed successfully"
rm -rf "${temp_dir}"
exit 0