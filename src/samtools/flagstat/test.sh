#!/bin/bash

test_dir="${meta_resources_dir}test_data"
echo ">>> Testing $meta_functionality_name"

"$meta_executable" \
  --bam "$test_dir/a.bam" \
  --bai "$test_dir/a.bam.bai" \
  --output "$test_dir/a.flagstat"

echo ">>> Checking whether output exists"
[ ! -f "$test_dir/a.flagstat" ] && echo "File 'a.flagstat' does not exist!" && exit 1

echo ">>> Checking whether output is non-empty"
[ ! -s "$test_dir/a.flagstat" ] && echo "File 'a.flagstat' is empty!" && exit 1

echo ">>> Checking whether output is correct"
diff "$test_dir/a.flagstat" "$test_dir/a_ref.flagstat" || \
    (echo "Output file chr19.flagstat does not match expected output" && exit 1)

rm "$test_dir/a.flagstat"

############################################################################################

echo ">>> Testing $meta_functionality_name with singletons in the input"

"$meta_executable" \
  --bam "$test_dir/test.paired_end.sorted.bam" \
  --bai "$test_dir/test.paired_end.sorted.bam.bai" \
  --output "$test_dir/test.paired_end.sorted.flagstat"

echo ">>> Checking whether output exists"
[ ! -f "$test_dir/test.paired_end.sorted.flagstat" ] && echo "File 'test.paired_end.sorted.flagstat' does not exist!" && exit 1

echo ">>> Checking whether output is non-empty"
[ ! -s "$test_dir/test.paired_end.sorted.flagstat" ] && echo "File 'test.paired_end.sorted.flagstat' is empty!" && exit 1

echo ">>> Checking whether output is correct"
diff "$test_dir/test.paired_end.sorted.flagstat" "$test_dir/test_ref.paired_end.sorted.flagstat" || \
    (echo "Output file chr19.flagstat does not match expected output" && exit 1)

rm "$test_dir/test.paired_end.sorted.flagstat"


echo "All tests succeeded!"
exit 0