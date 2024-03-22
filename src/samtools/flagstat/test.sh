#!/bin/bash

echo ">>> Testing $meta_functionality_name"

"$meta_executable" \
  --bam "$meta_resources_dir/test_data/chr19.bam" \
  --bai "$meta_resources_dir/test_data/chr19.bam.bai" \
  --output "$meta_resources_dir/test_data/chr19.flagstat"

echo ">>> Checking whether output exists"
[ ! -f "$meta_resources_dir/test_data/chr19.flagstat" ] && echo "File 'chr19.flagstat' does not exist!" && exit 1

echo ">>> Checking whether output is non-empty"
[ ! -s "$meta_resources_dir/test_data/chr19.flagstat" ] && echo "File 'chr19.flagstat' is empty!" && exit 1

echo ">>> Checking whether output is correct"
diff "$meta_resources_dir/test_data/chr19.flagstat" "$meta_resources_dir/test_data/chr19_ref.flagstat" || \
    (echo "Output file chr19.flagstat does not match expected output" && exit 1)

rm "$meta_resources_dir/test_data/chr19.flagstat"

echo "All tests succeeded!"
exit 0