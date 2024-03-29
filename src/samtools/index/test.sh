#!/bin/bash

test_dir="${meta_resources_dir}test_data"
echo ">>> Testing $meta_functionality_name"

echo ">>> Generating BAM index"
"$meta_executable" \
  --input "$test_dir/a.sorted.bam" \
  --bai \
  --output "$test_dir/a.sorted.bam.bai"

echo ">>> Check whether output exists"
[ ! -f "$test_dir/a.sorted.bam.bai" ] && echo "File 'a.sorted.bam.bai' does not exist!" && exit 1

echo ">>> Check whether output is empty"
[ ! -s "$test_dir/a.sorted.bam.bai" ] && echo "File 'a.sorted.bam.bai' is empty!" && exit 1

echo ">>> Check whether output is correct"
diff "$test_dir/a.sorted.bam.bai" "$test_dir/a_ref.sorted.bam.bai" || \
  (echo "File 'a.sorted.bam.bai' does not match expected output." && exit 1)

rm "$test_dir/a.sorted.bam.bai"

#################################################################################################

echo ">>> Generating CSI index"
"$meta_executable" \
  --input "$test_dir/a.sorted.bam" \
  --csi \
  --output "$test_dir/a.sorted.bam.csi"

echo ">>> Check whether output exists"
[ ! -f "$test_dir/a.sorted.bam.csi" ] && echo "File 'a.sorted.bam.csi' does not exist!" && exit 1

echo ">>> Check whether output is empty"
[ ! -s "$test_dir/a.sorted.bam.csi" ] && echo "File 'a.sorted.bam.csi' is empty!" && exit 1

echo ">>> Check whether output is correct"
diff "$test_dir/a.sorted.bam.csi" "$test_dir/a_ref.sorted.bam.csi" || \
  (echo "File 'a.sorted.bam.csi' does not match expected output." && exit 1)

rm "$test_dir/a.sorted.bam.csi"

#################################################################################################

echo ">>> Generating bam index with -M option"
"$meta_executable" \
  --input "$test_dir/a.sorted.bam" \
  --bai \
  --output "$test_dir/a_multiple.sorted.bam.bai" \
  --multiple

echo ">>> Check whether output exists"
[ ! -f "$test_dir/a_multiple.sorted.bam.bai" ] && echo "File 'a_multiple.sorted.bam.bai' does not exist!" && exit 1

echo ">>> Check whether output is empty"
[ ! -s "$test_dir/a_multiple.sorted.bam.bai" ] && echo "File 'a_multiple.sorted.bam.bai' is empty!" && exit 1

echo ">>> Check whether output is correct"
diff "$test_dir/a_multiple.sorted.bam.bai" "$test_dir/a_multiple_ref.sorted.bam.bai" || \
  (echo "File 'a_multiple.sorted.bam.bai' does not match expected output." && exit 1)


#################################################################################################

echo ">>> Generating BAM index with -m option"

"$meta_executable" \
  --input "$test_dir/a.sorted.bam" \
  --min_shift 4 \
  --bai \
  --output "$test_dir/a_4.sorted.bam.bai"

echo ">>> Check whether output exists"
[ ! -f "$test_dir/a_4.sorted.bam.bai" ] && echo "File 'a_4.sorted.bam.bai' does not exist!" && exit 1

echo ">>> Check whether output is empty"
[ ! -s "$test_dir/a_4.sorted.bam.bai" ] && echo "File 'a_4.sorted.bam.bai' is empty!" && exit 1

echo ">>> Check whether output is correct"
diff "$test_dir/a_4.sorted.bam.bai" "$test_dir/a_4_ref.sorted.bam.bai" || \
  (echo "File 'a_4.sorted.bam.bai' does not match expected output." && exit 1)

rm "$test_dir/a_4.sorted.bam.bai"

#################################################################################################


echo "All tests succeeded!"
exit 0