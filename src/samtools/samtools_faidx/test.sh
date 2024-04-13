#!/bin/bash

test_dir="${meta_resources_dir}/test_data"
echo ">>> Testing $meta_functionality_name"

"$meta_executable" \
  "$test_dir/test.fasta" \
  --output "$test_dir/test.fasta.fai" \
  --continue \
  --fai-idx

echo ">>> Checking whether output exists"
[ ! -f "$test_dir/test.fasta.fai" ] && echo "File 'test.fasta.fai' does not exist!" && exit 1
[ ! -f "$test_dir/test.fasta.fai.idx" ] && echo "File 'test.fasta.fai.idx' does not exist!" && exit 1
[ ! -f "$test_dir/test.fasta.gzi.idx" ] && echo "File 'test.fasta.gzi.idx' does not exist!" && exit 1

echo ">>> Checking whether output is non-empty"
[ ! -s "$test_dir/test.fasta.fai" ] && echo "File 'test.fasta.fai' is empty!" && exit 1
[ ! -s "$test_dir/test.fasta.fai.idx" ] && echo "File 'test.fasta.fai.idx' is empty!" && exit 1
[ ! -s "$test_dir/test.fasta.gzi.idx" ] && echo "File 'test.fasta.gzi.idx' is empty!" && exit 1

echo ">>> Checking whether output is correct"
diff "$test_dir/a.flagstat" "$test_dir/a_ref.flagstat" || \
    (echo "Output file a.flagstat does not match expected output" && exit 1)

rm "$test_dir/a.flagstat"

echo ">>> Test 2:"

"$meta_executable" \
  "$test_dir/test.fasta" \
  --output "$test_dir/test.fasta.fai" \
  --length 60 \
  --continue \
  --gzi-idx "$test_dir/test.fasta.gz.gzi" \


echo "All tests succeeded!"
exit 0