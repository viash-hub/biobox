#!/bin/bash

test_dir="${meta_resources_dir}/test_data"
echo ">>> Testing $meta_functionality_name"

"$meta_executable" \
  --input "$test_dir/test.fasta" \
  --output "$test_dir/test.fasta.fai"

echo "$meta_executable"
echo "$test_dir/test.fasta"

echo ">>> Checking whether output exists"
[ ! -f "$test_dir/test.fasta.fai" ] && echo "File 'test.fasta.fai' does not exist!" && exit 1

echo ">>> Checking whether output is non-empty"
[ ! -s "$test_dir/test.fasta.fai" ] && echo "File 'test.fasta.fai' is empty!" && exit 1

echo ">>> Checking whether output is correct"
diff "$test_dir/test.fasta.fai" "$test_dir/output/test.fasta.fai" || \
    (echo "Output file test.fasta.fai does not match expected output" && exit 1)

rm "$test_dir/test.fasta.fai"

####################################################################################################

echo ">>> Test 2: ${meta_functionality_name} with bgzipped input"

"$meta_executable" \
  --input "$test_dir/test.fasta.gz" \
  --output "$test_dir/test.fasta.gz.fai"

echo ">>> Checking whether output exists"1
[ ! -f "$test_dir/test.fasta.gz.fai" ] && echo "File 'test.fasta.gz.fai' does not exist!" && exit 1
[ ! -f "$test_dir/test.fasta.gz.gzi" ] && echo "File 'test.fasta.gz.gzi' does not exist!" && exit 1

echo ">>> Checking whether output is non-empty"
[ ! -s "$test_dir/test.fasta.gz.fai" ] && echo "File 'test.fasta.gz.fai' is empty!" && exit 1
[ ! -s "$test_dir/test.fasta.gz.gzi" ] && echo "File 'test.fasta.gz.gzi' is empty!" && exit 1

echo ">>> Checking whether output is correct"
diff "$test_dir/test.fasta.gz.fai" "$test_dir/output/test.fasta.gz.fai" || \
    (echo "Output file test_zip.fasta.gz.fai does not match expected output" && exit 1)
diff "$test_dir/test.fasta.gz.gzi" "$test_dir/output/test.fasta.gz.gzi" || \
    (echo "Output file test2.fasta.gz.gzi does not match expected output" && exit 1)

rm "$test_dir/test.fasta.gz.fai"
rm "$test_dir/test.fasta.gz.gzi"

####################################################################################################

echo ">>> Test 3: ${meta_functionality_name} with fastq input"

"$meta_executable" \
  --input "$test_dir/test.fastq" \
  --output "$test_dir/test.fastq.fai"

echo ">>> Checking whether output exists"
[ ! -f "$test_dir/test.fastq.fai" ] && echo "File 'test.fastq.fai' does not exist!" && exit 1

echo ">>> Checking whether output is non-empty"
[ ! -s "$test_dir/test.fastq.fai" ] && echo "File 'test.fastq.fai' is empty!" && exit 1

echo ">>> Checking whether output is correct"
diff "$test_dir/test.fastq.fai" "$test_dir/output/test.fastq.fai" || \
    (echo "Output file test.fastq.fai does not match expected output" && exit 1)

rm "$test_dir/test.fastq.fai"

####################################################################################################

echo ">>> Test 4: ${meta_functionality_name} with region file containing non-existent regions and 
      specific fasta line wrap length"

"$meta_executable" \
  --input "$test_dir/test.fasta" \
  --output "$test_dir/regions.fasta" \
  --length 10 \
  --continue \
  --region_file "$test_dir/test.regions" \
  --fai_idx "$test_dir/regions.fasta.fai"

echo ">>> Checking whether output exists"
[ ! -f "$test_dir/regions.fasta" ] && echo "File 'regions.fasta' does not exist!" && exit 1
[ ! -f "$test_dir/regions.fasta.fai" ] && echo "File 'regions.fasta.fai' does not exist!" && exit 1

echo ">>> Checking whether output is non-empty"
[ ! -s "$test_dir/regions.fasta" ] && echo "File 'regions.fasta' is empty!" && exit 1
[ ! -s "$test_dir/regions.fasta.fai" ] && echo "File 'regions.fasta.fai' is empty!" && exit 1

echo ">>> Checking whether output is correct"
diff "$test_dir/regions.fasta" "$test_dir/output/regions.fasta" || \
    (echo "Output file regions.fasta does not match expected output" && exit 1)
diff "$test_dir/regions.fasta.fai" "$test_dir/output/regions.fasta.fai" || \
    (echo "Output file regions.fasta.fai does not match expected output" && exit 1)

rm "$test_dir/regions.fasta"
rm "$test_dir/regions.fasta.fai"

####################################################################################################

echo "All tests succeeded!"
exit 0

