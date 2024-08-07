#!/bin/bash

## VIASH START
## VIASH END

test_dir="${meta_resources_dir}/test_data"

echo "> Run $meta_name with test data: identical files"
"$meta_executable" \
  --gff1 "$test_dir/1.gff" \
  --gff2 "$test_dir/1.gff" \
  --output_dir out_data

echo ">> Checking output"
[ ! -f "out_data/report.txt" ] && echo "Output file report.txt does not exist" && exit 1

echo ">> Check if output is empty"
[ ! -s "out_data/report.txt" ] && echo "Output file report.txt is empty" && exit 1

echo ">> Check if output matches expected output"
diff <(tail -n +2 "out_data/report.txt") <(tail -n +2 "$test_dir/agat_sp_compare_two_annotations_1.txt")
if [ $? -ne 0 ]; then
  echo "Output file report.txt does not match expected output"
  exit 1
fi

rm -rf out_data

echo "> Run $meta_name with test data: different files"
"$meta_executable" \
  --gff1 "$test_dir/file1.gff" \
  --gff2 "$test_dir/file2.gff" \
  --output_dir out_data

echo ">> Checking output"
[ ! -f "out_data/report.txt" ] && echo "Output file report.txt does not exist" && exit 1

echo ">> Check if output is empty"
[ ! -s "out_data/report.txt" ] && echo "Output file report.txt is empty" && exit 1

echo ">> Check if output matches expected output"
diff <(tail -n +2 "out_data/report.txt") <(tail -n +2 "$test_dir/agat_sp_compare_two_annotations_2.txt")
if [ $? -ne 0 ]; then
  echo "Output file report.txt does not match expected output"
  exit 1
fi

rm -rf out_data

echo "> Run $meta_name with test data: different files"
"$meta_executable" \
  --gff1 "$test_dir/file2.gff" \
  --gff2 "$test_dir/file1.gff" \
  --output_dir out_data

echo ">> Checking output"
[ ! -f "out_data/report.txt" ] && echo "Output file report.txt does not exist" && exit 1

echo ">> Check if output is empty"
[ ! -s "out_data/report.txt" ] && echo "Output file report.txt is empty" && exit 1

echo ">> Check if output matches expected output"
diff <(tail -n +2 "out_data/report.txt") <(tail -n +2 "$test_dir/agat_sp_compare_two_annotations_3.txt")
if [ $? -ne 0 ]; then
  echo "Output file report.txt does not match expected output"
  exit 1
fi

rm -rf out_data

echo "> Test successful"