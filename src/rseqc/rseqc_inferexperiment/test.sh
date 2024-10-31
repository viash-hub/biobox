#!/bin/bash

# define input and output for script
input_bam="$meta_resources_dir/test_data/sample.bam"
input_bed="$meta_resources_dir/test_data/test.bed12"
output="strandedness.txt"

echo ">>> Prepare test output data"

cat > "$meta_resources_dir/test_data/strandedness.txt" <<EOF


This is PairEnd Data
Fraction of reads failed to determine: 0.0000
Fraction of reads explained by "1++,1--,2+-,2-+": 1.0000
Fraction of reads explained by "1+-,1-+,2++,2--": 0.0000
EOF

cat > "$meta_resources_dir/test_data/strandedness2.txt" <<EOF
Unknown Data type
EOF

################################################################################
# run executable and tests

echo ">>> Test 1: Test with default parameters"

"$meta_executable" \
    --input_file "$input_bam" \
    --refgene "$input_bed" \
    --output "$output"

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

echo ">> Checking whether output can be found and has content"

[ ! -f "$output" ] && echo "$output is missing" && exit 1
[ ! -s "$output" ] && echo "$output is empty" && exit 1


echo ">> Checking whether output is correct"
diff "$output" "$meta_resources_dir/test_data/strandedness.txt" || { echo "Output is not correct"; exit 1; }

rm "$output"

################################################################################

echo ">>> Test 2: Test with non-default sample size and map quality"

"$meta_executable" \
    --input_file "$input_bam" \
    --refgene "$input_bed" \
    --output "$output" \
    --sample_size 150000 \
    --mapq 90

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

echo ">> Checking whether output can be found and has content"

[ ! -f "$output" ] && echo "$output is missing" && exit 1
[ ! -s "$output" ] && echo "$output is empty" && exit 1

echo ">> Checking whether output is correct"
diff "$output" "$meta_resources_dir/test_data/strandedness2.txt" || { echo "Output is not correct"; exit 1; }


echo "All tests passed"

exit 0