#!/bin/bash

# exit on error
set -e

## VIASH START
meta_executable="target/executable/seqtk/seqtk_subseq"
meta_resources_dir="src/seqtk"
## VIASH END

# TODO:
# - Fix Tab option test
# - Add strand aware test (create new fasta file with right configuration)

#########################################################################################
# Run basic test
mkdir test1
cd test1

echo "> Run seqtk_subseq on FASTA/Q file"
"$meta_executable" \
  --input "$meta_resources_dir/test_data/input.fasta" \
  --name_list "$meta_resources_dir/test_data/id.list" \
  --output "sub_sample.fq"

expected_output_basic=">KU562861.1
GGAGCAGGAGAGTGTTCGAGTTCAGAGATGTCCATGGCGCCGTACGAGAAGGTGATGGATGACCTGGCCAAGGGGCAGCAGTTCGCGACGCAGCTGCAGGGCCTCCTCCGGGACTCCCCCAAGGCCGGCCACATCATGGA
>MH150936.1
TAGAAGCTAATGAAAACTTTTCCTTTACTAAAAACCGTCAAACACGGTAAGAAACGCTTTTAATCATTTCAAAAGCAATCCCAATAGTGGTTACATCCAAACAAAACCCATTTCTTATATTTTCTCAAAAACAGTGAGAG"
output_basic=$(cat sub_sample.fq)

if [ "$output_basic" == "$expected_output_basic" ]; then
  echo "Test passed"
else
  echo "Test failed"
  echo "Expected:"
  echo "$expected_output_basic"
  echo "Got:"
  echo "$output_basic"
fi

#########################################################################################
# Run reg.bed as name list input test
cd ..
mkdir test2
cd test2

echo "> Run seqtk_subseq on FASTA/Q file with BED file as name list"
"$meta_executable" \
  --input "$meta_resources_dir/test_data/input.fasta" \
  --name_list "$meta_resources_dir/test_data/reg.bed" \
  --output "sub_sample.fq"

expected_output_basic=">KU562861.1:11-20
AGTGTTCGAG
>MH150936.1:11-20
TGAAAACTTT"
output_basic=$(cat sub_sample.fq)

if [ "$output_basic" == "$expected_output_basic" ]; then
  echo "Test passed"
else
  echo "Test failed"
  echo "Expected:"
  echo "$expected_output_basic"
  echo "Got:"
  echo "$output_basic"
fi

#########################################################################################
# Run tab option output test
cd ..
mkdir test3
cd test3

echo "> Run seqtk_subseq with TAB option"
"$meta_executable" \
  --tab \
  --input "$meta_resources_dir/test_data/input.fasta" \
  --name_list "$meta_resources_dir/test_data/reg.bed" \
  --output "sub_sample.fq"

expected_output_tabular=$'KU562861.1\t11\tAGTGTTCGAG\nMH150936.1\t11\tTGAAAACTTT'
output_tabular=$(cat sub_sample.fq)

if [ "$output_tabular" == "$expected_output_tabular" ]; then
  echo "Test passed"
else
  echo "Test failed"
  echo "Expected:"
  echo "$expected_output_tabular"
  echo "Got:"
  echo "$output_tabular"
fi

#########################################################################################
# Run line option output test
cd ..
mkdir test4
cd test4

echo "> Run seqtk_subseq with line length option"
"$meta_executable" \
  --sequence_line_length 5 \
  --input "$meta_resources_dir/test_data/input.fasta" \
  --name_list "$meta_resources_dir/test_data/reg.bed" \
  --output "sub_sample.fq"

expected_output_wrapped=">KU562861.1:11-20
AGTGT
TCGAG
>MH150936.1:11-20
TGAAA
ACTTT"
output_wrapped=$(cat sub_sample.fq)

if [ "$output_wrapped" == "$expected_output_wrapped" ]; then
  echo "Test passed"
else
  echo "Test failed"
  echo "Expected:"
  echo "$expected_output_wrapped"
  echo "Got:"
  echo "$output_wrapped"
fi

#########################################################################################
# Run Strand Aware option output test
cd ..
mkdir test5
cd test5

echo "> Run seqtk_subseq with strand aware option"
"$meta_executable" \
  --strand_aware \
  --input "$meta_resources_dir/test_data/input.fasta" \
  --name_list "$meta_resources_dir/test_data/reg.bed" \
  --output "sub_sample.fq"

# expected_output_wrapped=">KU562861.1:11-20
# AGTGT
# TCGAG
# >MH150936.1:11-20
# TGAAA
# ACTTT"
# output_wrapped=$(cat sub_sample.fq)

# if [ "$output_wrapped" == "$expected_output_wrapped" ]; then
#   echo "Line-wrapped output test passed"
# else
#   echo "Line-wrapped output test failed"
#   echo "Expected:"
#   echo "$expected_output_wrapped"
#   echo "Got:"
#   echo "$output_wrapped"
# fi