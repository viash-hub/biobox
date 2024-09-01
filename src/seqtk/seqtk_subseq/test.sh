#!/bin/bash

# exit on error
set -e

## VIASH START
meta_executable="target/executable/seqtk/seqtk_subseq"
meta_resources_dir="src/seqtk"
## VIASH END

# Create directories for tests
echo "Creating Test Data..."
mkdir test_data

# Create and populate input.fasta
cat > "test_data/input.fasta" <<EOL
>KU562861.1
GGAGCAGGAGAGTGTTCGAGTTCAGAGATGTCCATGGCGCCGTACGAGAAGGTGATGGATGACCTGGCCA
AGGGGCAGCAGTTCGCGACGCAGCTGCAGGGCCTCCTCCGGGACTCCCCCAAGGCCGGCCACATCATGGA
>GU056837.1
CTAATTTTATTTTTTTATAATAATTATTGGAGGAACTAAAACATTAATGAAATAATAATTATCATAATTA
TTAATTACATATTTATTAGGTATAATATTTAAGGAAAAATATATTTTATGTTAATTGTAATAATTAGAAC
>CP097510.1
CGATTTAGATCGGTGTAGTCAACACACATCCTCCACTTCCATTAGGCTTCTTGACGAGGACTACATTGAC
AGCCACCGAGGGAACCGACCTCCTCAATGAAGTCAGACGCCAAGAGCCTATCAACTTCCTTCTGCACAGC
>JAMFTS010000002.1
CCTAAACCCTAAACCCTAAACCCCCTACAAACCTTACCCTAAACCCTAAACCCTAAACCCTAAACCCTAA
ACCCGAAACCCTATACCCTAAACCCTAAACCCTAAACCCTAAACCCTAACCCAAACCTAATCCCTAAACC
>MH150936.1
TAGAAGCTAATGAAAACTTTTCCTTTACTAAAAACCGTCAAACACGGTAAGAAACGCTTTTAATCATTTC
AAAAGCAATCCCAATAGTGGTTACATCCAAACAAAACCCATTTCTTATATTTTCTCAAAAACAGTGAGAG
EOL

# Update id.list with new entries
cat > "test_data/id.list" <<EOL
KU562861.1
MH150936.1
EOL

# Create and populate reg.bed
cat > "test_data/reg.bed" <<EOL
KU562861.1$(echo -e "\t")10$(echo -e "\t")20$(echo -e "\t")region$(echo -e "\t")0$(echo -e "\t")+$(echo -e "\n")
MH150936.1$(echo -e "\t")10$(echo -e "\t")20$(echo -e "\t")region$(echo -e "\t")0$(echo -e "\t")-
EOL

#########################################################################################
# Run basic test
mkdir test1
cd test1

echo "> Run seqtk_subseq on FASTA/Q file"
"$meta_executable" \
  --input "../test_data/input.fasta" \
  --name_list "../test_data/id.list" \
  --output "sub_sample.fq"

expected_output_basic=">KU562861.1
GGAGCAGGAGAGTGTTCGAGTTCAGAGATGTCCATGGCGCCGTACGAGAAGGTGATGGATGACCTGGCCAAGGGGCAGCAGTTCGCGACGCAGCTGCAGGGCCTCCTCCGGGACTCCCCCAAGGCCGGCCACATCATGGA
>MH150936.1
TAGAAGCTAATGAAAACTTTTCCTTTACTAAAAACCGTCAAACACGGTAAGAAACGCTTTTAATCATTTCAAAAGCAATCCCAATAGTGGTTACATCCAAACAAAACCCATTTCTTATATTTTCTCAAAAACAGTGAGAG"
output_basic=$(cat sub_sample.fq)

if [ "$output_basic" != "$expected_output_basic" ]; then
  echo "Test failed"
  echo "Expected:"
  echo "$expected_output_basic"
  echo "Got:"
  echo "$output_basic"
  exit 1
fi

#########################################################################################
# Run reg.bed as name list input test
cd ..
mkdir test2
cd test2

echo "> Run seqtk_subseq on FASTA/Q file with BED file as name list"
"$meta_executable" \
  --input "../test_data/input.fasta" \
  --name_list "../test_data/reg.bed" \
  --output "sub_sample.fq"

expected_output_basic=">KU562861.1:11-20
AGTGTTCGAG
>MH150936.1:11-20
TGAAAACTTT"
output_basic=$(cat sub_sample.fq)

if [ "$output_basic" != "$expected_output_basic" ]; then
  echo "Test failed"
  echo "Expected:"
  echo "$expected_output_basic"
  echo "Got:"
  echo "$output_basic"
  exit 1
fi

#########################################################################################
# Run tab option output test
cd ..
mkdir test3
cd test3

echo "> Run seqtk_subseq with TAB option"
"$meta_executable" \
  --tab \
  --input "../test_data/input.fasta" \
  --name_list "../test_data/reg.bed" \
  --output "sub_sample.fq"

expected_output_tabular=$'KU562861.1\t11\tAGTGTTCGAG\nMH150936.1\t11\tTGAAAACTTT'
output_tabular=$(cat sub_sample.fq)

if [ "$output_tabular" != "$expected_output_tabular" ]; then
  echo "Test failed"
  echo "Expected:"
  echo "$expected_output_tabular"
  echo "Got:"
  echo "$output_tabular"
  exit 1
fi

#########################################################################################
# Run line option output test
cd ..
mkdir test4
cd test4

echo "> Run seqtk_subseq with line length option"
"$meta_executable" \
  --sequence_line_length 5 \
  --input "../test_data/input.fasta" \
  --name_list "../test_data/reg.bed" \
  --output "sub_sample.fq"

expected_output_wrapped=">KU562861.1:11-20
AGTGT
TCGAG
>MH150936.1:11-20
TGAAA
ACTTT"
output_wrapped=$(cat sub_sample.fq)

if [ "$output_wrapped" != "$expected_output_wrapped" ]; then
  echo "Test failed"
  echo "Expected:"
  echo "$expected_output_wrapped"
  echo "Got:"
  echo "$output_wrapped"
  exit 1
fi

#########################################################################################
# Run Strand Aware option output test
cd ..
mkdir test5
cd test5

echo "> Run seqtk_subseq with strand aware option"
"$meta_executable" \
  --strand_aware \
  --input "../test_data/input.fasta" \
  --name_list "../test_data/reg.bed" \
  --output "sub_sample.fq"

expected_output_wrapped=">KU562861.1:11-20
AGTGTTCGAG
>MH150936.1:11-20
AAAGTTTTCA"
output_wrapped=$(cat sub_sample.fq)

if [ "$output_wrapped" != "$expected_output_wrapped" ]; then
  echo "Test failed"
  echo "Expected:"
  echo "$expected_output_wrapped"
  echo "Got:"
  echo "$output_wrapped"
  exit 1
fi

echo "All tests succeeded!"
