#!/bin/bash

set -eo pipefail

## VIASH START
meta_executable="$PWD/target/executable/snpeff/snpeff"
meta_resources_dir="$PWD/src/snpeff"
## VIASH END

###########################################################################

# Test 1: Run SnpEff with only required parameters

mkdir test1
pushd test1 > /dev/null # cd test1 (stack)

echo "> Run Test 1: required parameters"
"$meta_executable" \
  -genome_version GRCh37.75 \
  -input "$meta_resources_dir/test_data/cancer.vcf" \
  -output out.vcf

# Check if output files are generated
output_files=("out.vcf" "snpEff_genes.txt" "snpEff_summary.html")

# Check if any of the files do not exist
for file in "${output_files[@]}"; do
    if [ ! -e "$file" ]; then
        echo "File $file does not exist."
    fi
done

# Check if files are empty
for file in "${output_files[@]}"; do
    if [ ! -s "$file" ]; then
        echo "File $file is empty."
    fi
done

popd > /dev/null # Remove directory from stack (LIFO)

echo "Test 1 succeeded."

###########################################################################

# Test 2: Run SnpEff with a different input + options

mkdir test2
pushd test2 > /dev/null

echo "> Run Test 2: different input + options"
"$meta_executable" \
  -genome_version GRCh37.75 \
  -input "$meta_resources_dir/test_data/test.vcf" \
  -interval "$meta_resources_dir/test_data/my_annotations.bed" \
  -no_stats \
  -output output.vcf

# Check if output.vcf exists
if [ ! -e "output.vcf" ]; then
    echo "File output.vcf does not exist."
fi

# These files should not exist
files=("snpEff_genes.txt" "snpEff_summary.html")
for file in "${files[@]}"; do
    if [ -e "$file" ]; then
        echo "Error: File $file exists."
    fi
done

# Check if output.vcf is empty
if [ ! -s "output.vcf" ]; then
    echo "File output.vcf is empty."
fi

popd > /dev/null

echo "Test 2 succeeded."

###########################################################################

# Test 3: Move the output files to other locations

mkdir test3
pushd test3 > /dev/null

mkdir temp

echo "> Run Test 3: move output files"
"$meta_executable" \
  -genome_version GRCh37.75 \
  -input "$meta_resources_dir/test_data/test.vcf" \
  -output output.vcf \
  -summary temp \
  -genes temp

# Check if output.vcf exists
if [ ! -e "output.vcf" ]; then
    echo "File output.vcf does not exist."
fi

# Check if the other output files have been moved to temp folder
output_files=("snpEff_genes.txt" "snpEff_summary.html")

# Check if any of the files do not exist
for file in "${output_files[@]}"; do
    if [ ! -e "temp/$file" ]; then
        echo "File $file does not exist in "temp" folder."
    fi
done

# Check if output.vcf is empty
if [ ! -s "output.vcf" ]; then
    echo "File output.vcf is empty."
fi

# Check if the other output files in temp folder are empty
for file in "${output_files[@]}"; do
    if [ ! -s "temp/$file" ]; then
        echo "File $file is empty."
    fi
done

popd > /dev/null

echo "Test 3 succeeded."

###########################################################################

echo "All tests successfully completed!"