#!/bin/bash

set -eo pipefail

## VIASH START
meta_executable="$PWD/target/executable/nanoplot/nanoplot"
meta_resources_dir="$PWD/src/nanoplot"
## VIASH END

# Files at runtime (.gz, .pickle and .feather)
wget "https://github.com/wdecoster/nanotest/archive/refs/heads/master.zip"
unzip master.zip

###########################################################################

# Test 1: Run NanoPlot with only input parameter (Fastq)

mkdir test1
pushd test1 > /dev/null # cd test1 (stack)

echo "> Run Test 1: one input (Fastq)"
"$meta_executable" \
  --fastq "$meta_resources_dir/test_data/test1.fastq" \
  --outdir output

# Check if output directory exists
if [[ ! -d output ]]; then
  echo "Output directory not found!"
  exit 1
fi

# Check if output files are generated
if [ "$(ls -1 "output" | wc -l)" -lt 1 ]; then # Apart from log file
  echo "Output files are not found!"
  exit 1
fi

# Check if files are empty
if find output -name "*.html" -type f -size 0 | grep -q .; then 
  echo "At least one HTML file is empty."
  exit 1
fi
if find output -name "*.png" -type f -size 0 | grep -q .; then 
  echo "At least one plot is empty."
  exit 1
fi
if find output -name "*.txt" -type f -size 0 | grep -q .; then 
  echo "NanoPlot summary file is empty."
  exit 1
fi

popd > /dev/null # Remove directory from stack (LIFO)

echo "Test 1 succeeded."

###########################################################################

# Test 2: Run NanoPlot with multiple inputs (Fastq)

mkdir test2
pushd test2 > /dev/null

echo "> Run Test 2: multiple inputs (Fastq)"
"$meta_executable" \
  --fastq "$meta_resources_dir/test_data/test1.fastq" "$meta_resources_dir/test_data/test2.fastq" \
  --outdir output

# Check if output directory exists
if [[ ! -d output ]]; then
  echo "Output directory not found!"
  exit 1
fi

# Check if output files are generated
if [ "$(ls -1 "output" | wc -l)" -lt 1 ]; then
  echo "Output files are not found!"
  exit 1
fi

# Check if files are empty
if find output -name "*.html" -type f -size 0 | grep -q .; then 
  echo "At least one HTML file is empty."
  exit 1
fi
if find output -name "*.png" -type f -size 0 | grep -q .; then 
  echo "At least one plot is empty."
  exit 1
fi
if find output -name "*.txt" -type f -size 0 | grep -q .; then 
  echo "NanoPlot summary file is empty."
  exit 1
fi

popd > /dev/null

echo "Test 2 succeeded."

###########################################################################

# Test 3: Run NanoPlot with multiple options-1

mkdir test3
pushd test3 > /dev/null

echo "> Run Test 3: multiple options-1"
"$meta_executable" \
  --fastq "$meta_resources_dir/test_data/test1.fastq" \
  --maxlength 40000 \
  --format jpg \
  --prefix biobox_ \
  --store \
  --color yellow \
  --info_in_report \
  --outdir output

# Check if output directory exists
if [[ ! -d output ]]; then
  echo "Output directory not found!"
  exit 1
fi

# Check if output files are generated
if [ "$(ls -1 "output" | wc -l)" -lt 1 ]; then
  echo "Output files are not found!"
  exit 1
fi

# Check if the extracted data exists (--store)
if ! ls output/*.pickle > /dev/null 2>&1; then
  echo "Extracted data is not found!"
  exit 1
fi

# Check if files are empty
if find output -name "*.html" -type f -size 0 | grep -q .; then 
  echo "At least one HTML file is empty."
  exit 1
fi
if find output -name "*.png" -type f -size 0 | grep -q .; then 
  echo "At least one plot is empty."
  exit 1
fi
if find output -name "*.txt" -type f -size 0 | grep -q .; then 
  echo "NanoPlot summary file is empty."
  exit 1
fi
if find output -name "*.pickle" -type f -size 0 | grep -q .; then 
  echo "Extracted data is empty."
  exit 1
fi

# Check if the output file starts with "biobox" prefix
if ! ls output/biobox* > /dev/null 2>&1; then
    echo "The prefix is not added to the output files."
    exit 1
fi

popd > /dev/null

echo "Test 3 succeeded."

###########################################################################

# Test 4: Run NanoPlot with multiple options-2

mkdir test4
pushd test4 > /dev/null

echo "> Run Test 4: multiple options-2"
"$meta_executable" \
  --fastq "$meta_resources_dir/test_data/test1.fastq" \
  --maxlength 40000 \
  --only_report \
  --raw \
  --outdir output

# Check if output directory exists
if [[ ! -d output ]]; then
  echo "Output directory not found!"
  exit 1
fi

# Check if output files are generated
if [ "$(ls -1 "output" | wc -l)" -ne 4 ]; then # 4 output files
  echo "Output files are not found!"
  exit 1
fi

# Check if the extracted data exists (--raw)
if ! ls output/*.tsv.gz > /dev/null 2>&1; then
  echo "Extracted data is not found!"
  exit 1
fi

# Check if files are empty
if find output -name "NanoPlot-report.html" -type f -size 0 | grep -q .; then 
  echo "NanoPlot report is empty."
  exit 1
fi
if find output -name "*.txt" -type f -size 0 | grep -q .; then 
  echo "NanoPlot summary file is empty."
  exit 1
fi
if find output -name "*.tsv.gz" -type f -size 0 | grep -q .; then 
  echo "Extracted data is empty."
  exit 1
fi

popd > /dev/null

echo "Test 4 succeeded."

###########################################################################

# Test 5: Run NanoPlot with different input (Fasta)

mkdir test5
pushd test5 > /dev/null

echo "> Run Test 5: Input Fasta"
"$meta_executable" \
  --fasta "$meta_resources_dir/test_data/test.fasta" \
  --outdir output

# Check if output directory exists
if [[ ! -d output ]]; then
  echo "Output directory not found!"
  exit 1
fi

# Check if output files are generated
if [ "$(ls -1 "output" | wc -l)" -lt 1 ]; then # Apart from log file
  echo "Output files are not found!"
  exit 1
fi

# Check if files are empty
if find output -name "*.html" -type f -size 0 | grep -q .; then 
  echo "At least one HTML file is empty."
  exit 1
fi
if find output -name "*.png" -type f -size 0 | grep -q .; then 
  echo "At least one plot is empty."
  exit 1
fi
if find output -name "*.txt" -type f -size 0 | grep -q .; then 
  echo "NanoPlot summary file is empty."
  exit 1
fi

popd > /dev/null

echo "Test 5 succeeded."

###########################################################################

# Test 6: Run NanoPlot with different input (Fastq_rich)

mkdir test6
pushd test6 > /dev/null

echo "> Run Test 6: Input Fastq_rich"
"$meta_executable" \
  --fastq_rich "$meta_resources_dir/test_data/test_rich.fastq" \
  --outdir output

# Check if output directory exists
if [[ ! -d output ]]; then
  echo "Output directory not found!"
  exit 1
fi

# Check if output files are generated
if [ "$(ls -1 "output" | wc -l)" -lt 1 ]; then # Apart from log file
  echo "Output files are not found!"
  exit 1
fi

# Check if files are empty
if find output -name "*.html" -type f -size 0 | grep -q .; then 
  echo "At least one HTML file is empty."
  exit 1
fi
if find output -name "*.png" -type f -size 0 | grep -q .; then 
  echo "At least one plot is empty."
  exit 1
fi
if find output -name "*.txt" -type f -size 0 | grep -q .; then 
  echo "NanoPlot summary file is empty."
  exit 1
fi

popd > /dev/null

echo "Test 6 succeeded."

###########################################################################

# Test 7: Run NanoPlot with different input (Fastq_minimal)

mkdir test7
pushd test7 > /dev/null

echo "> Run Test 7: Input Fasta"
"$meta_executable" \
  --fastq_minimal "../nanotest-master/reads.fastq.gz" \
  --outdir output

# Check if output directory exists
if [[ ! -d output ]]; then
  echo "Output directory not found!"
  exit 1
fi

# Check if output files are generated
if [ "$(ls -1 "output" | wc -l)" -lt 1 ]; then # Apart from log file
  echo "Output files are not found!"
  exit 1
fi

# Check if files are empty
if find output -name "*.html" -type f -size 0 | grep -q .; then 
  echo "At least one HTML file is empty."
  exit 1
fi
if find output -name "*.png" -type f -size 0 | grep -q .; then 
  echo "At least one plot is empty."
  exit 1
fi
if find output -name "*.txt" -type f -size 0 | grep -q .; then 
  echo "NanoPlot summary file is empty."
  exit 1
fi

popd > /dev/null

echo "Test 7 succeeded."

###########################################################################

# Test 8: Run NanoPlot with different input (Summary)

mkdir test8
pushd test8 > /dev/null

echo "> Run Test 8: Input Summary"
"$meta_executable" \
  --summary "$meta_resources_dir/test_data/test_summary.txt" \
  --outdir output

# Check if output directory exists
if [[ ! -d output ]]; then
  echo "Output directory not found!"
  exit 1
fi

# Check if output files are generated
if [ "$(ls -1 "output" | wc -l)" -lt 1 ]; then # Apart from log file
  echo "Output files are not found!"
  exit 1
fi

# Check if files are empty
if find output -name "*.html" -type f -size 0 | grep -q .; then 
  echo "At least one HTML file is empty."
  exit 1
fi
if find output -name "*.png" -type f -size 0 | grep -q .; then 
  echo "At least one plot is empty."
  exit 1
fi
if find output -name "*.txt" -type f -size 0 | grep -q .; then 
  echo "NanoPlot summary file is empty."
  exit 1
fi

popd > /dev/null

echo "Test 8 succeeded."

###########################################################################

# Test 9: Run NanoPlot with different input (BAM)

mkdir test9
pushd test9 > /dev/null

echo "> Run Test 9: Input BAM"
"$meta_executable" \
  --bam "$meta_resources_dir/test_data/test.bam" \
  --outdir output

# Check if output directory exists
if [[ ! -d output ]]; then
  echo "Output directory not found!"
  exit 1
fi

# Check if output files are generated
if [ "$(ls -1 "output" | wc -l)" -lt 1 ]; then # Apart from log file
  echo "Output files are not found!"
  exit 1
fi

# Check if files are empty
if find output -name "*.html" -type f -size 0 | grep -q .; then 
  echo "At least one HTML file is empty."
  exit 1
fi
if find output -name "*.png" -type f -size 0 | grep -q .; then 
  echo "At least one plot is empty."
  exit 1
fi
if find output -name "*.txt" -type f -size 0 | grep -q .; then 
  echo "NanoPlot summary file is empty."
  exit 1
fi

popd > /dev/null

echo "Test 9 succeeded."

###########################################################################

# Test 10: Run NanoPlot with different input (pickle)

mkdir test10
pushd test10 > /dev/null

echo "> Run Test 10: Input pickle"
"$meta_executable" \
  --pickle "../nanotest-master/alignment.pickle" \
  --outdir output

# Check if output directory exists
if [[ ! -d output ]]; then
  echo "Output directory not found!"
  exit 1
fi

# Check if output files are generated
if [ "$(ls -1 "output" | wc -l)" -lt 1 ]; then # Apart from log file
  echo "Output files are not found!"
  exit 1
fi

# Check if files are empty
if find output -name "*.html" -type f -size 0 | grep -q .; then 
  echo "At least one HTML file is empty."
  exit 1
fi
if find output -name "*.png" -type f -size 0 | grep -q .; then 
  echo "At least one plot is empty."
  exit 1
fi
if find output -name "*.txt" -type f -size 0 | grep -q .; then 
  echo "NanoPlot summary file is empty."
  exit 1
fi

popd > /dev/null

echo "Test 10 succeeded."

###########################################################################

# Test 11: Run NanoPlot with different input (feather)

mkdir test11
pushd test11 > /dev/null

echo "> Run Test 11: Input feather"
"$meta_executable" \
  --arrow "../nanotest-master/summary1.feather" \
  --outdir output

# Check if output directory exists
if [[ ! -d output ]]; then
  echo "Output directory not found!"
  exit 1
fi

# Check if output files are generated
if [ "$(ls -1 "output" | wc -l)" -lt 1 ]; then # Apart from log file
  echo "Output files are not found!"
  exit 1
fi

# Check if files are empty
if find output -name "*.html" -type f -size 0 | grep -q .; then 
  echo "At least one HTML file is empty."
  exit 1
fi
if find output -name "*.png" -type f -size 0 | grep -q .; then 
  echo "At least one plot is empty."
  exit 1
fi
if find output -name "*.txt" -type f -size 0 | grep -q .; then 
  echo "NanoPlot summary file is empty."
  exit 1
fi

popd > /dev/null

echo "Test 11 succeeded."

###########################################################################

# Test 12: Run NanoPlot with different output directory

mkdir test12
pushd test12 > /dev/null

echo "> Run Test 12: different output directory"
"$meta_executable" \
  --fastq "$meta_resources_dir/test_data/test1.fastq" \
  --outdir out

# Check if output directory exists
if [[ ! -d out ]]; then
  echo "Output directory not found!"
  exit 1
fi

# Check if output files are generated
if [ "$(ls -1 "out" | wc -l)" -lt 1 ]; then
  echo "Output files are not found!"
  exit 1
fi

# Check if files are empty
if find out -name "*.html" -type f -size 0 | grep -q .; then 
  echo "At least one HTML file is empty."
  exit 1
fi
if find out -name "*.png" -type f -size 0 | grep -q .; then 
  echo "At least one plot is empty."
  exit 1
fi
if find out -name "*.txt" -type f -size 0 | grep -q .; then 
  echo "NanoPlot summary file is empty."
  exit 1
fi

popd > /dev/null

echo "Test 12 succeeded."

###########################################################################

echo "All tests successfully completed!"