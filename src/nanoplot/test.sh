#!/bin/bash

set -eo pipefail

## VIASH START
meta_executable="$PWD/target/executable/nanoplot/nanoplot"
meta_resources_dir="$PWD/src/nanoplot"
## VIASH END

###########################################################################

# Test 1: Run NanoPlot with only input parameter

echo "> Run Test 1: one input"
"$meta_executable" \
  --fastq "$meta_resources_dir/test_data/test1.fastq"

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

rm -f output/* # Clear output

echo "Test 1 succeeded."

###########################################################################

# Test 2: Run NanoPlot with multiple inputs

echo "> Run Test 2: multiple inputs"
"$meta_executable" \
  --fastq "$meta_resources_dir/test_data/test1.fastq" "$meta_resources_dir/test_data/test2.fastq"

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

rm -f output/* # Clear output

echo "Test 2 succeeded."

###########################################################################

# Test 3: Run NanoPlot with multiple options-1

echo "> Run Test 3: multiple options-1"
"$meta_executable" \
  --fastq "$meta_resources_dir/test_data/test1.fastq" \
  --maxlength 40000 \
  --format jpg \
  --prefix biobox_ \
  --store \
  --color yellow \
  --info_in_report

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

rm -f output/* # Clear output

echo "Test 3 succeeded."

###########################################################################

# Test 4: Run NanoPlot with multiple options-2

echo "> Run Test 4: multiple options-2"
"$meta_executable" \
  --fastq "$meta_resources_dir/test_data/test1.fastq" \
  --maxlength 40000 \
  --only_report \
  --raw

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

rm -f output/* # Clear output

echo "Test 4 succeeded."

###########################################################################

# Test 5: Run NanoPlot with different input

echo "> Run Test 5: different input"
"$meta_executable" \
  --bam "$meta_resources_dir/test_data/test.bam"

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

rm -f output/* # Clear output

echo "Test 5 succeeded."

###########################################################################

# Test 6: Run NanoPlot with output options

echo "> Run Test 6: output options"
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

rm -rf out/ # Clear output

echo "Test 6 succeeded."

###########################################################################

echo "All tests successfully completed!"