#!/bin/bash

set -e

## VIASH START
meta_executable="$PWD/target/executable/nanoplot/nanoplot"
meta_resources_dir="$PWD/src/nanoplot"
## VIASH END

###########################################################################

# Test 1: Run NanoPlot with only input parameter

echo "> Run Test 1: one input"
"$meta_executable" \
  --fastq "$meta_resources_dir/test_data/test1.fastq"

# Check if HTML file is generated
if ! ls target/result/*.html 1> /dev/null 2>&1; then
  echo "Output HTML file not found"
  exit 1
fi

# Check if summary file is generated
if ! ls target/result/*.txt 1> /dev/null 2>&1; then
  echo "NanoPlot summary file not found"
  exit 1
fi

# Check if plots are generated
if ! ls target/result/*.png 1> /dev/null 2>&1; then
  echo "Plots are not found"
  exit 1
fi

# Check if files are empty
if find target/result -name "*.html" -type f -size 0 | grep -q .; then 
  echo "At least one HTML file is empty."
  exit 1
fi
if find target/result -name "*.png" -type f -size 0 | grep -q .; then 
  echo "At least one plot is empty."
  exit 1
fi
if find target/result -name "*.txt" -type f -size 0 | grep -q .; then 
  echo "NanoPlot summary file is empty."
  exit 1
fi

# rm -f target/result/* # Clear output

echo "Test 1 succeeded."

###########################################################################

# Test 2: Run NanoPlot with multiple inputs

echo "> Run Test 2: multiple inputs"
"$meta_executable" \
  --fastq "$meta_resources_dir/test_data/test1.fastq" "$meta_resources_dir/test_data/test2.fastq"

# Check if HTML file is generated
if ! ls target/result/*.html 1> /dev/null 2>&1; then
  echo "Output HTML file not found"
  exit 1
fi

# Check if summary file is generated
if ! ls target/result/*.txt 1> /dev/null 2>&1; then
  echo "NanoPlot summary file not found"
  exit 1
fi

# Check if plots are generated
if ! ls target/result/*.png 1> /dev/null 2>&1; then
  echo "Plots are not found"
  exit 1
fi

# Check if files are empty
if find target/result -name "*.html" -type f -size 0 | grep -q .; then 
  echo "At least one HTML file is empty."
  exit 1
fi
if find target/result -name "*.png" -type f -size 0 | grep -q .; then 
  echo "At least one plot is empty."
  exit 1
fi
if find target/result -name "*.txt" -type f -size 0 | grep -q .; then 
  echo "NanoPlot summary file is empty."
  exit 1
fi

# rm -f target/result/* # Clear output

echo "Test 2 succeeded."

###########################################################################

# Test 3: Run NanoPlot with multiple options

echo "> Run Test 3: multiple options"
"$meta_executable" \
  --fastq "$meta_resources_dir/test_data/test1.fastq" \
  --maxlength 40000 \
  --format jpg \
  --prefix biobox_ \
  --store \
  --color yellow \
  --info_in_report

# Check if HTML file is generated
if ! ls target/result/*.html 1> /dev/null 2>&1; then
  echo "Output HTML file not found"
  exit 1
fi

# Check if summary file is generated
if ! ls target/result/*.txt 1> /dev/null 2>&1; then
  echo "NanoPlot summary file not found"
  exit 1
fi

# Check if plots are generated (.jpg)
if ! ls target/result/*.jpg 1> /dev/null 2>&1; then
  echo "Plots are not found"
  exit 1
fi

# Check if the extracted data exists (--store)
if ! ls target/result/*.pickle 1> /dev/null 2>&1; then
  echo "Extracted data are not found"
  exit 1
fi

# Check if files are empty
if find target/result -name "*.html" -type f -size 0 | grep -q .; then 
  echo "At least one HTML file is empty."
  exit 1
fi
if find target/result -name "*.jpg" -type f -size 0 | grep -q .; then 
  echo "At least one plot is empty."
  exit 1
fi
if find target/result -name "*.txt" -type f -size 0 | grep -q .; then 
  echo "NanoPlot summary file is empty."
  exit 1
fi
if find target/result -name "*.pickle" -type f -size 0 | grep -q .; then 
  echo "Extracted data is empty."
  exit 1
fi

# Check if the output file starts with "biobox" prefic
if ! ls target/result/biobox* 1> /dev/null 2>&1; then
    echo "The prefix is not added to the output files."
else

# rm -f target/result/* # Clear output

echo "Test 3 succeeded."

###########################################################################

# Test 4: Run NanoPlot with output options

echo "> Run Test 4: output options"
"$meta_executable" \
  --fastq "$meta_resources_dir/test_data/test1.fastq" \
  --html "out/*.html" \
  --log "out/*.log" \
  --statsum "out/*.txt"

# Check if HTML file is generated
if ! ls out/*.html 1> /dev/null 2>&1; then
  echo "Output HTML file not found"
  exit 1
fi

# Check if summary file is generated
if ! ls out/*.txt 1> /dev/null 2>&1; then
  echo "NanoPlot summary file not found"
  exit 1
fi

# Check if plots are generated (.jpg)
if ! ls target/result/*.png 1> /dev/null 2>&1; then
  echo "Plots are not found"
  exit 1
fi

# Check if files are empty
if find out -name "*.html" -type f -size 0 | grep -q .; then 
  echo "At least one HTML file is empty."
  exit 1
fi
if find target/result -name "*.png" -type f -size 0 | grep -q .; then 
  echo "At least one plot is empty."
  exit 1
fi
if find out -name "*.txt" -type f -size 0 | grep -q .; then 
  echo "NanoPlot summary file is empty."
  exit 1
fi

# rm -f target/result/* # Clear output
# rm -rf out/ # Delete user-defined directory

echo "Test 4 succeeded."

###########################################################################

echo "All tests successfully completed!"