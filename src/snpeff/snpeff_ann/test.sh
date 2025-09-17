#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

###########################################################################

# create temporary directory and clean up on exit
TMPDIR=$(mktemp -d "$meta_temp_dir/$meta_name-XXXXXX")
echo "> Created $TMPDIR"
function clean_up {
  [[ -d "$TMPDIR" ]] && rm -rf "$TMPDIR"
}
trap clean_up EXIT
DATA_DIR="$TMPDIR/data"
mkdir "$DATA_DIR"
TEST_GENOME="$DATA_DIR/test"
mkdir "$TEST_GENOME"
echo "Downloading test data"
curl -o - https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz > "$TEST_GENOME/sequences.fa.gz"
curl -o - https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.chr.gtf.gz | zcat | sed -n '/^[#1]\s/p' | gzip > "$TEST_GENOME/genes.gtf.gz"
snpEff build -dataDir "$DATA_DIR" -noCheckCds -noCheckProtein -gtf22 -noLog -configOption "test.genome=GRCh38Chr1" -v test 

# Test 1: Run SnpEff with only required parameters

mkdir test1
pushd test1 > /dev/null # cd test1 (stack)

echo "> Run Test 1: required parameters"
"$meta_executable" \
  --genome_version test \
  --data_dir "$DATA_DIR" \
  --config_option "test.genome=GRCh38Chr1" \
  --input "$meta_resources_dir/test_data/cancer.vcf" \
  --output out.vcf

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
  --genome_version test \
  --data_dir "$DATA_DIR" \
  --config_option "test.genome=GRCh38Chr1" \
  --input "$meta_resources_dir/test_data/test.vcf" \
  --interval "$meta_resources_dir/test_data/my_annotations.bed" \
  --no_stats \
  --output output.vcf

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
  --genome_version test \
  --input "$meta_resources_dir/test_data/test.vcf" \
  --output output.vcf \
  --data_dir "$DATA_DIR" \
  --config_option "test.genome=GRCh38Chr1" \
  --summary temp \
  --genes temp

# Check if output.vcf exists
if [ ! -e "output.vcf" ]; then
    echo "File output.vcf does not exist."
fi

# Check if the other output files have been moved to temp folder
output_files=("snpEff_genes.txt" "snpEff_summary.html")

# Check if any of the files do not exist
for file in "${output_files[@]}"; do
    if [ ! -e "temp/$file" ]; then
        echo "File $file does not exist in 'temp' folder."
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