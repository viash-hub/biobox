#!/bin/bash

## VIASH START
## VIASH END

test_dir="${meta_resources_dir}/test_data"

echo "> Run $meta_name with test data"
"$meta_executable" \
  --gff "$test_dir/0_test.gff" \
  --output "output.gtf" 

echo ">> Checking output"
[ ! -f "output.gtf" ] && echo "Output file output.gtf does not exist" && exit 1

echo ">> Check if output is empty"
[ ! -s "output.gtf" ] && echo "Output file output.gtf is empty" && exit 1

echo ">> Check if the conversion resulted in the right GTF format"
idGFF=$(head -n 2 "$test_dir/0_test.gff" | grep -o 'ID=[^;]*' | cut -d '=' -f 2-)
expectedGTF="gene_id \"$idGFF\"; ID \"$idGFF\";"
extractedGTF=$(head -n 3 "output.gtf" | grep -o 'gene_id "[^"]*"; ID "[^"]*";')
[ "$extractedGTF" != "$expectedGTF" ] && echo "Output file output.gtf does not have the right format" && exit 1

rm output.gtf

echo "> Run $meta_name with test data and GTF version 2.5"
"$meta_executable" \
  --gff "$test_dir/0_test.gff" \
  --output "output.gtf" \
  --gtf_version "2.5"

echo ">> Check if the output file header display the right GTF version"
grep -q "##gtf-version 2.5" "output.gtf"
[ $? -ne 0 ] && echo "Output file output.gtf header does not display the right GTF version" && exit 1

echo "> Test successful"