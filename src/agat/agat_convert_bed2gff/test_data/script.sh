#!/bin/bash

# clone repo
if [ ! -d /tmp/agat_source ]; then
  git clone --depth 1 --single-branch --branch master https://github.com/NBISweden/AGAT /tmp/agat_source
fi

# copy test data
cp -r /tmp/agat_source/t/scripts_output/in/test.bed src/agat/agat_convert_bed2gff/test_data/test.bed
cp -r /tmp/agat_source/t/scripts_output/out/agat_convert_bed2gff_1.gff src/agat/agat_convert_bed2gff/test_data/agat_convert_bed2gff_1.gff