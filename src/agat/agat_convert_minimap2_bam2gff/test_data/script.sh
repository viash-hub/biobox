#!/bin/bash

# clone repo
if [ ! -d /tmp/agat_source ]; then
  git clone --depth 1 --single-branch --branch master https://github.com/NBISweden/AGAT /tmp/agat_source
fi

# copy test data
cp -r /tmp/agat_source/t/scripts_output/in/test_minimap2.sam src/agat/agat_convert_minimap2_bam2gff/test_data/
cp -r /tmp/agat_source/t/scripts_output/out/agat_convert_minimap2_bam2gff_1.gff src/agat/agat_convert_minimap2_bam2gff/test_data/