#!/bin/bash

# clone repo
if [ ! -d /tmp/agat_source ]; then
  git clone --depth 1 --single-branch --branch master https://github.com/NBISweden/AGAT /tmp/agat_source
fi

# copy test data
cp -r /tmp/agat_source/t/scripts_output/in/1.gff src/agat/agat_sp_add_intergenic_regions/test_data
cp -r /tmp/agat_source/t/scripts_output/out/agat_sp_add_intergenic_regions_1.gff src/agat/agat_sp_add_intergenic_regions/test_data
