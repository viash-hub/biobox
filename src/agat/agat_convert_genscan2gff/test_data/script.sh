#!/bin/bash

# clone repo
if [ ! -d /tmp/agat_source ]; then
  git clone --depth 1 --single-branch --branch master https://github.com/NBISweden/AGAT /tmp/agat_source
fi

# copy test data
cp -r /tmp/agat_source/t/scripts_output/in/test.genscan src/agat/agat_convert_genscan2gff/test_data/test.genscan
cp -r /tmp/agat_source/t/scripts_output/out/agat_convert_genscan2gff_1.gff src/agat/agat_convert_genscan2gff/test_data/agat_convert_genscan2gff_1.gff

