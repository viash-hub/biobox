#!/bin/bash

# clone repo
if [ ! -d /tmp/agat_source ]; then
  git clone --depth 1 --single-branch --branch master https://github.com/NBISweden/AGAT /tmp/agat_source
fi

# copy test data
cp -r /tmp/agat_source/t/scripts_output/in/1.gff src/agat/agat_sq_stat_basic/test_data/
cp -r /tmp/agat_source/t/scripts_output/out/agat_sq_stat_basic_1.gff src/agat/agat_sq_stat_basic/test_data/