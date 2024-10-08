#!/bin/bash

# clone repo
if [ ! -d /tmp/agat_source ]; then
  git clone --depth 1 --single-branch --branch master https://github.com/NBISweden/AGAT /tmp/agat_source
fi

# copy test data
cp -r /tmp/agat_source/t/gff_syntax/in/0_test.gff src/agat/agat_convert_sp_gxf2gxf/test_data
cp -r /tmp/agat_source/t/gff_syntax/out/0_correct_output.gff src/agat/agat_convert_sp_gxf2gxf/test_data
