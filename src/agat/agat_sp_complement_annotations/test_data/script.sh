#!/bin/bash

# clone repo
if [ ! -d /tmp/agat_source ]; then
  git clone --depth 1 --single-branch --branch master https://github.com/NBISweden/AGAT /tmp/agat_source
fi

# copy test data
cp -r /tmp/agat_source/t/gff_syntax/in/25_test.gff src/agat/agat_sp_complement_annotations/test_data
cp -r /tmp/agat_source/t/gff_syntax/in/9_test.gff src/agat/agat_sp_complement_annotations/test_data
cp -r /tmp/agat_source/t/scripts_output/in/agat_sp_complement_annotations/agat_sp_complement_annotations_ref.gff src/agat/agat_sp_complement_annotations/test_data
cp -r /tmp/agat_source/t/scripts_output/in/agat_sp_complement_annotations/agat_sp_complement_annotations_add.gff src/agat/agat_sp_complement_annotations/test_data
cp -r /tmp/agat_source/t/scripts_output/out/agat_sp_complement_annotations_1.gff src/agat/agat_sp_complement_annotations/test_data
cp -r /tmp/agat_source/t/scripts_output/out/agat_sp_complement_annotations_2.gff src/agat/agat_sp_complement_annotations/test_data