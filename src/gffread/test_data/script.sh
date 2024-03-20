#!/bin/bash

# clone repo
if [ ! -d /tmp/gffread_source ]; then
  git clone --depth 2 --single-branch --branch master https://github.com/gpertea/gffread.git /tmp/gffread_source
fi

# copy test data
cp -r /tmp/gffread_source/examples/* src/gffread/test_data
