#!/bin/bash

if [ ! -d /tmp/sortmerna_source ]; then
  git clone --depth 2 --single-branch --branch master https://github.com/snakemake/snakemake-wrappers.git /tmp/sortmerna_source
fi

# copy test data
cp -r /tmp/sortmerna_source/bio/sortmerna/test/* .
