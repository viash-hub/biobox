# multiqc test data

# Test data from https://github.com/snakemake/snakemake-wrappers/tree/master/bio/busco/test

if [ ! -d /tmp/snakemake-wrappers ]; then
  git clone --depth 1 --single-branch --branch master https://github.com/snakemake/snakemake-wrappers /tmp/snakemake-wrappers
fi

cp -r /tmp/snakemake-wrappers/bio/multiqc/test/samtools_stats/* src/multiqc/test_data
