# busco test data

# Test data was obtained from https://github.com/snakemake/snakemake-wrappers/tree/master/bio/busco/test

if [ ! -d /tmp/snakemake-wrappers ]; then
  git clone --depth 1 --single-branch --branch master https://github.com/snakemake/snakemake-wrappers /tmp/snakemake-wrappers
fi

cp -r /tmp/snakemake-wrappers/bio/busco/test/protein.fasta src/busco/test_data