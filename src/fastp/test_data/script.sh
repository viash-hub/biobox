# fastp test data

# Test data was obtained from https://github.com/snakemake/snakemake-wrappers/tree/master/bio/fastp/test

if [ ! -d /tmp/snakemake-wrappers ]; then
  git clone --depth 1 --single-branch --branch master https://github.com/snakemake/snakemake-wrappers /tmp/snakemake-wrappers
fi

cp -r /tmp/snakemake-wrappers/bio/fastp/test/reads/* src/fastp/test_data

