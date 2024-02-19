# pear test data

# Test data was obtained from https://github.com/snakemake/snakemake-wrappers/tree/master/bio/lofreq/call/test/data

if [ ! -d /tmp/snakemake-wrappers ]; then
  git clone --depth 1 --single-branch --branch master https://github.com/snakemake/snakemake-wrappers /tmp/snakemake-wrappers
fi

cp -r /tmp/snakemake-wrappers/bio/lofreq/call/test/data/* src/lofreq/call/test_data

