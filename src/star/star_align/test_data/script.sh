# arriba test data

# Test data was obtained from https://github.com/snakemake/snakemake-wrappers/tree/master/bio/arriba/test

if [ ! -d /tmp/snakemake-wrappers ]; then
  git clone --depth 1 --single-branch --branch master https://github.com/snakemake/snakemake-wrappers /tmp/snakemake-wrappers
fi

cp -r /tmp/snakemake-wrappers/bio/star/align/test/* src/star/star_align/test_data

# move reads
mv src/star/star_align/test_data/reads/* src/star/star_align/test_data
rm -r src/star/star_align/test_data/reads