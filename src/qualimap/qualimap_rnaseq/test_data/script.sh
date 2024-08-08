# qualimap test data

# Test data was obtained from https://github.com/snakemake/snakemake-wrappers/raw/master/bio/qualimap/rnaseq/test

if [ ! -d /tmp/snakemake-wrappers ]; then
  git clone --depth 1 --single-branch --branch master https://github.com/snakemake/snakemake-wrappers /tmp/snakemake-wrappers
fi

cp -r /tmp/snakemake-wrappers/bio/qualimap/rnaseq/test/mapped/a.bam src/qualimap/qualimap_rnaseq/test_data
cp -r /tmp/snakemake-wrappers/bio/qualimap/rnaseq/test/annotation.gtf src/qualimap/qualimap_rnaseq/test_data
