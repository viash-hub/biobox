# clone repo
if [ ! -d /tmp/snakemake-wrappers ]; then
  git clone --depth 1 --single-branch --branch master https://github.com/snakemake/snakemake-wrappers /tmp/snakemake-wrappers
fi

# copy test data
cp -r /tmp/snakemake-wrappers/bio/seqtk/test/* src/seqtk/seqtk_sample/test_data

rm src/seqtk/seqtk_sample/test_data/Snakefile
