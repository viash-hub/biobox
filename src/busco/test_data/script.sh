# busco test data

# Test data from https://github.com/snakemake/snakemake-wrappers/tree/master/bio/busco/test

if [ ! -d /tmp/snakemake-wrappers ]; then
  git clone --depth 1 --single-branch --branch master https://github.com/snakemake/snakemake-wrappers /tmp/snakemake-wrappers
fi

cp -r /tmp/snakemake-wrappers/bio/busco/test/protein.fasta src/busco/test_data

# Test data from busco test data at https://gitlab.com/ezlab/busco/-/tree/master/test_data?ref_type=heads
wget -O src/busco/test_data/genome.fna "https://gitlab.com/ezlab/busco/-/raw/master/test_data/eukaryota/genome.fna?ref_type=heads&inline=false"