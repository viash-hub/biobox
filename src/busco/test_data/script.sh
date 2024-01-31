# busco test data

# Test data from https://github.com/snakemake/snakemake-wrappers/tree/master/bio/busco/test

if [ ! -d /tmp/snakemake-wrappers ]; then
  git clone --depth 1 --single-branch --branch master https://github.com/snakemake/snakemake-wrappers /tmp/snakemake-wrappers
fi

cp -r /tmp/snakemake-wrappers/bio/busco/test/protein.fasta src/busco/test_data

# Test data from https://gitlab.com/ezlab/busco_protocol/-/raw/main/protocol1/Tglobosa_GCF_014133895.1_genome.fna
wget -O src/busco/test_data/genome.fna https://gitlab.com/ezlab/busco_protocol/-/raw/main/protocol1/Tglobosa_GCF_014133895.1_genome.fna