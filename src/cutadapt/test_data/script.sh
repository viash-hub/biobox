# cutadapt test data

# Test data was obtained from https://github.com/snakemake/snakemake-wrappers/tree/master/bio/cutadapt/test

if [ ! -d /tmp/snakemake-wrappers ]; then
  git clone --depth 1 --single-branch --branch master https://github.com/snakemake/snakemake-wrappers /tmp/snakemake-wrappers
fi

mkdir -p src/cutadapt/test_data/pe
mkdir src/cutadapt/test_data/se

cp -r /tmp/snakemake-wrappers/bio/cutadapt/se/test/reads/* src/cutadapt/test_data/se
cp -r /tmp/snakemake-wrappers/bio/cutadapt/pe/test/reads/* src/cutadapt/test_data/pe

rm -rf /tmp/snakemake-wrappers
