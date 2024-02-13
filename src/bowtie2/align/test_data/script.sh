# bowtie2 align test data

# Test data was obtained from https://github.com/snakemake/snakemake-wrappers/tree/master/bio/bowtie2/align/test

if [ ! -d /tmp/snakemake-wrappers ]; then
  git clone --depth 1 --single-branch --branch master https://github.com/snakemake/snakemake-wrappers /tmp/snakemake-wrappers
fi

cp -r /tmp/snakemake-wrappers/bio/bowtie2/align/test/genome.fasta src/bowtie2/align/test_data

INDEX_DIR="src/bowtie2/align/test_data/index"

if [ ! -d "$INDEX_DIR" ]; then
  mkdir -p "$INDEX_DIR"
fi

cp -r /tmp/snakemake-wrappers/bio/bowtie2/align/test/index/* $INDEX_DIR

READS_DIR="src/bowtie2/align/test_data/reads"

if [ ! -d "$READS_DIR" ]; then
  mkdir -p "$READS_DIR"
fi

cp -r /tmp/snakemake-wrappers/bio/bowtie2/align/test/reads/* $READS_DIR