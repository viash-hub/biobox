# clone repo
if [ ! -d /tmp/snakemake-wrappers ]; then
  git clone --depth 1 --single-branch --branch master https://github.com/snakemake/snakemake-wrappers /tmp/snakemake-wrappers
fi

# copy test data
cp -r /tmp/snakemake-wrappers/bio/seqtk/test/reads/* src/seqtk/seqtk_subseq/test_data

# remove a.fastq file
rm src/seqtk/seqtk_subseq/test_data/a.fastq

# unzip fastq files
gunzip src/seqtk/seqtk_subseq/test_data/*.gz
