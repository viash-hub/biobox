#!/bin/bash

## Fastq file ##
# Define the number of reads
NUM_READS=10
OUTPUT_FILE="./src/nanoplot/test_data/test1.fastq"

# Function to generate a random DNA sequence of given length
generate_sequence() {
    local length=$1 #assigns it the value of the first argument passed to the function 
    cat /dev/urandom | tr -dc 'ACGT' | fold -w $length | head -n 1
}

# Function to generate random quality scores of given length
generate_quality() {
    local length=$1
    local average_quality=$2
    local quality=""
    for ((i=0; i<length; i++)); do
        # Generate a quality score based on the average_quality
        quality+=$(awk -v avg=$average_quality 'BEGIN {printf "%c", int(rand()*10 + avg) + 33}')
    done
    echo $quality
}

echo -n "" > $OUTPUT_FILE #Create the fastq file
for i in $(seq 1 $NUM_READS); do
    # Randomly determine the read length (between 20 and 100 bases)
    read_length=$(shuf -i 20-100 -n 1)
    # Randomly determine the average quality (between 30 and 40)
    average_quality=$(shuf -i 0-40 -n 1)
    sequence=$(generate_sequence $read_length)
    quality=$(generate_quality $read_length $average_quality)
    echo "@read_$i" >> $OUTPUT_FILE
    echo $sequence >> $OUTPUT_FILE
    echo "+" >> $OUTPUT_FILE
    echo $quality >> $OUTPUT_FILE
    echo >> $OUTPUT_FILE  # Add a blank line between reads
done

NUM_READS=7
OUTPUT_FILE="./src/nanoplot/test_data/test2.fastq"
echo -n "" > $OUTPUT_FILE #Create another fastq file
for i in $(seq 1 $NUM_READS); do
    # Randomly determine the read length (between 20 and 100 bases)
    read_length=$(shuf -i 20-100 -n 1)
    # Randomly determine the average quality (between 30 and 40)
    average_quality=$(shuf -i 0-40 -n 1)
    sequence=$(generate_sequence $read_length)
    quality=$(generate_quality $read_length $average_quality)
    echo "@read_$i" >> $OUTPUT_FILE
    echo $sequence >> $OUTPUT_FILE
    echo "+" >> $OUTPUT_FILE
    echo $quality >> $OUTPUT_FILE
    echo >> $OUTPUT_FILE  # Add a blank line between reads
done

#########################################################################################

## Fasta file ##
wget -O src/nanoplot/test_data/test.fasta https://raw.githubusercontent.com/merenlab/reads-for-assembly/master/examples/files/fasta_01.fa

#########################################################################################

## Fastq_rich file ##
if [ ! -d nanotest ]; then
  git clone --depth 1 --single-branch --branch master https://github.com/wdecoster/nanotest/
fi

mv nanotest/reads-mixed-timestamp.fastq src/nanoplot/test_data/test_rich.fastq 

#########################################################################################

## Fastq_minimal file ##
mv nanotest/reads.fastq.gz src/nanoplot/test_data/test_minimal.fastq.gz

# Unzip file
gunzip -c src/nanoplot/test_data/test_minimal.fastq.gz > src/nanoplot/test_data/test_minimal.fastq

# Extract the first 5 reads (first 20 lines)
head -n 20 src/nanoplot/test_data/test_minimal.fastq > first_5_reads.fastq

mv first_5_reads.fastq src/nanoplot/test_data/test_minimal.fastq

# Zip the output file
gzip -f src/nanoplot/test_data/test_minimal.fastq

#########################################################################################

## Summary file ##
mv nanotest/sequencing_summary.txt src/nanoplot/test_data/test_summary.txt

#########################################################################################

## Bam file ##
if [ ! -d /tmp/snakemake-wrappers ]; then
  git clone --depth 1 --single-branch --branch master https://github.com/snakemake/snakemake-wrappers /tmp/snakemake-wrappers
fi

cp -r /tmp/snakemake-wrappers/bio/biobambam2/bamsormadup/test/mapped/a.bam src/nanoplot/test_data/test.bam

#########################################################################################

## Pickle file ##
mv nanotest/alignment.pickle src/nanoplot/test_data/test.pickle

#########################################################################################

## Feather/Arrow file ##
mv nanotest/summary1.feather src/nanoplot/test_data/test.feather

rm -rf nanotest