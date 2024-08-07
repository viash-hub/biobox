# echo -e "@SEQ_ID\nGATTTGGGGTTTTCCCAGTT\n+\n\!\!\!FFFEEEEDDDCCCGGGA@" > ./src/nanoplot/test_data/test.fastq

# echo -e ">SEQ_ID\nGATTTGGGGTTTTCCCAGTT" > ./src/nanoplot/test_data/test.fasta

# echo -e "filename\tread_id\tlength\tnum_events\tnum_events_insert\tnum_events_delete\nread.fastq\tread_1\t1000\t5000\t100\t50" > ./src/nanoplot/test_data/test_summary.txt

#!/bin/bash

# Define the number of reads
NUM_READS=10
OUTPUT_FILE="sample_varying.fastq"

# Function to generate a random DNA sequence of given length
generate_sequence() {
    local length=$1
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

# Create the FASTQ file
echo -n "" > $OUTPUT_FILE
for i in $(seq 1 $NUM_READS); do
    # Randomly determine the read length (between 20 and 100 bases)
    read_length=$(shuf -i 20-100 -n 1)
    # Randomly determine the average quality (between 30 and 40)
    average_quality=$(shuf -i 30-40 -n 1)
    sequence=$(generate_sequence $read_length)
    quality=$(generate_quality $read_length $average_quality)
    echo "@read_$i" >> $OUTPUT_FILE
    echo $sequence >> $OUTPUT_FILE
    echo "+" >> $OUTPUT_FILE
    echo $quality >> $OUTPUT_FILE
done

echo "FASTQ file '$OUTPUT_FILE' created with $NUM_READS reads of varying lengths and quality."