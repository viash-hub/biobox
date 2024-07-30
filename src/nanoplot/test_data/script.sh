echo -e "@SEQ_ID\nGATTTGGGGTTTTCCCAGTT\n+\n\!\!\!FFFEEEEDDDCCCGGGA@" > test.fastq

echo -e ">SEQ_ID\nGATTTGGGGTTTTCCCAGTT" > test.fasta

# echo -e ">ref\nGATTTGGGGTTTTCCCAGTT" > ref.fasta
# bwa index ref.fasta
# echo -e "@SEQ_ID\nGATTTGGGGTTTTCCCAGTT\n+\n'!!!FFFEEEEDDDCCCGGGA@'" > read.fastq
# bwa mem ref.fasta read.fastq > aln.sam
# samtools view -S -b aln.sam > test.bam
# samtools sort test.bam -o test_sorted.bam
# samtools index test_sorted.bam

echo -e "filename\tread_id\tlength\tnum_events\tnum_events_insert\tnum_events_delete\nread.fastq\tread_1\t1000\t5000\t100\t50" > test_summary.txt