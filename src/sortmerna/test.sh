#!/bin/bash

echo ">>> Testing $meta_functionality_name"

find $meta_resources_dir/test_data/rRNA -type f > test_data/rrna-db.txt

echo ">>> Testing for paired-end reads"
# out2 separates the read pairs into two files (one fwd and one rev)
# paired_in outputs both reads of a pair
# other is the output file for non-rRNA reads
"$meta_executable" \
    --output "rRNA_reads" \
    --other "non_rRNA_reads" \
    --input "$meta_resources_dir/test_data/reads_1.fq.gz;$meta_resources_dir/test_data/reads_2.fq.gz" \
    --ribo_database_manifest test_data/rrna-db.txt \
    --log test_log.log \
    --paired_in \
    --fastx \
    --out2
    

echo ">> Checking if the correct files are present"
[[ -f "rRNA_reads_fwd.fq.gz" ]] || [[ -f "rRNA_reads_rev.fq.gz" ]] || { echo "rRNA output fastq file is missing!"; exit 1; }
[[ -s "rRNA_reads_fwd.fq.gz" ]] && [[ -s "rRNA_reads_rev.fq.gz" ]] || { echo "rRNA output fastq file is empty!"; exit 1; }
[[ -f "non_rRNA_reads_fwd.fq.gz" ]] || [[ -f "non_rRNA_reads_rev.fq.gz" ]] || { echo "Non-rRNA output fastq file is missing!"; exit 1;}
gzip -dk non_rRNA_reads_fwd.fq.gz
gzip -dk non_rRNA_reads_rev.fq.gz
[[ ! -s "non_rRNA_reads_fwd.fq" ]] && [[ ! -s "non_rRNA_reads_rev.fq" ]] || { echo "Non-rRNA output fastq file is not empty!"; exit 1;}

rm -f rRNA_reads_fwd.fq.gz rRNA_reads_rev.fq.gz non_rRNA_reads_fwd.fq.gz non_rRNA_reads_rev.fq.gz test_log.log
rm -rf kvdb/


echo ">>> Testing for single-end reads"
"$meta_executable" \
    --aligned "rRNA_reads" \
    --other "non_rRNA_reads" \
    --input $meta_resources_dir/test_data/reads_1.fq.gz \
    --ref $meta_resources_dir/test_data/rRNA/database1.fa \
    --log test_log.log \
    --fastx

echo ">> Checking if the correct files are present"
[[ ! -f "rRNA_reads.fq.gz" ]] && echo "rRNA output fastq file is missing!" && exit 1
gzip -dk rRNA_reads.fq.gz
[[ -s "rRNA_reads.fq" ]] && echo "rRNA output fastq file is not empty!" && exit 1
[[ ! -f "non_rRNA_reads.fq.gz" ]] && echo "Non-rRNA output fastq file is missing!" && exit 1
[[ ! -s "non_rRNA_reads.fq.gz" ]] && echo "Non-rRNA output fastq file is empty!" && exit 1


echo ">>> All tests passed"
exit 0