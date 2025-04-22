#!/bin/bash

echo ">>> Testing $meta_name"

find $meta_resources_dir/test_data/rRNA -type f > test_data/rrna-db.txt

echo ">>> Testing for paired-end reads and database manifest"
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

################################################################################
echo ">>> Testing for paired-end reads and --ref and --paired_out arguments"
"$meta_executable" \
    --output "rRNA_reads" \
    --other "non_rRNA_reads" \
    --input "$meta_resources_dir/test_data/reads_1.fq.gz;$meta_resources_dir/test_data/reads_2.fq.gz" \
    --ref "$meta_resources_dir/test_data/rRNA/database1.fa;$meta_resources_dir/test_data/rRNA/database2.fa" \
    --log test_log.log \
    --paired_out \
    --fastx \
    --out2

echo ">> Checking if the correct files are present"
[[ -f "rRNA_reads_fwd.fq.gz" ]] || [[ -f "rRNA_reads_rev.fq.gz" ]] || { echo "rRNA output fastq file is missing!"; exit 1; }
gzip -dkf rRNA_reads_fwd.fq.gz
[[ ! -s "rRNA_reads_fwd.fq" ]] && [[ ! -s "rRNA_reads_rev.fq" ]] || { echo "rRNA output fastq file is not empty!"; exit 1; }
[[ -f "non_rRNA_reads_fwd.fq.gz" ]] || [[ -f "non_rRNA_reads_rev.fq.gz" ]] || { echo "Non-rRNA output fastq file is missing!"; exit 1;}
gzip -dkf non_rRNA_reads_fwd.fq.gz
gzip -dkf non_rRNA_reads_rev.fq.gz
[[ -s "non_rRNA_reads_fwd.fq" ]] && [[ -s "non_rRNA_reads_rev.fq" ]] || { echo "Non-rRNA output fastq file is empty!"; exit 1; }

rm -f rRNA_reads_fwd.fq.gz rRNA_reads_rev.fq.gz non_rRNA_reads_fwd.fq.gz non_rRNA_reads_rev.fq.gz test_log.log
rm -rf kvdb/

################################################################################

echo ">>> Testing for single-end reads and --ref argument"
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

rm -f rRNA_reads.fq.gz non_rRNA_reads.fq.gz test_log.log
rm -rf kvdb/

################################################################################

echo ">>> Testing for single-end reads with singleton output files"
"$meta_executable" \
    --aligned "rRNA_reads" \
    --other "non_rRNA_reads" \
    --input "$meta_resources_dir/test_data/reads_1.fq.gz;$meta_resources_dir/test_data/reads_2.fq.gz" \
    --ribo_database_manifest test_data/rrna-db.txt \
    --log test_log.log \
    --fastx \
    --sout

echo ">> Checking if the correct files are present"
[[ ! -f "rRNA_reads_paired.fq.gz" ]] && echo "Aligned paired fwd output fastq file is missing!" && exit 1
[[ ! -f "rRNA_reads_singleton.fq.gz" ]] && echo "Aligned singleton fwd output fastq file is missing!" && exit 1
[[ ! -f "non_rRNA_reads_fwd.fq" ]] && echo "Non-rRNA fwd output fastq file is missing!" && exit 1
[[ ! -f "non_rRNA_reads_rev.fq" ]] && echo "Non-rRNA rev output fastq file is missing!" && exit 1
[[ ! -f "non_rRNA_reads_singleton.fq.gz" ]] && echo "Non-rRNA singleton output fastq file is missing!" && exit 1
[[ ! -f "non_rRNA_reads_paired.fq.gz" ]] && echo "Non-rRNA paired output fastq file is missing!" && exit 1



echo ">>> All tests passed"
exit 0