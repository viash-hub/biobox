wget https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Pancreas/Visium_HD_Human_Pancreas_image.tif -O src/spaceranger/spaceranger_count/test_data/Visium_HD_Human_Pancreas_image.tif
wget https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Pancreas/Visium_HD_Human_Pancreas_probe_set.csv -O src/spaceranger/spaceranger_count/test_data/Visium_HD_Human_Pancreas_probe_set.csv
wget https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Pancreas/Visium_HD_Human_Pancreas_tissue_image.btf -O src/spaceranger/spaceranger_count/test_data/Visium_HD_Human_Pancreas_tissue_image.btf
wget https://s3-us-west-2.amazonaws.com/10x.files/samples/spatial-exp/3.0.0/Visium_HD_Human_Pancreas/Visium_HD_Human_Pancreas_fastqs.tar -O src/spaceranger/spaceranger_count/test_data/Visium_HD_Human_Pancreas_fastqs.tar

tar -xvf src/spaceranger/spaceranger_count/test_data/Visium_HD_Human_Pancreas_fastqs.tar -C src/spaceranger/spaceranger_count/test_data/

# seqtk sample -s100 src/spaceranger/spaceranger_count/test_data/Visium_HD_Human_Pancreas_fastqs/Visium_HD_Human_Pancreas_S1_L001_I1_001.fastq.gz 0.1 > src/spaceranger/spaceranger_count/test_data/subsampled_fastqs/Visium_HD_Human_Pancreas_S1_L001_I1_001.fastq.gz

# Input directory containing the original FASTQs
input_dir="src/spaceranger/spaceranger_count/test_data/Visium_HD_Human_Pancreas_fastqs"
# Output directory for the subsampled FASTQs
output_dir="src/spaceranger/spaceranger_count/test_data/subsampled_fastqs"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop over all FASTQ files in the input directory
for fastq_file in "$input_dir"/*.fastq.gz; do
  # Get the base filename (without directory)
  base_filename=$(basename "$fastq_file")
  
  # Define the output file path
  output_file="$output_dir/$base_filename"

  # Run seqtk on the current file (subsampling it)
  seqtk sample -s100 "$fastq_file" 1000 > "$output_file"
  
  # Optionally, print a message about the processed file
  echo "Subsampled $fastq_file to $output_file"
done


# # First find the line number where chr21 starts and chr22 starts
# chr21_start=$(grep -n '^>chr21' /home/jakubmajercik/Data_Intuitive/refdata-gex-GRCh38-2020-A/fasta/genome.fa | cut -d: -f1)
# chr22_start=$(grep -n '^>chr22' /home/jakubmajercik/Data_Intuitive/refdata-gex-GRCh38-2020-A/fasta/genome.fa | cut -d: -f1)

# # If chr22 doesn't exist (chr21 is last), get total line count instead
# if [ -z "$chr22_start" ]; then
#     chr22_start=$(wc -l < /home/jakubmajercik/Data_Intuitive/refdata-gex-GRCh38-2020-A/fasta/genome.fa)
# fi

# # Extract lines between chr21 and chr22
# # sed -n "${chr21_start},${chr22_start}p" /home/jakubmajercik/Data_Intuitive/refdata-gex-GRCh38-2020-A/fasta/genome.fa > chr21_subset.fa

# # Extract chr21 with header and sequence
# (             
# echo "# Subset of genome containing only chr21"
# awk '/^>chr21\s/{p=1;print;next} /^>/{p=0} p' /home/jakubmajercik/Data_Intuitive/refdata-gex-GRCh38-2020-A/fasta/genome.fa > src/spaceranger/spaceranger_count/test_data/chr21_subset.fa
# )    

# # Filter gene types
# spaceranger mkgtf src/spaceranger/spaceranger_count/test_data/chr21_subset.gtf src/spaceranger/spaceranger_count/test_data/filtered_chr21.gtf \
#     --attribute=gene_biotype:protein_coding \
#     --attribute=gene_biotype:lincRNA \
#     --attribute=gene_biotype:antisense

# # Create reference
# spaceranger mkref --genome=GRCh38_chr21 \
#     --fasta=src/spaceranger/spaceranger_count/test_data/chr21_subset.fa \
#     --genes=src/spaceranger/spaceranger_count/test_data/filtered_chr21.gtf \
#     --nthreads=12 \
#     --memgb=32

# # grep 'chr21' src/spaceranger/spaceranger_count/test_data/Visium_HD_Human_Pancreas_probe_set.csv > src/spaceranger/spaceranger_count/test_data/chr21_probe_set.csv

# # 1. Extract Ensembl gene IDs from chr21 GTF
# cat src/spaceranger/spaceranger_count/test_data/chr21_subset.gtf | grep -o 'gene_id "[^"]*"' | sort -u | cut -d'"' -f2 > src/spaceranger/spaceranger_count/test_data/chr21_genes.txt

# # 2. Create new probe set file with header and chr21 probes
# (
#     # Copy header lines (starting with #)
#     grep '^#' src/spaceranger/spaceranger_count/test_data/Visium_HD_Human_Pancreas_probe_set.csv > src/spaceranger/spaceranger_count/test_data/chr21_probe_set.csv
#     # Copy column header
#     head -n 6 src/spaceranger/spaceranger_count/test_data/Visium_HD_Human_Pancreas_probe_set.csv | tail -n 1 >> src/spaceranger/spaceranger_count/test_data/chr21_probe_set.csv
#     # Filter probes for chr21 genes
#     while read gene_id; do
#         grep "^$gene_id," src/spaceranger/spaceranger_count/test_data/Visium_HD_Human_Pancreas_probe_set.csv >> src/spaceranger/spaceranger_count/test_data/chr21_probe_set.csv
#     done < src/spaceranger/spaceranger_count/test_data/chr21_genes.txt
# )







# # 1. Extract all gene IDs from probe set (excluding header lines)
# tail -n +7 src/spaceranger/spaceranger_count/test_data/Visium_HD_Human_Pancreas_probe_set.csv | cut -d',' -f1 | sort -u > probeset_genes.txt

# # 2. Extract chr21 entries from original GTF that match these genes
# cat /home/jakubmajercik/Data_Intuitive/refdata-gex-GRCh38-2020-A/genes/genes.gtf | awk -F'\t' '$1=="chr21"' | grep -f probeset_genes.txt > new_chr21.gtf

# # 3. Create new filtered GTF with these genes
# spaceranger mkgtf new_chr21.gtf filtered_new_chr21.gtf \
#     --attribute=gene_biotype:protein_coding \
#     --attribute=gene_biotype:lincRNA \
#     --attribute=gene_biotype:antisense

# # 4. Create new reference with this GTF
# spaceranger mkref --genome=GRCh38_chr21_v2 \
#     --fasta=src/spaceranger/spaceranger_count/test_data/chr21_subset.fa \
#     --genes=filtered_new_chr21.gtf \
#     --nthreads=12 \
#     --memgb=32