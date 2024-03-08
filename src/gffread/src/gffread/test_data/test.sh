#!/bin/bash

## VIASH START
## VIASH END

echo "> Run gffread with a single genome"
"$meta_executable" \
    --input "$meta_resources_dir/test_data/genome.fasta" \



# write a test script in this format:
'''
#!/bin/bash

## VIASH START
## VIASH END

echo "> Run arriba with blacklist"
"$meta_executable" \
  --bam "$meta_resources_dir/test_data/A.bam" \
  --genome "$meta_resources_dir/test_data/genome.fasta" \
  --gene_annotation "$meta_resources_dir/test_data/annotation.gtf" \
  --blacklist "$meta_resources_dir/test_data/blacklist.tsv" \
  --fusions "fusions.tsv" \
  --fusions_discarded "fusions_discarded.tsv" \
  --interesting_contigs "1,2"

echo ">> Checking output"
[ ! -f "fusions.tsv" ] && echo "Output file fusions.tsv does not exist" && exit 1
[ ! -f "fusions_discarded.tsv" ] && echo "Output file fusions_discarded.tsv does not exist" && exit 1

echo ">> Check if output is empty"
[ ! -s "fusions.tsv" ] && echo "Output file fusions.tsv is empty" && exit 1
[ ! -s "fusions_discarded.tsv" ] && echo "Output file fusions_discarded.tsv is empty" && exit 1
'''

# with this information:

'''
rule test_gffread:
    input:
        fasta="genome.fasta",
        annotation="annotation.gtf",
        # ids="",  # Optional path to records to keep
        # nids="",  # Optional path to records to drop
        # seq_info="",  # Optional path to sequence information
        # sort_by="",  # Optional path to the ordered list of reference sequences
        # attr="",  # Optional annotation attributes to keep.
        # chr_replace="",  # Optional path to <original_ref_ID> <new_ref_ID>
    output:
        records="transcripts.fa",
        # dupinfo="",  # Optional path to clustering/merging information
    threads: 1
    log:
        "logs/gffread.log",
    params:
        extra="",
    wrapper:
        "master/bio/gffread"
'''






        