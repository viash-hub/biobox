# biobox x.x.x

## BREAKING CHANGES

* `star/star_align_reads`: Change all arguments from `--camelCase` to `--snake_case` (PR #62).

* `star/star_genome_generate`: Change all arguments from `--camelCase` to `--snake_case` (PR #62).

## NEW FUNCTIONALITY

* `star/star_align_reads`: Add star solo related arguments (PR #62).

* `bd_rhapsody/bd_rhapsody_make_reference`: Create a reference for the BD Rhapsody pipeline (PR #75).

* `umitools/umitools_dedup`: Deduplicate reads based on the mapping co-ordinate and the UMI attached to the read (PR #54).

* `seqtk`:
  - `seqtk/seqtk_sample`: Subsamples sequences from FASTA/Q files (PR #68).
  - `seqtk/seqtk_subseq`: Extract the sequences (complete or subsequence) from the FASTA/FASTQ files
                based on a provided sequence IDs or region coordinates file (PR #85).

* `agat`:
  - `agat_convert_sp_gff2gtf`: convert any GTF/GFF file into a proper GTF file (PR #76).
  - `/agat_convert_bed2gff`: convert bed file to gff format (PR #97).

* `bedtools`:
  - `bedtools/bedtools_intersect`: Allows one to screen for overlaps between two sets of genomic features (PR #94).
  - `bedtools/bedtools_sort`: Sorts a feature file (bed/gff/vcf) by chromosome and other criteria (PR #98).


* `agat/agat_convert_embl2gff`: convert an EMBL file into GFF format (PR #99).

## MINOR CHANGES

* `busco` components: update BUSCO to `5.7.1` (PR #72).

* Update CI to reusable workflow in `viash-io/viash-actions` (PR #86).

* Update several components in order to avoid duplicate code when using `unset` on boolean arguments (PR #133).

* Bump viash to `0.9.0-RC7` (PR #134)

## DOCUMENTATION

* Extend the contributing guidelines (PR #82):

  - Update format to Viash 0.9.

  - Descriptions should be formatted in markdown.

  - Add defaults to descriptions, not as a default of the argument.

  - Explain parameter expansion.

  - Mention that the contents of the output of components in tests should be checked.

* Add authorship to existing components (PR #88).

## BUG FIXES

* `pear`: fix component not exiting with the correct exitcode when PEAR fails (PR #70).

* `cutadapt`: fix `--par_quality_cutoff_r2` argument (PR #69).

* `cutadapt`: demultiplexing is now disabled by default. It can be re-enabled by using `demultiplex_mode` (PR #69).

* `multiqc`: update multiple separator to `;` (PR #81).


# biobox 0.1.0

## NEW FEATURES

* `arriba`: Detect gene fusions from RNA-seq data (PR #1).

* `fastp`: An ultra-fast all-in-one FASTQ preprocessor (PR #3).

* `busco`: 
    - `busco/busco_run`: Assess genome assembly and annotation completeness with single copy orthologs (PR #6).
    - `busco/busco_list_datasets`: Lists available busco datasets (PR #18).
    - `busco/busco_download_datasets`: Download busco datasets (PR #19).

* `cutadapt`: Remove adapter sequences from high-throughput sequencing reads (PR #7).

* `featurecounts`: Assign sequence reads to genomic features (PR #11).

* `bgzip`: Add bgzip functionality to compress and decompress files (PR #13).

* `pear`: Paired-end read merger (PR #10).

* `lofreq/call`: Call variants from a BAM file (PR #17).

* `lofreq/indelqual`: Insert indel qualities into BAM file (PR #17).

* `multiqc`: Aggregate results from bioinformatics analyses across many samples into a single report (PR #42).

* `star`:
    - `star/star_align_reads`: Align reads to a reference genome (PR #22).
    - `star/star_genome_generate`: Generate a genome index for STAR alignment (PR #58).

* `gffread`: Validate, filter, convert and perform other operations on GFF files (PR #29).  

* `salmon`:
    - `salmon/salmon_index`: Create a salmon index for the transcriptome to use Salmon in the mapping-based mode (PR #24).
    - `salmon/salmon_quant`: Transcript quantification from RNA-seq data (PR #24).

* `samtools`:
    - `samtools/samtools_flagstat`: Counts the number of alignments in SAM/BAM/CRAM files for each FLAG type (PR #31).
    - `samtools/samtools_idxstats`: Reports alignment summary statistics for a SAM/BAM/CRAM file (PR #32).
    - `samtools/samtools_index`: Index SAM/BAM/CRAM files (PR #35).
    - `samtools/samtools_sort`: Sort SAM/BAM/CRAM files (PR #36).
    - `samtools/samtools_stats`: Reports alignment summary statistics for a BAM file (PR #39).
    - `samtools/samtools_faidx`: Indexes FASTA files to enable random access to fasta and fastq files (PR #41).
    - `samtools/samtools_collate`: Shuffles and groups reads in SAM/BAM/CRAM files together by their names (PR #42).
    - `samtools/samtools_view`: Views and converts SAM/BAM/CRAM files (PR #48).
    - `samtools/samtools_fastq`: Converts a SAM/BAM/CRAM file to FASTQ (PR #52).
    - `samtools/samtools_fastq`: Converts a SAM/BAM/CRAM file to FASTA (PR #53).

* `umi_tools`:
    -`umi_tools/umi_tools_extract`: Flexible removal of UMI sequences from fastq reads (PR #71).

* `falco`: A C++ drop-in replacement of FastQC to assess the quality of sequence read data (PR #43).

* `bedtools`:
    - `bedtools_getfasta`: extract sequences from a FASTA file for each of the
                           intervals defined in a BED/GFF/VCF file (PR #59).

## MINOR CHANGES

* Uniformize component metadata (PR #23).

* Update to Viash 0.8.5 (PR #25).

* Update to Viash 0.9.0-RC3 (PR #51).

* Update to Viash 0.9.0-RC6 (PR #63).

* Switch to viash-hub/toolbox actions (PR #64).

## DOCUMENTATION

* Update README (PR #64).

## BUG FIXES

* Add escaping character before leading hashtag in the description field of the config file (PR #50).

* Format URL in biobase/bcl_convert description (PR #55).
