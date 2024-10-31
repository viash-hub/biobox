# biobox x.x.x

## NEW FUNCTIONALITY

* `agat`:
  - `agat/agat_convert_genscan2gff`: convert a genscan file into a GFF file (PR #100).
  - `agat/agat_sp_add_introns`: add intron features to gtf/gff file without intron features (PR #104).
  - `agat/agat_sp_filter_feature_from_kill_list`: remove features in a GFF file based on a kill list (PR #105).
  - `agat/agat_sp_merge_annotations`: merge different gff annotation files in one (PR #106).
  - `agat/agat_sp_statistics`: provides exhaustive statistics of a gft/gff file (PR #107).

* `bd_rhapsody/bd_rhapsody_sequence_analysis`: BD Rhapsody Sequence Analysis CWL pipeline (PR #96).

* `bedtools`:
   - `bedtools/bedtools_bamtobed`: Converts BAM alignments to BED6 or BEDPE format (PR #109).

* `rsem/rsem_calculate_expression`: Calculate expression levels (PR #93).

* `cellranger`:
  - `cellranger/cellranger_count`: Align fastq files using Cell Ranger count (PR #163).

## BREAKING CHANGES

* `rseqc`:
  - `rseqc/rseqc_inner_distance`: Calculate inner distance between read pairs (PR #159).
  - `rseqc/rseqc_inferexperiment`: Infer strandedness from sequencing reads (PR #158).
  - `rseqc/bam_stat`: Generate statistics from a bam file (PR #155).

* `nanoplot`: Plotting tool for long read sequencing data and alignments (PR #95).

## BUG FIXES

* `falco`: Fix a typo in the `--reverse_complement` argument (PR #157).

* `cutadapt`: Fix the the non-functional `action` parameter (PR #161).

## MINOR CHANGES

* `agat_convert_bed2gff`: change type of argument `inflate_off` from `boolean_false` to `boolean_true` (PR #160).

* `cutadapt`: change type of argument `no_indels` and `no_match_adapter_wildcards` from `boolean_false` to `boolean_true` (PR #160).

* Upgrade to Viash 0.9.0.

# biobox 0.2.0

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
  - `agat_convert_bed2gff`: convert bed file to gff format (PR #97).
  - `agat_convert_embl2gff`: convert an EMBL file into GFF format (PR #99).
  - `agat/agat_convert_sp_gff2gtf`: convert any GTF/GFF file into a proper GTF file (PR #76).
  - `agat/agat_convert_bed2gff`: convert bed file to gff format (PR #97).
  - `agat/agat_convert_embl2gff`: convert an EMBL file into GFF format (PR #99).
  - `agat/agat_convert_sp_gff2tsv`: convert gtf/gff file into tabulated file (PR #102).
  - `agat/agat_convert_sp_gxf2gxf`: fixes and/or standardizes any GTF/GFF file into full sorted GTF/GFF file (PR #103).


* `bedtools`:
  - `bedtools/bedtools_intersect`: Allows one to screen for overlaps between two sets of genomic features (PR #94).
  - `bedtools/bedtools_sort`: Sorts a feature file (bed/gff/vcf) by chromosome and other criteria (PR #98).
  - `bedtools/bedtools_genomecov`: Compute the coverage of a feature file (bed/gff/vcf/bam) among a genome (PR #128).
  - `bedtools/bedtools_groupby`: Summarizes a dataset column based upon common column groupings. Akin to the SQL "group by" command (PR #123).
  - `bedtools/bedtools_merge`: Merges overlapping BED/GFF/VCF entries into a single interval (PR #118).
  - `bedtools/bedtools_bamtofastq`: Convert BAM alignments to FASTQ files (PR #101).
  - `bedtools/bedtools_bedtobam`: Converts genomic feature records (bed/gff/vcf) to BAM format (PR #111).
  - `bedtools/bedtools_bed12tobed6`: Converts BED12 files to BED6 files (PR #140).
  - `bedtools/bedtools_links`: Creates an HTML file with links to an instance of the UCSC Genome Browser for all features / intervals in a (bed/gff/vcf) file (PR #137).
 
* `qualimap/qualimap_rnaseq`: RNA-seq QC analysis using qualimap (PR #74). 

* `rsem/rsem_prepare_reference`: Prepare transcript references for RSEM (PR #89).

* `bcftools`:
  - `bcftools/bcftools_concat`: Concatenate or combine VCF/BCF files (PR #145).
  - `bcftools/bcftools_norm`: Left-align and normalize indels, check if REF alleles match the reference, split multiallelic sites into multiple rows; recover multiallelics from multiple rows (PR #144).
  - `bcftools/bcftools_annotate`: Add or remove annotations from a VCF/BCF file (PR #143).
  - `bcftools/bcftools_stats`: Parses VCF or BCF and produces a txt stats file which can be plotted using plot-vcfstats (PR #142).
  - `bcftools/bcftools_sort`: Sorts BCF/VCF files by position and other criteria (PR #141).

* `fastqc`: High throughput sequence quality control analysis tool (PR #92).

* `sortmerna`: Local sequence alignment tool for mapping, clustering, and filtering rRNA from
  metatranscriptomic data (PR #146).

* `fq_subsample`: Sample a subset of records from single or paired FASTQ files (PR #147).

* `kallisto`:
    - `kallisto_index`: Create a kallisto index (PR #149).
    - `kallisto_quant`: Quantifying abundances of transcripts from RNA-Seq data, or more generally of target sequences using high-throughput sequencing reads (PR #152).

* `trimgalore`: Quality and adapter trimming for fastq files (PR #117). 


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
    - `umi_tools/umi_tools_extract`: Flexible removal of UMI sequences from fastq reads (PR #71).
    - `umi_tools/umi_tools_prepareforrsem`: Fix paired-end reads in name sorted BAM file to prepare for RSEM (PR #148).

* `falco`: A C++ drop-in replacement of FastQC to assess the quality of sequence read data (PR #43).

* `bedtools`:
    - `bedtools_getfasta`: extract sequences from a FASTA file for each of the
                           intervals defined in a BED/GFF/VCF file (PR #59).
                           
* `bbmap`:
    - `bbmap_bbsplit`: Split sequencing reads by mapping them to multiple references simultaneously (PR #138).


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
