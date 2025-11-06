# biobox unreleased

## NEW FUNCTIONALITY

* `ensembl_vep`: Added Ensembl Variant Effect Predictor (VEP) support for variant annotation and filtering (PR #205):
  * `ensembl_vep/vep_install`: Install VEP cache data, FASTA files, and plugins for variant effect prediction.
  * `ensembl_vep/filter_vep`: Filter and post-process VEP output by consequence type, phenotype, clinical significance, and more.
  * `ensembl_vep/vep`: Determine the effect of variants on genes, transcripts, and protein sequences with comprehensive annotation options.

# biobox 0.4.1

## MINOR CHANGES

* `falco`: Update falco to 1.2.5 (PR #201).

* `bases2fastq`: Bump from 2.2.0 to 2.2.1 (PR #202).

# biobox 0.4.0

## BREAKING CHANGES

* `fq_subsample` has been removed after its functionality was previously copied to `fq/fq_subsample`. Please use the latter instead. (PR #182).

* `snpeff` has been removed. Please use `snpeff/snpeff_ann` (which is a functional copy of `snpeff`) as this is the default subcommand when running this tool (PR #194)


## NEW FUNCTIONALITY

* `fq`: Added two new components for FASTQ file processing (PR #182):
  - `fq/fq_filter`: Filter FASTQ files based on record names or sequence patterns.
  - `fq/fq_generate`: Generate a random FASTQ file pair for testing and simulation purposes.

* `bwa`: Added BWA support for single-end and paired-end read alignment (PR #183).
  - `bwa/bwa_index`: Create BWA index files for reference genome alignment.
  - `bwa/bwa_mem`: BWA-MEM algorithm for sequence alignment supporting single-end and paired-end reads.
  - `bwa/bwa_aln`: BWA aln algorithm for aligning short sequence reads to a reference genome.
  - `bwa/bwa_samse`: BWA samse - generate single-end alignment in SAM format from BWA aln SAI files.
  - `bwa/bwa_sampe`: BWA sampe - generate paired-end alignment in SAM format from BWA aln SAI files.

* `bowtie2`: Add support for Bowtie2 alignment and indexing (PR #184).
  - `bowtie2/bowtie2_build`: Build Bowtie2 index files from reference sequences.
  - `bowtie2/bowtie2_align`: Align single-end and paired-end reads using Bowtie2.
  - `bowtie2/bowtie2_inspect`: Extract information from Bowtie2 index files.

* `bedtools`: Major expansion with 32 new components providing comprehensive genomic interval analysis (PR #188):
  - `bedtools/bedtools_annotate`: Annotate coverage based on overlaps with interval files
  - `bedtools/bedtools_bedpetobam`: Convert BEDPE to BAM format
  - `bedtools/bedtools_closest`: Find closest features between two interval files
  - `bedtools/bedtools_cluster`: Cluster nearby intervals
  - `bedtools/bedtools_complement`: Report intervals not covered by features
  - `bedtools/bedtools_coverage`: Compute coverage of features
  - `bedtools/bedtools_expand`: Expand blocked BED features
  - `bedtools/bedtools_fisher`: Compute Fisher's exact test for overlaps
  - `bedtools/bedtools_flank`: Create flanking intervals around features
  - `bedtools/bedtools_igv`: Create IGV batch scripts for visualization
  - `bedtools/bedtools_jaccard`: Compute Jaccard statistic between interval sets
  - `bedtools/bedtools_makewindows`: Make windows across genome or intervals
  - `bedtools/bedtools_map`: Map values from overlapping intervals
  - `bedtools/bedtools_maskfasta`: Mask FASTA sequences using intervals
  - `bedtools/bedtools_multicov`: Count coverage across multiple BAM files
  - `bedtools/bedtools_multiinter`: Identify common intervals across multiple files
  - `bedtools/bedtools_overlap`: Compute overlaps between paired-end reads and intervals
  - `bedtools/bedtools_pairtobed`: Find overlaps between paired-end reads and intervals
  - `bedtools/bedtools_pairtopair`: Find overlaps between paired-end read sets
  - `bedtools/bedtools_random`: Generate random intervals
  - `bedtools/bedtools_reldist`: Compute relative distances between features
  - `bedtools/bedtools_sample`: Sample random subsets of intervals
  - `bedtools/bedtools_shift`: Shift intervals by specified amounts
  - `bedtools/bedtools_shuffle`: Shuffle intervals while preserving size
  - `bedtools/bedtools_slop`: Extend intervals by specified amounts
  - `bedtools/bedtools_spacing`: Report spacing between intervals
  - `bedtools/bedtools_split`: Split BED12 features into individual intervals
  - `bedtools/bedtools_subtract`: Remove overlapping features
  - `bedtools/bedtools_summary`: Summarize interval statistics
  - `bedtools/bedtools_tag`: Tag BAM alignments with overlapping intervals
  - `bedtools/bedtools_unionbedg`: Combine multiple BEDGRAPH files
  - `bedtools/bedtools_window`: Find overlapping features within specified windows

* Developer tools: Added GitHub Copilot integration (PR #192):
  - `.github/copilot-instructions.md`: Complete coding assistant guide with biobox patterns, examples, and best practices
  - `.github/prompts/update-viash-component.prompt.md`: Step-by-step prompt for updating existing components
  - `.github/prompts/add-viash-component.prompt.md`: Comprehensive prompt for creating new components from scratch

## MAJOR CHANGES

* `bedtools`: Enhanced 11 existing bedtools components with improved functionality and standardized interfaces (PR #188):
  - `bedtools/bedtools_bamtobed`: Enhanced with additional output format options
  - `bedtools/bedtools_bamtofastq`: Improved paired-end read handling  
  - `bedtools/bedtools_bed12tobed6`: Standardized parameter handling
  - `bedtools/bedtools_bedtobam`: Enhanced genome file support
  - `bedtools/bedtools_genomecov`: Added scale and split options
  - `bedtools/bedtools_getfasta`: Improved FASTA extraction features
  - `bedtools/bedtools_groupby`: Enhanced grouping and operation options
  - `bedtools/bedtools_intersect`: Expanded intersection mode support
  - `bedtools/bedtools_links`: Improved link generation functionality
  - `bedtools/bedtools_merge`: Enhanced merging options and distance parameters
  - `bedtools/bedtools_sort`: Standardized sorting options

* `bcftools`: Updated components to version 1.22 with comprehensive improvements including enhanced argument coverage, improved script patterns, biobox standard compliance, and comprehensive testing overhaul (PR #193):
  * `bcftools_annotate`: Added `--verbosity` parameter; updated to use `meta_cpus` instead of `--threads` parameter
  * `bcftools_concat`: Renamed `--compact_PS` to `--compact_ps`, `--min_PQ` to `--min_pq`; added `--rm_dups`, `--drop_genotypes`, `--verbosity`, `--write_index` parameters; updated to use `meta_cpus` instead of `--threads` parameter
  * `bcftools_norm`: Renamed `--remove_duplicates` to `--rm_dup`, added `--remove_duplicates_flag` as boolean alias; added `--exclude`, `--include`, `--gff_annot`, `--multi_overlaps`, `--sort`, `--verbosity`, `--write_index` parameters; updated to use `meta_cpus` instead of `--threads` parameter
  * `bcftools_sort`: Removed `--max_mem` and `--temp_dir` parameters (now use `meta_memory_mb` and `meta_temp_dir` respectively); added `--verbosity`, `--write_index` parameters
  * `bcftools_stats`: Renamed `--allele_frequency_bins` to `--af_bins`, `--allele_frequency_bins_file` removed, `--allele_frequency_tag` to `--af_tag`, `--fasta_reference` to `--fasta_ref`, `--split_by_ID` to `--split_by_id`, `--targets_overlaps` to `--targets_overlap`

## MINOR CHANGES

* `bases2fastq`: Updated component with comprehensive argument support and latest practices (PR #190).

* `arriba`: Updated to v2.5.0 and refactored script and tests based on latest contributing guidelines (PR #187).

* `snpeff` has been updated to version `5.2f` (PR #194)

# BUG FIXES

* Fix the `commands` property from components being overwritten by the global configuration (which only included `ps`) (PR #196).

## DOCUMENTATION

* Major restructuring of the documentation pages (PR #185):
  - `CONTRIBUTING.md`: Streamlined guide with detailed sections moved to dedicated docs/ guides.
  - `README.md`: Streamlined content to guide people towards what they need.
  - `docs/COMPONENT_DEVELOPMENT.md`: New comprehensive guide covering component creation process.
  - `docs/SCRIPT_DEVELOPMENT.md`: New detailed guide for script development best practices.
  - `docs/TESTING.md`: New comprehensive testing guide.
  - `docs/DOCKER_GUIDE.md`: New Docker and engine best practices guide.

* `.github/PULL_REQUEST_TEMPLATE.md`: Fixed repository references to point to correct biobox repository instead of base template (PR #185).

# biobox 0.3.2

## NEW FUNCTIONALITY

* `fq`:
  - `fq/fq_lint`: Validate FASTQ files for common issues (PR #179).
  - `fq/fq_subsample`: Sample a subset of records from single or paired FASTQ files (PR #179).

## MAJOR CHANGES

* `fq_subsample`: This component has been deprecated in favour of `fq/fq_subsample`, and will be removed in biobox 0.4.0 (PR #179).

## MINOR CHANGES

* Update README (PR #177).

* Update author information (PR #180, PR #200).

* `fastqc`: add `--outdir` argument (PR #181).

# biobox 0.3.1

## NEW FUNCTIONALITY

* `bcl_convert`: add `force` argument (PR #171).
* `cellranger/cellranger_count`: Align fastq files using Cell Ranger count (PR #163).

## MINOR CHANGES

* Replace the deprecated use of the meta variable `functionality_name` by just `name` (PR #174).

* Bump viash to `0.9.4` (PR #175).

## DOCUMENTATION

* Update README (PR #176).

# biobox 0.3.0

## NEW FUNCTIONALITY

* `agat`:
  - `agat/agat_convert_genscan2gff`: convert a genscan file into a GFF file (PR #100).
  - `agat/agat_sp_add_introns`: add intron features to gtf/gff file without intron features (PR #104).
  - `agat/agat_sp_filter_feature_from_kill_list`: remove features in a GFF file based on a kill list (PR #105).
  - `agat/agat_sp_merge_annotations`: merge different gff annotation files in one (PR #106).
  - `agat/agat_sp_statistics`: provides exhaustive statistics of a gft/gff file (PR #107).
  - `agat/agat_sq_stat_basic`: provide basic statistics of a gtf/gff file (PR #110).

* `bd_rhapsody/bd_rhapsody_sequence_analysis`: BD Rhapsody Sequence Analysis CWL pipeline (PR #96).

* `bedtools`:
   - `bedtools/bedtools_bamtobed`: Converts BAM alignments to BED6 or BEDPE format (PR #109).

* `rsem/rsem_calculate_expression`: Calculate expression levels (PR #93).

* `cellranger`:
  - `cellranger/cellranger_mkref`: Build a Cell Ranger-compatible reference folder from user-supplied genome FASTA and gene GTF files (PR #164).

* `rseqc`:
  - `rseqc/rseqc_inner_distance`: Calculate inner distance between read pairs (PR #159).
  - `rseqc/rseqc_inferexperiment`: Infer strandedness from sequencing reads (PR #158).
  - `rseqc/bam_stat`: Generate statistics from a bam file (PR #155).

* `nanoplot`: Plotting tool for long read sequencing data and alignments (PR #95).

* `sgedemux`: demultiplexing sequencing data generated on Singular Genomics' sequencing instruments (PR #166).

* `bases2fasta`: demultiplexing sequencing data generated by Element Biosciences instruments (PR #167).

## BUG FIXES

* `falco`: Fix a typo in the `--reverse_complement` argument (PR #157).

* `cutadapt`: Fix the the non-functional `action` parameter (PR #161).

* `bbmap_bbsplit`: Change argument type of `build` to `file` and add output argument `index` (PR #162).

* `kallisto/kallisto_index`: Fix command script to use `--threads` option (PR #162).

* `kallisto/kallisto_quant`: Change type of argument `output_dir` to `file` and add output argument `log` (PR #162).

* `rsem/rsem_calculate_expression`: Fix output handling (PR #162).

* `sortmerna`: Change type pf argument `aligned` to `file`; update docker image; accept more than two reference files (PR #162).

* `umi_tools/umi_tools_extract`: Remove `umi_discard_reads` option and change `log2stderr` to input argument (PR #162).

* `star/star_genome_generate`: Fix passing of optional sjdb parameters (PR #170).

## MINOR CHANGES

* `agat_convert_bed2gff`: change type of argument `inflate_off` from `boolean_false` to `boolean_true` (PR #160).

* `cutadapt`: change type of argument `no_indels` and `no_match_adapter_wildcards` from `boolean_false` to `boolean_true` (PR #160).

* Upgrade to Viash 0.9.0.

* `bbmap_bbsplit`: Move to namespace `bbmap` (PR #162).

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
  - `agat/agat_convert_mfannot2gff`: convert MFannot "masterfile" annotation to gff format (PR #112).
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
