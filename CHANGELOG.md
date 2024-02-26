# base unreleased

## BREAKING CHANGES

* Change default `multiple_sep` to `;` (PR #25). This aligns with an upcoming breaking change in
  Viash 0.9.0 in order to avoid issues with the current default separator `:` unintentionally
  splitting up certain file paths.

## NEW FEATURES

* `arriba`: Detect gene fusions from RNA-seq data (PR #1).

* `fastp`: An ultra-fast all-in-one FASTQ preprocessor (PR #3).

* `busco`: 
    - `busco/busco_run`: Assess genome assembly and annotation completeness with single copy orthologs (PR #6).
    - `busco/busco_list_datasets`: Lists available busco datasets (PR #18).
    - `busco/busco_download_datasets`: Download busco datasets (PR #19).

* `featurecounts`: Assign sequence reads to genomic features (PR #11).

* `bgzip`: Add bgzip functionality to compress and decompress files (PR #13).

* `pear`: Paired-end read merger (PR #10).

* `lofreq/call`: Call variants from a BAM file (PR #17).

* `lofreq/indelqual`: Insert indel qualities into BAM file (PR #17).

* `salmon`:
    - `salmon/salmon_index`: Create a salmon index for the transcriptome to use Salmon in the mapping-based mode (PR #24).
    - `salmon/salmon_quant`: Transcript quantification from RNA-seq data (PR #24).

## MAJOR CHANGES

## MINOR CHANGES

* Uniformize component metadata (PR #23).

* Update to Viash 0.8.5 (PR #25).

## DOCUMENTATION

## BUG FIXES