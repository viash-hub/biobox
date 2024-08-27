name: bcftools_stats
namespace: bcftools
description: | 
  Parses VCF or BCF and produces a txt stats file which can be plotted using plot-vcfstats.
  When two files are given, the program generates separate stats for intersection
  and the complements. By default only sites are compared, -s/-S must given to include
  also sample columns.
keywords: [Stats, VCF, BCF]
links:
  homepage: https://samtools.github.io/bcftools/
  documentation: https://samtools.github.io/bcftools/bcftools.html#stats
  repository: https://github.com/samtools/bcftools
  issue_tracker: https://github.com/samtools/bcftools/issues
references:
  doi: https://doi.org/10.1093/gigascience/giab008
license: MIT/Expat, GNU
requirements:
  commands: [bcftools]
authors:
  - __merge__: /src/_authors/theodoro_gasperin.yaml
    roles: [ author, maintainer ]

argument_groups:
  - name: Inputs
    arguments:
      - name: --input
        alternatives: -i
        type: file
        multiple: true
        description: Input VCF/BCF file. Maximum of two files.
        required: true
    
  - name: Outputs
    arguments:
      - name: --output
        alternatives: -o
        direction: output
        type: file
        description: Output txt statistics file.
        required: true
         
  - name: Options
    arguments:
      
      - name: --allele_frequency_bins
        alternatives: --af_bins
        type: string
        description: | 
          Allele frequency bins, a list (0.1,0.5,1) or a file (0.1\n0.5\n1).

      - name: --allele_frequency_tag
        alternatives: --af_tag
        type: string
        description: | 
          Allele frequency tag to use, by default estimated from AN,AC or GT.

      - name: --first_allele_only
        alternatives: --first_only
        type: boolean_true
        description: | 
          Include only 1st allele at multiallelic sites

      - name: --collapse
        alternatives: --c
        type: string
        choices: [ snps, indels, both, all, some, none ]
        description: | 
          Treat as identical records with <snps|indels|both|all|some|none>.
          See https://samtools.github.io/bcftools/bcftools.html#common_options for details.

      - name: --depth
        alternatives: --d
        type: string
        description: | 
          Depth distribution: min,max,bin size [0,500,1]
        example: 0,500,1

      - name: --exclude
        alternatives: --e
        type: string
        description: | 
          Exclude sites for which the expression is true.
          See https://samtools.github.io/bcftools/bcftools.html#expressions for details.

      - name: --exons
        alternatives: --E
        type: file
        description: | 
          tab-delimited file with exons for indel frameshifts statistics. 
          The columns of the file are CHR, FROM, TO, with 1-based, inclusive, positions. 
          The file is BGZF-compressed and indexed with tabix
          e.g. 
            tabix -s1 -b2 -e3 file.gz

      - name: --apply_filters
        alternatives: --f
        type: string
        description: | 
          Require at least one of the listed FILTER strings (e.g. "PASS,.")

      - name: --fasta_reference
        alternatives: --F
        type: file
        description: | 
          Faidx indexed reference sequence file to determine INDEL context

      - name: --include
        alternatives: --i
        type: string
        description: | 
          Select sites for which the expression is true.
          See https://samtools.github.io/bcftools/bcftools.html#expressions for details.
      
      - name: --split_by_ID
        alternatives: --I
        type: boolean_true
        description: | 
          Collect stats for sites with ID separately (known vs novel)

      - name: --regions
        alternatives: --r
        type: string
        description: | 
          Restrict to comma-separated list of regions

      - name: --regions_file
        alternatives: --R
        type: file
        description: | 
          Restrict to regions listed in a file.

      - name: --regions_overlap
        type: string
        choices: ['pos', 'record', 'variant', '0', '1', '2']
        description: | 
          This option controls how overlapping records are determined: 
          set to 'pos' or '0' if the VCF record has to have POS inside a region (this corresponds to the default behavior of -t/-T); 
          set to 'record' or '1' if also overlapping records with POS outside a region should be included (this is the default behavior of -r/-R, 
          and includes indels with POS at the end of a region, which are technically outside the region); 
          or set to 'variant' or '2' to include only true overlapping variation (compare the full VCF representation "TA>T-" vs the true sequence variation "A>-").

      - name: --samples
        alternatives: --s
        type: string
        description: | 
          List of samples for sample stats, "-" to include all samples.

      - name: --samples_file
        alternatives: --S
        type: file
        description: | 
          File of samples to include.

      - name: --targets
        alternatives: --t
        type: string
        description: | 
          Similar as -r, --regions, but the next position is accessed by streaming the whole VCF/BCF 
          rather than using the tbi/csi index. Both -r and -t options can be applied simultaneously: -r uses the 
          index to jump to a region and -t discards positions which are not in the targets. Unlike -r, targets 
          can be prefixed with "^" to request logical complement. For example, "^X,Y,MT" indicates that 
          sequences X, Y and MT should be skipped. Yet another difference between the -t/-T and -r/-R is 
          that -r/-R checks for proper overlaps and considers both POS and the end position of an indel, 
          while -t/-T considers the POS coordinate only (by default; see also --regions-overlap and --targets-overlap). 
          Note that -t cannot be used in combination with -T.
      
      - name: --targets_file
        alternatives: --T
        type: file
        description: | 
          Similar to -R but streams rather than index-jumps.

      - name: --targets_overlaps
        type: string
        choices: ['pos', 'record', 'variant', '0', '1', '2']
        description: | 
          Include if POS in the region (0), record overlaps (1), variant overlaps (2).

      - name: --user_tstv
        alternatives: --u
        type: string
        description: | 
          Collect Ts/Tv stats for any tag using the given binning [0:1:100].
          A subfield can be selected as e.g. 'PV4[0]', here the first value of the PV4 tag.
      
      - name: --verbose 
        alternatives: --v
        type: boolean_true
        description: | 
          Produce verbose per-site and per-sample output.

resources:
  - type: bash_script
    path: script.sh

test_resources:
  - type: bash_script
    path: test.sh

engines:
  - type: docker
    image: debian:stable-slim
    setup:
      - type: apt
        packages: [bcftools, procps]
      - type: docker
        run: |
          echo "bcftools: \"$(bcftools --version | grep 'bcftools' | sed -n 's/^bcftools //p')\"" > /var/software_versions.txt
    test_setup:  
      - type: apt  
        packages: [tabix]

runners:
  - type: executable
  - type: nextflow


