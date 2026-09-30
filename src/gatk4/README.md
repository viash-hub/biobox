# Developing GATK4 components

This guide gives extra context and information for developing `gatk4` components.
It extends on the general guides in [`docs/`](../../docs/), you should read those first.

## Shared files

All GATK4 components use the shared files in [`utils/`](utils/).

| File | Use |
|---|---|
| `common_engines.yaml` | Shared engine definition merged into component configs. Change the container here to update all components. |
| `common_argument_groups.yaml` | Common arguments shared by several components. Merged into the config for components that support those options. |
| `script_helpers.sh` | Staging and argument helpers for `script.sh` (see below). |
| `test_helpers.sh` | Synthetic test file fixtures and other shared helpers for `test.sh`. |

Use the helpers instead of writing the same logic again in a component.
If a new component needs new logic that another component could use, add it to the helpers.

## Arguments

- Name arguments after the GATK flags (with underscores), and give the GATK short forms as `alternatives`
- Arguments should behave the same as the GATK tool.
  Do not copy defaults or behavior from reference workflows (the GATK WDL, `nf-core/sarek`, etc.).
- Put the GATK defaults in the descriptions (see [Description Formatting Guidelines](../../docs/COMPONENT_DEVELOPMENT.md#description-formatting-guidelines)).
  Get them from the generated `help.txt`.
- For tools with a few options, include all of them.
  For tools with many options, decide on a curated list of arguments.
- Examples of options that may be skipped:
  - DRAGEN, DRAGstr and PDHMM settings
  - Flow-based (Ultima) and Permutect options
  - Debug outputs
  - Deprecated options
  - General engine options that are not in `common_argument_groups.yaml`
- Include the settings for tool-default read filters when users need to change them
- GATK lists `--input` (reads) for many tools, but a tool does not always use it.
  Check each tool before you add it.
- An argument that GATK accepts as a string or a file needs two arguments so that Nextflow can stage the file.
  Check to see if these need to be mutually exclusive or if the tool accepts a mixture of strings and files.
- When you generate `choices` from `help.txt` carefully check the list.
  Possible values may be split over lines or only be listed in the description.

## Companion files

For some file types, GATK finds companion files (like indexes) by name next to the main file, but Viash and Nextflow stage each file separately.
Each companion file needs its own argument, and the script links it next to the main file in a temporary directory.

- **Reference**: `--reference`, `--reference_fai` and `--reference_dict`, staged with `stage_reference_trio`.
- **Reads**: `--input` and `--bai`, staged with `stage_bam_bai`, or with `stage_bams_bais` for `multiple: true` (matched by order).
  GATK matches samples (for example Mutect2 `--normal-sample`) by the `SM` tag, not by file name so it is ok for files to be renamed during staging.
- **VCF resources**: every VCF input has an optional `*_index` argument, staged with `stage_vcf_with_index`.
  Without an index, the helper copies the VCF, indexes it and writes a warning.
  This can be slow for large files such as gnomAD.
  - A plain `.vcf` works without an index, because GATK indexes it in memory
  - A `.vcf.gz` (or `.vcf.bgz`) without a `.tbi` fails ("An index is required but was not found")
  - A `.vcf.gz` made with `gzip` (not `bgzip`) cannot be indexed
- **Intervals**: `--intervals` and `--exclude_intervals` can be BED, interval_list or VCF files.
  Stage them with `stage_interval_files`, which matches the `*_index` files to the VCF files by order.

The staging helpers run in `$(...)`, where `set -e` does not apply.
So every step in them must exit explicitly on failure (`|| exit 1`).

## Outputs

**No implicit outputs.**
Every file that GATK writes must have an output argument.
The Nextflow runner keeps only declared outputs, so an undeclared file is kept by the executable runner but lost in Nextflow.

Some components write an additional file next to the main output by default.
These may need to be handled differently.

Examples:

| Side file | Components | Name GATK uses | How the component handles it |
|---|---|---|---|
| Mutect2 stats | `gatk4_mutect2` | `<output>.stats` | Required `--output_stats`, moved after the run |
| FilterMutectCalls filtering stats | `gatk4_filtermutectcalls` | `<output>.filteringStats.tsv` | Required `--filtering_stats`, passed to GATK |
| VCF index | `gatk4_mutect2`, `gatk4_filtermutectcalls` | `<output>.tbi` for `.vcf.gz`/`.vcf.bgz`, `<output>.idx` otherwise (`gatk_output_vcf_index_path`) | `--output_index`, required unless `--create_output_variant_index` is `false`, moved after the run |
| BAM index | `gatk4_mutect2` | `<output without .bam>.bai` (`gatk_output_bam_index_path`) | For example `--bam_output_index`, required unless `--create_output_bam_index` is `false` |

If two output arguments interact, the script should check the provided values are allowed before the run and stop with a clear error.

### Optional outputs

Some outputs are only written by GATK on request.
Sometimes there is an explicit boolean argument to control if an output is created or not but often it is controlled by whether a path has been provided.
The Viash Nextflow runner always gives every output argument a default path (`$id.$key.<name><ext>`) so the presence of a path cannot be used to control output creation in Nextflow.

This is handled by adding a matching `--create_*` boolean argument and setting `must_exist: false` on the output path argument.
The script should only create the output when `--create_*` is `true`.
If `--create_*` is `true` and no output path is provided the script should error.

## Resources

Set the JVM heap from `meta_memory_mb` with a fallback:

```bash
avail_mem_mb=$(( ${meta_memory_mb:-3072} * 8 / 10 ))
gatk --java-options "-Xmx${avail_mem_mb}M -XX:-UsePerfData" <Tool> "${cmd_args[@]}"
```

The `gatk` launcher has no default heap size of its own so a value must be set.
Without `-Xmx`, the JVM uses a part of the visible memory, which is not the memory given to the task.
The fallback value is taken from `nf-core` modules.

Give `--native-pair-hmm-threads` from `meta_cpus` for tools that have it.

## Testing

- The `test_helpers.sh` script has functions for creating synthetic test fixtures
- The GATK image has `samtools`, `bgzip` and `tabix` available for creating additional test files, as well as all GATK tools
- Confirm expected behavior in the real container, not only in the GATK documentation.
  For example, side file names, index requirements and tool-default read filters.
