# Biobox upstream fix plan

This file expands the two entries in [upstream-issues.md](upstream-issues.md) into a concrete fix plan.
It is written to be turned directly into issues/PRs against `viash-hub/biobox` (https://github.com/viash-hub/biobox).
Line/path references below were verified against the pinned `v0.4.2` tag for `snpeff_ann`, and against the `main` branch for `bwa_mem2_mem` (not yet tagged/released).

---

## 1. `snpeff/snpeff_ann`: `--stats` swallows the genome-version argument

**Repo state:** released, tagged `v0.4.2`.
**Files:** `src/snpeff/snpeff_ann/config.vsh.yaml`, `src/snpeff/snpeff_ann/script.sh`.

### Root cause

`--stats` (aliases `-s`, `--htmlStats`) is modeled as `type: boolean_true`, but snpEff's real `-stats` flag takes a mandatory filename argument.
`script.sh` emits it as a bare `${par_stats:+-stats}` token with nothing following it.
When `--stats` is set and no other optional flag happens to sit between it and the two trailing positionals, the command collapses to:

```
snpEff ann -stats "$par_genome_version" "$par_input" > "$par_output"
```

snpEff's parser grabs the token right after `-stats` as the output filename, so it consumes `$par_genome_version`, leaving `$par_input` (the VCF) misread as the genome version.
The `--csv_stats` argument right next to it already does this correctly (`type: file`, emitted as `${par_csv_stats:+-csvStats "$par_csv_stats"}`), and is the template to copy.

Note that snpEff writes `snpEff_summary.html` / `snpEff_genes.txt` by default whenever `-noStats` is absent, regardless of `-stats`.
The `mv` logic at the bottom of `script.sh` already only checks `$par_no_stats`, never `$par_stats`, confirming `--stats` is only meant to *rename/relocate* the default summary output, not to trigger it.

### Proposed config change

`config.vsh.yaml`, in the `Options` argument group:

```yaml
# before
- name: --stats
  alternatives: [-s, --htmlStats]
  type: boolean_true
  description: Create HTML summary file.

# after
- name: --stats
  alternatives: [-s, --htmlStats]
  type: file
  direction: output
  description: Create HTML summary file at this path (deprecated upstream in favor of --csv_stats).
```

### Proposed script change

`script.sh`, in the `unset_if_false` array: remove `par_stats` (no longer a boolean, nothing to unset).

`script.sh`, in the `snpEff ann` invocation:

```bash
# before
${par_stats:+-stats} \

# after
${par_stats:+-stats "$par_stats"} \
```

### Other files to update

- `help.txt`: regenerate or hand-edit the `--stats` entry to show it now takes a value.
- `test.sh` / `test_data`: add a case that passes `--stats <path>` and asserts the file is created at that path, to cover the regression directly (today's tests presumably don't exercise this flag, or the bug would have surfaced already).
- `CHANGELOG.md`: one line under a new release, e.g. `snpeff_ann: fix --stats/-htmlStats, was previously modeled as boolean_true and corrupted the command line when used (issue link)`.

### Risk / compatibility

Not a breaking change in practice: any existing caller that actually sets `stats: true` today gets a broken command line, so there is no working behavior to preserve.
Callers that leave `--stats` unset (the common case, matching this project's own workaround) are unaffected either way.

---

## 2. `bwa_mem2/bwa_mem2_mem`: `--index` does not co-stage sibling index files

**Repo state:** unreleased, present on `main` only (`src/bwa_mem2/bwa_mem2_mem`, `src/bwa_mem2/bwa_mem2_index`), not in any tagged version yet.
**Files:** `src/bwa_mem2/bwa_mem2_mem/config.vsh.yaml`, `src/bwa_mem2/bwa_mem2_mem/script.sh`.

### Root cause

A bwa-mem2 index is five files sharing one prefix (`.0123`, `.amb`, `.ann`, `.bwt.2bit.64`, `.pac`), produced as a set by `bwa_mem2_index`'s `--output` (already correctly typed as a directory-shaped `file` output).
`bwa_mem2_mem` declares `--index` as a single required `type: file` and passes it straight to `bwa-mem2 mem`.

Under the Nextflow (VDSL3) runner, each `type: file` argument becomes a Nextflow `path` input staged at the granularity of the exact value supplied: a single file path stages only that file into the task's work directory, a directory path stages the whole directory recursively.
Because `--index` points at one file, only that file is staged; the other four siblings, which live next to it in the original location, are not copied into the isolated task directory.
`bwa-mem2 mem` then fails to open `<staged-path>.amb`, `.ann`, etc. next to the one file that did get staged, since those paths do not exist in that directory.

This is specific to the Nextflow runner's per-task isolation; under the `executable` runner (also declared for this component) sibling files in the same source directory are often still reachable, which is presumably why this was not caught by casual/native testing.

### Proposed config change

`config.vsh.yaml`, `Input` argument group:

```yaml
# before
- name: "--index"
  type: file
  description: BWA index base name (prefix of .0123, .amb, .ann, .bwt.2bit.64, .pac files).
  required: true
  example: reference.fasta

# after
- name: "--index_directory"
  type: file
  description: Directory containing the BWA-MEM2 index files (.0123, .amb, .ann, .bwt.2bit.64, .pac).
  required: true
  example: bwa_index/
- name: "--index_prefix"
  type: string
  description: Filename of the original FASTA within the index directory, used as the index prefix. Auto-detected from the .0123 file if not given.
  example: genome.fasta
```

### Proposed script change

`script.sh`, replace the direct `"$par_index"` usage with a reconstructed base path, e.g.:

```bash
# before
    "$par_index"
    "$par_reads1"

# after
    index_prefix="${par_index_prefix:-$(basename "$(ls "$par_index_directory"/*.0123)" .0123)}"
    "${par_index_directory%/}/$index_prefix"
    "$par_reads1"
```

This mirrors this project's own `bwamem2/bwamem2_mem` component (`src/components/bwamem2/bwamem2_mem`), which already implements exactly this pattern and can serve as a working reference for the PR.

### Other files to update

- `test.sh`: update to build an index directory and pass `--index_directory` / `--index_prefix` instead of a single `--index` file.
- `CHANGELOG.md`: since this is unreleased, a plain entry under `NEW FUNCTIONALITY` (or amending the entry that introduced `bwa_mem2_mem`, if not yet released) is enough; no breaking-change note is needed because no tagged version has shipped the old `--index` shape yet.

### Risk / compatibility

None: the component has never been part of a tagged release, so changing its argument shape carries no breaking-change cost for existing consumers pinned to a released `biobox` version.
Anyone depending directly on `main` would need to update their call site, but that is an accepted cost of tracking an unreleased branch.

---

## Related, out of scope for now

The same architectural problem (single `--index` file/string standing in for a multi-file, shared-prefix index) also exists in already-released, tagged components:

- `bwa/bwa_mem`: `--index` is a single `type: file` for classic bwa's own 5-file index set (`.amb`, `.ann`, `.bwt`, `.pac`, `.sa`). Same fix shape as above, but this is a breaking API change since it is already released, and would need a changelog breaking-change note (there is precedent for this in biobox's own `CHANGELOG.md`, e.g. the v0.4.0 `snpeff` removal note).
- `bowtie2/bowtie2_align`: `--index` is `type: string` (a bare prefix, not even `type: file`), so under Nextflow none of the `.bt2` index files are staged at all. This is more broken than a missing-siblings problem: as declared, the component has no Nextflow-safe way to receive its index.

These are not part of the two tracked issues and are not blocking anything in this project today, but are worth filing as follow-up issues if/when this project or others start depending on those components under the Nextflow runner.
