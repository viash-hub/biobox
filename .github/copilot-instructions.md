## biobox – AI Coding Assistant Guide

**Purpose**: Enable AI assistants to efficiently contribute high-quality Viash components that follow biobox standards and pass all quality gates.

## Repository Overview

**biobox** is a curated collection of containerized bioinformatics tools built with the Viash framework. Each component:
- Lives under `src/<namespace>/<component>/` 
- Ships as both a containerized executable and Nextflow module
- Uses pinned biocontainer images for reproducibility
- Follows standardized patterns for arguments, testing, and error handling
- Must pass `viash test` validation in CI

**Key Files**:
- `_viash.yaml`: Top-level Viash project configuration
- `src/_utils/test_helpers.sh`: Centralized testing utilities
- `docs/`: Comprehensive development guides

## Component Structure (Reference: `src/bedtools/bedtools_annotate/`)

Each component consists of exactly 4 files:

### 1. `config.vsh.yaml` - Component Metadata & Interface
```yaml
# Core elements:
name: tool_name                    # snake_case, matches folder name
namespace: tool_namespace          # parent folder name
description: |                    # markdown with examples
  Brief tool description with [documentation](https://example.com) and usage examples
authors:
  - __merge__: /src/_authors/name.yaml  # reference author files
argument_groups:                  # organize by: inputs, outputs, options
  - name: inputs
    arguments:
      - name: --input_file        # snake_case, type-specific patterns
        type: file               # file, string, integer, boolean, double
        description: |           # markdown with examples and defaults
          Input description with example: `input.vcf`
        required: true
        multiple: false          # set true for arrays/lists
resources:
  - type: bash_script
    path: script.sh
test_resources:
  - type: bash_script 
    path: test.sh
  - path: /src/_utils/test_helpers.sh  # always include
engines:
  - type: docker
    image: quay.io/biocontainers/tool:version--build  # pinned versions
    setup:
      - type: docker
        run: |
          tool --version | head -1 | sed 's/tool /tool: /' > /var/software_versions.txt
runners:
  - type: executable
  - type: nextflow
```

### 2. `script.sh` - Bash Implementation
```bash
#!/bin/bash

## VIASH START
# Viash auto-injects par_* and meta_* variables here
## VIASH END

set -eo pipefail  # strict error handling

# Unset false boolean parameters (biobox standard)
[[ "$par_flag" == "false" ]] && unset par_flag
[[ "$par_option" == "false" ]] && unset par_option

# Handle multiple values (semicolon-separated from Viash)
if [[ -n "$par_files" ]]; then
  IFS=';' read -ra files_array <<< "$par_files"
fi

# Build command array (preferred pattern)
cmd_args=(
  tool_command
  --input "$par_input"
  --output "$par_output"
  ${par_flag:+--flag}                    # conditional inclusion
  ${par_option:+--option "$par_option"}  # conditional with value
  ${meta_cpus:+--threads "$meta_cpus"}   # use meta variables for resources
  ${meta_memory_gb:+--memory "${meta_memory_gb}G"}
)

# Handle multiple files in array
if [[ -n "$par_files" ]]; then
  for file in "${files_array[@]}"; do
    cmd_args+=(--file "$file")
  done
fi

# Execute command
"${cmd_args[@]}"
```

### 3. `test.sh` - Self-contained Tests
```bash
#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

# Always source centralized helpers
source "$meta_resources_dir/test_helpers.sh"

# Initialize test environment
setup_test_env

log "Starting tests for $meta_name"

# Create test data using helpers (never external downloads)
create_test_fasta "$meta_temp_dir/input.fasta" 5 100
create_test_bed "$meta_temp_dir/regions.bed" 3

# Test 1: Basic functionality
log "Starting TEST 1: Basic functionality"
"$meta_executable" \
  --input "$meta_temp_dir/input.fasta" \
  --output "$meta_temp_dir/output1.txt"

check_file_exists "$meta_temp_dir/output1.txt" "basic output"
check_file_not_empty "$meta_temp_dir/output1.txt" "basic output"
check_file_contains "$meta_temp_dir/output1.txt" "expected_pattern" "output content"
log "✅ TEST 1 completed successfully"

# Test 2: Advanced options
log "Starting TEST 2: Advanced options"
"$meta_executable" \
  --input "$meta_temp_dir/input.fasta" \
  --output "$meta_temp_dir/output2.txt" \
  --flag \
  --option "value"

check_file_exists "$meta_temp_dir/output2.txt" "advanced output"
log "✅ TEST 2 completed successfully"

# Always end with summary
print_test_summary "All tests completed successfully"
```

### 4. `help.txt` - Tool Help Reference
```
Raw output from: tool --help
(For development reference only)
```

## Essential Patterns for Success

### Script patterns that matter
- Place `## VIASH START` and `## VIASH END` at top; Viash injects `par_*` and `meta_*` vars there.
- Use arrays for commands; include options conditionally:
	- Unset false booleans: `[[ "$par_flag" == "false" ]] && unset par_flag`
	- Conditional include: `${par_flag:+--flag}` and `${meta_cpus:+--threads "$meta_cpus"}`
- Multiple values: Viash passes semicolon-separated strings; split to arrays:
	- `IFS=';' read -ra files_array <<< "$par_files"`
	- Support repeated flags in tests; Viash normalizes them.
- Resources: never add threads/memory as normal params. Use meta vars: `meta_cpus`, `meta_memory_gb`, `meta_temp_dir`, `meta_resources_dir`.

### Testing conventions (see `docs/TESTING.md` and `src/_utils/test_helpers.sh`)
- Tests must be self-contained. Generate data with helpers (`create_test_fasta/fastq/...`) and validate with `check_file_*`, `check_file_contains`, etc.
- Always add helpers as a test resource and source them: `source "$meta_resources_dir/test_helpers.sh"`.
- Use `$meta_executable` to run the built component; finish with `print_test_summary`.

### Docker/engine setup (see `docs/DOCKER_GUIDE.md`)
- Prefer biocontainers with pinned versions. Add version detection to `/var/software_versions.txt`:
	- Example (bedtools): `bedtools --version 2>&1 | head -1 | sed 's/.*bedtools v/bedtools: /' > /var/software_versions.txt`
- Avoid comments inside multiline `run: |` blocks; they become Dockerfile RUN lines.

### Authoring arguments and docs (see `docs/COMPONENT_DEVELOPMENT.md` and `docs/SCRIPT_DEVELOPMENT.md`)
- Use `--snake_case` names; write descriptions in markdown. Prefer file outputs over directories when possible.
- Exclude `--help/-h` and `--version`. Don't add CPU/memory flags—use meta vars.

## Core Workflows & Commands

### Local development
- Build: `viash build src/<ns>/<comp>/config.vsh.yaml --setup cachedbuild`
- Test single: `viash test src/<ns>/<comp>/config.vsh.yaml --keep true --verbose`
- Test all/namespace: `viash ns test --parallel` or `viash ns test -q <ns> --parallel`
- Run with resources: `viash run config.vsh.yaml --cpus 4 --memory 8GB -- --input x --output y`

### Reference examples
- `src/bedtools/bedtools_annotate/script.sh`: arrays, boolean unsetting, multi-value splitting, output redirection.
- `src/bedtools/bedtools_annotate/config.vsh.yaml`: argument groups, `multiple: true`, pinned image, version detection, runners.
- Central helpers: `src/_utils/test_helpers.sh`.

### Quality gates checklist
- Version detection in `/var/software_versions.txt`
- No CPU/memory args in config (use meta vars in script)
- Array-based arguments with conditional inclusion
- Boolean parameter unsetting pattern
- Multiple value splitting for arrays
- Two-space indentation in scripts
- Self-contained tests with centralized helpers
- Markdown descriptions with examples

## Prompt Files for Common Tasks

### Available prompts
- Update existing component: `.github/prompts/update-viash-component.prompt.md`
- Add new component: `.github/prompts/add-viash-component.prompt.md`

### How to run prompts
- In Chat, type `/` and select the prompt by name (e.g., `/update-viash-component`).
- Or open the prompt file in the editor and press the play button.
- Or run "Chat: Run Prompt" from the Command Palette and pick the prompt.

## AI Assistant Best Practices

### When working with this repository:
1. **Always read existing patterns first** - Use `src/bedtools/bedtools_annotate/` as the gold standard reference
2. **Follow the step-by-step prompts** - Use the provided prompt files for systematic component work
3. **Validate frequently** - Build and test components after each major change
4. **Use the helper functions** - Leverage `src/_utils/test_helpers.sh` for consistent test patterns
5. **Check quality gates** - Ensure all biobox standards are met before considering work complete

### Common Pitfalls to Avoid:
- Adding CPU/memory arguments to config (use meta variables instead)
- Forgetting to unset false boolean parameters
- Not using arrays for command building
- Including `--help` or `--version` arguments
- Writing tests that depend on external files
- Using inconsistent indentation (must be 2 spaces)
- Missing version detection in Docker setup

What to reference for examples
- `src/bedtools/bedtools_annotate/script.sh`: arrays, boolean unsetting, multi-value splitting, output redirection.
- `src/bedtools/bedtools_annotate/config.vsh.yaml`: argument groups, `multiple: true`, pinned image, version detection, runners.
- Central helpers: `src/_utils/test_helpers.sh`.

Questions to confirm/extend
- Any repo-wide defaults for Nextflow runner options beyond `runners: [executable,nextflow]`?
- Any additional meta variables or naming conventions to enforce (beyond current docs)?

Prompt files for common tasks
- Update existing component: `.github/prompts/update-viash-component.prompt.md`
- Add new component: `.github/prompts/add-viash-component.prompt.md`

How to run a prompt
- In Chat, type `/` and select the prompt by name (for example, `/update-viash-component`).
- Or open the prompt file in the editor and press the play button.
- Or run “Chat: Run Prompt” from the Command Palette and pick the prompt.

