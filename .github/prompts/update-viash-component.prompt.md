---
mode: agent
description: Update a Viash component following biobox conventions with step-by-step validation
model: Claude Sonnet 4
---

# Update Viash Component in biobox

## Context & Inputs
- **Component path**: `${input:componentPath}` (e.g., src/bcftools/bcftools_annotate)
- **Tool & version**: `${input:toolInfo}` (e.g., bcftools 1.20, quay.io/biocontainers/bcftools:1.20--h8b25389_0)
- **Changes needed**: `${input:changes}` (flags added/removed/modified)

## Reference Materials
- **Project patterns**: `.github/copilot-instructions.md`
- **Documentation**: `docs/COMPONENT_DEVELOPMENT.md`, `docs/SCRIPT_DEVELOPMENT.md`, `docs/TESTING.md`, `docs/DOCKER_GUIDE.md`
- **Example component**: `src/bedtools/bedtools_annotate/`
- **Test utilities**: `src/_utils/test_helpers.sh`

## Step-by-Step Update Process

### Step 1: Refresh Help Documentation
**Objective**: Update `help.txt` with latest tool help output

**Actions**:
1. If biocontainer exists: Run `docker run <image> <tool> --help` and save output to `help.txt`
2. If no container: Capture from authoritative documentation
3. Keep raw output for development reference

**Validation**: `help.txt` contains current command-line options

### Step 2: Update Configuration (`config.vsh.yaml`)
**Objective**: Synchronize arguments with current tool capabilities

**Actions**:
1. Compare `help.txt` against existing `argument_groups`
2. Add missing arguments using `--snake_case` naming
3. Remove deprecated arguments
4. Update argument types (`file`, `string`, `integer`, `boolean`, `double`)
5. Set `multiple: true` for array inputs
6. Write clear markdown descriptions with examples
7. Exclude `--help/-h` and `--version` (handled by Viash)
8. Never add CPU/memory parameters (use `meta_*` variables)

**Validation**: All current tool options represented, no deprecated options remain

### Step 3: Update Docker Engine & Version Detection
**Objective**: Pin specific biocontainer version with version detection

**Actions**:
1. Pin `quay.io/biocontainers/<tool>:<version>--<build>` in engines section
2. Add version detection command to write `tool: <version>` to `/var/software_versions.txt`
3. Avoid comments in multiline `run: |` blocks

**Example**:
```yaml
engines:
  - type: docker
    image: quay.io/biocontainers/bcftools:1.20--h8b25389_0
    setup:
      - type: docker
        run: |
          bcftools --version 2>&1 | head -1 | sed 's/bcftools /bcftools: /' > /var/software_versions.txt
```

**Validation**: Container builds successfully, version detection works

### Step 4: Refactor Script (`script.sh`)
**Objective**: Follow biobox scripting patterns

**Actions**:
1. Maintain `## VIASH START` / `## VIASH END` blocks
2. Use `set -eo pipefail` for error handling
3. Apply two-space indentation consistently
4. Unset false boolean parameters: `[[ "$par_flag" == "false" ]] && unset par_flag`
5. Build command as array with conditional parameter inclusion:
   - `${par_flag:+--flag}`
   - `${par_value:+--option "$par_value"}`
   - `${meta_cpus:+--threads "$meta_cpus"}`
6. Handle multiple values by splitting: `IFS=';' read -ra files_array <<< "$par_files"`
7. Quote file paths and create output directories as needed

**Validation**: Script follows all biobox patterns, handles all parameter types correctly

### Step 5: Rewrite Tests (`test.sh`)
**Objective**: Create comprehensive, self-contained tests

**Actions**:
1. Source test helpers: `source "$meta_resources_dir/test_helpers.sh"`
2. Initialize environment: `setup_test_env`
3. Generate test data using helpers (e.g., `create_test_fasta`, `create_test_bed`)
4. Create multiple test scenarios:
   - Basic functionality with required parameters
   - Advanced features with optional flags
   - Edge cases if applicable
5. Validate outputs with `check_file_*` functions
6. Use structured logging: `log "Starting TEST 1: Description"`
7. End with `print_test_summary "All tests completed successfully"`

**Validation**: Tests are self-contained, comprehensive, and pass reliably

### Step 6: Build & Test Validation
**Objective**: Verify component works correctly

**Actions**:
1. Build component: `viash build ${input:componentPath}/config.vsh.yaml --setup cachedbuild`
2. Run tests: `viash test ${input:componentPath}/config.vsh.yaml --keep true --verbose`
3. Check version detection in container
4. Optional: Run namespace tests: `viash ns test -q <namespace> --parallel`

**Validation**: All tests pass, no build errors, version detection works

### Step 7: Quality Assurance Review
**Objective**: Ensure all biobox standards are met

**Checklist**:
- [ ] `/var/software_versions.txt` contains correct version line
- [ ] No CPU/memory arguments in config (script uses `meta_*` variables)
- [ ] Script uses arrays for command building
- [ ] Boolean parameters are unset when false
- [ ] Multiple values are properly split and handled
- [ ] Two-space indentation throughout
- [ ] Argument descriptions are clear and include examples
- [ ] Tests are self-contained and comprehensive
- [ ] All files follow biobox conventions

## Expected Output

**Updated Files**:
- `help.txt`: Current tool help output
- `config.vsh.yaml`: Updated arguments, pinned container, version detection
- `script.sh`: Refactored following biobox patterns
- `test.sh`: Comprehensive self-contained tests

**Summary Report**:
- Brief rationale for each file change
- Test execution results (PASS/FAIL)
- Any issues encountered and resolutions
- Verification that all quality gates are met

## Success Criteria
- Component builds without errors
- All tests pass consistently
- Follows all biobox coding standards
- Version detection works correctly
- Ready for production use
