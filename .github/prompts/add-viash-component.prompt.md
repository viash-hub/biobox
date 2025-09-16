---
mode: agent
description: Add a new Viash component following biobox conventions with comprehensive scaffolding
model: Claude Sonnet 4
---

# Add New Viash Component in biobox

## Context & Inputs

- **Component path**: `${input:componentPath}` (e.g., src/namespace/tool_name)
- **Tool name**: `${input:toolName}` (e.g., bcftools annotate)
- **Docker image**: `${input:dockerImage}` (optional, e.g., quay.io/biocontainers/tool:version--build)
- **Help source**: `${input:helpSource}` (providedText | providedFile | container | docs | none)
- **Help content**: `${input:helpText}` (optional, raw --help output)
- **Help file**: `${input:helpFile}` (optional, path to existing help.txt)

## Reference Materials

- **Project patterns**: `.github/copilot-instructions.md`
- **Documentation**: `docs/COMPONENT_DEVELOPMENT.md`, `docs/SCRIPT_DEVELOPMENT.md`, `docs/TESTING.md`, `docs/DOCKER_GUIDE.md`
- **Example component**: `src/bedtools/bedtools_annotate/`
- **Test utilities**: `src/_utils/test_helpers.sh`

## Step-by-Step Creation Process

### Step 1: Create Component Structure

**Objective**: Set up the required 4-file component structure

**Actions**:

1. Create directory: `${input:componentPath}`
2. Create files: `config.vsh.yaml`, `script.sh`, `test.sh`, `help.txt`
3. Ensure proper directory permissions

**Validation**: All required files exist in correct location

### Step 2: Generate Help Documentation

**Objective**: Populate `help.txt` with tool help information

**Actions**:

- If `${input:helpText}` provided: Write directly to `help.txt`
- If `${input:helpFile}` provided: Copy content to `help.txt`
- If `${input:dockerImage}` provided: Run `docker run <image> <tool> --help > help.txt`
- If none provided: Create placeholder with TODO markers

**Validation**: `help.txt` contains useful tool documentation

### Step 3: Draft Configuration (`config.vsh.yaml`)

**Objective**: Create complete component configuration

**Template Structure**:

```yaml
name: {component_name}           # from folder name (snake_case)
namespace: {namespace}           # from parent folder
description: |
  {tool_description_with_links_and_examples}
keywords: [bioinformatics, {domain}, {tool}, {additional_keywords}] # additional keywords could be input/output file formats, methods, or other relevant terms
links:
  - repository: {tool_repo_url}
  - documentation: {docs_url}
authors:
  - __merge__: /src/_authors/{author_name}.yaml
argument_groups:
  - name: "Inputs"
    arguments:
      - name: --input
        type: file
        description: Input file description with example
        required: true
  - name: "Outputs"
    arguments:
      - name: --output
        type: file
        direction: output
        description: Output file description
        required: true
  - name: "Options"
    arguments:
      # Tool-specific options here
requirements:                            # add if tool needs specific CLI commands
  commands: [{cli_tool}]                 # list of required command-line tools
resources:
  - type: bash_script
    path: script.sh
test_resources:
  - type: bash_script
    path: test.sh
  - path: /src/_utils/test_helpers.sh
engines:
  - type: docker
    image: {pinned_biocontainer_image}
    setup:
      - type: docker
        run: |
          {tool} --version | head -1 | sed 's/{tool} /{tool}: /' > /var/software_versions.txt
runners:
  - type: executable
  - type: nextflow
```



**Actions**:

1. Parse component and namespace from path
2. Research tool for description, links, keywords
3. If `help.txt` available: Extract arguments and map to Viash structure
4. Use `--snake_case` naming throughout
5. Set appropriate argument types and requirements
6. Include comprehensive markdown descriptions with examples
7. Exclude `--help/-h` and `--version` arguments
8. Never add CPU/memory arguments (use meta variables)

**Validation**: Configuration is complete and follows biobox standards

### Step 4: Implement Script (`script.sh`)

**Objective**: Create working bash implementation

**Template Structure**:

```bash
#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# Unset false boolean parameters
[[ "$par_flag" == "false" ]] && unset par_flag

# Handle multiple values if needed
if [[ -n "$par_multiple_files" ]]; then
  IFS=';' read -ra files_array <<< "$par_multiple_files"
fi

# Build command array
cmd_args=(
  {tool_command}
  --input "$par_input"
  --output "$par_output"
  ${par_flag:+--flag}
  ${par_option:+--option "$par_option"}
  ${meta_cpus:+--threads "$meta_cpus"}
  ${meta_memory_gb:+--memory "${meta_memory_gb}G"}
)

# Add multiple files if present
if [[ -n "$par_multiple_files" ]]; then
  for file in "${files_array[@]}"; do
    cmd_args+=(--file "$file")
  done
fi

# Create output directory if needed
mkdir -p "$(dirname "$par_output")"

# Execute command
"${cmd_args[@]}"
```

**Actions**:

1. Start with proper Viash blocks and error handling
2. Implement boolean parameter unsetting pattern
3. Handle multiple values with array splitting
4. Build command as array with conditional inclusion
5. Use meta variables for resources (CPU, memory)
6. Quote all file paths properly
7. Create output directories as needed
8. Apply two-space indentation consistently

**Validation**: Script follows all biobox patterns and handles parameters correctly

### Step 5: Write Comprehensive Tests (`test.sh`)

**Objective**: Create self-contained, comprehensive test suite

**Template Structure**:

```bash
#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

# Source centralized test helpers
source "$meta_resources_dir/test_helpers.sh"

# Initialize test environment
setup_test_env

log "Starting tests for $meta_name"

# Generate test data (choose appropriate helpers)
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

# Test 2: Optional parameters
log "Starting TEST 2: Optional parameters"
"$meta_executable" \
  --input "$meta_temp_dir/input.fasta" \
  --output "$meta_temp_dir/output2.txt" \
  --flag \
  --option "value"

check_file_exists "$meta_temp_dir/output2.txt" "output with options"
log "✅ TEST 2 completed successfully"

# Test 3: Multiple inputs (if applicable)
if [[ -n "$par_multiple_support" ]]; then
  log "Starting TEST 3: Multiple inputs"
  "$meta_executable" \
    --multiple_files "$meta_temp_dir/file1.txt;$meta_temp_dir/file2.txt" \
    --output "$meta_temp_dir/output3.txt"
  
  check_file_exists "$meta_temp_dir/output3.txt" "multiple inputs output"
  log "✅ TEST 3 completed successfully"
fi

print_test_summary "All tests completed successfully"
```

**Actions**:

1. Source test helpers and initialize environment
2. Generate appropriate test data using helper functions
3. Create 2-3 test scenarios covering:
   - Basic required functionality
   - Optional parameters/flags
   - Multiple inputs if supported
   - Edge cases if relevant
4. Use structured logging for each test
5. Validate outputs with helper functions
6. End with summary

**Validation**: Tests are self-contained, comprehensive, and executable

### Step 6: Build & Validate Component

**Objective**: Verify the complete component works

**Actions**:

1. Build: `viash build ${input:componentPath}/config.vsh.yaml --setup cachedbuild`
2. Test: `viash test ${input:componentPath}/config.vsh.yaml --keep true --verbose`
3. Check Docker image builds correctly
4. Verify version detection works
5. Run basic functionality test manually

**Validation**: Component builds and tests pass successfully

### Step 7: Final Quality Review

**Objective**: Ensure component meets all biobox standards

**Quality Checklist**:

- [ ] All 4 files present and properly structured
- [ ] Configuration uses correct naming conventions
- [ ] Script follows biobox patterns (arrays, conditionals, meta vars)
- [ ] Tests are self-contained and comprehensive
- [ ] Docker image is pinned with version detection
- [ ] No CPU/memory arguments in config
- [ ] Markdown descriptions are clear with examples
- [ ] Component passes `viash test` successfully
- [ ] Follows two-space indentation
- [ ] Uses proper error handling

## Expected Output

**Created Files**:

- `${input:componentPath}/config.vsh.yaml`: Complete component configuration
- `${input:componentPath}/script.sh`: Working bash implementation
- `${input:componentPath}/test.sh`: Comprehensive test suite
- `${input:componentPath}/help.txt`: Tool help documentation

**Summary Report**:

- Brief description of component functionality
- Notes on any special features or requirements
- Test execution results
- Instructions for usage and further development

## Success Criteria

- Complete 4-file component structure
- Component builds without errors
- All tests pass reliably
- Follows all biobox coding standards
- Ready for integration and use
- Well-documented and maintainable

## Notes

- When uncertain about patterns, refer to `src/bedtools/bedtools_annotate/`
- Prefer pinned biocontainer images when available
- Keep components focused and single-purpose
- Ensure tests cover realistic usage scenarios
