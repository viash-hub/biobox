# Contributing Guidelines

We encourage contributions from the community! This guide will help you get started with creating new components for the biobox repository.

**Quick overview:** Fork → Develop → Test → Submit PR

## Quick Start

### Essential Config Template

```yaml
name: your_tool
namespace: category
description: Brief description of what the tool does
keywords: [tag1, tag2]
links:
  homepage: https://tool-homepage.com
  documentation: https://tool-docs.com
  repository: https://github.com/user/repo
references:
  doi: 10.1000/journal.12345
license: MIT/Apache-2.0/GPL-3.0
requirements:
  commands: [your-tool, dependency-tool]
authors:
  - __merge__: /src/_authors/your_name.yaml
    roles: [author, maintainer]
argument_groups:
  - name: Input
    arguments: [...]
  - name: Output
    arguments: [...]
  - name: Arguments
    arguments: [...]
resources:
  - type: bash_script
    path: script.sh
test_resources:
  - type: bash_script
    path: test.sh
  - path: /src/_utils/test_helpers.sh
engines:
  - type: docker
    image: quay.io/biocontainers/tool:version--build_string
    setup:
      - type: docker
        run:
          - tool --version 2>&1 | head -1 | sed 's/.*version /tool: /' > /var/software_versions.txt
runners:
  - type: executable
  - type: nextflow
```

### Essential Commands

```bash
# Create component structure
mkdir -p src/namespace/tool_name
touch src/namespace/tool_name/{script.sh,test.sh,config.vsh.yaml}

# Generate help file
docker run container tool --help > src/namespace/tool_name/help.txt

# Test your component
viash test src/namespace/tool_name/config.vsh.yaml

# Build for testing
viash build src/namespace/tool_name/config.vsh.yaml --setup cachedbuild

# Lint / test the whole namespace
viash ns test --parallel
```

## Conventions

These conventions are **mandatory** for all new components and are actively being back-filled across the repository (see `TODO.md`). When in doubt, look at a recently merged component in the same namespace.

### Required config fields

Every `config.vsh.yaml` MUST declare:

- `name`
- `description`
- `license`
- `authors` (at least one entry, using `__merge__: /src/_authors/<you>.yaml`)
- `argument_groups`
- `resources`
- `test_resources`
- `engines`
- `runners`

Plus `namespace:` when the component is part of a multi-component tool — see below.

If you are a new author, add a YAML file under `src/_authors/<firstname_lastname>.yaml` in the same PR.

### Single-component vs. multi-component tools

Two repo layouts are supported:

- **Single-component tool** (one wrapper per tool, e.g. `fastqc`, `cutadapt`, `multiqc`): place the component at the top level, `src/<tool>/{config.vsh.yaml,script.sh,test.sh,...}`. Do **not** declare a `namespace:` field.
- **Multi-component tool** (tool exposes multiple subcommands, e.g. `samtools sort`, `samtools index`): place each sub-component at `src/<tool>/<tool>_<subcommand>/`, and set `namespace: <tool>` in every sub-component's config.

Rule of thumb: start top-level. Only move to the `namespace:` layout when a second sub-component shows up.

### Argument group names

Use the **singular** form and the standard ordering:

1. `Input`
2. `Output`
3. `Arguments` (tool options that are neither input nor output)

Do not use `Inputs`, `Outputs`, `Options`, or tool-specific group names at the top level.

### Parameter naming

**Mirror the wrapped tool.** Use the same long-option name the underlying tool exposes. If the tool calls it `--reads`, `--bcl_input_directory`, or `-1`, reflect that in the Viash argument name (with leading `--`). This keeps the component's CLI discoverable for users already familiar with the tool, and it keeps the script trivial (flag name matches `par_*` stem).

Guidelines:

- Prefer the tool's **long** option as the Viash `name:`; add short forms under `alternatives:` (e.g. `alternatives: [-1]`).
- If the tool only exposes a short flag, use a descriptive long name and put the short flag in `alternatives:`.
- Substitute `-` with `_` where Viash requires it (Viash maps `--my-flag` ↔ `par_my_flag`).
- Resource-related options (threads, CPUs, memory) are the exception: **do not expose** them as arguments. Use `meta_cpus` and `meta_memory_*` inside the script instead. See the `meta_*` section of [Script Development Guide](docs/SCRIPT_DEVELOPMENT.md).
- Booleans: pick the variant that matches the tool's semantics. They are **not** interchangeable.

#### Boolean variants

| Type             | CLI behaviour                              | Use when                                                                                  |
| ---------------- | ------------------------------------------ | ----------------------------------------------------------------------------------------- |
| `boolean_true`   | `--flag` present ⇒ true, absent ⇒ false    | Default is false; passing the flag turns the option on. Most common case.                 |
| `boolean_false`  | `--flag` present ⇒ false, absent ⇒ true    | Default is true; passing the flag turns the option off (e.g. `--no-cache`).               |
| `boolean`        | User must pass `--flag=true` / `--flag=false` | Tri-state intent: the user explicitly chooses, or the wrapped tool itself takes `true`/`false` as a value (not just a presence toggle). |

Scripting implication: with `boolean_true` / `boolean_false` you typically `unset` the "false" case and use `${par_x:+--flag}`. With bare `boolean` you must forward the actual value (e.g. `--flag="$par_x"`) because both `true` and `false` are meaningful inputs to the tool.

| Intent            | Source of truth           | Example                                  |
| ----------------- | ------------------------- | ---------------------------------------- |
| Input / output    | The wrapped tool's flag   | `--reads`, `--bcl_input_directory`, `-1` |
| Threads / memory  | Viash `meta_*` variables  | `${meta_cpus:+--threads "$meta_cpus"}`   |
| Presence-toggle flag | `boolean_true` / `boolean_false` | `${par_verbose:+--verbose}`        |
| Explicit true/false value | `boolean`        | `--my-opt="$par_my_opt"`                 |

When the underlying flag is cryptic, document it in the argument's `info:` block:

```yaml
- name: --input_1
  alternatives: [-1]
  type: file
  required: true
  info:
    orig_arg: -1
```

### Script conventions

- Language: `bash` (use `type: bash_script`) unless the tool genuinely requires Python.
- Shebang: `#!/bin/bash`.
- Error handling: `set -eo pipefail` on a dedicated line (no `-x`, no trailing whitespace).
- Build arguments with a `cmd_args=( ... )` array; invoke the tool with `"${cmd_args[@]}"`.
- Optional booleans: unset false values first, then use `${par_x:+--flag}` expansion.
- Many booleans (5+): use the `unset_if_false` loop pattern from `docs/SCRIPT_DEVELOPMENT.md`.
- Resource hints: `${meta_cpus:+--threads "$meta_cpus"}`, `${meta_memory_gb:+--memory "${meta_memory_gb}G"}`.
- Multi-value params (`multiple: true`): split with `IFS=';' read -ra arr <<< "$par_x"`.

See [Script Development Guide](docs/SCRIPT_DEVELOPMENT.md) for the full pattern catalogue.

### Test conventions

- Source the shared helpers: `source "$meta_resources_dir/test_helpers.sh"`.
- Generate test data inside the test script whenever possible (keeps the repo lean).
- Assert on file existence **and** content where feasible.
- Each component ships a `test.sh` (or multiple test scripts) under `test_resources`.

See [Testing Guide](docs/TESTING.md).

## Development Guides

### 🔧 [Component Development Guide](docs/COMPONENT_DEVELOPMENT.md)
How to create components: config templates, metadata, arguments, containers, help files, and Docker setup.

### 📝 [Script Development Guide](docs/SCRIPT_DEVELOPMENT.md)
Writing good bash scripts: array-based commands, error handling, conditional parameters, boolean flags, and parameter patterns.

### 🐍 [Python Script Development Guide](docs/PYTHON_SCRIPT_DEVELOPMENT.md)
Writing Python-based components: `par` / `meta` dicts, subprocess handling, when to pick Python over bash.

### ✅ [Testing Guide](docs/TESTING.md)
Testing your components: self-contained tests, generating test data, output validation, and testing multiple scenarios.

### 🐳 [Docker Guide](docs/DOCKER_GUIDE.md)
Working with containers: choosing biocontainers, version pinning, detecting software versions, and container best practices.

## Contribution Process

### Submitting Your Component

1. **Test thoroughly**: Ensure your component passes all tests.

   ```bash
   viash test src/namespace/tool_name/config.vsh.yaml
   # or, for the whole namespace:
   viash ns test --query "^namespace"
   ```

2. **Add changelog entry**: Document your change in `CHANGELOG.md` under the `## Unreleased` section. Use the existing style, for example:

   ```markdown
   ## Unreleased

   ### New functionality

   * `namespace/tool_name`: Add component wrapping `tool` vX.Y.Z (PR #123).

   ### Bug fixes

   * `namespace/tool_name`: Fix handling of empty `--input` (PR #124).
   ```

3. **Self-review your changes** against this checklist:
   - Required config fields present (see [Conventions](#conventions)).
   - Argument group names are singular (`Input`, `Output`, `Arguments`).
   - Argument names mirror the wrapped tool's flags.
   - No `par_threads` / `par_cores` / `par_memory` arguments (use `meta_*`).
   - Boolean arguments use the variant matching the tool's semantics (`boolean_true`, `boolean_false`, or bare `boolean`).
   - Script starts with `#!/bin/bash` + `set -eo pipefail`.
   - `viash test` passes locally.
   - `CHANGELOG.md` updated.

4. **Create a pull request**:
   - Target `main`.
   - Include a clear description of the changes.
   - Link to any relevant issues or discussions.
   - Keep PRs focused; one component per PR is ideal.

### Review Process

- All contributions go through code review.
- Components must pass automated CI tests (`viash ns test`).
- Docker containers must be pinned to an explicit version (no `latest`).
- Documentation must be complete and accurate.
- Reviewers will flag deviations from the [Conventions](#conventions) section.

## Getting Help

### Resources

- **[Viash Documentation](https://viash.io/)**
- **[GitHub Discussions](https://github.com/viash-hub/biobox/discussions)**
- **[Issue Tracker](https://github.com/viash-hub/biobox/issues)**

### Common Questions

**Q: How do I find the right Docker container?**
A: Search for "biocontainer [tool_name]" or check [quay.io/biocontainers](https://quay.io/organization/biocontainers).

**Q: My component fails to build. What should I check?**
A: Verify the Docker image exists, check syntax in `config.vsh.yaml`, and ensure all required commands are available.

**Q: How do I handle tools with complex argument patterns?**
A: Check existing similar components for patterns, or ask in GitHub Discussions.

**Q: Can I create custom Docker containers?**
A: Yes, but biocontainers are preferred when available. See the [Docker Guide](docs/DOCKER_GUIDE.md) for details.

**Q: My tool is written in Python/R, can I still contribute?**
A: Yes. Prefer `bash_script` when the tool is a CLI. For non-trivial logic, use `python_script`. Follow the same `par_*` / `meta_*` conventions; see existing Python components under `src/bd_rhapsody/` and `src/star/` for examples.

---

Happy contributing!
