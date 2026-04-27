# Python Script Development Guide

This guide covers best practices for writing Python-based component scripts in biobox. Most biobox components are bash wrappers (see [Script Development Guide](SCRIPT_DEVELOPMENT.md)); Python is used when the wrapper needs real logic and/or complex data structures — parsing structured output, conditionals beyond a few branches, or data transformations that would be painful in bash.

## Table of Contents
- [When to Choose Python over Bash](#when-to-choose-python-over-bash)
- [Script Structure and Template](#script-structure-and-template)
- [The `par` and `meta` Dicts](#the-par-and-meta-dicts)
- [Argument Handling](#argument-handling)
- [Calling External Tools](#calling-external-tools)
- [Resource Usage](#resource-usage)
- [Error Handling and Logging](#error-handling-and-logging)
- [File and Path Handling](#file-and-path-handling)
- [Dependencies and Containers](#dependencies-and-containers)
- [Testing](#testing)
- [Real-World Example](#real-world-example)
- [Common Pitfalls](#common-pitfalls)

## When to Choose Python over Bash

Default to bash. Reach for Python when:

- The wrapper parses structured data (JSON, YAML, TSV with complex rules) and feeds it to the tool.
- Argument construction has non-trivial branching (more than a handful of if/else).
- The wrapper has to coordinate multiple tool invocations and shuffle their outputs.
- Input validation is substantive (e.g. matched paired-end file lists, schema checks).

If the job is "arrange CLI flags and execute one tool", stay in bash. It is shorter and the surrounding ecosystem in this repo assumes bash.

In the config, use a Python resource:

```yaml
resources:
  - type: python_script
    path: script.py
```

## Script Structure and Template

```python
import subprocess
from pathlib import Path

## VIASH START
par = {
    "input": "test.fastq",
    "output": "out.bam",
    "verbose": False,
}
meta = {
    "cpus": 4,
    "memory_gb": 8,
    "temp_dir": "/tmp",
    "resources_dir": ".",
}
## VIASH END

# ... your logic ...

cmd = [
    "your_tool",
    "--input", par["input"],
    "--output", par["output"],
]
if par.get("verbose"):
    cmd.append("--verbose")
if meta.get("cpus"):
    cmd += ["--threads", str(meta["cpus"])]

subprocess.run(cmd, check=True)
```

### The `## VIASH START` / `## VIASH END` Block

Viash replaces everything between these markers at build time. Put dev-time stubs in that block so you can run `python script.py` directly during iteration:

- `par`: dict with one entry per argument declared in `config.vsh.yaml`.
- `meta`: dict with Viash runtime metadata.

When the component runs under Viash, `par` and `meta` are populated from the actual invocation — your stubs are discarded.

## The `par` and `meta` Dicts

In Python components, arguments are delivered as a `par` dict rather than environment variables:

- `par["<arg_name>"]` — the value of an argument (without the leading `--`, and with dashes replaced by underscores).
- Missing-or-not-provided optional arguments are `None`, **not** absent from the dict.
- Booleans arrive as Python `bool` values (`True` / `False`), not the strings `"true"` / `"false"`.
- `multiple: true` arguments arrive as a Python `list` of the appropriate type.
- File arguments arrive as `str` paths; wrap in `pathlib.Path` if you need path operations.

`meta` mirrors the bash `meta_*` variables. The common keys:

| Key | Type | Purpose |
| --- | --- | --- |
| `meta["cpus"]` | `int` or `None` | CPU core budget |
| `meta["memory_gb"]`, `meta["memory_mb"]`, ... | `int` or `None` | Memory budgets in various units |
| `meta["temp_dir"]` | `str` | Scratch directory |
| `meta["resources_dir"]` | `str` | Path to component resource files |
| `meta["name"]` | `str` | Component name |
| `meta["config"]` | `str` | Path to the compiled config (useful for introspection) |

## Argument Handling

```python
# Required argument — Viash guarantees this is not None
cmd = ["tool", "--input", par["input"]]

# Optional argument
if par.get("reference") is not None:
    cmd += ["--reference", par["reference"]]

# Boolean flag (presence toggle — boolean_true / boolean_false in the config)
if par.get("verbose"):
    cmd.append("--verbose")

# Bare boolean (value forwarded — user passed --verbose=true/false)
if par.get("verbose") is not None:
    cmd += ["--verbose", str(par["verbose"]).lower()]

# multiple: true — already a list
for fq in par["reads"]:
    cmd += ["--reads", fq]
```

Note the boolean split mirrors the three Viash types (`boolean_true`, `boolean_false`, bare `boolean`) described in [CONTRIBUTING.md](../CONTRIBUTING.md#boolean-variants). With presence-toggle types, a `True` value means "emit the flag"; with bare `boolean`, forward the value.

## Calling External Tools

Always `subprocess.run(cmd, check=True)` with `cmd` as a list:

```python
subprocess.run(
    ["samtools", "view", "-b", par["input"], "-o", par["output"]],
    check=True,
)
```

- Never `shell=True` unless you have a very good reason (piping through shell operators, and even then prefer composing in Python).
- `check=True` turns non-zero exits into `CalledProcessError`, so failures stop the script.
- Capture output with `capture_output=True, text=True` when you need to parse it; otherwise let it stream to the component's stderr/stdout.

For commands that need the shell (pipelines, redirections), prefer Python-level plumbing:

```python
# Instead of: subprocess.run("tool | grep -v '^#' > out.tsv", shell=True, check=True)
with open(par["output"], "w") as out:
    p1 = subprocess.Popen(["tool"], stdout=subprocess.PIPE)
    subprocess.run(["grep", "-v", "^#"], stdin=p1.stdout, stdout=out, check=True)
    p1.stdout.close()
    p1.wait()
```

## Resource Usage

```python
cmd = ["tool", "--input", par["input"], "--output", par["output"]]

if meta.get("cpus"):
    cmd += ["--threads", str(meta["cpus"])]
if meta.get("memory_gb"):
    cmd += ["--memory", f"{meta['memory_gb']}G"]
```

Never expose `par_threads` / `par_cores` / `par_memory` in the config. Resource budgets come from Viash (`meta_*`), so they stay consistent with Nextflow / other runners.

## Error Handling and Logging

- Let exceptions propagate. Viash treats a non-zero exit as a failure; a stack trace is the best diagnostic you can leave.
- Do not wrap the whole script in `try/except Exception: pass` or similar. Silent failures are a CLAUDE.md no-go.
- Use `print(..., flush=True)` or the `logging` module for progress messages. `flush=True` matters for streamed logs in CI.
- Validate inputs at boundaries and raise with a clear message:

```python
if par["input_r2"] and len(par["input"]) != len(par["input_r2"]):
    raise ValueError(
        f"Paired-end R1/R2 count mismatch: {len(par['input'])} vs {len(par['input_r2'])}"
    )
```

## File and Path Handling

- Prefer `pathlib.Path` over manual string manipulation.
- Use `meta["temp_dir"]` for scratch work; wrap with `tempfile.TemporaryDirectory(dir=meta["temp_dir"])` when you want automatic cleanup.
- Create output parent dirs explicitly when the tool does not:

```python
Path(par["output"]).parent.mkdir(parents=True, exist_ok=True)
```

## Dependencies and Containers

- Keep imports minimal; stick to the Python standard library when possible.
- Any third-party package must be installed in the Docker image via the `engines:` setup block:

```yaml
engines:
  - type: docker
    image: python:3.12-slim
    setup:
      - type: python
        packages:
          - pyyaml
          - pandas
```

- Do not add `requirements.txt` or `pyproject.toml` at the component level — the Viash setup block is the source of truth.

## Testing

Component tests still drive the component via a `test.sh` (not `test.py`), because the test runs the compiled executable, not the raw Python. Inside that test script the component looks like any other:

```bash
"$meta_executable" \
  --input sample.fastq \
  --output out.bam
```

If your Python has pure-Python helper logic worth unit-testing in isolation, factor it into a separate module under the component directory and invoke it from both `script.py` and a small pytest file. Keep it rare; the component's integration test is the primary safety net.

See the [Testing Guide](TESTING.md) for the general test contract.

## Real-World Example

`src/star/star_align_reads/script.py` is the reference implementation. It:

1. Reads the compiled config via `meta["config"]` to introspect argument metadata.
2. Assembles STAR's unusual comma-separated `readFilesIn` argument from Viash's list-valued `par["input"]`.
3. Detects gzipped / bzipped inputs and sets `readFilesCommand` accordingly.
4. Rewrites output-path arguments so STAR writes to a scratch directory, then moves the produced files to the user-specified locations.

It is longer than a typical bash wrapper because it needs to do real work — exactly the situation where Python earns its place.

## Common Pitfalls

### 1. Forgetting `check=True`

```python
# Wrong — failures slip through silently
subprocess.run(cmd)

# Correct
subprocess.run(cmd, check=True)
```

### 2. Passing `cmd` as a string without `shell=True`

```python
# Wrong — raises FileNotFoundError
subprocess.run("tool --input foo")

# Correct — list form
subprocess.run(["tool", "--input", "foo"], check=True)
```

### 3. Treating `par["flag"]` as a string

```python
# Wrong — par["flag"] is already a bool
if par["flag"] == "true":
    ...

# Correct
if par["flag"]:
    ...
```

### 4. Hardcoding thread/memory values

```python
# Wrong
cmd += ["--threads", "8"]

# Correct
if meta.get("cpus"):
    cmd += ["--threads", str(meta["cpus"])]
```

### 5. Missing `flush=True` on progress logs

Streamed logs in CI can look frozen if Python buffers stdout. Use `print(..., flush=True)` or `logging.basicConfig(...)` with a stream handler.
