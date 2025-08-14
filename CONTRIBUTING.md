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
  - name: Inputs
    arguments: [...]
  - name: Outputs  
    arguments: [...]
resources:
  - type: bash_script
    path: script.sh
test_resources:
  - type: bash_script
    path: test.sh
engines:
  - type: docker
    image: quay.io/biocontainers/tool:version--build_string
    setup:
      - type: docker
        run: |
          tool --version 2>&1 | head -1 | sed 's/.*version /tool: /' > /var/software_versions.txt
runners:
  - type: executable
  - type: nextflow
```

### Essential Commands

```bash
# Create component structure
mkdir -p src/namespace/tool_name

# Generate help file
docker run container tool --help > src/namespace/tool_name/help.txt

# Test your component
viash test src/namespace/tool_name/config.vsh.yaml

# Build for testing
viash build src/namespace/tool_name/config.vsh.yaml --setup cachedbuild
```

### Key Best Practices

- **Generate test data** in test scripts (don't store static files)
- **Use array-based arguments**: `cmd_args=(...)` instead of repetitive `cmd_args+=()`
- **Use `__merge__`** for author info: `__merge__: /src/_authors/name.yaml`
- **Specify requirements**: List all commands in `requirements.commands`
- **Pin Docker versions**: Use specific versions with build strings
- **Add version detection**: Create `/var/software_versions.txt`

## Detailed Development Guides

For comprehensive instructions on each aspect of component development, see our detailed guides:

<details>
<summary><strong>Component Development Process</strong></summary>

See: **[Component Development Guide](docs/COMPONENT_DEVELOPMENT.md)**

This guide covers:
- Creating config templates
- Adding metadata and arguments
- Finding suitable containers
- Generating help files
- Setting up inputs/outputs
- Docker engine configuration

</details>

<details>
<summary><strong>Script Development Best Practices</strong></summary>

See: **[Script Development Guide](docs/SCRIPT_DEVELOPMENT.md)**

This guide covers:
- Array-based command building
- Error handling with `set -eo pipefail`
- Conditional parameter inclusion
- Boolean flag management
- Parameter expansion patterns

</details>

<details>
<summary><strong>Testing Guidelines</strong></summary>

See: **[Testing Guide](docs/TESTING.md)**

This guide covers:
- Generating test data programmatically
- Self-contained test design
- Output validation techniques
- Helper function patterns
- Multi-scenario testing

</details>

<details>
<summary><strong>Docker Configuration</strong></summary>

See: **[Docker Guide](docs/DOCKER_GUIDE.md)**

This guide covers:
- Biocontainer selection
- Version pinning strategies
- Software version detection
- Custom container creation
- Container best practices

</details>

## Contribution Process

### Submitting Your Component

1. **Test thoroughly**: Ensure your component passes all tests
   ```bash
   viash test src/namespace/tool_name/config.vsh.yaml
   ```

2. **Build successfully**: Verify the component builds without errors
   ```bash
   viash build src/namespace/tool_name/config.vsh.yaml --setup cachedbuild
   ```

3. **Follow naming conventions**: Use `snake_case` for component names and arguments

4. **Update documentation**: Add your component to relevant documentation if needed

5. **Create a pull request**: Submit your changes with a clear description

### Review Process

- All contributions go through code review
- Components must pass automated tests
- Docker containers must be properly versioned
- Documentation must be complete and accurate

## Getting Help

### Resources

- **[Viash Documentation](https://viash.io/)**
- **[GitHub Discussions](https://github.com/viash-io/biobox/discussions)**
- **[Issue Tracker](https://github.com/viash-io/biobox/issues)**

### Common Questions

**Q: How do I find the right Docker container?**  
A: Search for "biocontainer [tool_name]" or check [quay.io/biocontainers](https://quay.io/organization/biocontainers)

**Q: My component fails to build. What should I check?**  
A: Verify the Docker image exists, check syntax in config.vsh.yaml, and ensure all required commands are available

**Q: How do I handle tools with complex argument patterns?**  
A: Check existing similar components for patterns, or ask in GitHub Discussions

**Q: Can I create custom Docker containers?**  
A: Yes, but prefer biocontainers when available. See the [Docker Guide](docs/DOCKER_GUIDE.md) for details.

---

Happy contributing!
