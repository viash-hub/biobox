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

- **Follow modern standards**: Use current coding patterns and component structure
- **Ensure reproducibility**: Pin versions and document dependencies clearly
- **Generate test data**: Create self-contained tests that don't rely on external files
- **Write clean code**: Use consistent naming and clear, maintainable scripts

For detailed implementation guidelines, check out our development guides:

## Development Guides

### 🔧 [Component Development Guide](docs/COMPONENT_DEVELOPMENT.md)
How to create components: config templates, metadata, arguments, containers, help files, and Docker setup.

### 📝 [Script Development Guide](docs/SCRIPT_DEVELOPMENT.md) 
Writing good scripts: array-based commands, error handling, conditional parameters, boolean flags, and parameter patterns.

### ✅ [Testing Guide](docs/TESTING.md)
Testing your components: self-contained tests, generating test data, output validation, and testing multiple scenarios.

### 🐳 [Docker Guide](docs/DOCKER_GUIDE.md)
Working with containers: choosing biocontainers, version pinning, detecting software versions, and container best practices.

## Contribution Process

### Submitting Your Component

1. **Test thoroughly**: Ensure your component passes all tests
   ```bash
   viash test src/namespace/tool_name/config.vsh.yaml
   ```

2. **Add changelog entry**: Document your changes in `CHANGELOG.md` under the "Unreleased" section

3. **Review your changes**: Check your code for:
   - Consistent naming and coding conventions
   - Clear, maintainable code structure
   - Proper error handling and edge cases
   - Complete documentation and helpful comments

4. **Create a pull request**: Submit your changes.
  - Include a clear description of the changes you've made
  - Link to any relevant issues or discussions
  - Review the changes critically before submitting the PR

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
