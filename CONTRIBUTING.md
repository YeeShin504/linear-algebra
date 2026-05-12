# Contributing to ma1522-linear-algebra


Thank you for your interest in contributing! All types of contributions are welcome, including code, documentation, tests, and bug reports.

## Getting Started

- This project requires **Python 3.10 or higher**.
- Dependencies are managed with [uv](https://github.com/astral-sh/uv) and specified in `pyproject.toml` and `uv.lock`.
- To set up the development environment:

```sh
uv venv
uv pip install -e .[dev,docs]
```

## Project Structure

Key files and directories:
- `src/ma1522/symbolic.py` — Main `Matrix` class with linear algebra operations.
- `src/ma1522/custom_types.py` — Type definitions for decompositions (SVD, QR, LU, etc.) and solutions.
- `tests/` — Test suite with unit tests and doctests.
- `docs/` — Markdown documentation and tutorial notebooks.

## Linting & Formatting

Use [ruff](https://docs.astral.sh/ruff/) for linting and formatting:

```sh
ruff check src/ tests/
ruff format src/ tests/
```

## Pre-commit Hooks

This project uses [pre-commit](https://pre-commit.com/) to automate code quality checks on every commit.

### Setup

Install the git hooks:

```sh
pre-commit install
```

### Running Manually

Run all hooks on all files:

```sh
pre-commit run --all-files
```

Run a specific hook:

```sh
pre-commit run ruff --all-files
pre-commit run check-yaml --all-files
```

### What Hooks Are Configured

The `.pre-commit-config.yaml` file defines:
- **ruff**: Linting and code style checks (enforces docstrings, naming conventions, complexity rules, etc.)
- **check-yaml**: Validates YAML syntax (with `--unsafe` flag for Python-specific YAML tags in `mkdocs.yml`)
- **check-merge-conflict**: Detects merge conflict markers
- **end-of-file-fixer**: Ensures files end with a newline

These hooks run automatically before each commit to catch issues early. If a hook fails, the commit is blocked. Fix the issues and try again.

## Running Tests

Run tests using [pytest](https://docs.pytest.org/):

```sh
# Run all tests
pytest

# Run specific test file
pytest tests/test_decompositions.py

# Run with verbose output
pytest -v

# Run doctests
pytest --doctest-modules src/ma1522/symbolic.py
```

### Note on Performance

Some decomposition tests (especially SVD on complex symbolic matrices) may take several seconds. This is expected behavior. The `singular_value_decomposition` method includes a configurable timeout (default: 30 seconds) with automatic fallback to numerical QR-based completion for long-running computations.

## Building Documentation

- Documentation is built with [MkDocs](https://www.mkdocs.org/):

```sh
mkdocs serve
```

## Submitting Contributions

- Open an issue or pull request for any change.
- Please include clear commit messages and relevant tests/docs if applicable.
- Add tests for new features and ensure all tests pass locally before submitting.
- Update docstrings with examples and type hints.
- There are no strict rules for branch naming or PR structure—just be clear and respectful.


---

Thank you for helping improve this project!
