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

## Linting & Formatting

- Use [ruff](https://docs.astral.sh/ruff/) for linting and formatting:

```sh
ruff check src/ tests/
ruff format src/ tests/
```

## Running Tests

- Run tests using [pytest](https://docs.pytest.org/):

```sh
pytest
```

## Building Documentation

- Documentation is built with [MkDocs](https://www.mkdocs.org/):

```sh
mkdocs serve
```


## Submitting Contributions

- Open an issue or pull request for any change.
- Please include clear commit messages and relevant tests/docs if applicable.
- There are no strict rules for branch naming or PR structure—just be clear and respectful.


---

Thank you for helping improve this project!