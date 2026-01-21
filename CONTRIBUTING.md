# Contributing to NRXN1 Structural Variant Analysis Pipeline

Thank you for your interest in contributing! This document provides guidelines for contributing to this project.

## Getting Started

1. Fork the repository
2. Clone your fork:
   ```bash
   git clone https://github.com/yourusername/nrxn1-sv-analysis.git
   cd nrxn1-sv-analysis
   ```
3. Create a virtual environment:
   ```bash
   python -m venv venv
   source venv/bin/activate
   pip install -r requirements.txt
   ```
4. Install development dependencies:
   ```bash
   pip install -e ".[dev]"
   ```

## Development Workflow

### Branch Naming

- `feature/` - New features
- `fix/` - Bug fixes
- `docs/` - Documentation updates
- `refactor/` - Code refactoring

### Code Style

We use the following tools for code quality:

- **Black** for code formatting
- **Flake8** for linting
- **MyPy** for type checking

Run before committing:
```bash
black src/ tests/
flake8 src/ tests/
mypy src/
```

### Testing

Run tests with pytest:
```bash
pytest tests/ -v
pytest tests/ -v --cov=src --cov-report=html
```

All new features should include tests.

## Pull Request Process

1. Create a feature branch from `main`
2. Make your changes
3. Run tests and linting
4. Update documentation if needed
5. Submit a pull request

### PR Checklist

- [ ] Code follows project style guidelines
- [ ] Tests pass locally
- [ ] New features include tests
- [ ] Documentation updated
- [ ] Commit messages are clear

## Reporting Issues

When reporting issues, please include:

- Python version
- Operating system
- Steps to reproduce
- Expected vs actual behavior
- Error messages/logs

## Feature Requests

We welcome feature requests! Please:

1. Check existing issues first
2. Describe the use case
3. Explain expected behavior

## Code of Conduct

- Be respectful and inclusive
- Provide constructive feedback
- Focus on the code, not the person

## License

By contributing, you agree that your contributions will be licensed under the MIT License.
