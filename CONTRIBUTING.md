# Contributing to QC Library

Thank you for your interest in contributing to the QC Library! This document outlines the standards and workflows for contributing to the project.

## Code Style

We adhere to the **C11** standard. Please ensure your compiler supports C11.

- **Indentation**: 4 spaces. No tabs.
- **Braces**: Opening brace on the same line (K&R style).
- **Naming**: `snake_case` for variables and functions. `PascalCase` for types/structs (e.g., `QuantumRegister`).
- **Comments**: Use Doxygen-style comments (`/** ... */`) for public API functions.

## Thread Safety

- **RNG**: The library uses `_Thread_local` for its internal random number generator state. Do not use global mutable state that is not thread-safe.
- **OpenMP**: Core loops should be parallelized using OpenMP directives where appropriate.

## Testing

All changes must pass the test suite.

```bash
# Run standard tests
make test
```

If you are adding new gates or features, please add corresponding tests in `tests/test_suite.c`.

## Pull Request Process

1. Fork the repository.
2. Create a feature branch (`git checkout -b feature/amazing-feature`).
3. Commit your changes.
4. Push to the branch.
5. Open a Pull Request.

## Performance

- Use `QC_USE_DOUBLE` (default) for development.
- performance critical code should be verified with `-O3` and release build flags.
