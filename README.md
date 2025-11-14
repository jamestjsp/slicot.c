# slicot.c

C11 translation of SLICOT (Subroutine Library In Control Theory) from Fortran77.

[![Tests](https://img.shields.io/badge/tests-pytest-blue.svg)](tests/python/)
[![Python](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/)

---

## Quick Start

### 1. Setup Python Environment

```bash
# Create isolated virtual environment
./scripts/setup_venv.sh
source venv/bin/activate
```

### 2. Build

```bash
# Configure & build C library
cmake --preset macos-arm64-debug
cmake --build --preset macos-arm64-debug-build

# Install Python bindings
pip install -e .
```

### 3. Test

```bash
# Run all tests
pytest tests/python/ -v

# Or via CTest
ctest --preset macos-arm64-debug-test
```

---

## Documentation

- **[CLAUDE.md](CLAUDE.md)** - Development workflow (TDD, translation patterns)
- **[PYTHON_SETUP.md](PYTHON_SETUP.md)** - Virtual environment guide
- **[TESTING.md](TESTING.md)** - Testing infrastructure
- **[SLICUTLET_ANALYSIS.md](SLICUTLET_ANALYSIS.md)** - Reference implementation analysis

---

## Key Features

✅ **Column-major storage** - Fortran-compatible, BLAS/LAPACK ready
✅ **Python bindings** - NumPy arrays with auto-computed dimensions
✅ **pytest tests** - Modern testing with SciPy references
✅ **Virtual environment** - Isolated Python dependencies
✅ **TDD workflow** - RED → GREEN → REFACTOR → VERIFY

---

## Translation Status

| Category | Status | Routines |
|----------|--------|----------|
| **MB01** | ✅ 1/X | mb01qd |
| **MB03** | ✅ 1/X | mb03oy |
| **Level 0** | 📋 2/297 | (no SLICOT dependencies) |

---

## CMake Presets

Available presets:
- `macos-x64-debug` / `macos-x64-release`
- `macos-arm64-debug` / `macos-arm64-release`

Build options:
- `BUILD_SHARED_LIBS`: Build shared libraries (default: ON)
- `BUILD_TESTING`: Build tests (default: ON)
- `BUILD_PYTHON_BINDINGS`: Build Python bindings (default: ON)
