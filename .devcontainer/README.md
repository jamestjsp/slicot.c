# Dev Container

Ubuntu-based C++ dev environment with:
- GCC, CMake, gfortran
- Python 3.11, pytest
- BLAS/LAPACK libraries
- GitHub CLI

## Usage

VSCode: Reopen in Container
CLI: `devcontainer open .`

## First Run

Container runs `postCreateCommand`:
1. Installs BLAS/LAPACK, gfortran, ninja
2. Sets up Python venv
3. Builds C library (linux-x64-debug)
4. Installs Python package in editable mode

## Build & Test

```bash
source venv/bin/activate
cmake --preset linux-x64-debug
cmake --build --preset linux-x64-debug-build
pytest tests/python/ -v
```
