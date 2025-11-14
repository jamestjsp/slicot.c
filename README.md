# slicot.c

C11 translation of SLICOT (Subroutine Library In Control Theory) from Fortran77.

## Quick Start

```bash
# Setup & build
./scripts/setup_venv.sh && source venv/bin/activate
cmake --preset macos-arm64-debug
cmake --build --preset macos-arm64-debug-build
pip install -e .

# Test
pytest tests/python/ -v
```

**Presets:** `macos-{x64,arm64}-{debug,release}`

## Translation Status

| Category | Status | Routines |
|----------|--------|----------|
| **MB01** | ✅ 1/X | mb01qd |
| **MB03** | ✅ 1/X | mb03oy |
| **TG01** | ✅ 1/X | tg01fd |
| **Level 0** | 📋 3/297 | (no SLICOT deps) |

**Tests:** 21/21 passing (mb01qd: 7, mb03oy: 9, tg01fd: 5)

## Features

- Column-major storage (Fortran-compatible)
- Python bindings (NumPy arrays)
- TDD workflow (RED→GREEN→REFACTOR)

## Docs

- **[CLAUDE.md](CLAUDE.md)** - Development workflow & translation patterns
