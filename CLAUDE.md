# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview
C11 translation of SLICOT (Subroutine Library In Control Theory) from Fortran77.

**Reference:** Fortran77 source in `SLICOT-Reference/src/`

## Build Commands

```bash
# Configure and build
cmake --preset macos-arm64-debug
cmake --build --preset macos-arm64-debug-build

# Test
ctest --preset macos-arm64-debug-test              # All tests
ctest --preset macos-arm64-debug-test -R AB01MD    # Specific test
ctest --preset macos-arm64-debug-test -V           # Verbose

# Clean rebuild
rm -rf build/macos-arm64-debug && cmake --preset macos-arm64-debug && cmake --build --preset macos-arm64-debug-build
```

**Presets:** `macos-x64-{debug,release}`, `macos-arm64-{debug,release}`

## Architecture

### Directory Structure

- `src/` - C11 sources, organized by SLICOT chapter (AB/, MB/, MC/, etc.)
- `include/slicot.h` - Public API declarations
- `tests/unit/` - GTest unit tests
- `SLICOT-Reference/src/` - Fortran source (XXYYZZ.f format)
- `SLICOT-Reference/doc/` - HTML documentation with examples
- `tools/extract_dependencies.py` - Dependency analyzer
- `.claude/skills/slicot-knowledge/` - Translation knowledge base

### SLICOT Naming: `XXYYZZ`
- `XX` - Chapter (AB, MB, MC, etc.)
- `YY` - Section number
- `ZZ` - Variant

**Translation:** `AB01MD.f` → `src/AB/ab01md.c` (lowercase)

## Translation Workflow: TDD Cycle

**🚨 MANDATORY: RED → GREEN → REFACTOR → VERIFY**

### Phase 1: RED - Write Tests First

1. **Check dependencies:**
   ```bash
   python tools/extract_dependencies.py SLICOT-Reference/src/ AB01MD
   ```

2. **Create test** `tests/unit/test_ab01md.cpp` with 3+ tests:
   - Basic case (from `SLICOT-Reference/doc/AB01MD.html`)
   - Edge case (N=0)
   - Error case (invalid parameters)

3. **Add to** `tests/unit/CMakeLists.txt`:
   ```cmake
   set(UNIT_TEST_SOURCES
       test_ab01md.cpp
   )
   ```

4. **Verify failure:**
   ```bash
   cmake --build --preset macos-arm64-debug-build  # Should fail - function doesn't exist
   ```

5. **Commit:**
   ```bash
   git commit -m "RED: Add tests for AB01MD"
   ```

### Phase 2: GREEN - Implement

1. **Create** `src/AB/ab01md.c`:
   ```c
   #include "slicot.h"
   #include <stdint.h>

   void ab01md(int32_t n, double* a, int32_t* info) {
       *info = 0;
       if (n < 0) { *info = -1; return; }
       if (n == 0) return;
       // ... translate from Fortran ...
   }
   ```

2. **Update build files:**
   - Add to `src/CMakeLists.txt` → `SLICOT_SOURCES`
   - Add declaration to `include/slicot.h`

3. **Verify tests pass:**
   ```bash
   cmake --build --preset macos-arm64-debug-build
   ctest --preset macos-arm64-debug-test -R AB01MD -V
   ```

4. **Commit:**
   ```bash
   git commit -m "GREEN: Implement AB01MD"
   ```

### Phase 3: REFACTOR - Clean Up

1. Review: minimal comments, verify BLAS/LAPACK usage
2. Test: `ctest --preset macos-arm64-debug-test -R AB01MD -V`
3. Commit: `git commit -m "REFACTOR: Clean up AB01MD"`

### Phase 4: VERIFY - Full Suite

```bash
ctest --preset macos-arm64-debug-test  # All tests must pass
```

## Translation Strategy

**Bottom-up by dependency level:**

```bash
# Find Level 0 (leaves - only BLAS/LAPACK deps)
python tools/extract_dependencies.py SLICOT-Reference/src/ | grep "Level 0"
```

- Level 0 first (297 routines, can parallelize)
- Then Level 1, 2, etc.
- Always check dependencies before translating

## Key Patterns

**Types:**
- `INTEGER` → `int32_t`
- `DOUBLE PRECISION` → `double`

**Arrays:**
- 0-based indexing (Fortran uses 1-based)
- Column-major layout preserved
- Leading dimensions (lda, ldb) preserved

**Error codes (`info` parameter):**
- `0` = success
- `< 0` = parameter `-info` invalid
- `> 0` = algorithm error

**CRITICAL:** Never reimplement LAPACK routines - use existing macros/implementations

## Quality Checklist

- [ ] Tests written FIRST (RED)
- [ ] Tests initially failed
- [ ] Min 3 tests (basic, edge, error)
- [ ] All tests pass (GREEN)
- [ ] Code cleaned (REFACTOR)
- [ ] Full suite passes (VERIFY)
- [ ] TDD commits (RED/GREEN/REFACTOR)
- [ ] Test data from SLICOT docs
- [ ] BLAS/LAPACK used correctly

## Reference

- `tools/README.md` - Detailed workflow, examples
- `.claude/skills/slicot-knowledge/SKILL.md` - Comprehensive knowledge base
- `SLICOT-Reference/doc/*.html` - Official documentation
