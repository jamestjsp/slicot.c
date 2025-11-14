# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview
C11 translation of SLICOT (Subroutine Library In Control Theory) from Fortran77.

**Reference:** Fortran77 source in `SLICOT-Reference/src/`

## Build System

**CMake with presets** for cross-platform builds. Uses Ninja generator.

### Common Commands

**Configure:**
```bash
cmake --preset macos-arm64-debug    # or macos-x64-debug, macos-arm64-release, macos-x64-release
```

**Build:**
```bash
cmake --build --preset macos-arm64-debug-build
```

**Test:**
```bash
# All tests
ctest --preset macos-arm64-debug-test

# Specific test pattern
ctest --preset macos-arm64-debug-test -R AB01MD

# Verbose output
ctest --preset macos-arm64-debug-test -V
```

**Clean rebuild:**
```bash
rm -rf build/macos-arm64-debug
cmake --preset macos-arm64-debug
cmake --build --preset macos-arm64-debug-build
```

### Build Options

Set in `CMakeLists.txt` or via `-D` flag:
- `BUILD_SHARED_LIBS` (default: ON) - Build shared vs static library
- `BUILD_TESTING` (default: ON) - Build GTest test suite
- `BUILD_EXAMPLES` (default: ON) - Build example programs

## Architecture

### Directory Structure

```
slicot.c/
├── src/                      # C11 library sources (empty, to be populated)
│   └── CMakeLists.txt        # Add sources to SLICOT_SOURCES list
├── include/
│   └── slicot.h              # Public API header
├── tests/
│   ├── unit/                 # GTest unit tests
│   └── integration/          # GTest integration tests
├── examples/                 # Example programs
├── SLICOT-Reference/         # Fortran77 reference implementation
│   ├── src/                  # Fortran source files (XXYYZZ.f format)
│   └── doc/                  # HTML documentation
├── tools/
│   ├── extract_dependencies.py  # Analyze Fortran dependencies
│   └── README.md             # Translation workflow guide
└── .claude/skills/slicot-knowledge/  # SLICOT translation knowledge base
```

### SLICOT Organization

Routines named `XXYYZZ`:
- `XX` - Chapter (AB=Analysis/Benchmark, MB=Matrix ops, MC=Polynomial, etc.)
- `YY` - Section number
- `ZZ` - Variant

Translation convention: `AB01MD.f` → `src/AB/ab01md.c` (lowercase)

## Translation Workflow: TDD Cycle

**🚨 MANDATORY: Follow strict Test-Driven Development (TDD)**

All translations MUST follow the **RED → GREEN → REFACTOR → VERIFY** cycle.

### TDD Phases (Non-Negotiable)

#### Phase 1: RED - Write Tests First

**Tests MUST fail initially** (function doesn't exist yet).

1. **Check dependencies:**
   ```bash
   python tools/extract_dependencies.py SLICOT-Reference/src/ AB01MD
   ```

2. **Create test file** `tests/unit/test_ab01md.cpp`:
   ```cpp
   #include <gtest/gtest.h>
   #include "slicot.h"

   TEST(AB01MD, BasicCase) {
       // Parse test data from SLICOT-Reference/doc/AB01MD.html
       int n = 3;
       double a[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};  // Column-major
       int info;

       ab01md(n, a, &info);  // Function doesn't exist - will fail compilation

       EXPECT_EQ(info, 0);
       EXPECT_NEAR(a[0], 1.5, 1e-14);
   }

   TEST(AB01MD, EdgeCaseZeroN) {
       int n = 0, info;
       double a[] = {0.0};
       ab01md(n, a, &info);
       EXPECT_EQ(info, 0);  // Quick return
   }

   TEST(AB01MD, InvalidParameter) {
       int n = -1, info;
       double a[] = {0.0};
       ab01md(n, a, &info);
       EXPECT_LT(info, 0);  // Negative info = parameter error
   }
   ```

3. **Add to `tests/unit/CMakeLists.txt`:**
   ```cmake
   set(UNIT_TEST_SOURCES
       test_ab01md.cpp    # Add here
   )
   ```

4. **Verify tests FAIL:**
   ```bash
   cmake --build --preset macos-arm64-debug-build  # Compilation error expected
   ```

5. **Commit RED phase:**
   ```bash
   git add tests/unit/test_ab01md.cpp tests/unit/CMakeLists.txt
   git commit -m "RED: Add tests for AB01MD"
   ```

#### Phase 2: GREEN - Minimal Implementation

**Write minimal code to make tests pass.**

1. **Create source** `src/AB/ab01md.c`:
   ```c
   #include "slicot.h"
   #include <stdint.h>

   void ab01md(int32_t n, double* a, int32_t* info) {
       *info = 0;

       // Parameter validation
       if (n < 0) {
           *info = -1;
           return;
       }

       // Quick return
       if (n == 0) return;

       // Minimal implementation from Fortran
       // ... translate algorithm ...
   }
   ```

2. **Add to `src/CMakeLists.txt`:**
   ```cmake
   set(SLICOT_SOURCES
       AB/ab01md.c    # Add here
   )
   ```

3. **Add declaration to `include/slicot.h`:**
   ```c
   void ab01md(int32_t n, double* a, int32_t* info);
   ```

4. **Build and verify tests PASS:**
   ```bash
   cmake --build --preset macos-arm64-debug-build
   ctest --preset macos-arm64-debug-test -R AB01MD -V
   ```

5. **Commit GREEN phase:**
   ```bash
   git add src/AB/ab01md.c src/CMakeLists.txt include/slicot.h
   git commit -m "GREEN: Implement AB01MD"
   ```

#### Phase 3: REFACTOR - Clean Up

**Improve code quality, ensure no regressions.**

1. **Review implementation:**
   - Add minimal comments for complex logic
   - Verify numerical correctness
   - Check BLAS/LAPACK usage

2. **Verify tests still pass:**
   ```bash
   ctest --preset macos-arm64-debug-test -R AB01MD -V
   ```

3. **Commit REFACTOR phase:**
   ```bash
   git add src/AB/ab01md.c
   git commit -m "REFACTOR: Clean up AB01MD implementation"
   ```

#### Phase 4: VERIFY - Final Validation

**Ensure no regressions in full test suite.**

```bash
# Full test suite
ctest --preset macos-arm64-debug-test

# Confirm all tests pass
echo $?  # Should be 0
```

### Translation Strategy

**Dependency-driven bottom-up approach:**

1. **Find leaf routines (Level 0):**
   ```bash
   python tools/extract_dependencies.py SLICOT-Reference/src/ | grep "Level 0"
   ```
   - Level 0 routines only depend on BLAS/LAPACK
   - Can translate in parallel
   - Foundation for higher levels

2. **Translate by levels:**
   - Level 0 first (297 routines - highest priority)
   - Level 1 next (depend only on Level 0)
   - Continue until all levels complete

3. **For each routine:**
   - Run dependency analysis
   - Ensure dependencies exist
   - Follow TDD cycle (RED → GREEN → REFACTOR → VERIFY)
   - Minimum 3 tests per routine

### Key Translation Patterns

**Types:**
- `INTEGER` → `int32_t`
- `DOUBLE PRECISION` → `double`
- `COMPLEX*16` → TBD (when needed)

**Arrays:**
- 0-based indexing (adjust from Fortran's 1-based)
- Column-major layout (Fortran-compatible)
- Leading dimensions (lda, ldb, etc.) preserved

**BLAS/LAPACK:**
- Use macros when available
- Check reference implementation
- Never reimplement LAPACK routines

**Error handling:**
- `info` parameter (by pointer)
- `info = 0` success
- `info < 0` parameter -info is invalid
- `info > 0` algorithm-specific error

## Dependency Analysis Tool

`tools/extract_dependencies.py` parses Fortran source to identify:
- Call graph between SLICOT routines
- LAPACK/BLAS dependencies
- Translation order (by level)

**Usage:**
```bash
# Full analysis
python tools/extract_dependencies.py SLICOT-Reference/src/

# Specific routine
python tools/extract_dependencies.py SLICOT-Reference/src/ AB01MD

# Find leaf routines (Level 0)
python tools/extract_dependencies.py SLICOT-Reference/src/ | grep "Level 0"
```

**Output symbols:**
- `✓` Routine exists in source
- `✗` External dependency

Level 0 routines can be translated in parallel as they only depend on BLAS/LAPACK.

## Testing Strategy

- **Framework:** Google Test (GTest)
- **Test data:** Extract from `SLICOT-Reference/doc/*.html` Example sections
- **Tolerances:** 1e-14 for EXPECT_NEAR (tight numerical validation)
- **Coverage:** Basic cases, edge cases (N=0, singular matrices), error conditions

## Code Conventions

- **C Standard:** C11
- **Naming:** Lowercase for functions (`ab01md`), match Fortran routine names
- **Arrays:** Column-major (Fortran-compatible), 0-based indexing
- **Comments:** Minimal - only for bug fixes, TODO, known issues

## Reference Documentation

- `tools/README.md` - Complete translation workflow, patterns, examples
- `.claude/skills/slicot-knowledge/SKILL.md` - Comprehensive SLICOT knowledge base
- `SLICOT-Reference/doc/*.html` - Official routine documentation with examples
