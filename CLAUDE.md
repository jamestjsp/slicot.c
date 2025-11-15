# CLAUDE.md

C11 translation of SLICOT (Subroutine Library In Control Theory) from Fortran77.

**Reference:** `SLICOT-Reference/src/` (Fortran77), `SLICOT-Reference/doc/` (HTML docs with examples)

## Quick Commands

```bash
# Setup (first time)
./scripts/setup_venv.sh && source venv/bin/activate

# Build & test
source venv/bin/activate  # Always activate first!
cmake --preset macos-arm64-debug
cmake --build --preset macos-arm64-debug-build
pip install -e .
pytest tests/python/ -v

# Clean rebuild
rm -rf build/macos-arm64-debug && cmake --preset macos-arm64-debug && cmake --build --preset macos-arm64-debug-build && pip install -e .

# Memory/pointer bug detection (sanitizers - recommended before PR)
source venv/bin/activate
cmake --preset macos-arm64-debug-sanitizers
cmake --build --preset macos-arm64-debug-sanitizers-build
pip install -e .
pytest tests/python/ -v
```

**Presets:**
- Standard: `macos-{x64,arm64}-{debug,release}`
- Sanitizers: `macos-arm64-debug-{asan,ubsan,sanitizers}`

## Directory Structure

```
src/XX/routine.c          # C11 implementation (XX=AB,MB,MC...)
include/slicot.h          # Public API + Doxygen docs
include/slicot_types.h    # Type aliases (i32, f64)
include/slicot_blas.h     # BLAS/LAPACK wrappers (SLC_DGEMM, etc.)
python/slicot_module.c    # Python/C extension (manual, not generated)
python/slicot/__init__.py # Python exports
tests/python/test_*.py    # pytest tests
tools/extract_dependencies.py  # Dependency analyzer
```

**Naming:** `AB01MD.f` (Fortran) → `src/AB/ab01md.c` (C, lowercase)

## TDD Workflow: RED → GREEN → REFACTOR → VERIFY

### RED: Write Test First

1. **Use skill**: `/skill slicot-knowledge` for parsing docs & test data
2. Check dependencies: `python3 tools/extract_dependencies.py SLICOT-Reference/src/ AB01MD`
3. Find test data (priority order):
   - `SLICOT-Reference/doc/AB01MD.html` (parse with skill)
   - `SLICOT-Reference/examples/TAB01MD.f` + `data/AB01MD.dat` + `results/AB01MD.res`
   - `SLICOT-Reference/benchmark_data/BB01*.dat`
   - NumPy/SciPy synthetic data (get approval first)
4. Create `tests/python/test_ab01md.py` with 3+ tests (basic, edge, error)
5. Verify failure: `pytest tests/python/test_ab01md.py -v`
6. Commit: `git commit -m "RED: Add tests for AB01MD"`

**Test template:**
```python
import pytest
import numpy as np
from slicot import ab01md

def test_ab01md_basic():
    a = np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float64, order='F')
    result, info = ab01md(a)
    assert info == 0
    np.testing.assert_allclose(result, expected, rtol=1e-14)
```

### GREEN: Implement

1. Create `src/AB/ab01md.c`
2. Update:
   - `src/CMakeLists.txt` → add to `SLICOT_SOURCES`
   - `include/slicot.h` → add declaration + Doxygen
   - `python/slicot_module.c` → add wrapper (see existing for pattern)
   - `python/slicot/__init__.py` → export function
3. Build & test: `cmake --build --preset macos-arm64-debug-build && pip install -e . && pytest tests/python/test_ab01md.py -v`
4. Commit: `git commit -m "GREEN: Implement AB01MD"`

### REFACTOR: Clean Up

1. Minimal comments (only for bugs/TODOs/known issues)
2. Verify BLAS/LAPACK usage
3. Test: `pytest tests/python/test_ab01md.py -v`
4. Commit: `git commit -m "REFACTOR: Clean up AB01MD"`

### VERIFY: Full Suite

```bash
pytest tests/python/ -v  # All tests must pass

# RECOMMENDED: Run with sanitizers to catch memory/pointer bugs
cmake --preset macos-arm64-debug-sanitizers
cmake --build --preset macos-arm64-debug-sanitizers-build
pip install -e .
pytest tests/python/ -v
```

**Sanitizer benefits**: Catches index conversion bugs (out-of-bounds), use-after-free, double-free, integer overflow automatically

## Key Patterns

### Types
- `INTEGER` → `i32` (int32_t)
- `DOUBLE PRECISION` → `f64` (double)
- `LOGICAL` → `bool`

### Arrays (Column-Major)
- **Index:** `a[i + j*lda]` (row i, col j)
- **Storage:** Column-major (Fortran-compatible, zero BLAS/LAPACK overhead)
- **NumPy:** `order='F'` mandatory in tests
- **Example:**
  ```c
  // Math: [1 2]  Memory: [1, 3, 2, 4]  (col0, col1)
  //       [3 4]
  double a[] = {1.0, 3.0, 2.0, 4.0};  // NOT {1, 2, 3, 4}
  ```

### Error Codes
- `info = 0`: success
- `info < 0`: parameter `-info` invalid
- `info > 0`: algorithm error

### BLAS/LAPACK
- Use `SLC_DGEMM()` etc. from `include/slicot_blas.h`
- Pass scalars by pointer: `SLC_DGEMM("N", "N", &m, &n, &k, &alpha, a, &lda, ...)`
- Never reimplement LAPACK routines

### Index Conversion (Fortran ↔ C)
- **SLICOT routines return 1-based indices** (e.g., `mb03oy` pivot arrays)
- **Always convert to 0-based before C array access:** `k = iwork[j] - 1;`
- **CRITICAL: Validate bounds after conversion to prevent out-of-bounds access**
- Common bugs:
  - Using 1-based value directly as C index causes wrong array access
  - Using converted index without bounds check risks buffer overflow

**Correct pattern:**
```c
k = iwork[j] - 1;           // Convert 1-based to 0-based
if (k < 0 || k >= n) {      // REQUIRED: Validate bounds
    break;                  // Handle invalid index
}
// Safe to use k as array index
if (iwork[k] < 0) {
    // ...
}
```

**Vulnerable pattern (DO NOT USE):**
```c
k = iwork[j] - 1;
if (iwork[k] < 0) {         // BUG: No bounds check before iwork[k]!
    // ...
}
```

### Python Wrapper Pattern

**Input arrays (modify in-place):**
```c
// Use NPY_ARRAY_FARRAY to preserve Fortran order
a_array = (PyArrayObject*)PyArray_FROM_OTF(a_obj, NPY_DOUBLE,
                                           NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
```

**Output arrays (allocated in C):**
```c
// Use PyArray_New with explicit Fortran strides
npy_intp dims[2] = {m, n};
npy_intp strides[2] = {sizeof(f64), m * sizeof(f64)};  // Column-major
q_array = PyArray_New(&PyArray_Type, 2, dims, NPY_DOUBLE, strides, q, 0, NPY_ARRAY_FARRAY, NULL);
PyArray_ENABLEFLAGS((PyArrayObject*)q_array, NPY_ARRAY_OWNDATA);
```

**CRITICAL: In-place modification - return input array directly:**
```c
// When routine modifies input array in-place (e.g., B becomes U):
// CORRECT: Return the modified input array
PyObject *result = Py_BuildValue("Odi", b_array, scale, info);
Py_DECREF(a_array);
Py_DECREF(b_array);

// WRONG: Creating wrapper around same memory causes double-free
// PyObject *u_array = PyArray_New(..., b_data, ...);  // DON'T DO THIS!
// PyArray_ENABLEFLAGS(u_array, NPY_ARRAY_OWNDATA);    // b_data already owned by b_array
// Result: segfault from double-free on cleanup
```

**Register in method table:**
```c
static PyMethodDef SlicotMethods[] = {
    {"routine", py_routine, METH_VARARGS, "Docstring..."},
    // ...
};
```

**Export in `python/slicot/__init__.py`:**
```python
from ._slicot import routine
__all__ = ['routine', ...]
```

## Quality Checklist

- [ ] Tests from `SLICOT-Reference/doc/*.html`
- [ ] NumPy arrays use `order='F'`
- [ ] Python wrapper uses `NPY_ARRAY_FARRAY`
- [ ] Function exported in `__init__.py`
- [ ] Doxygen docs in `include/slicot.h`
- [ ] Min 3 tests (basic, edge, error)
- [ ] TDD commits (RED/GREEN/REFACTOR)

## Translation Strategy

Bottom-up by dependency level (use `tools/extract_dependencies.py`):

```bash
# Find Level 0 (297 routines, only BLAS/LAPACK deps, parallelizable)
python3 tools/extract_dependencies.py SLICOT-Reference/src/ | grep "Level 0"
```

Translate Level 0 first, then Level 1, 2, etc.

## Debugging: Fortran Diagnostic Workflow

When C implementation produces wrong results, compare against Fortran reference:

```bash
# Enable diagnostic tools
cmake --preset macos-arm64-debug -DBUILD_FORTRAN_DIAG=ON

# Run side-by-side comparison (builds, executes, compares)
cmake --build --preset macos-arm64-debug-build --target diag_all

# View results
grep -A3 "OUTPUT" build/macos-arm64-debug/fortran_diag/*_trace.txt
```

**When to use:**
- Test fails but unsure where C diverges from Fortran
- Need to isolate bug to specific computation step
- Debugging complex multi-step algorithms (Schur, Lyapunov, etc.)

**Output:** High-precision traces (16 digits) in `build/macos-arm64-debug/fortran_diag/`

**See:** `fortran_diag/README.md` for adding new routines

## Reference Docs

- `fortran_diag/README.md` - Diagnostic workflow for debugging C vs Fortran
- `tools/README.md` - Detailed workflow examples
- `.claude/skills/slicot-knowledge/SKILL.md` - Translation knowledge base
- `scripts/setup_venv.sh` - Virtual environment setup
