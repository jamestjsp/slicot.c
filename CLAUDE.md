# CLAUDE.md

C11 translation of SLICOT (Subroutine Library In Control Theory) from Fortran77.

**Reference:** `SLICOT-Reference/src/` (Fortran77), `SLICOT-Reference/doc/` (HTML docs)

## Quick Commands

```bash
# Setup (first time)
./scripts/setup_venv.sh && source venv/bin/activate

# Build & test
source venv/bin/activate
cmake --preset linux-x64-debug
cmake --build --preset linux-x64-debug-build
pip install -e .
pytest tests/python/ -v

# Clean rebuild
rm -rf build/linux-x64-debug && cmake --preset linux-x64-debug && cmake --build --preset linux-x64-debug-build && pip install -e .

# Sanitizers (run before PR to catch memory bugs)
cmake --preset linux-x64-debug-sanitizers
cmake --build --preset linux-x64-debug-sanitizers-build
pip install -e . --force-reinstall --no-deps
# NOTE: Command substitution doesn't work in Claude Code - use explicit path
# Linux devcontainer (aarch64): /usr/lib/gcc/aarch64-linux-gnu/13/libasan.so
# Find path with: gcc -print-file-name=libasan.so
LD_PRELOAD=/usr/lib/gcc/aarch64-linux-gnu/13/libasan.so ASAN_OPTIONS=detect_leaks=0 pytest tests/python/ -v

# Valgrind (alternative memory checker - SLOW, ~50x overhead)
# Use for specific tests, not full suite. Focus on "definitely lost" = 0
cmake --preset linux-x64-debug
cmake --build --preset linux-x64-debug-build
pip install -e . --force-reinstall --no-deps
# Quick single-file test (recommended):
timeout 120 valgrind --leak-check=full --suppressions=python.supp python -m pytest tests/python/test_specific.py -v
# Or quick inline test:
timeout 60 valgrind --leak-check=full --suppressions=python.supp python -c "import slicot; print('ok')"

# Debug C vs Fortran
cmake --preset linux-x64-debug -DBUILD_FORTRAN_DIAG=ON
cmake --build --preset linux-x64-debug-build --target diag_all
grep -A3 "OUTPUT" build/linux-x64-debug/fortran_diag/*_trace.txt
```

**Presets:** `linux-x64-{debug,release}`, `macos-arm64-{debug,release}`, `linux-x64-debug-{asan,ubsan,sanitizers}`

## Memory Debugging Tools

**Sanitizers (ASAN/UBSAN):** Fast, catches errors at runtime
- Use for CI/CD and regular development
- Best for: buffer overflows, use-after-free, undefined behavior
- Overhead: ~2x slowdown
- **Required before PR:** Run full test suite with ASAN

**Valgrind:** Comprehensive, slower but more thorough
- Use for specific test files only (too slow for full suite)
- Best for: memory leaks, uninitialized reads, invalid frees
- Overhead: ~20-50x slowdown (use `timeout` to avoid hanging)
- Suppressions file: `python.supp` filters Python false positives
- **Focus on:** "definitely lost" bytes (must be 0)
- **Ignore:** "possibly lost" and "still reachable" (Python internals)

## Directory Structure

```
src/XX/routine.c          # C11 (XX=AB,MB,MC...)
include/slicot.h          # Public API + Doxygen
include/slicot_types.h    # Type aliases (i32, f64)
include/slicot_blas.h     # BLAS/LAPACK wrappers
python/slicot_module.c    # Python/C extension
python/slicot/__init__.py # Python exports
tests/python/test_*.py    # pytest tests
tools/extract_dependencies.py  # Dependency analyzer
```

**Naming:** `AB01MD.f` → `src/AB/ab01md.c`

## Translation Workflow

**Parallel Work**: Use `parallel-workflow-orchestrator` droid for managing concurrent development streams.
**Use agent:** `slicot-fortran-translator` and `/skill slicot-knowledge` handles full TDD workflow (RED/GREEN/REFACTOR/VERIFY)
**Check deps:** `python3 tools/extract_dependencies.py SLICOT-Reference/src/ ROUTINE_NAME`

## Critical Patterns

### Types
- `INTEGER` → `i32`, `DOUBLE PRECISION` → `f64`, `LOGICAL` → `bool`
- **Exception:** LAPACK callbacks (DGEES/DGGES SELECT) MUST use `int`, not `bool` (ABI: FORTRAN LOGICAL=4 bytes, C bool=1 byte)

### Column-Major Arrays
```c
// Index: a[i + j*lda] (row i, col j)
// Memory layout: column-by-column
double a[] = {1.0, 3.0, 2.0, 4.0};  // [[1,2], [3,4]]
```
**NumPy tests:** Always use `order='F'`

### Index Conversion (CRITICAL - Security)
SLICOT returns 1-based indices. Always convert & validate:

```c
// CORRECT
k = iwork[j] - 1;
if (k < 0 || k >= n) break;  // REQUIRED bounds check
// Now safe: iwork[k]

// WRONG - Buffer overflow risk!
k = iwork[j] - 1;
if (iwork[k] < 0) { }  // Missing bounds check before access
```

### Python Wrapper (CRITICAL - Memory)
**In-place modification:**
```c
// CORRECT: Return modified input array
PyObject *result = Py_BuildValue("Odi", b_array, scale, info);
Py_DECREF(a_array);
Py_DECREF(b_array);

// WRONG: Double-free crash!
// PyObject *u_array = PyArray_New(..., b_data, ...);
// PyArray_ENABLEFLAGS(u_array, NPY_ARRAY_OWNDATA);
```

**Input arrays:**
```c
a_array = PyArray_FROM_OTF(a_obj, NPY_DOUBLE, NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
```

**Output arrays:**
```c
npy_intp dims[2] = {m, n};
npy_intp strides[2] = {sizeof(f64), m * sizeof(f64)};
q_array = PyArray_New(&PyArray_Type, 2, dims, NPY_DOUBLE, strides, q, 0, NPY_ARRAY_FARRAY, NULL);
PyArray_ENABLEFLAGS(q_array, NPY_ARRAY_OWNDATA);
```

### BLAS/LAPACK
- Use `SLC_DGEMM()` etc. from `slicot_blas.h`
- Pass scalars by pointer: `SLC_DGEMM("N", "N", &m, &n, &k, &alpha, a, &lda, ...)`

### Error Codes
`info = 0` success, `info < 0` param error, `info > 0` algorithm error

## Test Data Strategy

**Threshold rule**: Use NPZ files for datasets with ≥50 values or >10 lines of data.

| Scenario | Strategy | Example |
|----------|----------|---------|
| Small data (<50 values) | Embed inline | 3x3 matrix, short vector |
| Large data (≥50 values) | NPZ file in `tests/python/data/` | 1000-sample time series |
| Shared between tests | ALWAYS use NPZ | Same I/O data for IB01AD/IB01BD |

**NPZ file pattern:**
```python
# Creating test data file (one-time)
np.savez('tests/python/data/routine_test_data.npz', u=u, y=y, expected_a=a)

# Loading in test
def load_test_data():
    data_path = os.path.join(os.path.dirname(__file__), 'data', 'routine_test_data.npz')
    data = np.load(data_path)
    return data['u'], data['y'], data['expected_a']
```

**Why**: Keeps test files readable (<400 lines), prevents data duplication, enables data sharing between related tests.

## Docs

- `fortran_diag/README.md` - C vs Fortran debugging
- `tools/README.md` - Workflow examples
- `.claude/skills/slicot-knowledge/SKILL.md` - Translation knowledge
