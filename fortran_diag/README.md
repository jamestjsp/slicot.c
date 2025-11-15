# Fortran Diagnostic Tools

Side-by-side comparison of SLICOT Fortran reference vs C implementation to debug divergence.

## Quick Start

```bash
# Enable diagnostic tools
cmake --preset macos-arm64-debug -DBUILD_FORTRAN_DIAG=ON

# Build and run full workflow
cmake --build --preset macos-arm64-debug-build --target diag_all

# View results
grep -A3 "OUTPUT U" build/macos-arm64-debug/fortran_diag/*_trace.txt
```

Output: `build/macos-arm64-debug/fortran_diag/`
- `fortran_trace.txt` - Fortran execution
- `c_trace.txt` - C execution
- `output/diff.txt` - Comparison report

## Available Targets

```bash
cmake --build --preset macos-arm64-debug-build --target <target>
```

| Target | Description |
|--------|-------------|
| `sg03bd_diag_c` | Build C diagnostic only |
| `sg03bd_diag_fortran` | Build Fortran diagnostic only |
| `run_fortran_diag` | Run Fortran |
| `run_c_diag` | Run C |
| `compare_diag` | Compare outputs |
| `diag_all` | Full workflow (recommended) |

## Adding a New Routine

To debug a different routine (e.g., `AB01MD`):

### 1. Create Fortran Diagnostic

Copy and modify `fortran/sg03bd_diag.f`:

```bash
cd fortran_diag/fortran
cp sg03bd_diag.f ab01md_diag.f
```

Edit `ab01md_diag.f`:
- Change routine call: `CALL SG03BD(...)` → `CALL AB01MD(...)`
- Update parameters to match AB01MD signature
- Adjust matrix dimensions and I/O

### 2. Create C Diagnostic

Copy and modify `c/sg03bd_diag.c`:

```bash
cd fortran_diag/c
cp sg03bd_diag.c ab01md_diag.c
```

Edit `ab01md_diag.c`:
- Change function call: `sg03bd(...)` → `ab01md(...)`
- Update parameters to match AB01MD signature
- Adjust matrix dimensions and I/O

### 3. Update CMakeLists.txt

Add targets in `fortran_diag/CMakeLists.txt`:

```cmake
# Test data paths
set(AB01MD_TEST_DATA "${SLICOT_FORTRAN_DIR}/examples/data/AB01MD.dat")

# C diagnostic
add_executable(ab01md_diag_c c/ab01md_diag.c)
target_link_libraries(ab01md_diag_c PRIVATE slicot ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} m)
target_include_directories(ab01md_diag_c PRIVATE ${PROJECT_SOURCE_DIR}/include ${PROJECT_BINARY_DIR}/include)

# Fortran diagnostic
add_executable(ab01md_diag_fortran fortran/ab01md_diag.f)
target_link_libraries(ab01md_diag_fortran PRIVATE slicot_fortran ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

# Run targets
add_custom_target(run_ab01md_fortran
    COMMAND ab01md_diag_fortran < ${AB01MD_TEST_DATA} > ab01md_fortran.txt 2>&1 || true
    DEPENDS ab01md_diag_fortran
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_custom_target(run_ab01md_c
    COMMAND ab01md_diag_c < ${AB01MD_TEST_DATA} > ab01md_c.txt 2>&1 || true
    DEPENDS ab01md_diag_c
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_custom_target(compare_ab01md
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/scripts/compare.py ab01md_fortran.txt ab01md_c.txt || true
    DEPENDS run_ab01md_fortran run_ab01md_c
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_custom_target(ab01md_diag_all DEPENDS compare_ab01md)
```

### 4. Build and Run

```bash
# Reconfigure
cmake --preset macos-arm64-debug -DBUILD_FORTRAN_DIAG=ON

# Run new diagnostic
cmake --build --preset macos-arm64-debug-build --target ab01md_diag_all

# View traces
cat build/macos-arm64-debug/fortran_diag/ab01md_fortran.txt
cat build/macos-arm64-debug/fortran_diag/ab01md_c.txt
```

## Test Data Sources

SLICOT reference examples:
- Data: `SLICOT-Reference/examples/data/<ROUTINE>.dat`
- Expected: `SLICOT-Reference/examples/results/<ROUTINE>.res`
- Example program: `SLICOT-Reference/examples/T<ROUTINE>.f`

## Output Format

High precision (16 digits, E23.16 format):
- Parameters (N, M, etc.)
- Input matrices
- Output matrices
- Scalars (INFO, SCALE, etc.)
- Eigenvalues/vectors

## SG03BD Bug Confirmed

Fortran (correct): `U = [1.600, 0.680, 0.204]` (diagonal)
C (wrong): `U = [0.437, 0.581, 0.080]` (diagonal)
Error: 60-73% off

Bug isolated to Lyapunov solver (SG03BV) - all other outputs match.

## References

- C impl: `src/SG/sg03bd.c`
- Fortran ref: `SLICOT-Reference/src/SG03BD.f`
- Example: `SLICOT-Reference/examples/TSG03BD.f`
- Docs: `SLICOT-Reference/doc/SG03BD.html`
