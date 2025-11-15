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

**NEW (Automated):** Use the generator script

```bash
# Generate skeleton files and update CMake
./tools/generate_diagnostic.sh AB01MD
```

This creates:
- `fortran_diag/fortran/ab01md_diag.f` (from template)
- `fortran_diag/c/ab01md_diag.c` (from template)
- Updates `fortran_diag/CMakeLists.txt` automatically

### Fill in the Parsing Logic

The generated skeletons have `TODO` comments marking what needs to be filled in:

**Reference:** Check `SLICOT-Reference/examples/TAB01MD.f` for READ statements

**Fortran (`fortran_diag/fortran/ab01md_diag.f`):**
1. Update parameter list (N, M, etc.)
2. Add READ statements matching data file format
3. Add matrix/vector declarations
4. Update CALL statement with routine signature
5. Print output matrices

**C (`fortran_diag/c/ab01md_diag.c`):**
1. Update variable declarations
2. Add scanf() for parameters
3. Add matrix reading loops (remember column-major storage!)
4. Update function call with routine signature
5. Print output matrices

**Validation helpers included:**
- `frobenius_norm()` - Frobenius norm (confirms matrix magnitude)
- `checksum()` - Sum of all elements (quick sanity check)
- `max_abs()` - Maximum absolute value

### Build and Run

```bash
# Reconfigure (if needed)
cmake --preset macos-arm64-debug -DBUILD_FORTRAN_DIAG=ON

# Run diagnostic
cmake --build --preset macos-arm64-debug-build --target ab01md_diag_all

# View results
cat build/macos-arm64-debug/fortran_diag/ab01md_fortran.txt
cat build/macos-arm64-debug/fortran_diag/ab01md_c.txt
cat build/macos-arm64-debug/fortran_diag/output/ab01md_diff.txt
```

### Manual Approach (Old)

<details>
<summary>Click to expand manual steps (not recommended)</summary>

Copy and modify existing files manually, update CMakeLists.txt by hand.
See git history for original instructions. The automated approach above is faster and less error-prone.

</details>

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

## SG03BD Bug Confirmed (Fixed)

**Historical bug (now fixed):**
- Fortran (correct): `U = [1.600, 0.680, 0.204]` (diagonal)
- C (wrong): `U = [0.437, 0.581, 0.080]` (diagonal)
- Error: 60-73% off

Bug was isolated to Lyapunov solver (SG03BV) using this diagnostic framework.

**New validation features:**
- Input/output validation with Frobenius norms, checksums, max abs values
- Helps confirm parsing correctness before comparing results
- Example output:
  ```
  === INPUT VALIDATION ===
  A Frobenius norm: 1.2345678901234567E+01
  A checksum:       3.4567890123456789E+00
  ```

## References

- C impl: `src/SG/sg03bd.c`
- Fortran ref: `SLICOT-Reference/src/SG03BD.f`
- Example: `SLICOT-Reference/examples/TSG03BD.f`
- Docs: `SLICOT-Reference/doc/SG03BD.html`
