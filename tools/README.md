# SLICOT Dependency Analysis

## extract_dependencies.py

Analyzes Fortran 77 dependencies to determine translation order.

### Usage

```bash
# Overall summary
python3 tools/extract_dependencies.py SLICOT-Reference/src/

# Specific routine
python3 tools/extract_dependencies.py SLICOT-Reference/src/ MB01QD
```

### Output

- **Level 0**: Leaf routines (no SLICOT deps) - translate first
- **Level 1+**: Depend on lower levels - translate after dependencies
- **LAPACK/BLAS deps**: External calls needed

### Example

```bash
$ python3 tools/extract_dependencies.py SLICOT-Reference/src/ MB01QD

ROUTINE: MB01QD
Dependency Level: 0

No SLICOT dependencies (LEAF ROUTINE)

LAPACK/BLAS Dependencies (2):
  LSAME
  DLAMCH

Used By (2 routines):
  MB01PD
  MB03DD
```

**Interpretation:**
- Level 0 = translate immediately
- Needs LSAME, DLAMCH (inline or use from LAPACK)
- Unblocks MB01PD and MB03DD

