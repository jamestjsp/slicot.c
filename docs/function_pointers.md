# Function Pointer Patterns in SLICOT C Translation

## Overview
MD03BD and related routines require function pointer callbacks for:
- **FCN**: Error function and Jacobian evaluation
- **QRFACT**: QR factorization with pivoting
- **LMPARM**: Levenberg-Marquardt parameter computation

## C Pattern

### Function Pointer Type Definitions
```c
// FCN callback type
typedef void (*slicot_fcn_t)(
    i32 *iflag, i32 *m, i32 *n,
    i32 *ipar, i32 *lipar,
    f64 *dpar1, i32 *ldpar1,
    f64 *dpar2, i32 *ldpar2,
    f64 *x, i32 *nfevl,
    f64 *e, f64 *j, i32 *ldj,
    f64 *dwork, i32 *ldwork,
    i32 *info
);

// QRFACT callback type
typedef void (*slicot_qrfact_t)(
    i32 *n, i32 *ipar, i32 *lipar,
    f64 *fnorm, f64 *j, i32 *ldj,
    f64 *e, f64 *jnorms, f64 *gnorm,
    i32 *ipvt, f64 *dwork, i32 *ldwork,
    i32 *info
);

// LMPARM callback type
typedef void (*slicot_lmparm_t)(
    const char *cond, i32 *n,
    i32 *ipar, i32 *lipar,
    f64 *r, i32 *ldr, i32 *ipvt,
    f64 *diag, f64 *qtb, f64 *delta,
    f64 *par, i32 *ranks, f64 *x, f64 *rx,
    f64 *tol, f64 *dwork, i32 *ldwork,
    i32 *info
);
```

### MD03BD Function Signature
```c
void md03bd(
    const char *xinit, const char *scale, const char *cond,
    slicot_fcn_t fcn,        // Function pointer
    slicot_qrfact_t qrfact,  // Function pointer
    slicot_lmparm_t lmparm,  // Function pointer
    i32 *m, i32 *n, i32 *itmax,
    f64 *factor, i32 *nprint,
    i32 *ipar, i32 *lipar,
    f64 *dpar1, i32 *ldpar1,
    f64 *dpar2, i32 *ldpar2,
    f64 *x, f64 *diag,
    i32 *nfev, i32 *njev,
    f64 *ftol, f64 *xtol, f64 *gtol, f64 *tol,
    i32 *iwork, f64 *dwork, i32 *ldwork,
    i32 *iwarn, i32 *info
);
```

### Usage Example
```c
// Implementation of FCN callback
void my_fcn(i32 *iflag, i32 *m, i32 *n, ...) {
    if (*iflag == 1) {
        // Compute error functions
    } else if (*iflag == 2) {
        // Compute Jacobian
    }
}

// Call MD03BD
slicot_fcn_t fcn = &my_fcn;
slicot_qrfact_t qrfact = &nf01bs;
slicot_lmparm_t lmparm = &nf01bp;

md03bd(
    "G", "I", "E",
    fcn, qrfact, lmparm,
    &m, &n, &itmax,
    ...
);
```

## Python Wrapper Strategy

### Approach 1: Predefined Callbacks
Provide standard implementations (NF01BE, NF01BF, NF01BP, NF01BS) as named options:
```python
# User specifies callback by name
slicot.md03bd(..., fcn='nf01be', qrfact='nf01bs', lmparm='nf01bp')
```

### Approach 2: Python Callable (Advanced)
Support Python callables via ctypes callbacks:
```python
def my_fcn(iflag, m, n, ipar, ...):
    if iflag == 1:
        return error_vector
    elif iflag == 2:
        return jacobian_matrix

slicot.md03bd(..., fcn=my_fcn, qrfact='nf01bs', lmparm='nf01bp')
```

**Recommended**: Start with Approach 1 (predefined), add Approach 2 in later phase if needed.

## Implementation Notes

1. **Type safety**: Use typedefs in `slicot_types.h`
2. **Null checks**: Validate function pointers before calling
3. **IFLAG handling**: FCN must handle all IFLAG values (0,1,2,3)
4. **Error propagation**: Check INFO from callbacks
5. **Python wrapping**: Map string names to C function pointers

## References
- MD03BD: `SLICOT-Reference/src/MD03BD.f` lines 52-315
- NF01BE, NF01BF: Example FCN implementations
- NF01BS: Example QRFACT implementation
- NF01BP: Example LMPARM implementation
