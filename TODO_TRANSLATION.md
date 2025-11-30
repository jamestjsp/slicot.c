# SLICOT Translation TODO

Target routines from user list, excluding already-translated ones.

**Status:** 75/627 routines translated (12%)

---

## Already Translated (Skip)

- IB01AD, IB01BD, IB01CD, IB01RD, IB03BD
- SB02MT, SB03OD, SG03BD
- TB01VD, TG01FD
- MA02AD, MA02ED, MA02FD, MB03OY, MB03UD, MB01SD (+ 60 more helpers)

---

## Tier 1: Leaf Routines (Level 0) - No SLICOT Dependencies

Can be translated in parallel. Only depend on BLAS/LAPACK.

| Routine | Purpose | LAPACK Deps |
|---------|---------|-------------|
| **AB04MD** | Bilinear transform | DGEMM, DGETRF, DGETRI, DGETRS, DLASCL, DSCAL, DTRSM |
| **AB05MD** | Series connection | DGEMM, DLACPY, DLASET |
| **AB05ND** | Feedback connection | DCOPY, DGEMM, DGEMV, DGETRF, DGETRS, DLACPY, DLASET |
| **AB05OD** | State-space parallel | DLACPY, DLASET |
| **AB05PD** | Feedback matrices | DAXPY, DLACPY, DLASET |
| **AB05QD** | Rowwise concatenation | DLACPY, DLASET |
| **AB07ND** | Dual system | DCOPY, DGECON, DGEMM, DGEMV, DGETRF, DGETRI, DLACPY |
| **AB13MD** | Hankel norm bound | DCOPY, DGEMV, DLACPY, DLASET, DSCAL, DSYCON, DSYSV, DSYTRF, DSYTRS |
| **DK01MD** | Data smoothing | (pure Fortran) |
| **MB03VY** | Accumulate transformations | DLASET, DORGHR, DORGQR |
| **MB05ND** | Matrix exponential integral | DAXPY, DCOPY, DGEMM, DGEMV, DGESV, DLACPY, DLASET, DSCAL |
| **MC01TD** | Polynomial GCD | DCOPY, DRSCL |
| **SB10JD** | Hankel singular values | DGEMM, DGESVD, DLACPY, DLASET, DSCAL |
| **TB01ID** | Matrix scaling | DSCAL |
| **TC01OD** | Polynomial reverse | DCOPY, DSWAP |
| **TF01MD** | Transfer response | DCOPY, DGEMM, DGEMV, DLASET |
| **TF01RD** | Transfer matrix | DGEMM, DLACPY |
| **TG01AD** | Descriptor balance | DAXPY, DCOPY, DSCAL |

---

## Tier 2: Level 1 Dependencies

| Routine | Level | SLICOT Deps | Purpose |
|---------|-------|-------------|---------|
| **AB08NZ** | 1 | TB01IZ | Complex normal rank |
| **DE01OD** | 1 | DG01MD | Convolution via FFT |
| **DG01ND** | 1 | DG01MD | Inverse FFT |
| **MB03VD** | 1 | MB04PY | Hessenberg decomp |
| **MB03WD** | 1 | MB04PY | Skew-Hamiltonian Schur |
| **MB05MD** | 1 | MB05MY | Matrix exponential |
| **SB02MD** | 1 | SB02MU | Algebraic Riccati (discrete) |
| **SB02OD** | 1 | SB02OY | Algebraic Riccati (continuous) |
| **SB04ND** | 1 | SB04NV, SB04NW, SB04NX, SB04NY | Sylvester equation |
| **TB05AD** | 1 | MB02RZ, MB02SZ, MB02TZ | Complex transfer matrix |
| **TC04AD** | 1 | AB07MD, TC01OD | Polynomial to state-space |

---

## Tier 3: Level 2 Dependencies

| Routine | Level | SLICOT Deps | Purpose |
|---------|-------|-------------|---------|
| **AB01MD** | 2 | MB01PD | Staircase form |
| **AB01ND** | 2 | MB01PD, MB03OY | Controllability staircase |
| **AB08MD** | 2 | AB08NX, TB01ID | Normal rank |
| **AB08ND** | 2 | AB08NX, TB01ID | Kronecker indices |
| **AG08BD** | 2 | AG08BY, MA02BD, MA02CD, TB01XD, TG01AD, TG01FD | Descriptor normal rank |
| **MB02ED** | 2 | MB02CX, MB02CY | Toeplitz factorization |
| **MB03RD** | 2 | MA02AD, MB03QX, MB03RX, MB03RY | Real Schur form |
| **SB01BD** | 2 | MB03QD, MB03QY, SB01BX, SB01BY | Pole assignment |
| **SB03MD** | 2 | MB01RD, SB03MX, SB03MY | Lyapunov equation |
| **SB04MD** | 2 | SB04MU, SB04MY | Sylvester (Hessenberg) |
| **SB04QD** | 2 | SB04QU, SB04QY | Sylvester (Schur) |
| **SG02AD** | 2 | MB01SD, MB02PD, MB02VD, SB02OY | Generalized Riccati |
| **SG03AD** | 2 | MB01RD, MB01RW, SG03AX, SG03AY | Generalized Lyapunov |

---

## Tier 4: Level 3 Dependencies

| Routine | Level | SLICOT Deps | Purpose |
|---------|-------|-------------|---------|
| **AB01OD** | 3 | AB01ND | Controllable part |
| **TB01PD** | 3 | AB07MD, TB01ID, TB01UD, TB01XD | Minimal realization |
| **TB03AD** | 3 | AB07MD, MA02GD, TB01ID, TB01UD, TB01YD, TB03AY, TC01OD | Transfer to polynomial |

---

## Tier 5: Level 4 Dependencies

| Routine | Level | SLICOT Deps | Purpose |
|---------|-------|-------------|---------|
| **AB09AX** | 4 | MA02AD, MA02DD, MB03UD, SB03OU | Balanced truncation helper |
| **AB13BD** | 4 | SB03OU, SB08DD | H-infinity norm |
| **AB13DD** | 4 | MA02AD, MB01SD, MB03XD, TB01ID, TG01BD | Hankel norm (descriptor) |
| **SB10DD** | 4 | MA02AD, MB01RU, MB01RX, SB02OD, SB02SD | H2 control |
| **TB04AD** | 4 | AB07MD, TB01XD, TB04AY | State-space to transfer |
| **TD04AD** | 4 | TB01PD, TB01XD, TD03AY | Polynomial to state-space |
| **AB13ED** | 1 | MA02ED, MB04ZD | Passivity index |
| **AB13FD** | 1 | MA02ED, MB04ZD | L-infinity norm |

---

## Tier 6: Level 5+ Dependencies (Complex)

| Routine | Level | SLICOT Deps | Purpose |
|---------|-------|-------------|---------|
| **AB09AD** | 5 | AB09AX, TB01ID, TB01WD | Balanced truncation |
| **AB09BD** | 5 | AB09BX, TB01ID, TB01WD | Singular perturbation |
| **AB09MD** | 5 | AB09AX, TB01ID, TB01KD | Balanced truncation (discrete) |
| **AB09ND** | 5 | AB09BX, TB01ID, TB01KD | Singular perturbation (discrete) |
| **SB10AD** | 6 | SB10LD, SB10PD, SB10QD, SB10RD | H-infinity synthesis |
| **SB10FD** | 6 | SB10PD, SB10QD, SB10RD | Suboptimal H-infinity |
| **SB10HD** | 6 | SB10UD, SB10VD, SB10WD | Loop-shaping H-infinity |
| **SB10YD** | 6 | AB04MD, DG01MD, SB10ZP | Frequency weighting |
| **IB03AD** | 6 | IB01AD✓, IB01BD✓, IB01CD✓, MD03AD, TB01VD✓, TB01VY✓, TF01MX✓ | Wiener system (data) |

---

## Common Dependencies (Shared Infrastructure)

These are needed by multiple target routines - prioritize early:

| Routine | Used By | Status |
|---------|---------|--------|
| **MB03OY** | AB01ND, AB08NX, AB08NY | ✅ Done |
| **MA02AD** | AB09AX, AB13DD, MB03RD, SB10DD | ✅ Done |
| **TB01ID** | AB08MD, AB08ND, AB09AD, AB09BD, AB09MD, AB09ND, AB13DD, TB01PD, TB03AD | ❌ TODO |
| **MB01PD** | AB01MD, AB01ND | ❌ TODO |
| **AB07MD** | TB01PD, TB03AD, TB04AD, TC04AD | ❌ TODO |
| **TB01XD** | AG08BD, TB01PD, TB04AD, TD04AD | ❌ TODO |
| **SB02OY** | SB02OD, SG02AD | ❌ TODO |
| **MB04ZD** | AB13ED, AB13FD | ❌ TODO |
| **DG01MD** | DE01OD, DG01ND, SB10YD | ❌ TODO |

---

## Suggested Translation Order

### Phase A: Foundation Leaves (17 routines)
```
TB01ID, AB04MD, AB05MD, AB05ND, AB05OD, AB05PD, AB05QD,
AB07ND, AB13MD, DK01MD, MC01TD, MB05ND, SB10JD, TC01OD,
TF01MD, TF01RD, TG01AD
```

### Phase B: Shared Infrastructure (9 routines)
```
MB01PD, AB07MD, TB01XD, SB02OY, MB04ZD, DG01MD,
MB04PY, SB02MU, MB05MY
```

### Phase C: Level 1-2 Drivers (15 routines)
```
AB01MD, AB01ND, DE01OD, DG01ND, MB03VD, MB03WD, MB05MD,
SB02MD, SB02OD, SB04ND, TB05AD, TC04AD, AB08MD, AB08ND
```

### Phase D: Level 2-3 Solvers (8 routines)
```
SB01BD, SB03MD, SB04MD, SB04QD, SG02AD, SG03AD, AB01OD, TB01PD
```

### Phase E: Level 3-4 Transformations (6 routines)
```
TB03AD, TB04AD, TD04AD, AB13BD, AB13DD, SB10DD
```

### Phase F: Model Reduction (4 routines)
```
AB09AD, AB09BD, AB09MD, AB09ND
```

### Phase G: H-infinity Synthesis (5 routines)
```
SB10AD, SB10FD, SB10HD, SB10YD, IB03AD
```

### Phase H: Analysis (3 routines)
```
AB13ED, AB13FD, AG08BD
```

---

## Notes

- Routines marked ✅ are already translated
- LAPACK/BLAS routines wrapped via `slicot_blas.h`
- Use `tools/extract_dependencies.py` for full trees
- Common deps can be grouped across phases
