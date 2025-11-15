# IB03BD Translation Order

**Issue**: #17
**Strategy**: Bottom-up (Level 0 → Level 6)

## Phase 1: Level 0 Utilities (4 routines)

| Routine | LOC | Purpose | Files |
|---------|-----|---------|-------|
| MA02AD | ~50 | Matrix copy | `src/MA/ma02ad.c` |
| TB01WD | ~200 | Similarity transformation | `src/TB/tb01wd.c` |
| TB01VY | ~300 | Inverse normal form | `src/TB/tb01vy.c` |
| TF01MX | ~440 | State-space simulation | `src/TF/tf01mx.c` |

**Dependencies**: BLAS/LAPACK only
**Branch**: `feat/phase-1-level-0-utils`
**PR Checklist**:
- [ ] C implementations + tests
- [ ] Python wrappers
- [ ] Tests pass with sanitizers

---

## Phase 2: MD03BD Core (1 routine + deps)

| Routine | LOC | Purpose | Dependencies |
|---------|-----|---------|--------------|
| MD03BD | 1191 | Levenberg-Marquardt optimizer | Function pointers |

**Function Pointer Support**: FCN, QRFACT, LMPARM callbacks
**Branch**: `feat/phase-2-md03bd-core`
**PR Checklist**:
- [ ] Function pointer typedefs in `slicot_types.h`
- [ ] MD03BD implementation
- [ ] Template callbacks (NF01BE, NF01BF, NF01BS, NF01BP)
- [ ] Tests with predefined callbacks
- [ ] Python wrapper (string-based callback selection)

---

## Phase 3: NF01 Support (6 routines)

| Routine | LOC | Purpose |
|---------|-----|---------|
| NF01BE | ~150 | FCN: error function |
| NF01BF | ~150 | FCN: Wiener system error |
| NF01BP | ~200 | LMPARM: parameter computation |
| NF01BS | ~200 | QRFACT: QR factorization |
| MD03BA | ~50 | QRFACT init support |
| MD03BB | ~50 | LMPARM init support |

**Dependencies**: MD03BD (Phase 2)
**Branch**: `feat/phase-3-nf01-support`

---

## Phase 4: Level 1-2 Matrix Algorithms

### IB Family (3 routines)
| Routine | Level | Purpose |
|---------|-------|---------|
| IB01OD | 1 | Hankel matrix processing |
| IB01ND | 2 | Covariance/state sequence |
| IB01QD | 2 | System matrices application |
| IB01RD | 2 | Initial state estimation |

### SB Family (2 routines)
| Routine | Level | Purpose |
|---------|-------|---------|
| SB02MT | 2 | Riccati preprocessing |
| SB02ND | 2 | Newton's Riccati method |

**Branch**: `feat/phase-4-level-1-2-algos`
**Dependencies**: Level 0 utilities

---

## Phase 5: State Estimation (complete from Phase 4)
_Merged into Phase 4_

---

## Phase 6: SB03OD Lyapunov Solver

| Routine | LOC | Purpose | Level |
|---------|-----|---------|-------|
| SB03OD | ~600 | Lyapunov equation solver | 3 |

**Branch**: `feat/phase-6-sb03od`
**Dependencies**: Check SB03OD deps (may need additional SB/MB routines)

---

## Phase 7: Large Components

| Routine | LOC | Level | Purpose |
|---------|-----|-------|---------|
| IB01MD | 1459 | 3 | R factor computation |
| IB01PD | ~800 | 3 | System matrices from R |

**Branch**: `feat/phase-7-large-components`
**Dependencies**: IB01MY, MB04OD, IB01PX, IB01PY, MB02QY (check recursively)
**Note**: May require splitting into sub-phases

---

## Phase 8: SB02RD Riccati Solver

| Routine | LOC | Level |
|---------|-----|-------|
| SB02RD | ~1200 | 4 |

**Branch**: `feat/phase-8-sb02rd`
**Dependencies**: Full Riccati dependency chain (complex)
**Note**: High complexity, extensive testing required

---

## Phase 9: Integration Routines

| Routine | LOC | Level | Purpose |
|---------|-----|-------|---------|
| TB01VD | 488 | 4 | Output normal form conversion |
| IB01CD | 808 | 3 | Initial state computation |

**Branch**: `feat/phase-9-integration`
**Dependencies**: MA02AD, SB03OD (Phase 6), IB01QD, IB01RD, TB01WD

---

## Phase 10: IB01AD (MOESP)

| Routine | LOC | Level |
|---------|-----|-------|
| IB01AD | 698 | 4 |

**Branch**: `feat/phase-10-ib01ad`
**Dependencies**: IB01MD, IB01ND, IB01OD (Phases 4, 7)
**Test**: Standalone MOESP algorithm validation

---

## Phase 11: IB01BD (N4SID)

| Routine | LOC | Level |
|---------|-----|-------|
| IB01BD | 776 | 5 |

**Branch**: `feat/phase-11-ib01bd`
**Dependencies**: IB01PD, MA02AD, SB02MT, SB02ND, SB02RD (Phases 1, 4, 8)
**Test**: Standalone N4SID algorithm validation

---

## Phase 12: IB03BD Main Entry Point

| Routine | LOC | Level |
|---------|-----|-------|
| IB03BD | 1072 | 6 |

**Branch**: `feat/phase-12-ib03bd-main`
**Dependencies**: All previous phases
**PR Checklist**:
- [ ] Main entry point implementation
- [ ] Complex workspace management (LDWORK calculation)
- [ ] 4 initialization modes (INIT='L','S','B','N')
- [ ] Python wrapper
- [ ] Test with `IB03BD.dat` → validate `IB03BD.res`
  - Expected: IWARN=12, residual=0.2995840
  - Iterations=42, fcn_evals=898, jac_evals=295
- [ ] Full integration test suite
- [ ] Update main README

---

## Workflow Per Phase

1. **Branch**: `git checkout -b feat/phase-N-description`
2. **Code**: TDD with slicot-fortran-translator agent
3. **Test**: `pytest tests/python/ -v`
4. **Sanitize**: `cmake --preset macos-arm64-debug-sanitizers && cmake --build --preset macos-arm64-debug-sanitizers-build && pytest tests/python/ -v`
5. **PR**: Create with phase checklist
6. **Review**: Address comments
7. **Merge**: `gh pr merge`
8. **Sync**: `git checkout main && git pull`
9. **Next**: Phase N+1

---

## Critical Notes

1. **Dependencies**: Run `extract_dependencies.py` before each phase to verify
2. **Column-major**: All arrays use Fortran order (`order='F'`)
3. **Index conversion**: 1-based → 0-based with bounds checks
4. **Function pointers**: Established in Phase 2, used in Phases 3, 12
5. **Workspace**: LDWORK calculations critical (IB03BD lines 236-287)
6. **Error codes**: INFO=0 success, INFO<0 param error, INFO>0 algorithm error

---

## Open Questions

1. **SB03OD dependencies**: Need full recursive check
2. **IB01MD subtree**: May require additional MB/IB routines
3. **SB02RD complexity**: Consider sub-phases?
4. **Phase splitting**: Phases 7-8 may overflow context, monitor closely

---

## References

- GitHub Issue: #17
- Function Pointers: `docs/function_pointers.md`
- Full Dependencies: `docs/ib03bd_full_dependencies.txt`
- Dependency Chain: `docs/ib03bd_dependency_chain.txt`
