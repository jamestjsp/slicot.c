# Test Suite Review: Phase 1 Level 0 Routines

**Routines Reviewed**: MA02AD, TB01WD, TB01VY, TF01MX
**Date**: 2025-11-15
**Status**: Phase 1 (feat/phase-1-level-0-utils)

## Summary

| Routine | Test Source | Quality | Coverage | Action Required |
|---------|------------|---------|----------|-----------------|
| MA02AD  | **Synthetic** | Good | Excellent | Enhance with property tests |
| TB01WD  | **SLICOT Example** (HTML doc) | Excellent | Good | Add error/edge cases |
| TB01VY  | **Synthetic** | Fair | Basic | Generate with control package |
| TF01MX  | **Synthetic** | Good | Good | Generate with control package |

---

## 1. MA02AD - Matrix Transposition

### Test Data Source
**100% Synthetic** - No SLICOT examples available
- HTML doc: "None" for Program Text/Data/Results
- No example files in `SLICOT-Reference/examples/`
- No benchmark data

### Current Test Quality: **Good**
**File**: `tests/python/test_ma02ad.py`

**Coverage Analysis**:
- ✓ Full transpose (JOB='F')
- ✓ Upper triangular (JOB='U')
- ✓ Lower triangular (JOB='L')
- ✓ Rectangular matrices (both orientations)
- ✓ Upper trapezoid (M > N)
- ✓ Edge cases: zero rows, zero columns, single element
- ✓ Fortran order (`order='F'`)
- ✓ Tight tolerance (rtol=1e-14)

**Strengths**:
1. Comprehensive parameter coverage (all JOB modes)
2. Multiple matrix shapes (square, rectangular, trapezoidal)
3. Excellent edge case handling
4. Correct Fortran ordering

**Gaps**:
1. No error handling tests (invalid JOB parameter)
2. No large matrix stress tests
3. Missing property-based tests (e.g., (A^T)^T = A)

### Recommendations

**Priority: Medium** (already good coverage)

1. **Add property-based tests**:
```python
def test_ma02ad_transpose_inverse():
    """Test (A^T)^T = A"""
    a = np.random.randn(5, 3).astype(float, order='F')
    b = ma02ad('F', a)
    c = ma02ad('F', b)
    np.testing.assert_allclose(c, a, rtol=1e-14)

def test_ma02ad_triangular_symmetry():
    """Test symmetric matrix U and L transposes"""
    a = np.random.randn(4, 4).astype(float, order='F')
    a = a + a.T  # Make symmetric
    u = ma02ad('U', a)
    l = ma02ad('L', a)
    # Upper transpose should match lower
    np.testing.assert_allclose(u, l, rtol=1e-14)
```

2. **Add error handling**:
```python
def test_ma02ad_invalid_job():
    """Test invalid JOB parameter"""
    a = np.zeros((2, 2), order='F')
    with pytest.raises(ValueError):
        ma02ad('X', a)  # Invalid JOB
```

3. **Use Python control package** (optional):
Not applicable - MA02AD is pure linear algebra, control package won't help.

**Verdict**: Keep synthetic tests. Already excellent coverage. Add properties/error tests only.

---

## 2. TB01WD - Schur Form Transformation

### Test Data Source
**SLICOT Example (HTML doc)** - Authoritative source
- HTML doc: Complete example with 5×5 system
- Example file: `SLICOT-Reference/examples/TTB01WD.f`
- Data file: `SLICOT-Reference/examples/data/TB01WD.dat`

### Current Test Quality: **Excellent**
**File**: `tests/python/test_tb01wd.py`

**Data Extraction Correctness**:
✓ Correctly parsed from HTML doc (lines 26-90)
✓ Proper READ statement interpretation:
  - `READ ((A(I,J), J=1,N), I=1,N)` → Row-wise reading
  - `READ ((B(I,J), J=1,M), I=1,N)` → Row-wise reading
  - `READ ((C(I,J), J=1,N), I=1,P)` → Row-wise reading
✓ Fortran order (`order='F'`)
✓ Correct tolerance for HTML precision (rtol=1e-3, atol=1e-4)

**Coverage Analysis**:
- ✓ Main functionality (HTML example, N=5, M=2, P=3)
- ✓ Eigenvalue validation (WR, WI)
- ✓ Matrix transformations (A, B, C, U)
- ✓ Zero dimension quick return (N=0)
- ✓ Parameter validation (N<0, M<0, P<0)
- ✗ Missing: Multiple test cases with different system sizes
- ✗ Missing: Edge case with N>0 but M=0 or P=0
- ✗ Missing: Numerical conditioning tests

**Strengths**:
1. Uses authoritative SLICOT example data
2. Proper HTML tolerance handling
3. Comprehensive output validation
4. Good parameter error testing

**Gaps**:
1. Only one functional test (single system size)
2. No tests with M=0 or P=0 (valid cases)
3. No ill-conditioned matrix tests
4. Could leverage existing `.dat` file more explicitly

### Recommendations

**Priority: Low** (already excellent, just missing variety)

1. **Extract data from TTB01WD.f/.dat directly**:
```python
def test_tb01wd_from_dat_file():
    """Test using data directly from TB01WD.dat"""
    # Read SLICOT-Reference/examples/data/TB01WD.dat
    # Parse first line as dimensions
    # Parse subsequent matrices
    # Validate against known results from TTB01WD.f output
```

2. **Add dimension edge cases**:
```python
def test_tb01wd_no_inputs():
    """Test system with M=0 (no inputs)"""
    n, m, p = 3, 0, 2
    a = np.random.randn(n, n).astype(float, order='F')
    b = np.zeros((n, max(1, m)), order='F')
    c = np.random.randn(p, n).astype(float, order='F')
    # Should compute Schur form even with M=0

def test_tb01wd_no_outputs():
    """Test system with P=0 (no outputs)"""
    n, m, p = 3, 2, 0
    # Similar test
```

3. **Use Python control package** for additional test generation:
```python
import control

def test_tb01wd_random_stable_system():
    """Generate stable system using control package"""
    # Create random stable state-space system
    sys = control.rss(4, 2, 3)  # 4 states, 2 inputs, 3 outputs
    a, b, c, d = control.ssdata(sys)

    # Ensure stable (eigenvalues have negative real parts)
    eigs = np.linalg.eigvals(a)
    assert np.all(np.real(eigs) < 0)

    # Run TB01WD
    a_out, b_out, c_out, u, wr, wi, info = slicot.tb01wd(
        4, 2, 3, a.astype(float, order='F'),
        b.astype(float, order='F'),
        c.astype(float, order='F')
    )

    # Verify Schur form properties
    # Verify eigenvalues preserved
    eigs_schur = wr + 1j*wi
    np.testing.assert_allclose(sorted(eigs.real), sorted(wr), rtol=1e-10)
```

**Verdict**: Tests are excellent. Minor enhancements only for broader coverage.

---

## 3. TB01VY - Output Normal Form to State-Space

### Test Data Source
**100% Synthetic** - No SLICOT examples
- HTML doc: "None" for Program Text/Data/Results
- No example files in `SLICOT-Reference/examples/`
- No benchmark data

### Current Test Quality: **Fair**
**File**: `tests/python/test_tb01vy.py`

**Coverage Analysis**:
- ✓ Basic APPLY='N' test
- ✓ Basic APPLY='A' test
- ✓ Zero dimensions (N=0, M=0, L=0)
- ✓ Error handling (LTHETA, invalid APPLY, negative dims)
- ✓ Larger system (N=3, M=2, L=2)
- ✓ Fortran order
- ✗ Missing: Validation against known transformation
- ✗ Missing: Inverse transformation test (state-space → THETA → state-space)
- ✗ Missing: Numerical stability tests
- ✗ Missing: Comparison with reference implementation

**Weaknesses**:
1. Tests only verify output shapes and non-NaN values
2. No validation of actual transformation correctness
3. Random data in larger test (line 147) - not reproducible
4. No tests verifying mathematical properties

**Critical Issue**: Tests don't validate correctness, only that code runs without crashing.

### Recommendations

**Priority: HIGH** - Tests need significant improvement

1. **Generate reference data using Python control/slycot**:
```python
import control
import slycot  # Reference implementation

def test_tb01vy_vs_slycot():
    """Validate against slycot reference implementation"""
    n, m, l = 2, 1, 2

    # Create THETA vector
    theta = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
                     dtype=float, order='F')

    # Run slycot reference
    a_ref, b_ref, c_ref, d_ref, x0_ref = slycot.tb01vy(
        n, m, l, theta, apply='N'
    )

    # Run our implementation
    a, b, c, d, x0, info = slicot.tb01vy(n, m, l, theta, apply='N')

    assert info == 0
    np.testing.assert_allclose(a, a_ref, rtol=1e-14)
    np.testing.assert_allclose(b, b_ref, rtol=1e-14)
    np.testing.assert_allclose(c, c_ref, rtol=1e-14)
    np.testing.assert_allclose(d, d_ref, rtol=1e-14)
    np.testing.assert_allclose(x0, x0_ref, rtol=1e-14)
```

2. **Add round-trip tests**:
```python
def test_tb01vy_roundtrip():
    """Test state-space → THETA → state-space preserves system"""
    # Start with known state-space system
    a_orig = np.array([[0.5, 0.1], [0.0, 0.8]], order='F')
    b_orig = np.array([[1.0], [0.5]], order='F')
    # ... (need inverse transformation routine)
```

3. **Generate systems with Python control**:
```python
def test_tb01vy_control_system():
    """Generate test using control package"""
    sys = control.rss(3, 2, 2, strictly_proper=False)
    a, b, c, d = control.ssdata(sys)

    # Convert to output normal form parameters (via slycot)
    # Then test TB01VY conversion back to state-space
    # Verify system equivalence
```

4. **Mathematical property tests**:
```python
def test_tb01vy_observability():
    """Verify output normal form has full-rank observability matrix"""
    # Output normal form should have specific structure
    # Test observability matrix rank
```

**Verdict**: Replace most synthetic tests with slycot-validated tests or control package generation.

---

## 4. TF01MX - Matrix Product Time Response

### Test Data Source
**100% Synthetic** - No SLICOT examples
- HTML doc: "None" for Program Text/Data/Results
- No example files in `SLICOT-Reference/examples/`
- No benchmark data

### Current Test Quality: **Good**
**File**: `tests/python/test_tf01mx.py`

**Coverage Analysis**:
- ✓ Basic 2×2 system with hand-calculated results (lines 10-59)
- ✓ No inputs (M=0) case
- ✓ Zero states (N=0) case (feedthrough only)
- ✓ Zero outputs (NY=0) case
- ✓ Parameter validation (N<0)
- ✓ Workspace validation test
- ✓ Fortran order
- ✓ Tight tolerance (rtol=1e-14)

**Strengths**:
1. Hand-verified expected outputs (lines 36-53)
2. Multiple edge cases (M=0, N=0, NY=0)
3. Clear documentation of manual calculations
4. Good parameter validation

**Weaknesses**:
1. Only simple test cases (2×2 system)
2. No comparison with reference implementation
3. No tests with complex dynamics (oscillatory, unstable)
4. Missing tests for different parameter combinations

### Recommendations

**Priority: MEDIUM** - Good foundation, needs validation

1. **Validate against Python control package**:
```python
import control

def test_tf01mx_vs_control_lsim():
    """Compare with control.lsim (linear simulation)"""
    n, m, p = 2, 1, 1

    # System matrices
    a = np.array([[0.5, 0.1], [0.0, 0.8]], order='F')
    b = np.array([[1.0], [0.5]], order='F')
    c = np.array([[1.0, 0.0]], order='F')
    d = np.array([[0.0]], order='F')

    # Create state-space system
    sys = control.ss(a, b, c, d, dt=1.0)  # Discrete-time

    # Input sequence
    ny = 3
    u = np.array([[1.0], [0.5], [0.0]], order='F')
    x0 = np.array([1.0, 0.5], order='F')

    # control package simulation
    t = np.arange(ny)
    _, y_ref, _ = control.forced_response(sys, t, u.flatten(), x0)

    # Our implementation
    s = np.vstack([np.hstack([a, b]), np.hstack([c, d])])
    y, x_final, info = slicot.tf01mx(n, m, p, ny, s, u, x0)

    assert info == 0
    np.testing.assert_allclose(y.flatten(), y_ref, rtol=1e-13)
```

2. **Add slycot validation**:
```python
import slycot

def test_tf01mx_vs_slycot():
    """Validate against slycot.tf01mx"""
    n, m, p, ny = 2, 1, 1, 3
    s = np.array([[0.5, 0.1, 1.0],
                  [0.0, 0.8, 0.5],
                  [1.0, 0.0, 0.0]], order='F')
    u = np.array([[1.0], [0.5], [0.0]], order='F')
    x = np.array([1.0, 0.5], order='F')

    # slycot reference
    y_ref, x_ref = slycot.tf01mx(n, m, p, ny, s, u, x)

    # Our implementation
    y, x_final, info = slicot.tf01mx(n, m, p, ny, s, u, x)

    assert info == 0
    np.testing.assert_allclose(y, y_ref, rtol=1e-14)
    np.testing.assert_allclose(x_final, x_ref, rtol=1e-14)
```

3. **Generate random systems with control package**:
```python
def test_tf01mx_random_stable():
    """Test with randomly generated stable system"""
    sys = control.drss(4, 2, 3)  # Discrete-time random state-space
    a, b, c, d = control.ssdata(sys)

    # Generate random input
    ny = 10
    u = np.random.randn(ny, 2)
    x0 = np.random.randn(4)

    # Simulate with control package
    t = np.arange(ny)
    _, y_ref, x_ref = control.forced_response(sys, t, u.T, x0)

    # Our implementation
    s = np.vstack([np.hstack([a, b]), np.hstack([c, d])]).astype(float, order='F')
    y, x_final, info = slicot.tf01mx(4, 2, 3, ny, s, u.astype(float, order='F'),
                                     x0.astype(float, order='F'))

    assert info == 0
    np.testing.assert_allclose(y, y_ref.T, rtol=1e-12)
```

**Verdict**: Good synthetic tests. Add control package validation for confidence.

---

## Overall Recommendations

### Immediate Actions (Before Next PR)

1. **TB01VY (HIGH PRIORITY)**:
   - Add slycot validation tests
   - Replace shape-only tests with numerical validation
   - Add control package random system tests

2. **MA02AD (LOW PRIORITY)**:
   - Add property-based tests (transpose inverse)
   - Add error handling tests

3. **TF01MX (MEDIUM PRIORITY)**:
   - Add control.forced_response validation
   - Add slycot comparison tests
   - Keep existing hand-calculated tests

4. **TB01WD (LOW PRIORITY)**:
   - Already excellent
   - Optionally add M=0/P=0 edge cases
   - Optionally add control.rss generated tests

### Python Control Package Strategy

**Available capabilities**:
- `control.ss()` - Create state-space systems
- `control.rss()` - Random continuous-time stable system
- `control.drss()` - Random discrete-time stable system
- `control.forced_response()` - Simulate with arbitrary inputs
- `control.ssdata()` - Extract A, B, C, D matrices

**Use for**:
1. Generating diverse test systems (stable, unstable, MIMO)
2. Cross-validation with established implementations
3. Property testing (controllability, observability, stability)
4. Stress testing with large random systems

**Don't use for**:
1. Primary validation (prefer slycot reference or SLICOT examples)
2. Simple utility routines (like MA02AD)

### Long-Term Testing Strategy

1. **Test Hierarchy** (priority order):
   - SLICOT example data (gold standard)
   - slycot validation (reference implementation)
   - Python control package (property/stress tests)
   - Hand-calculated synthetic (simple cases only)

2. **Coverage Goals**:
   - Every routine: ≥1 authoritative test (SLICOT or slycot)
   - Every routine: Error handling tests
   - Every routine: Edge cases (zero dimensions, single element)
   - Complex routines: Property tests, stress tests

3. **Test Quality Metrics**:
   - Source traceability (SLICOT example > slycot > synthetic)
   - Numerical validation (actual values, not just shapes)
   - Coverage breadth (parameters, dimensions, edge cases)
   - Reproducibility (no unexplained random data)

---

## Action Items

- [ ] Enhance TB01VY with slycot validation
- [ ] Add control package tests to TF01MX
- [ ] Add property tests to MA02AD
- [ ] Document test data sources in test file headers
- [ ] Create test generation utilities using control package
- [ ] Add slycot as test dependency (already in requirements.txt)

**Estimated Effort**: 4-6 hours total
- TB01VY: 2-3 hours (high priority)
- TF01MX: 1-2 hours
- MA02AD: 1 hour
- Documentation: 30 min
