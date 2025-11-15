"""
Tests for TB01VY - Output normal form to state-space conversion
"""
import numpy as np
import pytest

slicot = pytest.importorskip("slicot")


def test_tb01vy_basic_apply_n():
    """Test basic functionality with APPLY='N'"""
    n = 2
    m = 1
    l = 2

    # Create THETA parameter vector
    # THETA layout: [A,C params (N*L)], [B params (N*M)], [D params (L*M)], [x0 (N)]
    # Total: N*(L+M+1) + L*M = 2*(2+1+1) + 2*1 = 8 + 2 = 10
    theta = np.array([
        # A,C parameters (N*L = 4)
        0.1, 0.2, 0.3, 0.4,
        # B parameters (N*M = 2)
        0.5, 0.6,
        # D parameters (L*M = 2)
        0.7, 0.8,
        # x0 (N = 2)
        0.9, 1.0
    ], dtype=float, order='F')

    a, b, c, d, x0, info = slicot.tb01vy(n, m, l, theta, apply='N')

    assert info == 0
    assert a.shape == (n, n)
    assert b.shape == (n, m)
    assert c.shape == (l, n)
    assert d.shape == (l, m)
    assert x0.shape == (n,)

    # Check D matrix (direct copy from THETA)
    np.testing.assert_allclose(d, [[0.7], [0.8]], rtol=1e-14)

    # Check B matrix (direct copy from THETA)
    np.testing.assert_allclose(b, [[0.5], [0.6]], rtol=1e-14)

    # Check x0 (direct copy from THETA)
    np.testing.assert_allclose(x0, [0.9, 1.0], rtol=1e-14)

    # A and C are computed via orthogonal transformations
    # Just verify dimensions and no NaNs
    assert not np.any(np.isnan(a))
    assert not np.any(np.isnan(c))


def test_tb01vy_basic_apply_a():
    """Test basic functionality with APPLY='A' (bijective mapping)"""
    n = 2
    m = 1
    l = 2

    theta = np.array([
        0.1, 0.2, 0.3, 0.4,
        0.5, 0.6,
        0.7, 0.8,
        0.9, 1.0
    ], dtype=float, order='F')

    a, b, c, d, x0, info = slicot.tb01vy(n, m, l, theta, apply='A')

    assert info == 0
    assert a.shape == (n, n)
    assert b.shape == (n, m)
    assert c.shape == (l, n)
    assert d.shape == (l, m)
    assert x0.shape == (n,)

    # D, B, x0 should be same regardless of APPLY
    np.testing.assert_allclose(d, [[0.7], [0.8]], rtol=1e-14)
    np.testing.assert_allclose(b, [[0.5], [0.6]], rtol=1e-14)
    np.testing.assert_allclose(x0, [0.9, 1.0], rtol=1e-14)

    # A and C will differ from APPLY='N'
    assert not np.any(np.isnan(a))
    assert not np.any(np.isnan(c))


def test_tb01vy_zero_dimensions():
    """Test edge cases with zero dimensions"""
    # N=0: Should return quickly
    n, m, l = 0, 1, 1
    theta = np.array([0.5], dtype=float, order='F')  # D only

    a, b, c, d, x0, info = slicot.tb01vy(n, m, l, theta, apply='N')

    assert info == 0
    assert d.shape == (l, m)
    np.testing.assert_allclose(d, [[0.5]], rtol=1e-14)

    # M=0: No inputs
    n, m, l = 2, 0, 1
    theta = np.array([0.1, 0.2, 0.3, 0.4], dtype=float, order='F')  # A,C params + x0

    a, b, c, d, x0, info = slicot.tb01vy(n, m, l, theta, apply='N')

    assert info == 0
    assert a.shape == (n, n)
    assert c.shape == (l, n)

    # L=0: No outputs (special case returns x0 only)
    n, m, l = 2, 1, 0
    theta = np.array([0.5, 0.6, 0.9, 1.0], dtype=float, order='F')  # B + x0

    a, b, c, d, x0, info = slicot.tb01vy(n, m, l, theta, apply='N')

    assert info == 0
    np.testing.assert_allclose(x0, [0.9, 1.0], rtol=1e-14)


def test_tb01vy_error_handling():
    """Test parameter validation"""
    n, m, l = 2, 1, 2

    # LTHETA too small
    theta_short = np.array([0.1, 0.2], dtype=float, order='F')

    with pytest.raises(ValueError, match="ltheta"):
        slicot.tb01vy(n, m, l, theta_short, apply='N')

    # Invalid APPLY
    theta = np.zeros(10, dtype=float, order='F')

    with pytest.raises(ValueError, match="apply"):
        slicot.tb01vy(n, m, l, theta, apply='X')

    # Negative dimensions
    with pytest.raises(ValueError, match="n"):
        slicot.tb01vy(-1, m, l, theta, apply='N')


def test_tb01vy_larger_system():
    """Test with larger system dimensions"""
    n = 3
    m = 2
    l = 2

    # Total: N*(L+M+1) + L*M = 3*(2+2+1) + 2*2 = 15 + 4 = 19
    # Use APPLY='A' for random data since norm(theta_i) may exceed 1
    theta = np.random.RandomState(42).randn(19)
    theta = theta.astype(float, order='F')

    a, b, c, d, x0, info = slicot.tb01vy(n, m, l, theta, apply='A')

    assert info == 0
    assert a.shape == (n, n)
    assert b.shape == (n, m)
    assert c.shape == (l, n)
    assert d.shape == (l, m)
    assert x0.shape == (n,)

    # All outputs should be finite
    assert np.all(np.isfinite(a))
    assert np.all(np.isfinite(b))
    assert np.all(np.isfinite(c))
    assert np.all(np.isfinite(d))
    assert np.all(np.isfinite(x0))
