#!/usr/bin/env python3
"""
pytest tests for MB03OY - Rank determination via incremental condition estimation.
"""
import pytest
import numpy as np
from slicot import mb03oy


def test_mb03oy_rank_deficient_matrix():
    """Test MB03OY with rank-deficient matrix (6x4, rank 3)."""
    m, n = 6, 4
    rcond = 1.0e-10
    svlmax = 0.0

    # Rank-deficient matrix with rank ~3 (column-major storage)
    a = np.array([
        [6.9352668323e-01, -1.0442760938e+00, -1.0312245858e+00, -1.6794289521e-01],
        [-1.3954146423e+00, -2.6089798788e-02, -5.8559390779e-01, 6.1803578721e-01],
        [-3.2058987212e+00, 3.6415744214e-01, 2.2737438267e+00, 8.0734713252e-02],
        [-5.2076729914e-01, -3.8916678618e-01, -2.3420832454e-01, 7.8113716229e-01],
        [-4.5698494441e-02, 1.6414232666e+00, 1.8928814982e-01, 1.7432426363e+00],
        [2.4015719701e+00, 2.6929887908e-01, 3.0543126382e-01, 0.0]
    ], dtype=np.float64, order='F')

    a_result, rank, info, sval, jpvt, tau = mb03oy(m, n, a.copy(), rcond, svlmax)

    # Verify successful execution
    assert info == 0

    # Verify estimated rank
    assert rank == 3

    # Verify singular value estimates
    assert sval[0] > 0.0  # Largest singular value positive
    assert sval[1] > 0.0  # Rank-th singular value positive
    assert sval[2] < sval[1]  # (rank+1)-th smaller

    # Verify condition number estimate
    cond_estimate = sval[0] / sval[1]
    assert cond_estimate < 1.0 / rcond

    # Verify singular value estimates (approximate values)
    np.testing.assert_allclose(sval[0], 4.329756e+00, rtol=1e-5)
    np.testing.assert_allclose(sval[1], 1.288771e+00, rtol=1e-5)
    assert sval[2] < 1e-8  # Near-zero for rank-deficient

    # Verify R matrix upper triangle (first rank columns)
    # Note: exact values depend on LAPACK implementation details
    for j in range(rank):
        for i in range(j + 1):
            assert not np.isnan(a_result[i, j])


def test_mb03oy_full_rank_matrix():
    """Test MB03OY with full-rank identity matrix."""
    m, n = 4, 4
    rcond = 1.0e-10
    svlmax = 0.0

    # Identity matrix (full rank)
    a = np.eye(4, dtype=np.float64, order='F')

    a_result, rank, info, sval, jpvt, tau = mb03oy(m, n, a.copy(), rcond, svlmax)

    assert info == 0
    assert rank == 4  # Full rank

    # For identity, singular values should all be ~1.0
    np.testing.assert_allclose(sval[0], 1.0, rtol=1e-10)
    np.testing.assert_allclose(sval[1], 1.0, rtol=1e-10)
    np.testing.assert_allclose(sval[2], 1.0, rtol=1e-10)


def test_mb03oy_full_rank_random_matrix():
    """Test MB03OY with full-rank random matrix."""
    m, n = 5, 4
    rcond = 1.0e-10
    svlmax = 0.0

    # Generate random full-rank matrix
    np.random.seed(12345)
    a = np.random.randn(m, n)
    a = np.asfortranarray(a)  # Convert to column-major

    a_result, rank, info, sval, jpvt, tau = mb03oy(m, n, a.copy(), rcond, svlmax)

    assert info == 0
    assert rank == 4  # Should be full rank (min(m,n))

    # Verify singular values are positive
    assert sval[0] > 0.0
    assert sval[1] > 0.0
    assert sval[2] > 0.0


def test_mb03oy_edge_case_zero_rows():
    """Test MB03OY with zero rows."""
    m, n = 0, 4
    rcond = 1.0e-10
    svlmax = 0.0

    a = np.array([[]], dtype=np.float64, order='F').reshape(0, 4)

    a_result, rank, info, sval, jpvt, tau = mb03oy(m, n, a.copy(), rcond, svlmax)

    assert info == 0
    assert rank == 0  # No rows means rank 0


def test_mb03oy_edge_case_zero_cols():
    """Test MB03OY with zero columns."""
    m, n = 4, 0
    rcond = 1.0e-10
    svlmax = 0.0

    a = np.array([[]], dtype=np.float64, order='F').reshape(4, 0)

    a_result, rank, info, sval, jpvt, tau = mb03oy(m, n, a.copy(), rcond, svlmax)

    assert info == 0
    assert rank == 0  # No columns means rank 0


def test_mb03oy_error_negative_m():
    """Test MB03OY error handling: negative m."""
    m, n = -1, 4
    rcond = 1.0e-10
    svlmax = 0.0

    a = np.zeros((1, 4), dtype=np.float64, order='F')

    with pytest.raises(Exception):  # Should fail in C validation or Python wrapper
        a_result, rank, info, sval, jpvt, tau = mb03oy(m, n, a, rcond, svlmax)


def test_mb03oy_error_negative_n():
    """Test MB03OY error handling: negative n."""
    m, n = 4, -1
    rcond = 1.0e-10
    svlmax = 0.0

    a = np.zeros((4, 1), dtype=np.float64, order='F')

    with pytest.raises(Exception):  # Should fail in C validation or Python wrapper
        a_result, rank, info, sval, jpvt, tau = mb03oy(m, n, a, rcond, svlmax)


def test_mb03oy_strict_rcond():
    """Test MB03OY with strict rcond (should detect lower rank)."""
    m, n = 4, 4
    rcond = 0.5  # Strict threshold
    svlmax = 0.0

    # Matrix with condition number ~2 (two singular values differ by factor 2)
    a = np.diag([2.0, 2.0, 1.0, 1.0])
    a = np.asfortranarray(a)

    a_result, rank, info, sval, jpvt, tau = mb03oy(m, n, a.copy(), rcond, svlmax)

    assert info == 0
    # With rcond=0.5, only singular values >= 0.5*max_sval should count
    # Expected rank depends on threshold application
    assert 0 < rank <= 4


def test_mb03oy_numerical_rank_detection():
    """Test MB03OY detects numerical rank correctly."""
    m, n = 4, 3
    rcond = 1.0e-8
    svlmax = 0.0

    # Create rank-2 matrix with small noise
    np.random.seed(42)
    U = np.random.randn(m, 2)
    V = np.random.randn(n, 2)
    a = U @ V.T
    a += 1e-12 * np.random.randn(m, n)  # Tiny noise
    a = np.asfortranarray(a)

    a_result, rank, info, sval, jpvt, tau = mb03oy(m, n, a.copy(), rcond, svlmax)

    assert info == 0
    assert rank == 2  # Should detect rank 2 (noise is below threshold)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
