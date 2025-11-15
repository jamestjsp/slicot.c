"""Tests for MA02AD - Matrix transposition"""
import numpy as np
import pytest
from slicot import ma02ad


def test_ma02ad_full_transpose():
    """Test full matrix transpose (JOB='F')"""
    # 3x2 matrix
    a = np.array([[1.0, 4.0],
                  [2.0, 5.0],
                  [3.0, 6.0]], order='F')

    b = ma02ad('F', a)

    # Expected: 2x3 transpose
    expected = np.array([[1.0, 2.0, 3.0],
                        [4.0, 5.0, 6.0]], order='F')

    np.testing.assert_allclose(b, expected, rtol=1e-14)


def test_ma02ad_upper_triangular():
    """Test upper triangular transpose (JOB='U')"""
    # 3x3 matrix - only upper triangle used
    a = np.array([[1.0, 2.0, 3.0],
                  [4.0, 5.0, 6.0],
                  [7.0, 8.0, 9.0]], order='F')

    b = ma02ad('U', a)

    # Expected: B(j,i) = A(i,j) for i <= j (upper triangle)
    # B is transpose, so upper triangle of A becomes upper triangle of B
    expected = np.array([[1.0, 0.0, 0.0],
                        [2.0, 5.0, 0.0],
                        [3.0, 6.0, 9.0]], order='F')

    np.testing.assert_allclose(b, expected, rtol=1e-14)


def test_ma02ad_lower_triangular():
    """Test lower triangular transpose (JOB='L')"""
    # 3x3 matrix - only lower triangle used
    a = np.array([[1.0, 2.0, 3.0],
                  [4.0, 5.0, 6.0],
                  [7.0, 8.0, 9.0]], order='F')

    b = ma02ad('L', a)

    # Expected: B(j,i) = A(i,j) for i >= j (lower triangle)
    # Lower triangle of A becomes lower triangle of B after transpose
    expected = np.array([[1.0, 4.0, 7.0],
                        [0.0, 5.0, 8.0],
                        [0.0, 0.0, 9.0]], order='F')

    np.testing.assert_allclose(b, expected, rtol=1e-14)


def test_ma02ad_rectangular():
    """Test rectangular matrix transpose"""
    # 2x4 matrix
    a = np.array([[1.0, 3.0, 5.0, 7.0],
                  [2.0, 4.0, 6.0, 8.0]], order='F')

    b = ma02ad('F', a)

    # Expected: 4x2 transpose
    expected = np.array([[1.0, 2.0],
                        [3.0, 4.0],
                        [5.0, 6.0],
                        [7.0, 8.0]], order='F')

    np.testing.assert_allclose(b, expected, rtol=1e-14)


def test_ma02ad_upper_trapezoid():
    """Test upper trapezoidal transpose (M > N)"""
    # 4x2 matrix - upper trapezoid
    a = np.array([[1.0, 3.0],
                  [2.0, 4.0],
                  [5.0, 6.0],
                  [7.0, 8.0]], order='F')

    b = ma02ad('U', a)

    # Expected: B(j,i) = A(i,j) for i <= min(j,m)
    # Only elements where row <= col in original
    expected = np.array([[1.0, 0.0, 0.0, 0.0],
                        [3.0, 4.0, 0.0, 0.0]], order='F')

    np.testing.assert_allclose(b, expected, rtol=1e-14)


def test_ma02ad_zero_rows():
    """Test edge case: zero rows"""
    a = np.array([], dtype=np.float64).reshape(0, 3, order='F')

    b = ma02ad('F', a)

    assert b.shape == (3, 0)


def test_ma02ad_zero_cols():
    """Test edge case: zero columns"""
    a = np.array([], dtype=np.float64).reshape(3, 0, order='F')

    b = ma02ad('F', a)

    assert b.shape == (0, 3)


def test_ma02ad_single_element():
    """Test edge case: single element"""
    a = np.array([[5.0]], order='F')

    b = ma02ad('F', a)

    np.testing.assert_allclose(b, np.array([[5.0]], order='F'), rtol=1e-14)
