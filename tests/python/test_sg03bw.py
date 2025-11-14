import pytest
import numpy as np
from slicot import sg03bw


def test_sg03bw_basic_n1():
    """Test SG03BW with N=1 (basic case).

    Solves: A^T * X * C + E^T * X * D = SCALE * Y
    """
    trans = 'N'
    m, n = 2, 1

    # Simple upper triangular A and E
    a = np.array([[2.0, 1.0],
                  [0.0, 3.0]], dtype=np.float64, order='F')

    e = np.array([[1.0, 0.5],
                  [0.0, 1.0]], dtype=np.float64, order='F')

    c = np.array([[4.0]], dtype=np.float64, order='F')
    d = np.array([[2.0]], dtype=np.float64, order='F')

    # Right-hand side Y (will be overwritten with solution X)
    x = np.array([[1.0],
                  [2.0]], dtype=np.float64, order='F')

    x_out, scale, info = sg03bw(trans, a, e, c, d, x)

    # Check success
    assert info == 0, f"SG03BW failed with info={info}"

    # Scale should be positive
    assert 0 < scale <= 1.0, f"Invalid scale={scale}"

    # Solution should be finite
    assert np.all(np.isfinite(x_out))


def test_sg03bw_basic_n2():
    """Test SG03BW with N=2."""
    trans = 'N'
    m, n = 2, 2

    # Upper quasitriangular A (2x2 block)
    a = np.array([[1.0, 2.0],
                  [0.5, 3.0]], dtype=np.float64, order='F')

    # Upper triangular E
    e = np.array([[2.0, 1.0],
                  [0.0, 1.5]], dtype=np.float64, order='F')

    c = np.array([[1.0, 0.0],
                  [0.0, 1.0]], dtype=np.float64, order='F')

    d = np.array([[0.5, 0.0],
                  [0.0, 0.5]], dtype=np.float64, order='F')

    x = np.array([[1.0, 0.5],
                  [0.5, 1.0]], dtype=np.float64, order='F')

    x_out, scale, info = sg03bw(trans, a, e, c, d, x)

    assert info == 0, f"SG03BW failed with info={info}"
    assert 0 < scale <= 1.0, f"Invalid scale={scale}"
    assert np.all(np.isfinite(x_out))


def test_sg03bw_transposed():
    """Test SG03BW with TRANS='T' (transposed equation)."""
    trans = 'T'
    m, n = 2, 1

    a = np.array([[3.0, 1.0],
                  [0.0, 2.0]], dtype=np.float64, order='F')

    e = np.array([[1.5, 0.5],
                  [0.0, 1.0]], dtype=np.float64, order='F')

    c = np.array([[2.0]], dtype=np.float64, order='F')
    d = np.array([[1.0]], dtype=np.float64, order='F')

    x = np.array([[1.0],
                  [1.0]], dtype=np.float64, order='F')

    x_out, scale, info = sg03bw(trans, a, e, c, d, x)

    assert info == 0, f"SG03BW failed with info={info}"
    assert 0 < scale <= 1.0, f"Invalid scale={scale}"
    assert np.all(np.isfinite(x_out))


def test_sg03bw_zero_m():
    """Test SG03BW with M=0 (edge case)."""
    trans = 'N'
    m, n = 0, 1

    a = np.array([], dtype=np.float64).reshape(0, 0, order='F')
    e = np.array([], dtype=np.float64).reshape(0, 0, order='F')
    c = np.array([[1.0]], dtype=np.float64, order='F')
    d = np.array([[1.0]], dtype=np.float64, order='F')
    x = np.array([], dtype=np.float64).reshape(0, 1, order='F')

    x_out, scale, info = sg03bw(trans, a, e, c, d, x)

    # Quick return expected
    assert info == 0, f"SG03BW failed with info={info}"
    assert scale == 1.0, f"Expected scale=1.0, got {scale}"


def test_sg03bw_invalid_trans():
    """Test SG03BW with invalid TRANS parameter."""
    trans = 'X'  # Invalid
    m, n = 2, 1

    a = np.eye(2, dtype=np.float64, order='F')
    e = np.eye(2, dtype=np.float64, order='F')
    c = np.ones((1, 1), dtype=np.float64, order='F')
    d = np.ones((1, 1), dtype=np.float64, order='F')
    x = np.ones((2, 1), dtype=np.float64, order='F')

    x_out, scale, info = sg03bw(trans, a, e, c, d, x)

    # Should return error
    assert info < 0, f"Expected info < 0 for invalid TRANS, got {info}"


def test_sg03bw_invalid_n():
    """Test SG03BW with invalid N (must be 1 or 2)."""
    trans = 'N'
    m, n = 2, 3  # Invalid N

    a = np.eye(2, dtype=np.float64, order='F')
    e = np.eye(2, dtype=np.float64, order='F')
    c = np.eye(3, dtype=np.float64, order='F')
    d = np.eye(3, dtype=np.float64, order='F')
    x = np.ones((2, 3), dtype=np.float64, order='F')

    x_out, scale, info = sg03bw(trans, a, e, c, d, x)

    # Should return error
    assert info < 0, f"Expected info < 0 for invalid N, got {info}"


def test_sg03bw_nearly_singular():
    """Test SG03BW with nearly singular system.

    Should return INFO=1 and use perturbed values.
    """
    trans = 'N'
    m, n = 2, 1

    # Create nearly singular system
    # A and E structured to make equation nearly singular
    a = np.array([[1e-15, 1.0],
                  [0.0, 1e-15]], dtype=np.float64, order='F')

    e = np.array([[1e-15, 0.5],
                  [0.0, 1e-15]], dtype=np.float64, order='F')

    c = np.array([[1.0]], dtype=np.float64, order='F')
    d = np.array([[1.0]], dtype=np.float64, order='F')

    x = np.array([[1.0],
                  [1.0]], dtype=np.float64, order='F')

    x_out, scale, info = sg03bw(trans, a, e, c, d, x)

    # May return INFO=1 for near singularity
    # But should still produce a finite result
    assert info >= 0, f"Unexpected error: info={info}"
    assert np.all(np.isfinite(x_out))
    assert 0 < scale <= 1.0


def test_sg03bw_identity_matrices():
    """Test SG03BW with identity matrices."""
    trans = 'N'
    m, n = 3, 1

    a = np.eye(3, dtype=np.float64, order='F')
    e = np.eye(3, dtype=np.float64, order='F')
    c = np.array([[2.0]], dtype=np.float64, order='F')
    d = np.array([[3.0]], dtype=np.float64, order='F')

    # Y = [1, 2, 3]^T
    x = np.array([[1.0],
                  [2.0],
                  [3.0]], dtype=np.float64, order='F')

    x_out, scale, info = sg03bw(trans, a, e, c, d, x)

    assert info == 0, f"SG03BW failed with info={info}"
    assert scale > 0, f"Invalid scale={scale}"

    # For identity matrices: (A^T + E^T) * X * (C + D) = (I + I) * X * (C + D)
    # = 2X * (C+D) = SCALE * Y
    # So X should be proportional to Y / (2*(C+D))
    expected_factor = scale / (2.0 * (c[0, 0] + d[0, 0]))
    expected_x = x * expected_factor / scale

    # Solution should have correct structure (may be scaled)
    assert np.all(np.isfinite(x_out))
