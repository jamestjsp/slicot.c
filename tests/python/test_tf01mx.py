import numpy as np
import pytest

try:
    from slicot import tf01mx
except ImportError:
    pytest.skip("tf01mx not available", allow_module_level=True)


def test_tf01mx_basic():
    """Test TF01MX with simple 2x2 system."""
    # System: x(k+1) = A*x(k) + B*u(k), y(k) = C*x(k) + D*u(k)
    # A = [0.5, 0.1; 0.0, 0.8], B = [1.0; 0.5], C = [1.0, 0.0], D = [0.0]
    n = 2
    m = 1
    p = 1
    ny = 3

    # System matrix S = [A B; C D]
    s = np.array([
        [0.5, 0.1, 1.0],
        [0.0, 0.8, 0.5],
        [1.0, 0.0, 0.0]
    ], dtype=float, order='F')

    # Input sequence u(k) for k=1,2,3
    u = np.array([
        [1.0],
        [0.5],
        [0.0]
    ], dtype=float, order='F')

    # Initial state
    x = np.array([1.0, 0.5], dtype=float, order='F')

    # Expected outputs computed manually:
    # k=1: y(1) = C*x(1) + D*u(1) = 1.0*1.0 + 0.0*0.0 + 0.0*1.0 = 1.0
    #      x(2) = A*x(1) + B*u(1) = [0.5*1.0+0.1*0.5; 0.0*1.0+0.8*0.5] + [1.0*1.0; 0.5*1.0]
    #           = [0.55; 0.4] + [1.0; 0.5] = [1.55; 0.9]
    # k=2: y(2) = C*x(2) + D*u(2) = 1.0*1.55 + 0.0*0.9 + 0.0*0.5 = 1.55
    #      x(3) = A*x(2) + B*u(2) = [0.5*1.55+0.1*0.9; 0.0*1.55+0.8*0.9] + [1.0*0.5; 0.5*0.5]
    #           = [0.865; 0.72] + [0.5; 0.25] = [1.365; 0.97]
    # k=3: y(3) = C*x(3) + D*u(3) = 1.0*1.365 + 0.0*0.97 + 0.0*0.0 = 1.365
    #      x(4) = A*x(3) + B*u(3) = [0.5*1.365+0.1*0.97; 0.0*1.365+0.8*0.97] + [0.0; 0.0]
    #           = [0.7795; 0.776]

    y_expected = np.array([
        [1.0],
        [1.55],
        [1.365]
    ], dtype=float, order='F')

    x_final_expected = np.array([0.7795, 0.776], dtype=float, order='F')

    y, x_final, info = tf01mx(n, m, p, ny, s, u, x)

    assert info == 0
    np.testing.assert_allclose(y, y_expected, rtol=1e-14)
    np.testing.assert_allclose(x_final, x_final_expected, rtol=1e-14)


def test_tf01mx_no_inputs():
    """Test TF01MX with system having no inputs (M=0)."""
    # System: x(k+1) = A*x(k), y(k) = C*x(k)
    # A = [0.9, 0.0; 0.0, 0.8], C = [1.0, 1.0]
    n = 2
    m = 0
    p = 1
    ny = 2

    # System matrix S = [A; C]
    s = np.array([
        [0.9, 0.0],
        [0.0, 0.8],
        [1.0, 1.0]
    ], dtype=float, order='F')

    # No inputs
    u = np.zeros((ny, 0), dtype=float, order='F')

    # Initial state
    x = np.array([2.0, 1.0], dtype=float, order='F')

    # Expected outputs:
    # k=1: y(1) = C*x(1) = 1.0*2.0 + 1.0*1.0 = 3.0
    #      x(2) = A*x(1) = [0.9*2.0; 0.8*1.0] = [1.8; 0.8]
    # k=2: y(2) = C*x(2) = 1.0*1.8 + 1.0*0.8 = 2.6
    #      x(3) = A*x(2) = [0.9*1.8; 0.8*0.8] = [1.62; 0.64]

    y_expected = np.array([
        [3.0],
        [2.6]
    ], dtype=float, order='F')

    x_final_expected = np.array([1.62, 0.64], dtype=float, order='F')

    y, x_final, info = tf01mx(n, m, p, ny, s, u, x)

    assert info == 0
    np.testing.assert_allclose(y, y_expected, rtol=1e-14)
    np.testing.assert_allclose(x_final, x_final_expected, rtol=1e-14)


def test_tf01mx_zero_states():
    """Test TF01MX with zero states (N=0)."""
    # Non-dynamic system: y(k) = D*u(k)
    n = 0
    m = 2
    p = 1
    ny = 2

    # System matrix S = [D]
    s = np.array([
        [1.5, 0.5]
    ], dtype=float, order='F')

    # Input sequence
    u = np.array([
        [1.0, 2.0],
        [0.5, 1.0]
    ], dtype=float, order='F')

    # No state
    x = np.zeros(0, dtype=float, order='F')

    # Expected outputs: y(k) = D*u(k)
    # k=1: y(1) = 1.5*1.0 + 0.5*2.0 = 2.5
    # k=2: y(2) = 1.5*0.5 + 0.5*1.0 = 1.25

    y_expected = np.array([
        [2.5],
        [1.25]
    ], dtype=float, order='F')

    y, x_final, info = tf01mx(n, m, p, ny, s, u, x)

    assert info == 0
    np.testing.assert_allclose(y, y_expected, rtol=1e-14)
    assert len(x_final) == 0


def test_tf01mx_zero_outputs():
    """Test TF01MX with zero output steps (NY=0)."""
    n = 2
    m = 1
    p = 1
    ny = 0

    s = np.array([[0.5, 0.1, 1.0], [0.0, 0.8, 0.5], [1.0, 0.0, 0.0]],
                 dtype=float, order='F')
    u = np.zeros((0, m), dtype=float, order='F')
    x = np.array([1.0, 0.5], dtype=float, order='F')

    y, x_final, info = tf01mx(n, m, p, ny, s, u, x)

    assert info == 0
    assert y.shape == (0, p)
    np.testing.assert_allclose(x_final, [1.0, 0.5], rtol=1e-14)


def test_tf01mx_invalid_n():
    """Test TF01MX with invalid N."""
    with pytest.raises(ValueError, match="N must be >= 0"):
        n = -1
        m = 1
        p = 1
        ny = 1
        s = np.zeros((2, 2), dtype=float, order='F')
        u = np.zeros((1, 1), dtype=float, order='F')
        x = np.zeros(1, dtype=float, order='F')
        tf01mx(n, m, p, ny, s, u, x)


def test_tf01mx_invalid_ldwork():
    """Test TF01MX workspace validation."""
    # This should trigger LDWORK validation
    # For M>0, LDWORK >= 2*N+M+P
    n = 2
    m = 1
    p = 1
    ny = 1

    s = np.array([[0.5, 0.1, 1.0], [0.0, 0.8, 0.5], [1.0, 0.0, 0.0]],
                 dtype=float, order='F')
    u = np.array([[1.0]], dtype=float, order='F')
    x = np.array([1.0, 0.5], dtype=float, order='F')

    # This should work fine (wrapper allocates workspace)
    y, x_final, info = tf01mx(n, m, p, ny, s, u, x)
    assert info == 0
