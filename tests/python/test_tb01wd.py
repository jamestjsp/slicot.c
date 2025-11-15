"""
Tests for TB01WD - Orthogonal similarity transformation to real Schur form

Test data from SLICOT-Reference/doc/TB01WD.html example
"""
import numpy as np
import pytest

try:
    import slicot
    HAS_SLICOT = True
except ImportError:
    HAS_SLICOT = False


@pytest.mark.skipif(not HAS_SLICOT, reason="slicot module not available")
class TestTB01WD:
    def test_example_from_html(self):
        """Test TB01WD with example from HTML documentation

        From HTML Program Data - READ statements show row-wise reading:
        READ ( NIN, FMT = * ) ( ( A(I,J), J = 1,N ), I = 1,N )
        READ ( NIN, FMT = * ) ( ( B(I,J), J = 1,M ), I = 1, N )
        READ ( NIN, FMT = * ) ( ( C(I,J), J = 1,N ), I = 1,P )
        """
        n, m, p = 5, 2, 3

        # A matrix (5x5) - row-wise reading
        a = np.array([
            [-0.04165,    4.9200,   -4.9200,         0,         0],
            [-1.387944,   -3.3300,         0,         0,         0],
            [   0.5450,         0,         0,   -0.5450,         0],
            [        0,         0,    4.9200,  -0.04165,    4.9200],
            [        0,         0,         0, -1.387944,   -3.3300]
        ], dtype=np.float64, order='F')

        # B matrix (5x2) - row-wise reading
        b = np.array([
            [       0,         0],
            [  3.3300,         0],
            [       0,         0],
            [       0,         0],
            [       0,    3.3300]
        ], dtype=np.float64, order='F')

        # C matrix (3x5) - row-wise reading
        c = np.array([
            [1, 0, 0, 0, 0],
            [0, 0, 1, 0, 0],
            [0, 0, 0, 1, 0]
        ], dtype=np.float64, order='F')

        # Expected eigenvalues (real, imag pairs)
        expected_wr = np.array([-0.7483, -0.7483, -1.6858, -1.6858, -1.8751])
        expected_wi = np.array([2.9940, -2.9940, 2.0311, -2.0311, 0.0])

        # Expected transformed A (U'*A*U in real Schur form)
        expected_a = np.array([
            [-0.7483,  -8.6406,   0.0000,   0.0000,   1.1745],
            [ 1.0374,  -0.7483,   0.0000,   0.0000,  -2.1164],
            [ 0.0000,   0.0000,  -1.6858,   5.5669,   0.0000],
            [ 0.0000,   0.0000,  -0.7411,  -1.6858,   0.0000],
            [ 0.0000,   0.0000,   0.0000,   0.0000,  -1.8751]
        ], dtype=np.float64, order='F')

        # Expected transformed B (U'*B)
        expected_b = np.array([
            [-0.5543,   0.5543],
            [-1.6786,   1.6786],
            [-0.8621,  -0.8621],
            [ 2.1912,   2.1912],
            [-1.5555,   1.5555]
        ], dtype=np.float64, order='F')

        # Expected transformed C (C*U)
        expected_c = np.array([
            [ 0.6864,  -0.0987,   0.6580,   0.2589,  -0.1381],
            [-0.0471,   0.6873,   0.0000,   0.0000,  -0.7249],
            [-0.6864,   0.0987,   0.6580,   0.2589,   0.1381]
        ], dtype=np.float64, order='F')

        # Expected transformation matrix U
        expected_u = np.array([
            [ 0.6864,  -0.0987,   0.6580,   0.2589,  -0.1381],
            [-0.1665,  -0.5041,  -0.2589,   0.6580,  -0.4671],
            [-0.0471,   0.6873,   0.0000,   0.0000,  -0.7249],
            [-0.6864,   0.0987,   0.6580,   0.2589,   0.1381],
            [ 0.1665,   0.5041,  -0.2589,   0.6580,   0.4671]
        ], dtype=np.float64, order='F')

        # Call tb01wd
        a_out, b_out, c_out, u, wr, wi, info = slicot.tb01wd(n, m, p, a, b, c)

        # Check success
        assert info == 0, f"tb01wd failed with info={info}"

        # Check eigenvalues (tolerance for HTML 4-decimal display)
        np.testing.assert_allclose(wr, expected_wr, rtol=1e-3, atol=1e-4)
        np.testing.assert_allclose(wi, expected_wi, rtol=1e-3, atol=1e-4)

        # Check transformed matrices (tolerance for HTML 4-decimal display)
        np.testing.assert_allclose(a_out, expected_a, rtol=1e-3, atol=1e-4)
        np.testing.assert_allclose(b_out, expected_b, rtol=1e-3, atol=1e-4)
        np.testing.assert_allclose(c_out, expected_c, rtol=1e-3, atol=1e-4)
        np.testing.assert_allclose(u, expected_u, rtol=1e-3, atol=1e-4)

    def test_zero_dimension(self):
        """Test quick return for N=0"""
        n, m, p = 0, 0, 0
        a = np.zeros((1, 1), dtype=np.float64, order='F')
        b = np.zeros((1, 1), dtype=np.float64, order='F')
        c = np.zeros((1, 1), dtype=np.float64, order='F')

        a_out, b_out, c_out, u, wr, wi, info = slicot.tb01wd(n, m, p, a, b, c)

        assert info == 0

    def test_invalid_parameters(self):
        """Test error handling for invalid parameters"""
        # Invalid N < 0
        with pytest.raises((ValueError, RuntimeError)):
            slicot.tb01wd(-1, 1, 1, np.zeros((1,1), order='F'),
                         np.zeros((1,1), order='F'), np.zeros((1,1), order='F'))

        # Invalid M < 0
        with pytest.raises((ValueError, RuntimeError)):
            slicot.tb01wd(1, -1, 1, np.zeros((1,1), order='F'),
                         np.zeros((1,1), order='F'), np.zeros((1,1), order='F'))

        # Invalid P < 0
        with pytest.raises((ValueError, RuntimeError)):
            slicot.tb01wd(1, 1, -1, np.zeros((1,1), order='F'),
                         np.zeros((1,1), order='F'), np.zeros((1,1), order='F'))
