"""
Tests for SB10JD - Convert descriptor state-space to regular state-space.

SB10JD converts the descriptor state-space system:
    E*dx/dt = A*x + B*u
         y = C*x + D*u
into regular state-space form:
    dx/dt = Ad*x + Bd*u
        y = Cd*x + Dd*u

Uses SVD decomposition for descriptor elimination.
"""
import numpy as np
import pytest
from slicot import sb10jd


class TestSB10JDBasic:
    """Basic functionality tests."""

    def test_identity_descriptor(self):
        """
        Test with E = I (identity). Output should match input.

        When E is identity, the descriptor system is already in regular form.
        """
        n, m, np_ = 3, 2, 2

        a = np.array([
            [-1.0, 0.5, 0.2],
            [0.0, -2.0, 0.3],
            [0.0, 0.0, -3.0]
        ], order='F', dtype=np.float64)

        b = np.array([
            [1.0, 0.0],
            [0.5, 1.0],
            [0.0, 0.5]
        ], order='F', dtype=np.float64)

        c = np.array([
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0]
        ], order='F', dtype=np.float64)

        d = np.array([
            [0.0, 0.0],
            [0.0, 0.0]
        ], order='F', dtype=np.float64)

        e = np.eye(n, order='F', dtype=np.float64)

        a_orig = a.copy()
        b_orig = b.copy()
        c_orig = c.copy()
        d_orig = d.copy()

        ad, bd, cd, dd, nsys, info = sb10jd(a, b, c, d, e)

        assert info == 0
        assert nsys == n
        np.testing.assert_allclose(ad[:nsys, :nsys], a_orig, rtol=1e-12)
        np.testing.assert_allclose(bd[:nsys, :], b_orig, rtol=1e-12)
        np.testing.assert_allclose(cd[:, :nsys], c_orig, rtol=1e-12)
        np.testing.assert_allclose(dd, d_orig, rtol=1e-12)

    def test_scaled_identity_descriptor(self):
        """
        Test with E = 2*I (scaled identity).

        When E = 2*I, the regular form is:
            dx/dt = (1/2)*A*x + (1/2)*B*u
        Random seed: 42 (for reproducibility)
        """
        np.random.seed(42)
        n, m, np_ = 3, 2, 2

        a = np.array([
            [-2.0, 1.0, 0.0],
            [0.0, -4.0, 1.0],
            [0.0, 0.0, -6.0]
        ], order='F', dtype=np.float64)

        b = np.array([
            [2.0, 0.0],
            [0.0, 2.0],
            [1.0, 1.0]
        ], order='F', dtype=np.float64)

        c = np.array([
            [1.0, 0.0, 1.0],
            [0.0, 1.0, 0.0]
        ], order='F', dtype=np.float64)

        d = np.array([
            [0.1, 0.0],
            [0.0, 0.1]
        ], order='F', dtype=np.float64)

        e = 2.0 * np.eye(n, order='F', dtype=np.float64)

        ad, bd, cd, dd, nsys, info = sb10jd(a, b, c, d, e)

        assert info == 0
        assert nsys == n

        a_expected = 0.5 * a
        b_expected = 0.5 * b

        np.testing.assert_allclose(ad[:nsys, :nsys], a_expected, rtol=1e-10)
        np.testing.assert_allclose(bd[:nsys, :], b_expected, rtol=1e-10)


class TestSB10JDRankDeficient:
    """Tests for rank-deficient E matrices (descriptor elimination)."""

    def test_rank_deficient_e_simple(self):
        """
        Test with rank-deficient E (nsys < n).

        E = [[1, 0], [0, 0]] has rank 1, so nsys should be 1.
        """
        n, m, np_ = 2, 1, 1

        a = np.array([
            [1.0, 2.0],
            [3.0, 4.0]
        ], order='F', dtype=np.float64)

        b = np.array([
            [1.0],
            [1.0]
        ], order='F', dtype=np.float64)

        c = np.array([
            [1.0, 0.0]
        ], order='F', dtype=np.float64)

        d = np.array([
            [0.0]
        ], order='F', dtype=np.float64)

        e = np.array([
            [1.0, 0.0],
            [0.0, 0.0]
        ], order='F', dtype=np.float64)

        ad, bd, cd, dd, nsys, info = sb10jd(a, b, c, d, e)

        assert info == 0
        assert nsys == 1

    def test_rank_deficient_e_2x3(self):
        """
        Test rank-2 E with n=3.

        E matrix with rank 2 should produce nsys=2.
        """
        n, m, np_ = 3, 2, 2

        a = np.array([
            [-1.0, 1.0, 0.0],
            [0.0, -2.0, 1.0],
            [1.0, 0.0, -3.0]
        ], order='F', dtype=np.float64)

        b = np.array([
            [1.0, 0.0],
            [0.0, 1.0],
            [0.5, 0.5]
        ], order='F', dtype=np.float64)

        c = np.array([
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0]
        ], order='F', dtype=np.float64)

        d = np.zeros((np_, m), order='F', dtype=np.float64)

        e = np.array([
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0]
        ], order='F', dtype=np.float64)

        ad, bd, cd, dd, nsys, info = sb10jd(a, b, c, d, e)

        assert info == 0
        assert nsys == 2


class TestSB10JDTransferFunction:
    """Test transfer function preservation property."""

    def test_transfer_function_preservation(self):
        """
        Verify transfer function is preserved.

        G(s) = C*(s*E - A)^(-1)*B + D should equal
        Gd(s) = Cd*(s*I - Ad)^(-1)*Bd + Dd at test frequencies.

        Random seed: 100 (for reproducibility)
        """
        np.random.seed(100)
        n, m, np_ = 3, 2, 2

        a = np.array([
            [-2.0, 0.5, 0.0],
            [0.0, -3.0, 0.5],
            [0.0, 0.0, -4.0]
        ], order='F', dtype=np.float64)

        b = np.array([
            [1.0, 0.0],
            [0.5, 1.0],
            [0.0, 0.5]
        ], order='F', dtype=np.float64)

        c = np.array([
            [1.0, 0.0, 0.5],
            [0.0, 1.0, 0.0]
        ], order='F', dtype=np.float64)

        d = np.array([
            [0.1, 0.0],
            [0.0, 0.1]
        ], order='F', dtype=np.float64)

        e = np.eye(n, order='F', dtype=np.float64)

        a_orig = a.copy()
        b_orig = b.copy()
        c_orig = c.copy()
        d_orig = d.copy()
        e_orig = e.copy()

        ad, bd, cd, dd, nsys, info = sb10jd(a, b, c, d, e)

        assert info == 0

        test_freqs = [0.1, 1.0, 10.0, 100.0]
        for omega in test_freqs:
            s = 1j * omega

            g_desc = c_orig @ np.linalg.solve(s * e_orig - a_orig, b_orig) + d_orig

            ad_sub = ad[:nsys, :nsys]
            bd_sub = bd[:nsys, :]
            cd_sub = cd[:, :nsys]
            g_regular = cd_sub @ np.linalg.solve(s * np.eye(nsys) - ad_sub, bd_sub) + dd

            np.testing.assert_allclose(g_regular, g_desc, rtol=1e-10)


class TestSB10JDZeroE:
    """Test with zero or near-zero E matrix."""

    def test_zero_e_matrix(self):
        """
        Test with E = 0 (all singular values below tolerance).

        When E = 0, A becomes -1/eps*I on diagonal, B and C become zero.
        """
        n, m, np_ = 2, 1, 1

        a = np.array([
            [1.0, 2.0],
            [3.0, 4.0]
        ], order='F', dtype=np.float64)

        b = np.array([
            [1.0],
            [1.0]
        ], order='F', dtype=np.float64)

        c = np.array([
            [1.0, 0.0]
        ], order='F', dtype=np.float64)

        d = np.array([
            [0.5]
        ], order='F', dtype=np.float64)

        e = np.zeros((n, n), order='F', dtype=np.float64)

        ad, bd, cd, dd, nsys, info = sb10jd(a, b, c, d, e)

        assert info == 0
        assert nsys == n

        np.testing.assert_allclose(bd[:nsys, :], 0.0, atol=1e-10)
        np.testing.assert_allclose(cd[:, :nsys], 0.0, atol=1e-10)


class TestSB10JDEdgeCases:
    """Edge cases and boundary conditions."""

    def test_n_zero(self):
        """Test N=0 quick return."""
        n, m, np_ = 0, 2, 2

        a = np.zeros((1, 1), order='F', dtype=np.float64)
        b = np.zeros((1, m), order='F', dtype=np.float64)
        c = np.zeros((np_, 1), order='F', dtype=np.float64)
        d = np.zeros((np_, m), order='F', dtype=np.float64)
        e = np.zeros((1, 1), order='F', dtype=np.float64)

        ad, bd, cd, dd, nsys, info = sb10jd(a, b, c, d, e, n=0, m=m, np=np_)

        assert info == 0
        assert nsys == 0

    def test_m_zero(self):
        """Test M=0 (no inputs)."""
        n, m, np_ = 2, 0, 1

        a = np.array([
            [-1.0, 0.5],
            [0.0, -2.0]
        ], order='F', dtype=np.float64)

        b = np.zeros((n, 1), order='F', dtype=np.float64)
        c = np.array([[1.0, 0.0]], order='F', dtype=np.float64)
        d = np.zeros((np_, 1), order='F', dtype=np.float64)
        e = np.eye(n, order='F', dtype=np.float64)

        ad, bd, cd, dd, nsys, info = sb10jd(a, b, c, d, e, m=0)

        assert info == 0
        assert nsys == n

    def test_np_zero(self):
        """Test NP=0 (no outputs)."""
        n, m, np_ = 2, 1, 0

        a = np.array([
            [-1.0, 0.5],
            [0.0, -2.0]
        ], order='F', dtype=np.float64)

        b = np.array([[1.0], [0.5]], order='F', dtype=np.float64)
        c = np.zeros((1, n), order='F', dtype=np.float64)
        d = np.zeros((1, m), order='F', dtype=np.float64)
        e = np.eye(n, order='F', dtype=np.float64)

        ad, bd, cd, dd, nsys, info = sb10jd(a, b, c, d, e, np=0)

        assert info == 0
        assert nsys == n


class TestSB10JDErrorHandling:
    """Test error handling."""

    def test_negative_n(self):
        """Test negative N returns error."""
        a = np.zeros((1, 1), order='F', dtype=np.float64)
        b = np.zeros((1, 1), order='F', dtype=np.float64)
        c = np.zeros((1, 1), order='F', dtype=np.float64)
        d = np.zeros((1, 1), order='F', dtype=np.float64)
        e = np.zeros((1, 1), order='F', dtype=np.float64)

        with pytest.raises(ValueError):
            sb10jd(a, b, c, d, e, n=-1)

    def test_negative_m(self):
        """Test negative M returns error."""
        a = np.zeros((2, 2), order='F', dtype=np.float64)
        b = np.zeros((2, 1), order='F', dtype=np.float64)
        c = np.zeros((1, 2), order='F', dtype=np.float64)
        d = np.zeros((1, 1), order='F', dtype=np.float64)
        e = np.zeros((2, 2), order='F', dtype=np.float64)

        with pytest.raises(ValueError):
            sb10jd(a, b, c, d, e, m=-1)


class TestSB10JDLargerSystems:
    """Test larger systems for robustness."""

    def test_5x5_full_rank(self):
        """
        Test 5x5 system with full rank E.

        Random seed: 200 (for reproducibility)
        """
        np.random.seed(200)
        n, m, np_ = 5, 3, 2

        a = -np.eye(n) + 0.1 * np.random.randn(n, n)
        a = np.asfortranarray(a)

        b = np.random.randn(n, m)
        b = np.asfortranarray(b)

        c = np.random.randn(np_, n)
        c = np.asfortranarray(c)

        d = 0.1 * np.random.randn(np_, m)
        d = np.asfortranarray(d)

        e = np.eye(n) + 0.01 * np.random.randn(n, n)
        e = np.asfortranarray(e)

        a_orig = a.copy()
        b_orig = b.copy()
        c_orig = c.copy()
        d_orig = d.copy()
        e_orig = e.copy()

        ad, bd, cd, dd, nsys, info = sb10jd(a, b, c, d, e)

        assert info == 0
        assert nsys == n

        s = 1j * 1.0
        g_desc = c_orig @ np.linalg.solve(s * e_orig - a_orig, b_orig) + d_orig
        ad_sub = ad[:nsys, :nsys]
        bd_sub = bd[:nsys, :]
        cd_sub = cd[:, :nsys]
        g_regular = cd_sub @ np.linalg.solve(s * np.eye(nsys) - ad_sub, bd_sub) + dd

        np.testing.assert_allclose(g_regular, g_desc, rtol=1e-9)

    def test_6x6_rank_deficient(self):
        """
        Test 6x6 system with rank-4 E.

        Random seed: 300 (for reproducibility)
        """
        np.random.seed(300)
        n, m, np_ = 6, 2, 3

        a = np.random.randn(n, n)
        a = np.asfortranarray(a)

        b = np.random.randn(n, m)
        b = np.asfortranarray(b)

        c = np.random.randn(np_, n)
        c = np.asfortranarray(c)

        d = np.random.randn(np_, m)
        d = np.asfortranarray(d)

        u_e = np.linalg.qr(np.random.randn(n, n))[0]
        s_e = np.diag([2.0, 1.5, 1.0, 0.5, 0.0, 0.0])
        e = u_e @ s_e @ u_e.T
        e = np.asfortranarray(e)

        ad, bd, cd, dd, nsys, info = sb10jd(a, b, c, d, e)

        assert info == 0
        assert nsys == 4
