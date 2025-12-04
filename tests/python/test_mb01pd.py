"""
Tests for MB01PD: Matrix scaling to safe numerical range.

MB01PD scales a matrix so its norm is in [SMLNUM, BIGNUM] range, or undoes
such scaling.
"""

import numpy as np
import pytest


class TestMB01PDBasic:
    """Basic functionality tests for MB01PD matrix scaling."""

    def test_scale_full_matrix_no_change_needed(self):
        """
        Test scaling a full matrix whose norm is already in safe range.

        When ANRM is in [SMLNUM, BIGNUM], no scaling should occur.
        Random seed: 42 (for reproducibility)
        """
        from slicot import mb01pd

        np.random.seed(42)
        m, n = 4, 3
        a = np.random.randn(m, n).astype(float, order='F')
        a_original = a.copy()
        anrm = np.linalg.norm(a, ord='fro')

        a_out, info = mb01pd('S', 'G', m, n, 0, 0, anrm, 0, None, a)

        assert info == 0
        np.testing.assert_allclose(a_out, a_original, rtol=1e-14)

    def test_scale_full_matrix_very_small_norm(self):
        """
        Test scaling a full matrix with very small norm up to SMLNUM.

        When ANRM < SMLNUM, matrix is scaled by SMLNUM/ANRM.
        Random seed: 123 (for reproducibility)
        """
        from slicot import mb01pd

        np.random.seed(123)
        m, n = 3, 3
        a = np.random.randn(m, n).astype(float, order='F')
        a = a * 1e-310
        a_original = a.copy()
        anrm = np.linalg.norm(a, ord='fro')

        a_out, info = mb01pd('S', 'G', m, n, 0, 0, anrm, 0, None, a)

        assert info == 0
        a_norm_after = np.linalg.norm(a_out, ord='fro')
        assert a_norm_after > anrm or anrm == 0.0

    def test_scale_full_matrix_very_large_norm(self):
        """
        Test scaling a full matrix with very large norm down to BIGNUM.

        When ANRM > BIGNUM, matrix is scaled by BIGNUM/ANRM.
        Random seed: 456 (for reproducibility)
        """
        from slicot import mb01pd

        np.random.seed(456)
        m, n = 3, 3
        a = np.random.randn(m, n).astype(float, order='F')
        a = a * 1e300
        a_original = a.copy()
        anrm = np.linalg.norm(a, ord='fro')

        a_out, info = mb01pd('S', 'G', m, n, 0, 0, anrm, 0, None, a)

        assert info == 0
        a_norm_after = np.linalg.norm(a_out, ord='fro')
        assert a_norm_after < anrm

    def test_undo_scaling(self):
        """
        Test undoing scaling: scale then undo should restore original.

        Mathematical property: (scale then undo) is identity.
        Random seed: 789 (for reproducibility)
        """
        from slicot import mb01pd

        np.random.seed(789)
        m, n = 4, 4
        a = np.random.randn(m, n).astype(float, order='F')
        a = a * 1e-310
        a_original = a.copy()
        anrm = np.linalg.norm(a, ord='fro')

        if anrm == 0.0:
            pytest.skip("Zero norm matrix - no scaling needed")

        a_scaled, info = mb01pd('S', 'G', m, n, 0, 0, anrm, 0, None, a)
        assert info == 0

        a_restored, info = mb01pd('U', 'G', m, n, 0, 0, anrm, 0, None, a_scaled)
        assert info == 0

        np.testing.assert_allclose(a_restored, a_original, rtol=1e-14)


class TestMB01PDTriangular:
    """Tests for triangular matrix types."""

    def test_lower_triangular(self):
        """
        Test scaling a lower triangular matrix.

        Only lower triangle is scaled; upper part should remain unchanged.
        Random seed: 111 (for reproducibility)
        """
        from slicot import mb01pd

        np.random.seed(111)
        n = 4
        a = np.tril(np.random.randn(n, n)).astype(float, order='F')
        a = a * 1e-310
        a_original = a.copy()
        anrm = np.linalg.norm(a, ord='fro')

        if anrm == 0.0:
            pytest.skip("Zero norm matrix")

        a_scaled, info = mb01pd('S', 'L', n, n, 0, 0, anrm, 0, None, a)
        assert info == 0

        a_restored, info = mb01pd('U', 'L', n, n, 0, 0, anrm, 0, None, a_scaled)
        assert info == 0

        np.testing.assert_allclose(a_restored, a_original, rtol=1e-14)

    def test_upper_triangular(self):
        """
        Test scaling an upper triangular matrix.

        Only upper triangle is scaled.
        Random seed: 222 (for reproducibility)
        """
        from slicot import mb01pd

        np.random.seed(222)
        n = 4
        a = np.triu(np.random.randn(n, n)).astype(float, order='F')
        a = a * 1e-310
        a_original = a.copy()
        anrm = np.linalg.norm(a, ord='fro')

        if anrm == 0.0:
            pytest.skip("Zero norm matrix")

        a_scaled, info = mb01pd('S', 'U', n, n, 0, 0, anrm, 0, None, a)
        assert info == 0

        a_restored, info = mb01pd('U', 'U', n, n, 0, 0, anrm, 0, None, a_scaled)
        assert info == 0

        np.testing.assert_allclose(a_restored, a_original, rtol=1e-14)

    def test_upper_hessenberg(self):
        """
        Test scaling an upper Hessenberg matrix.

        Upper Hessenberg = upper triangular + one subdiagonal.
        Random seed: 333 (for reproducibility)
        """
        from slicot import mb01pd

        np.random.seed(333)
        n = 4
        a = np.triu(np.random.randn(n, n), -1).astype(float, order='F')
        a = a * 1e-310
        a_original = a.copy()
        anrm = np.linalg.norm(a, ord='fro')

        if anrm == 0.0:
            pytest.skip("Zero norm matrix")

        a_scaled, info = mb01pd('S', 'H', n, n, 0, 0, anrm, 0, None, a)
        assert info == 0

        a_restored, info = mb01pd('U', 'H', n, n, 0, 0, anrm, 0, None, a_scaled)
        assert info == 0

        np.testing.assert_allclose(a_restored, a_original, rtol=1e-14)


class TestMB01PDBlock:
    """Tests for block structured matrices."""

    def test_block_lower_triangular(self):
        """
        Test scaling a block lower triangular matrix.

        Random seed: 444 (for reproducibility)
        """
        from slicot import mb01pd

        np.random.seed(444)
        n = 6
        nbl = 2
        nrows = np.array([3, 3], dtype=np.int32)
        a = np.zeros((n, n), dtype=float, order='F')
        a[:3, :3] = np.random.randn(3, 3)
        a[3:, :3] = np.random.randn(3, 3)
        a[3:, 3:] = np.random.randn(3, 3)
        a = a * 1e-310
        a_original = a.copy()
        anrm = np.linalg.norm(a, ord='fro')

        if anrm == 0.0:
            pytest.skip("Zero norm matrix")

        a_scaled, info = mb01pd('S', 'L', n, n, 0, 0, anrm, nbl, nrows, a)
        assert info == 0

        a_restored, info = mb01pd('U', 'L', n, n, 0, 0, anrm, nbl, nrows, a_scaled)
        assert info == 0

        np.testing.assert_allclose(a_restored, a_original, rtol=1e-14)

    def test_block_upper_triangular(self):
        """
        Test scaling a block upper triangular matrix.

        Random seed: 555 (for reproducibility)
        """
        from slicot import mb01pd

        np.random.seed(555)
        n = 6
        nbl = 2
        nrows = np.array([3, 3], dtype=np.int32)
        a = np.zeros((n, n), dtype=float, order='F')
        a[:3, :3] = np.random.randn(3, 3)
        a[:3, 3:] = np.random.randn(3, 3)
        a[3:, 3:] = np.random.randn(3, 3)
        a = a * 1e-310
        a_original = a.copy()
        anrm = np.linalg.norm(a, ord='fro')

        if anrm == 0.0:
            pytest.skip("Zero norm matrix")

        a_scaled, info = mb01pd('S', 'U', n, n, 0, 0, anrm, nbl, nrows, a)
        assert info == 0

        a_restored, info = mb01pd('U', 'U', n, n, 0, 0, anrm, nbl, nrows, a_scaled)
        assert info == 0

        np.testing.assert_allclose(a_restored, a_original, rtol=1e-14)


class TestMB01PDEdgeCases:
    """Edge case tests for MB01PD."""

    def test_zero_dimensions(self):
        """Test with zero dimensions - quick return."""
        from slicot import mb01pd

        a = np.array([[1.0]], dtype=float, order='F')
        a_out, info = mb01pd('S', 'G', 0, 0, 0, 0, 0.0, 0, None, a)
        assert info == 0

    def test_zero_norm(self):
        """Test with zero norm - quick return."""
        from slicot import mb01pd

        a = np.zeros((3, 3), dtype=float, order='F')
        a_original = a.copy()

        a_out, info = mb01pd('S', 'G', 3, 3, 0, 0, 0.0, 0, None, a)
        assert info == 0
        np.testing.assert_array_equal(a_out, a_original)

    def test_m_not_equal_n_full(self):
        """
        Test with non-square full matrix.

        Random seed: 666 (for reproducibility)
        """
        from slicot import mb01pd

        np.random.seed(666)
        m, n = 5, 3
        a = np.random.randn(m, n).astype(float, order='F')
        a = a * 1e-310
        a_original = a.copy()
        anrm = np.linalg.norm(a, ord='fro')

        if anrm == 0.0:
            pytest.skip("Zero norm matrix")

        a_scaled, info = mb01pd('S', 'G', m, n, 0, 0, anrm, 0, None, a)
        assert info == 0

        a_restored, info = mb01pd('U', 'G', m, n, 0, 0, anrm, 0, None, a_scaled)
        assert info == 0

        np.testing.assert_allclose(a_restored, a_original, rtol=1e-14)


class TestMB01PDErrors:
    """Error handling tests for MB01PD."""

    def test_invalid_scun(self):
        """Test error for invalid SCUN parameter."""
        from slicot import mb01pd

        a = np.zeros((3, 3), dtype=float, order='F')
        a_out, info = mb01pd('X', 'G', 3, 3, 0, 0, 1.0, 0, None, a)
        assert info == -1

    def test_invalid_type(self):
        """Test error for invalid TYPE parameter."""
        from slicot import mb01pd

        a = np.zeros((3, 3), dtype=float, order='F')
        a_out, info = mb01pd('S', 'X', 3, 3, 0, 0, 1.0, 0, None, a)
        assert info == -2

    def test_negative_m(self):
        """Test error for negative M."""
        from slicot import mb01pd

        a = np.zeros((3, 3), dtype=float, order='F')
        a_out, info = mb01pd('S', 'G', -1, 3, 0, 0, 1.0, 0, None, a)
        assert info == -3

    def test_negative_n(self):
        """Test error for negative N."""
        from slicot import mb01pd

        a = np.zeros((3, 3), dtype=float, order='F')
        a_out, info = mb01pd('S', 'G', 3, -1, 0, 0, 1.0, 0, None, a)
        assert info == -4

    def test_negative_anrm(self):
        """Test error for negative ANRM."""
        from slicot import mb01pd

        a = np.zeros((3, 3), dtype=float, order='F')
        a_out, info = mb01pd('S', 'G', 3, 3, 0, 0, -1.0, 0, None, a)
        assert info == -7

    def test_negative_nbl(self):
        """Test error for negative NBL."""
        from slicot import mb01pd

        a = np.zeros((3, 3), dtype=float, order='F')
        a_out, info = mb01pd('S', 'G', 3, 3, 0, 0, 1.0, -1, None, a)
        assert info == -8


class TestMB01PDMathematicalProperties:
    """Mathematical property validation tests."""

    def test_scaling_preserves_relative_structure(self):
        """
        Test that scaling preserves relative element magnitudes.

        Property: A_scaled / norm(A_scaled) should equal A / norm(A).
        Random seed: 888 (for reproducibility)
        """
        from slicot import mb01pd

        np.random.seed(888)
        m, n = 4, 4
        a = np.random.randn(m, n).astype(float, order='F')
        a = a * 1e-310
        a_original = a.copy()
        anrm = np.linalg.norm(a, ord='fro')

        if anrm == 0.0:
            pytest.skip("Zero norm matrix")

        a_scaled, info = mb01pd('S', 'G', m, n, 0, 0, anrm, 0, None, a)
        assert info == 0

        a_scaled_norm = np.linalg.norm(a_scaled, ord='fro')
        if a_scaled_norm > 0:
            np.testing.assert_allclose(
                a_scaled / a_scaled_norm,
                a_original / anrm,
                rtol=1e-14
            )

    def test_scaling_is_linear(self):
        """
        Test that scaling is a linear operation (multiplication by constant).

        Property: A_scaled = c * A for some constant c > 0.
        Random seed: 999 (for reproducibility)
        """
        from slicot import mb01pd

        np.random.seed(999)
        m, n = 3, 3
        a = np.random.randn(m, n).astype(float, order='F')
        a = a * 1e-310
        a_original = a.copy()
        anrm = np.linalg.norm(a, ord='fro')

        if anrm == 0.0:
            pytest.skip("Zero norm matrix")

        a_scaled, info = mb01pd('S', 'G', m, n, 0, 0, anrm, 0, None, a)
        assert info == 0

        nz_mask = a_original != 0
        if np.any(nz_mask):
            ratios = a_scaled[nz_mask] / a_original[nz_mask]
            np.testing.assert_allclose(ratios, ratios[0], rtol=1e-14)
            assert ratios[0] > 0

    def test_idempotence_when_in_range(self):
        """
        Test that scaling twice has same effect as scaling once when in range.

        Property: If A is already in range, scaling is identity.
        Random seed: 1010 (for reproducibility)
        """
        from slicot import mb01pd

        np.random.seed(1010)
        m, n = 4, 4
        a = np.random.randn(m, n).astype(float, order='F')
        a_original = a.copy()
        anrm = np.linalg.norm(a, ord='fro')

        a_scaled1, info = mb01pd('S', 'G', m, n, 0, 0, anrm, 0, None, a)
        assert info == 0

        np.testing.assert_allclose(a_scaled1, a_original, rtol=1e-14)
