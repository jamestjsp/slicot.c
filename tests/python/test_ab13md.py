"""
Tests for AB13MD - Upper bound on structured singular value (mu).

Computes an upper bound on the structured singular value for a given
square complex matrix and a given block structure of the uncertainty.

The structured singular value (mu) is used for robustness analysis in
control systems with mixed parametric uncertainty and unmodeled dynamics.
"""

import numpy as np
import pytest
from numpy.testing import assert_allclose


def test_ab13md_basic():
    """
    Test basic functionality using SLICOT HTML doc example.

    Input:
    - N = 6 (matrix order)
    - M = 5 (number of blocks)
    - NBLOCK = [1, 1, 2, 1, 1] (block sizes)
    - ITYPE = [1, 1, 2, 2, 2] (block types: 1=real, 2=complex)
    - Z = 6x6 complex matrix from example data

    Expected output:
    - BOUND = 41.74753408 (upper bound on mu)
    - INFO = 0 (success)
    """
    from slicot import ab13md

    z = np.array([
        [-1.0+6.0j,   2.0-3.0j,   3.0+8.0j,  -1.0+6.0j,   2.0-3.0j,   3.0+8.0j],
        [ 3.0+8.0j,  -5.0-9.0j,  -6.0+2.0j,   3.0+8.0j,  -5.0-9.0j,  -6.0+2.0j],
        [ 4.0+2.0j,  -2.0+5.0j,  -6.0-7.0j,   4.0+2.0j,  -2.0+5.0j,  -6.0-7.0j],
        [-4.0+11.0j,  8.0-7.0j,  12.0-1.0j,  -4.0+11.0j,  8.0-7.0j,  12.0-1.0j],
        [ 5.0-4.0j,  -4.0-8.0j,   1.0-3.0j,   5.0-4.0j,  -4.0-8.0j,   1.0-3.0j],
        [-6.0+14.0j,  2.0-5.0j,   4.0+16.0j, -6.0+14.0j,  2.0-5.0j,   4.0+16.0j]
    ], order='F', dtype=np.complex128)

    nblock = np.array([1, 1, 2, 1, 1], dtype=np.int32)
    itype = np.array([1, 1, 2, 2, 2], dtype=np.int32)

    bound, d, g, x, info = ab13md(z, nblock, itype)

    assert info == 0
    assert_allclose(bound, 41.74753408, rtol=1e-3, atol=1e-4)
    assert d.shape == (6,)
    assert g.shape == (6,)
    assert all(d > 0)


def test_ab13md_single_complex_block():
    """
    Test with single N-by-N complex block (exact: largest singular value).

    When NBLOCK = [N] and ITYPE = [2] (single complex block),
    mu(Z) = sigma_max(Z), the largest singular value.

    Random seed: 42 (for reproducibility)
    """
    from slicot import ab13md

    np.random.seed(42)
    n = 4

    z_real = np.random.randn(n, n)
    z_imag = np.random.randn(n, n)
    z = (z_real + 1j * z_imag).astype(np.complex128, order='F')

    nblock = np.array([n], dtype=np.int32)
    itype = np.array([2], dtype=np.int32)

    bound, d, g, x, info = ab13md(z, nblock, itype)

    assert info == 0

    sigma_max = np.linalg.svd(z, compute_uv=False)[0]
    assert_allclose(bound, sigma_max, rtol=1e-10, atol=1e-12)


def test_ab13md_with_fact():
    """
    Test FACT = 'F' mode (reuse X from previous call).

    Calling with similar matrices should work with FACT='F'.
    Random seed: 123 (for reproducibility)
    """
    from slicot import ab13md

    np.random.seed(123)
    n = 4
    m = 2

    z_real = np.random.randn(n, n)
    z_imag = np.random.randn(n, n)
    z = (z_real + 1j * z_imag).astype(np.complex128, order='F')

    nblock = np.array([2, 2], dtype=np.int32)
    itype = np.array([2, 2], dtype=np.int32)

    bound1, d1, g1, x1, info1 = ab13md(z, nblock, itype, fact='N')
    assert info1 == 0

    z_perturbed = z + 0.01 * (np.random.randn(n, n) + 1j * np.random.randn(n, n))
    z_perturbed = z_perturbed.astype(np.complex128, order='F')

    bound2, d2, g2, x2, info2 = ab13md(z_perturbed, nblock, itype, fact='F', x=x1)
    assert info2 == 0

    assert bound2 > 0


def test_ab13md_mixed_blocks():
    """
    Test with mixed real and complex blocks.

    Random seed: 456 (for reproducibility)
    """
    from slicot import ab13md

    np.random.seed(456)
    n = 5
    nblock = np.array([1, 1, 1, 2], dtype=np.int32)
    itype = np.array([1, 2, 1, 2], dtype=np.int32)

    z_real = np.random.randn(n, n)
    z_imag = np.random.randn(n, n)
    z = (z_real + 1j * z_imag).astype(np.complex128, order='F')

    bound, d, g, x, info = ab13md(z, nblock, itype)

    assert info == 0
    assert bound > 0
    assert d.shape == (n,)
    assert g.shape == (n,)


def test_ab13md_block_size_mismatch_error():
    """
    Test error handling: sum of block sizes != N.

    Should return INFO = 2.
    """
    from slicot import ab13md

    n = 6
    z = np.zeros((n, n), order='F', dtype=np.complex128)
    nblock = np.array([1, 2, 2], dtype=np.int32)
    itype = np.array([2, 2, 2], dtype=np.int32)

    bound, d, g, x, info = ab13md(z, nblock, itype)
    assert info == 2


def test_ab13md_real_block_size_error():
    """
    Test error handling: real block with size > 1.

    Real blocks (ITYPE=1) must have size 1. Should return INFO = 3.
    """
    from slicot import ab13md

    n = 4
    z = np.zeros((n, n), order='F', dtype=np.complex128)
    nblock = np.array([2, 2], dtype=np.int32)
    itype = np.array([1, 2], dtype=np.int32)

    bound, d, g, x, info = ab13md(z, nblock, itype)
    assert info == 3


def test_ab13md_invalid_block_type_error():
    """
    Test error handling: invalid block type.

    Block type must be 1 or 2. Should return INFO = 4.
    """
    from slicot import ab13md

    n = 4
    z = np.zeros((n, n), order='F', dtype=np.complex128)
    nblock = np.array([2, 2], dtype=np.int32)
    itype = np.array([2, 3], dtype=np.int32)

    bound, d, g, x, info = ab13md(z, nblock, itype)
    assert info == 4


def test_ab13md_zero_block_size_error():
    """
    Test error handling: zero or negative block size.

    Block sizes must be positive. Should return INFO = 1.
    """
    from slicot import ab13md

    n = 4
    z = np.zeros((n, n), order='F', dtype=np.complex128)
    nblock = np.array([0, 2, 2], dtype=np.int32)
    itype = np.array([2, 2, 2], dtype=np.int32)

    bound, d, g, x, info = ab13md(z, nblock, itype)
    assert info == 1


def test_ab13md_zero_matrix():
    """
    Test edge case: zero matrix.

    mu(0) = 0 for any block structure.
    """
    from slicot import ab13md

    n = 4
    z = np.zeros((n, n), order='F', dtype=np.complex128)
    nblock = np.array([2, 2], dtype=np.int32)
    itype = np.array([2, 2], dtype=np.int32)

    bound, d, g, x, info = ab13md(z, nblock, itype)

    assert info == 0
    assert bound == 0.0


def test_ab13md_hermitian_matrix():
    """
    Test with Hermitian matrix.

    Random seed: 789 (for reproducibility)
    """
    from slicot import ab13md

    np.random.seed(789)
    n = 4

    z_temp = np.random.randn(n, n) + 1j * np.random.randn(n, n)
    z = ((z_temp + z_temp.conj().T) / 2).astype(np.complex128, order='F')

    nblock = np.array([2, 2], dtype=np.int32)
    itype = np.array([2, 2], dtype=np.int32)

    bound, d, g, x, info = ab13md(z, nblock, itype)

    assert info == 0
    assert bound > 0


def test_ab13md_bound_nonnegative():
    """
    Test mathematical property: mu bound is always non-negative.

    Random seed: 888 (for reproducibility)
    """
    from slicot import ab13md

    np.random.seed(888)

    for _ in range(5):
        n = np.random.randint(2, 6)
        m = np.random.randint(1, n)
        sizes = []
        remaining = n
        for i in range(m - 1):
            size = np.random.randint(1, remaining - (m - i - 1) + 1)
            sizes.append(size)
            remaining -= size
        sizes.append(remaining)

        nblock = np.array(sizes, dtype=np.int32)
        itype = np.random.choice([1, 2], size=m).astype(np.int32)
        for i in range(m):
            if itype[i] == 1 and nblock[i] > 1:
                itype[i] = 2

        z = (np.random.randn(n, n) + 1j * np.random.randn(n, n)).astype(np.complex128, order='F')

        bound, d, g, x, info = ab13md(z, nblock, itype)

        if info == 0:
            assert bound >= 0


def test_ab13md_scaling_property():
    """
    Test mathematical property: mu(alpha * Z) = |alpha| * mu(Z).

    Random seed: 999 (for reproducibility)
    """
    from slicot import ab13md

    np.random.seed(999)
    n = 4

    z = (np.random.randn(n, n) + 1j * np.random.randn(n, n)).astype(np.complex128, order='F')
    nblock = np.array([2, 2], dtype=np.int32)
    itype = np.array([2, 2], dtype=np.int32)

    bound1, d1, g1, x1, info1 = ab13md(z, nblock, itype)
    assert info1 == 0

    alpha = 2.5 + 1.5j
    z_scaled = (alpha * z).astype(np.complex128, order='F')

    bound2, d2, g2, x2, info2 = ab13md(z_scaled, nblock, itype)
    assert info2 == 0

    assert_allclose(bound2, abs(alpha) * bound1, rtol=1e-3, atol=1e-6)
