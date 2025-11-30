"""
Tests for IB03BD - Wiener system identification using Levenberg-Marquardt algorithm.

IB03BD computes a set of parameters for approximating a Wiener system in a
least-squares sense using a neural network approach. The Wiener system consists
of a linear part (state-space model) and a static nonlinearity.

System structure:
    x(t+1) = A*x(t) + B*u(t)      (linear state-space)
    z(t)   = C*x(t) + D*u(t)
    y(t)   = f(z(t), wb(1:L))     (nonlinear function)

The parameter vector X is partitioned as X = (wb(1),...,wb(L), theta) where:
- wb(i): neural network parameters for nonlinear part
- theta: linear part parameters in output normal form

Mathematical property: The algorithm minimizes sum(||e(t)||^2) where e(t) = y(t) - Y(t).
"""

import numpy as np
import pytest
import os

# Test data from SLICOT reference example
# Parameters from IB03BD.dat:
# NOBR=10, M=1, L=1, NSMP=1024, N=4, NN=12, ITMAX1=500, ITMAX2=1000
# TOL1=1e-5, TOL2=1e-5, INIT='B'

# Expected results from IB03BD.res:
# IWARN = 12, residual = 0.2995840
# Iterations = 42, fcn_evals = 898, jac_evals = 295


def load_test_data():
    """Load test data from npz file."""
    data_path = os.path.join(os.path.dirname(__file__), 'data', 'ib03bd_test_data.npz')
    data = np.load(data_path)
    return data['u'], data['y'], data['seed']


def test_ib03bd_html_example():
    """
    Validate basic functionality using SLICOT HTML doc example.

    From IB03BD.html and IB03BD.dat:
    - Input: 1024 samples of input-output data for Wiener system
    - Uses INIT='B' to initialize both linear and nonlinear parts
    - NOBR=10 block rows, N=4 state order, NN=12 neurons

    Expected from IB03BD.res:
    - IWARN = 12 (convergence with warnings)
    - Residual norm = 0.2995840
    - 42 iterations, 898 function evaluations, 295 Jacobian evaluations

    Test validates numerical correctness of the optimization result.
    """
    from slicot import ib03bd

    u, y, seed = load_test_data()

    # Reshape to column vectors (NSMP x M and NSMP x L)
    u = u.reshape(-1, 1).astype(float, order='F')
    y = y.reshape(-1, 1).astype(float, order='F')

    # Parameters from example
    nobr = 10
    m = 1
    l = 1
    nsmp = 1024
    n = 4
    nn = 12
    itmax1 = 500
    itmax2 = 1000
    tol1 = 1e-5
    tol2 = 1e-5
    init = 'B'

    # BSN = NN*(L+2) + 1 = 12*3 + 1 = 37
    # LTHS = N*(L+M+1) + L*M = 4*3 + 1 = 13
    # NX = BSN*L + LTHS = 37 + 13 = 50
    bsn = nn * (l + 2) + 1
    lths = n * (l + m + 1) + l * m
    nx = bsn * l + lths

    # Seed for random initialization (from data file)
    dwork_seed = seed.astype(float, order='F')

    x, iwarn, info, dwork_out = ib03bd(
        init, nobr, m, l, nsmp, n, nn, itmax1, itmax2,
        u, y, tol1, tol2, dwork_seed
    )

    # Check success (INFO should be 0)
    assert info == 0, f"IB03BD failed with INFO={info}"

    # IWARN = 12 means:
    # - k=0 (no warnings from IB01AD/IB01BD/IB01CD)
    # - j=1 (initialization step converged with both reductions at most TOL1)
    # - i=2 (main optimization: relative error between iterates at most TOL2)
    assert iwarn == 12, f"Expected IWARN=12, got {iwarn}"

    # Check residual norm (sum of squares)
    residual = dwork_out[1]  # DWORK(2) contains residual
    expected_residual = 0.2995840
    assert abs(residual - expected_residual) < 1e-3, \
        f"Residual mismatch: got {residual}, expected {expected_residual}"

    # Check number of iterations
    iterations = int(dwork_out[2])  # DWORK(3) contains iterations
    assert iterations == 42, f"Expected 42 iterations, got {iterations}"

    # Check parameter vector length
    assert len(x) == nx, f"Expected X length {nx}, got {len(x)}"


def test_ib03bd_expected_solution():
    """
    Validate that computed solution matches expected values from IB03BD.res.

    Expected final solution (first 10 values):
    14.1294, 1.1232, 6.4322, -11.2418, 7.6380, -33.4730, -64.7203, 747.1515, -0.4623, -92.6092

    Random seed: Uses SLICOT seed [1998, 1999, 2000, 2001] for reproducibility.
    """
    from slicot import ib03bd

    u, y, seed = load_test_data()

    u = u.reshape(-1, 1).astype(float, order='F')
    y = y.reshape(-1, 1).astype(float, order='F')

    nobr = 10
    m = 1
    l = 1
    nsmp = 1024
    n = 4
    nn = 12
    itmax1 = 500
    itmax2 = 1000
    tol1 = 1e-5
    tol2 = 1e-5
    init = 'B'

    dwork_seed = seed.astype(float, order='F')

    x, iwarn, info, dwork_out = ib03bd(
        init, nobr, m, l, nsmp, n, nn, itmax1, itmax2,
        u, y, tol1, tol2, dwork_seed
    )

    assert info == 0

    # Expected solution from IB03BD.res
    expected_x = np.array([
        14.1294, 1.1232, 6.4322, -11.2418, 7.6380, -33.4730, -64.7203, 747.1515, -0.4623, -92.6092,
        6.1682, -0.7672, 0.1194, 0.3558, 0.9091, 0.2948, 1.3465, 0.0093, 0.0560, -0.0035,
        -0.4179, -0.0455, -2.0871, -0.9196, 1.0777, 0.9213, 0.5373, 1.0412, -0.3978, 7.6832,
        -6.8614, -31.6119, -0.1092, -9.8984, 0.1257, 0.4056, 0.0472, 7.5819, -13.3969, 2.4869,
        -66.0727, -0.8411, -0.7040, 1.9641, 1.3059, -0.2046, -0.9326, 0.0040, 0.4032, 0.1479
    ])

    # Compare computed solution to expected
    # Use rtol=1e-3 due to precision in .res file (4 decimal places)
    np.testing.assert_allclose(x[:50], expected_x, rtol=1e-3, atol=1e-3)


def test_ib03bd_init_l_linear_only():
    """
    Test INIT='L' mode - initialize linear part only.

    When INIT='L', X must already contain initial parameters for the
    nonlinear part wb(1:L). The routine only initializes the linear part.

    Random seed: 42 (for reproducibility)
    """
    from slicot import ib03bd

    np.random.seed(42)

    u, y, _ = load_test_data()
    u = u.reshape(-1, 1).astype(float, order='F')
    y = y.reshape(-1, 1).astype(float, order='F')

    nobr = 10
    m = 1
    l = 1
    nsmp = 1024
    n = 4
    nn = 12
    itmax1 = 500
    itmax2 = 100  # Fewer iterations for faster test
    tol1 = 1e-4
    tol2 = 1e-4
    init = 'L'

    bsn = nn * (l + 2) + 1  # 37
    nths = bsn * l  # 37
    lths = n * (l + m + 1) + l * m  # 13
    nx = nths + lths  # 50

    # Initialize nonlinear part parameters (random small values)
    x_init = np.random.randn(nths).astype(float, order='F') * 0.1

    x, iwarn, info, dwork_out = ib03bd(
        init, nobr, m, l, nsmp, n, nn, itmax1, itmax2,
        u, y, tol1, tol2, x_init=x_init
    )

    # Should succeed or give a warning
    assert info == 0 or info > 0, f"Unexpected error INFO={info}"

    # Result should have correct size
    assert len(x) == nx, f"Expected X length {nx}, got {len(x)}"


def test_ib03bd_init_n_no_init():
    """
    Test INIT='N' mode - no initialization, use provided initial approximation.

    When INIT='N', X must already contain the complete initial parameters
    for both nonlinear and linear parts.

    Random seed: 123 (for reproducibility)
    """
    from slicot import ib03bd

    np.random.seed(123)

    u, y, seed = load_test_data()
    u = u.reshape(-1, 1).astype(float, order='F')
    y = y.reshape(-1, 1).astype(float, order='F')

    nobr = 10  # Not used when INIT='N', but needed for signature
    m = 1
    l = 1
    nsmp = 1024
    n = 4
    nn = 12
    itmax1 = 0  # Not used when INIT='N'
    itmax2 = 50  # Fewer iterations for faster test
    tol1 = 1e-4
    tol2 = 1e-4
    init = 'N'

    bsn = nn * (l + 2) + 1
    nths = bsn * l
    lths = n * (l + m + 1) + l * m
    nx = nths + lths

    # Initialize with random values
    x_init = np.random.randn(nx).astype(float, order='F') * 0.1

    x, iwarn, info, dwork_out = ib03bd(
        init, nobr, m, l, nsmp, n, nn, itmax1, itmax2,
        u, y, tol1, tol2, x_init=x_init
    )

    # Should succeed (might give warning for max iterations)
    assert info == 0, f"IB03BD failed with INFO={info}"

    assert len(x) == nx


def test_ib03bd_parameter_errors():
    """
    Test that invalid parameters are detected.

    IB03BD should return INFO < 0 for invalid arguments:
    - NOBR <= 0 when INIT='L' or 'B'
    - L <= 0 when INIT='L' or 'B'
    - N out of range
    - NN < 0
    """
    from slicot import ib03bd

    u, y, seed = load_test_data()
    u = u.reshape(-1, 1).astype(float, order='F')
    y = y.reshape(-1, 1).astype(float, order='F')

    dwork_seed = seed.astype(float, order='F')

    # Test NOBR <= 0
    with pytest.raises(Exception):
        x, iwarn, info, _ = ib03bd(
            'B', 0, 1, 1, 1024, 4, 12, 500, 1000,
            u, y, 1e-5, 1e-5, dwork_seed
        )


def test_ib03bd_convergence_property():
    """
    Validate convergence property: residual decreases through optimization.

    The Levenberg-Marquardt algorithm should decrease the sum of squared errors.

    Mathematical property:
    - Final residual should be less than initial residual
    - Residual = sum(||y(t) - Y(t)||^2)

    Random seed: Uses SLICOT seed for reproducibility.
    """
    from slicot import ib03bd

    u, y, seed = load_test_data()
    u = u.reshape(-1, 1).astype(float, order='F')
    y = y.reshape(-1, 1).astype(float, order='F')

    nobr = 10
    m = 1
    l = 1
    nsmp = 1024
    n = 4
    nn = 12
    itmax1 = 500
    itmax2 = 1000
    tol1 = 1e-5
    tol2 = 1e-5
    init = 'B'

    dwork_seed = seed.astype(float, order='F')

    x, iwarn, info, dwork_out = ib03bd(
        init, nobr, m, l, nsmp, n, nn, itmax1, itmax2,
        u, y, tol1, tol2, dwork_seed
    )

    assert info == 0

    # Final residual should be positive and finite
    residual = dwork_out[1]
    assert residual > 0, "Residual should be positive"
    assert np.isfinite(residual), "Residual should be finite"

    # Should be less than initial naive residual (sum of squared outputs)
    naive_residual = np.sum(y**2)
    assert residual < naive_residual, \
        f"Final residual {residual} should be less than naive {naive_residual}"


def test_ib03bd_small_system():
    """
    Test with a smaller system for faster execution.

    Uses fewer samples and neurons to verify basic functionality.
    Random seed: 456 (for reproducibility)
    """
    from slicot import ib03bd

    np.random.seed(456)

    # Generate small synthetic data
    nsmp = 100
    m = 1
    l = 1
    nobr = 5
    n = 2
    nn = 3

    # Create simple input-output data
    u = np.random.randn(nsmp, m).astype(float, order='F')
    y = np.random.randn(nsmp, l).astype(float, order='F')

    itmax1 = 50
    itmax2 = 50
    tol1 = 1e-3
    tol2 = 1e-3
    init = 'B'

    # Random seed for initialization
    dwork_seed = np.array([42.0, 43.0, 44.0, 45.0], dtype=float, order='F')

    x, iwarn, info, dwork_out = ib03bd(
        init, nobr, m, l, nsmp, n, nn, itmax1, itmax2,
        u, y, tol1, tol2, dwork_seed
    )

    # Should complete without error
    assert info == 0, f"IB03BD failed with INFO={info}"

    # Check output dimensions
    bsn = nn * (l + 2) + 1
    lths = n * (l + m + 1) + l * m
    nx = bsn * l + lths
    assert len(x) == nx
