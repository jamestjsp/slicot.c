import unittest
import numpy as np
from slicot import mb04ow
from scipy.linalg import qr

class TestMB04OW(unittest.TestCase):
    def test_mb04ow_basic(self):
        """
        Validate MB04OW against scipy.linalg.qr
        
        MB04OW updates the QR factorization of:
             ( U  )   = Q * ( R )
             ( x' )         ( 0 )
        
        where U = [ U1 U2 ]
                  [ 0  T  ]
        
        It also updates [B; C; d'] by applying Q'.
        """
        # Reproducibility
        np.random.seed(42)
        
        m, n, p = 4, 3, 2
        lda = m + 2
        ldt = n + 2
        ldb = m + 2
        ldc = n + 2
        
        # Generate random matrices (Fortran order)
        # U1: mxm upper triangular
        # U2: mxn full
        # T: nxn upper triangular
        u1 = np.triu(np.random.rand(m, m))
        u2 = np.random.rand(m, n)
        t = np.triu(np.random.rand(n, n))
        
        # Store U1, U2 in A
        a = np.zeros((lda, n + m), order='F')
        a[:m, :m] = u1
        a[:m, m:m+n] = u2
        
        # Store T in T_arr
        t_arr = np.zeros((ldt, n), order='F')
        t_arr[:n, :n] = t
        
        # B, C
        b = np.random.rand(m, p)
        b_arr = np.zeros((ldb, p), order='F')
        b_arr[:m, :p] = b
        
        c = np.random.rand(n, p)
        c_arr = np.zeros((ldc, p), order='F')
        c_arr[:n, :p] = c
        
        # x (m+n vector)
        x = np.random.rand(m + n)
        x_arr = x.copy() # MB04OW modifies X
        
        # d (p vector)
        d = np.random.rand(p)
        d_arr = d.copy() # MB04OW modifies D
        
        # Construct the big matrix K for verification
        # K = [ U1  U2  B ]
        #     [ 0   T   C ]
        #     [ x'      d']
        #
        # But carefully map x: x is length m+n
        # x[0:m] corresponds to U1 columns
        # x[m:m+n] corresponds to U2/T columns
        
        # Dimensions of K: (m + n + 1) x (m + n + p)
        K = np.zeros((m + n + 1, m + n + p))
        
        # Row 0..m-1
        K[:m, :m] = u1
        K[:m, m:m+n] = u2
        K[:m, m+n:] = b
        
        # Row m..m+n-1
        K[m:m+n, m:m+n] = t
        K[m:m+n, m+n:] = c
        
        # Row m+n (the update row)
        K[m+n, :m+n] = x
        K[m+n, m+n:] = d
        
        # Perform QR on K using scipy
        # We want R such that K = Q * R
        # Since K has more columns than rows, scipy qr with mode='economic' returns (M, M) R.
        # But here rows=m+n+1, cols=m+n+p. Typically p is small, so cols > rows is possible.
        # MB04OW assumes we zero out the LAST row (index m+n).
        # And it returns the updated triangular factor R.
        
        # Actually, MB04OW computes Q such that Q^T * [ U; x'] = [ R; 0 ]
        # So Q^T * K = [ Updated_Upper_Part ]
        #              [ 0 ... 0  Resid ]
        
        # Let's compute QR of K
        Q_full, R_full = qr(K, mode='full')
        
        # Note: QR decomposition is not unique (signs of diagonals).
        # However, R^T R = K^T K.
        # Also, if we enforce positive diagonal for R, it is unique.
        # MB04OW uses DLARTG which might not enforce positive diagonal.
        # So we should compare absolute values of diagonal, or better:
        # Reconstruct K from output of MB04OW and check error.
        
        # Call MB04OW
        # mb04ow(m, n, p, a, t, x, b, c, d, incx=1, incd=1)
        # Note: python wrapper should handle array dimensions.
        # We assume the wrapper signature:
        # a_out, t_out, x_out, b_out, c_out, d_out = mb04ow(a, t, x, b, c, d)
        # OR inplace. Let's assume standard SLICOT wrapper pattern (return modified arrays).
        
        a_out, t_out, x_out, b_out, c_out, d_out = mb04ow(m, n, p, a, t_arr, x_arr, b_arr, c_arr, d_arr)
        
        # Reconstruct the updated matrix from outputs
        # The top m+n rows should contain the updated R and transformed B,C
        # The last row is implicitly zero in the first m+n cols, and contains transformed d in last p cols.
        
        # But wait, MB04OW doesn't return the full last row of the transformed matrix?
        # "The transformations performed are also applied to ... ( B' C' d )' "
        # Yes, d is updated.
        # x is updated (to zeroes? or contains information about rotations? 
        # "On exit, the content of X is changed." (Likely garbage or rotation info)
        # The documentation says "Q(j) being chosen to annihilate the jth element of x."
        # So mathematically, after rotations, x should be 0.
        # In practice, X is workspace.
        
        # Construct K_out from outputs
        K_out = np.zeros_like(K)
        
        # Extract R1, R2 from a_out
        r1 = a_out[:m, :m] # Upper triangular
        r2 = a_out[:m, m:m+n]
        
        # Extract R3 from t_out
        r3 = t_out[:n, :n] # Upper triangular
        
        # Extract B_out, C_out
        b_new = b_out[:m, :p]
        c_new = c_out[:n, :p]
        
        # Extract d_new
        d_new = d_out
        
        # Fill K_out
        K_out[:m, :m] = np.triu(r1)
        K_out[:m, m:m+n] = r2
        K_out[:m, m+n:] = b_new
        
        K_out[m:m+n, m:m+n] = np.triu(r3)
        K_out[m:m+n, m+n:] = c_new
        
        K_out[m+n, m+n:] = d_new
        
        # Now, K_out should be Q^T * K for some orthogonal Q.
        # Thus K_out^T * K_out should equal K^T * K.
        
        gram_in = K.T @ K
        gram_out = K_out.T @ K_out
        
        np.testing.assert_allclose(gram_out, gram_in, rtol=1e-12, atol=1e-14)
        
        # Also check that R1, R3 are upper triangular
        np.testing.assert_allclose(np.tril(r1, -1), 0, atol=1e-15)
        np.testing.assert_allclose(np.tril(r3, -1), 0, atol=1e-15)

if __name__ == '__main__':
    unittest.main()
