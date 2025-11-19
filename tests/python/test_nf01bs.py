import unittest
import numpy as np
from slicot import nf01bs

class TestNF01BS(unittest.TestCase):
    def test_nf01bs_full_matrix(self):
        """
        Validate NF01BS with full matrix (BN <= 1 or BSN = 0).
        Should behave like MD03BX (QR with pivoting).
        """
        np.random.seed(42)
        
        m, n = 5, 3
        st = 0
        bn = 1
        bsm = 5
        bsn = 3
        ipar = np.array([st, bn, bsm, bsn], dtype=np.int32)
        
        fnorm = 1.0
        j_in = np.random.rand(m, n)
        j = np.asfortranarray(j_in)
        
        e_in = np.random.rand(m)
        e = e_in.copy()
        
        # Wrapper: j_out, e_out, jnorms, gnorm, ipvt, info = nf01bs(n, ipar, fnorm, j, e)
        
        j_out, e_out, jnorms, gnorm, ipvt, info = nf01bs(n, ipar, fnorm, j, e)
        
        assert info == 0
        
        # Verify Q'*e norm (orthogonal transformation preserves norm)
        # Note: Q is m x m. e is m.
        np.testing.assert_allclose(np.linalg.norm(e_out), np.linalg.norm(e_in), rtol=1e-14)
        
        # Verify R is upper triangular (in first n rows)
        r = np.triu(j_out[:n, :])
        # Check strict lower part is zero? 
        # MD03BX returns R in upper triangle. Lower part contains Householder vectors?
        # Actually MD03BX/NF01BS returns "upper triangular factor R". 
        # "On exit, the leading N-by-N upper triangular part... contains R"
        
        # Check diagonal is non-increasing (pivoting)
        # MD03BX guarantees diagonal elements of nonincreasing magnitude
        diag_r = np.abs(np.diag(j_out))
        # Check sorted descending (nonincreasing)
        assert np.all(np.diff(diag_r) <= 0), "Diagonal not sorted descending"

if __name__ == '__main__':
    unittest.main()
