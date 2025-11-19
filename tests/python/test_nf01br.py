import unittest
import numpy as np
from slicot import nf01br

class TestNF01BR(unittest.TestCase):
    def test_nf01br_full_rank(self):
        """
        Validate NF01BR with full rank matrix (BN <= 1 or BSN = 0).
        Solves R*x = b.
        """
        np.random.seed(42)
        
        # Case 1: Full upper triangular matrix (BN <= 1)
        # n = 4, ipar structure: st, bn, bsm, bsn
        # Let's use BN=1, BSM=4, BSN=4, ST=0 => N=4.
        
        n = 4
        st = 0
        bn = 1
        bsm = 4
        bsn = 4
        ipar = np.array([st, bn, bsm, bsn], dtype=np.int32)
        
        ldr = n
        r = np.triu(np.random.rand(n, n))
        # Ensure non-singular
        for i in range(n):
            r[i, i] += 2.0
            
        r_arr = np.asfortranarray(r)
        
        b = np.random.rand(n)
        b_in = b.copy()
        
        # sdiag and s are not used for UPLO='U'
        sdiag = np.zeros(n)
        s = np.zeros((1, 1))
        ranks = np.zeros(bn+1, dtype=np.int32)
        
        # Call nf01br(cond, uplo, trans, n, ipar, r, sdiag, s, b, ranks, tol)
        # cond='N', uplo='U', trans='N'
        
        b_out, ranks_out, info = nf01br('N', 'U', 'N', n, ipar, r_arr, sdiag, s, b_in, ranks, 0.0)
        
        assert info == 0
        
        # Verify R*x = b
        x = b_out
        b_rec = r @ x
        np.testing.assert_allclose(b_rec, b, rtol=1e-14)

    def test_nf01br_transpose(self):
        """
        Validate NF01BR with full rank matrix, solving R'*x = b.
        """
        np.random.seed(42)
        n = 4
        st = 0
        bn = 1
        bsm = 4
        bsn = 4
        ipar = np.array([st, bn, bsm, bsn], dtype=np.int32)
        
        r = np.triu(np.random.rand(n, n))
        for i in range(n): r[i, i] += 2.0
        r_arr = np.asfortranarray(r)
        
        b = np.random.rand(n)
        b_in = b.copy()
        
        sdiag = np.zeros(n)
        s = np.zeros((1, 1))
        ranks = np.zeros(bn+1, dtype=np.int32)
        
        # Solve R'*x = b
        b_out, ranks_out, info = nf01br('N', 'U', 'T', n, ipar, r_arr, sdiag, s, b_in, ranks, 0.0)
        
        assert info == 0
        
        x = b_out
        b_rec = r.T @ x
        np.testing.assert_allclose(b_rec, b, rtol=1e-14)

if __name__ == '__main__':
    unittest.main()
