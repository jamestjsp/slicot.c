import unittest
import numpy as np
try:
    from slicot import md03ba
except ImportError:
    pass

class TestMD03BA(unittest.TestCase):
    def test_md03ba_basic(self):
        try:
            from slicot import md03ba
        except ImportError:
            self.fail("Could not import md03ba")

        np.random.seed(42)
        m, n = 5, 3
        ipar = np.array([m], dtype=np.int32)
        
        j_in = np.random.rand(m, n)
        j = np.asfortranarray(j_in)
        
        e_in = np.random.rand(m)
        e = e_in.copy()
        fnorm = np.linalg.norm(e)
        
        # Wrapper signature:
        # j_out, e_out, jnorms, gnorm, ipvt, info = md03ba(n, ipar, fnorm, j, e)
        
        j_out, e_out, jnorms, gnorm, ipvt, info = md03ba(n, ipar, fnorm, j, e)
        
        assert info == 0
        
        # Verify E norm preservation (Q is orthogonal)
        norm_e_out = np.linalg.norm(e_out)
        np.testing.assert_allclose(norm_e_out, fnorm, rtol=1e-14)
        
        # Verify dimensions
        assert j_out.shape == (m, n)
        assert e_out.shape == (m,)
        assert jnorms.shape == (n,)
        assert ipvt.shape == (n,)

if __name__ == '__main__':
    unittest.main()
