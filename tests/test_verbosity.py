import pytest
import sympy as sym
from ma1522 import Matrix

class TestVerbosityCoverage:
    """Tests designed to hit printing/verbosity lines in symbolic.py for coverage."""
    
    def test_ref_verbosity(self):
        A = Matrix.from_str("1 2; 3 4")
        # REF verbosity
        A.ref(verbosity=1)
        A.ref(verbosity=2)

    def test_rref_standard(self):
        A = Matrix.from_str("1 2; 3 4")
        # RREF does not take verbosity, but evaluate_cases does
        A.rref()
        
    def test_solve_verbosity(self):
        A = Matrix.from_str("1 2; 3 4")
        b = Matrix.from_str("5; 6")
        A.solve(b, verbosity=1)
        A.solve(b, verbosity=2)
        
    def test_inverse_verbosity(self):
        A = Matrix.from_str("1 2; 3 4")
        A.inverse(verbosity=1)
        A.inverse(verbosity=2)
        
    def test_ref_plu_verbosity(self):
        A = Matrix.from_str("1 2; 3 4")
        A.ref(verbosity=1)
        
    def test_det_standard(self):
        A = Matrix.from_str("1 2; 3 4")
        A.det()
        
    def test_diagonalize_verbosity(self):
        A = Matrix.from_str("1 0; 0 1")
        A.diagonalize(verbosity=1)
        A.is_diagonalizable(verbosity=1)
        
    def test_svd_verbosity(self):
        A = Matrix.from_str("1 2; 3 4")
        A.singular_value_decomposition(verbosity=1)
        
    def test_gram_schmidt_verbosity(self):
        A = Matrix.from_str("1 0; 0 1").T
        A.gram_schmidt(verbosity=1)
        
    def test_orthogonality_verbosity(self):
        A = Matrix.from_str("1 0; 0 1").T
        A.is_vec_orthogonal(verbosity=1)
        A.is_mat_orthogonal(verbosity=1)
        
    def test_subspace_verbosity(self):
        A = Matrix.from_str("1 0; 0 1").T
        B = Matrix.from_str("1 0; 0 1").T
        A.is_subspace_of(B, verbosity=1)
        A.is_subspace_of(B, verbosity=2)
        A.is_same_subspace(B, verbosity=1)
        
    def test_evaluate_cases_verbosity(self):
        a = sym.Symbol("a")
        A = Matrix([[a, 1], [1, 1]])
        rhs = Matrix([1, 1])
        A.evaluate_cases(rhs=rhs, verbosity=1)
        
    def test_stochastic_verbosity(self):
        A = Matrix([[sym.S.Half, sym.S.Half], [sym.S.Half, sym.S.Half]])
        A.is_stochastic(verbosity=1)
