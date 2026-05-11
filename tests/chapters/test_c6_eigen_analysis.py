"""Include the following methods in the tests
- cpoly
- is_diagonalizable
- eigenvects_associated
- diagonalize
- is_orthogonally_diagonalizable
- orthogonally_diagonalize
- is_stochastic
- equilibrium_vectors
- singular_value_decomposition
- fast_svd
"""

import sympy as sym
import pytest

from ma1522 import Matrix, SVD


class TestChapter6:
    def test_cpoly(self):
        mat = Matrix([[1, 2], [3, 4]])
        poly = mat.cpoly()
        poly = sym.Poly(poly)
        assert poly.all_coeffs() == [1, -5, -2]

    def test_cpoly_2(self):
        res = Matrix.diag(-1, 4, sym.Rational(1, 2)).cpoly()
        roots = set()
        for expr in res.args:  # type: ignore
            roots.add(sym.solve(expr)[0])
        assert roots == set((4, sym.Rational(1, 2), -1))

    def test_is_diagonalizable(self):
        mat = Matrix([[1, 1], [0, 1]])
        assert mat.is_diagonalizable() is False
        mat = Matrix([[1, 0], [0, 2]])
        assert mat.is_diagonalizable() is True

    def test_diagonalize(self):
        mat = Matrix([[1, 2], [0, 3]])
        pdp = mat.diagonalize()
        assert (pdp.P @ pdp.D @ pdp.P.inv() - mat).norm() < 1e-10

    def test_diagonalize_verbose_filters_complex_roots(self):
        mat = Matrix([[0, -1], [1, 0]])
        with pytest.raises(Exception):
            mat.diagonalize(reals_only=True, verbosity=1)

    def test_eigenvects_associated(self):
        mat = Matrix([[2, 0], [0, 3]])
        assert mat.eigenvects_associated(2) == [Matrix([[1], [0]])]

    def test_is_orthogonally_diagonalizable(self):
        mat = Matrix([[1, 2], [2, 1]])
        assert mat.is_orthogonally_diagonalizable() is True
        mat = Matrix([[1, 2], [3, 4]])
        assert mat.is_orthogonally_diagonalizable() is False

    def test_orthogonally_diagonalize(self):
        mat = Matrix([[1, 2], [2, 1]])
        pdp = mat.orthogonally_diagonalize(verbosity=0)
        assert (pdp.P @ pdp.D @ pdp.P.T - mat).norm() < 1e-10

    def test_orthogonally_diagonalize_without_factoring(self):
        mat = Matrix.eye(2)
        pdp = mat.orthogonally_diagonalize(factor=False, verbosity=0)
        assert (pdp.P @ pdp.D @ pdp.P.T - mat).norm() < 1e-10

    def test_is_stochastic(self):
        assert Matrix([[sym.Rational(1, 2), sym.Rational(1, 3)], [sym.Rational(1, 2), sym.Rational(2, 3)]]).is_stochastic(verbosity=0)
        assert not Matrix([[1, 1], [-1, 0]]).is_stochastic(verbosity=0)
        assert not Matrix([[1], [0]]).is_stochastic(verbosity=0)

    def test_equilibrium_vectors(self):
        mat = Matrix([[0.8, 0.3], [0.2, 0.7]])
        eq_vects = mat.equilibrium_vectors()
        assert mat @ eq_vects == eq_vects

    def test_fast_svd(self):
        mat = Matrix([[1, 2], [3, 4]])
        svd = mat.fast_svd(option="np", identify=False)
        U, S, V = svd.U, svd.S, svd.V
        assert (Matrix(U) @ Matrix(S) @ Matrix(V).T - mat).norm() < 1e-10

    def test_fast_svd_invalid_option_falls_back_to_numpy_result(self):
        mat = Matrix([[1, 2], [3, 4]])
        with pytest.warns(SyntaxWarning, match="Invalid option"):
            svd = mat.fast_svd(option="bad")  # type: ignore[arg-type]
        assert not isinstance(svd, SVD)

    def test_singular_value_decomposition(self):
        mat = Matrix([[1, 2], [3, 4]])
        svd = mat.singular_value_decomposition(verbosity=2)
        assert (svd.U @ svd.S @ svd.V.T - mat).norm() < 1e-9
