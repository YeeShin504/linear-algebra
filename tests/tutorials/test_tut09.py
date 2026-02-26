from numpy import char
import pytest
import sympy as sym
from sympy.matrices.matrixbase import MatrixError

from ma1522 import Matrix


class TestTutorial09:
    def test_question_1a(self):
        mat = Matrix.from_str("1 2 0; 0 1 1; 1 0 -2")
        aug = Matrix.from_str("300; 300; 300")
        assert not aug.is_subspace_of(mat)

    def test_question_1b(self):
        mat = Matrix.from_str("1 2 0; 0 1 1; 1 0 -2")
        aug = Matrix.from_str("300; 300; 300")
        sol = mat.solve_least_squares(aug)
        part, gen = sol.sep_part_gen()
        assert part == Matrix([200, 100, 0])
        assert len(gen.free_symbols) == 1
        free = gen.free_symbols.pop()
        assert gen == free * Matrix([2, -1, 1])

    def test_question_2b(self):
        A = Matrix.from_str("1 1 0; 1 1 0; 1 1 1; 0 1 1")
        qr = A.QRdecomposition(full=True)
        Q, R = qr.Q, qr.R
        assert Q.shape == (4, 4)
        gram = (Q.T @ Q).applyfunc(lambda x: x.simplify())
        assert gram == Matrix.eye(4)
        Q_thin = Q.select_cols(0, 1, 2)
        R_top = R.select_rows(0, 1, 2)
        assert (Q_thin @ R_top).applyfunc(sym.simplify) == A

    def test_question_3a(self):
        X = Matrix.from_str("1 1 2; 1 2 1; 2 1 1")
        p_X = X**3 - 4 * X**2 - X + 4 * Matrix.eye(3)
        assert p_X.applyfunc(sym.simplify) == Matrix.zeros(3, 3)

    def test_question_3b(self):
        X = Matrix.from_str("1 1 2; 1 2 1; 2 1 1")
        # X is symmetric — all eigenvalues must be real
        char_poly = X.cpoly()
        assert isinstance(char_poly, sym.Mul)

        x = sym.Symbol("x", real=True)
        expected = (x - 4) * (x + 1) * (x - 1)
        assert sym.simplify(char_poly - expected) == 0

    def test_question_4a(self):
        A = Matrix.from_str("1 -3 3; 3 -5 3; 6 -6 4")
        assert A.is_diagonalizable()
        PDP = A.diagonalize(verbosity=1)
        assert (PDP.P.inv() @ A @ PDP.P).applyfunc(sym.simplify) == PDP.D

    def test_question_4b(self):
        A = Matrix.from_str("9 8 6 3; 0 -1 3 -4; 0 0 2 0; 0 0 0 3")
        PDP = A.diagonalize(reals_only=True)
        assert (PDP.P.inv() @ A @ PDP.P).applyfunc(sym.simplify) == PDP.D
        # Diagonal entries are the eigenvalues
        eigenvals = set(PDP.D.diagonal())
        assert eigenvals == {9, -1, 2, 3}

    def test_question_4c(self):
        A = Matrix.from_str("1 0 0; 1 1 0; 0 1 1")
        # Not diagonalizable over reals → raises MatrixError
        with pytest.raises(MatrixError):
            A.diagonalize(reals_only=True)

    def test_question_4d(self):
        A = Matrix.from_str("0 0 1 0; 0 0 0 1; 1 0 0 0; 0 1 0 0")
        PDP = A.diagonalize(reals_only=True)
        assert (PDP.P.inv() @ A @ PDP.P).applyfunc(sym.simplify) == PDP.D

    def test_question_4e(self):
        A = Matrix.from_str("-1 1 1; 1 1 -1; -4 2 3")
        # Has complex eigenvalues → not diagonalizable over reals
        assert not A.is_diagonalizable(reals_only=True)
        # but diagonlizable over complex numbers
        assert A.is_diagonalizable(reals_only=False)
        PDP = A.diagonalize(reals_only=False)
        assert (PDP.P.inv() @ A @ PDP.P).applyfunc(sym.simplify) == PDP.D
