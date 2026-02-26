import pytest
import sympy as sym

from ma1522 import Matrix
from ma1522.custom_types import ScalarFactor


class TestTutorial07:
    def test_question_1b(self):
        mat = Matrix.from_str("1 3 -2 0; 2 6 -5 -2; 0 0 5 10")
        aug = Matrix.zeros(rows=3, cols=1)

        sol = mat.solve(aug)[0]
        assert len(sol.free_symbols) == 2

        x = Matrix([-3, 1, 0, 0])
        z = Matrix([-4, 0, -2, 1])

        sol_dict = sol.sep_unk()
        assert x in sol_dict.values()
        assert z in sol_dict.values()

    def test_question_1c(self):
        A = Matrix.from_str("1 3 -2 0; 2 6 -5 -2; 0 0 5 10").T
        v = Matrix.create_unk_matrix(4, 1, "x")
        sol = sym.solve(A.T @ v, v)
        v_sub = v.subs(sol)
        assert len(v_sub.free_symbols) == 2

        W_perp = A.orthogonal_complement()
        assert W_perp.cols == A.rows - A.rank()

        # Every column of W_perp must be orthogonal to every column of A
        for j in range(W_perp.cols):
            v = W_perp.select_cols(j)
            assert (A.T @ v).applyfunc(sym.simplify) == Matrix.zeros(A.cols, 1)

    def test_question_3a(self):
        v1 = Matrix.from_str("1; 2; -1")
        v2 = Matrix.from_str("1; 0; 1")
        assert v1.dot(v1) == 6
        assert v1.dot(v2) == 0
        assert v2.dot(v1) == 0
        assert v2.dot(v2) == 2

    def test_question_3b(self):
        v1 = Matrix.from_str("1; 2; -1")
        v2 = Matrix.from_str("1; 0; 1")
        V = v1.row_join(v2)
        gram = V.T @ V
        assert gram == Matrix([[6, 0], [0, 2]])

    def test_question_4a(self):
        S = Matrix.from_str("1 1 1 1 1; 1 2 -1 -2 0; 1 -1 1 -1 0").T
        assert S.is_linearly_independent(verbosity=2)

    def test_question_4b(self):
        S = Matrix.from_str("1 1 1 1 1; 1 2 -1 -2 0; 1 -1 1 -1 0").T
        # All pairs are orthogonal
        assert S.is_vec_orthogonal()

    def test_question_4c(self):
        S = Matrix.from_str("1 1 1 1 1; 1 2 -1 -2 0; 1 -1 1 -1 0").T
        W_perp = S.orthogonal_complement()
        assert W_perp.shape == (5, 2)
        assert W_perp.rank() == 2
        for j in range(W_perp.cols):
            v = W_perp.select_cols(j)
            assert (S.T @ v).applyfunc(sym.simplify) == Matrix.zeros(3, 1)

    def test_question_4d(self):
        S = Matrix.from_str("1 1 1 1 1; 1 2 -1 -2 0; 1 -1 1 -1 0").T
        T_norm = S.normalized(factor=True)
        assert T_norm.full.shape == S.shape
        assert T_norm.full.is_vec_orthogonal()
        assert T_norm.diag.is_diagonal()
        assert T_norm.order == "FD"

    def test_question_4ef(self):
        S = Matrix.from_str("1 1 1 1 1; 1 2 -1 -2 0; 1 -1 1 -1 0").T
        v = Matrix.from_str("2; 0; 1; 1; -1")
        W = S
        v_proj = v.proj_comp(W, verbosity=1)

        # v_proj must be in W (check consistency)
        aug = W.row_join(v_proj)
        assert aug.rank() == W.rank()

        # v - v_proj must be orthogonal to W
        W_perp = S.orthogonal_complement()
        v_norm = v - v_proj
        aug2 = W_perp.row_join(v_norm)
        assert aug2.rank() == W_perp.rank()

    def test_question_5a(self):
        S = Matrix.from_str("1 2 2 -1; 1 1 -1 1; -1 1 -1 -1; -2 1 1 2").T
        assert S.is_vec_orthogonal()
        assert S.rank() == 4

    def test_question_5b(self):
        S = Matrix.from_str("1 2 2 -1; 1 1 -1 1; -1 1 -1 -1; -2 1 1 2").T
        with pytest.warns(UserWarning):  # trivial orthogonal complement
            w = S.orthogonal_complement(verbosity=1)
            assert len(w) == 0

    def test_question_5d(self):
        S = Matrix.from_str("1 2 2 -1; 1 1 -1 1; -1 1 -1 -1; -2 1 1 2").T
        v = Matrix.from_str("0; 1; 2; 3")
        v_S = v.coords_relative(S)
        assert (S @ v_S).applyfunc(sym.simplify) == v
        expected = Matrix(
            [
                sym.Rational(3, 10),
                sym.Rational(1, 2),
                -1,
                sym.Rational(9, 10),
            ]
        )
        assert v_S.applyfunc(sym.simplify) == expected

    def test_question_5e(self):
        S = Matrix.from_str("1 2 2 -1; 1 1 -1 1; -1 1 -1 -1; -2 1 1 2").T
        T = S.copy().normalized(factor=True)
        assert isinstance(T, ScalarFactor)
        assert T.full.is_vec_orthogonal()
        T_eval = T.eval()
        assert isinstance(T_eval, Matrix)
        assert T_eval.shape == S.shape
        assert T_eval.is_mat_orthogonal()

        w_S = Matrix.from_str("1; 2; 1; 1")
        w = S @ w_S
        w_T = w.coords_relative(T_eval)
        assert (T_eval @ w_T).applyfunc(sym.simplify) == w
