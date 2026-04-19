import sympy as sym
import pytest

from ma1522 import Matrix


class TestTutorial11:
    def test_question_1a(self):
        in_vec = Matrix.from_str("x; y")
        out_vec = Matrix.from_str("x+y; y-x")
        mat = in_vec.standard_matrix(out_vec, matrices=1)
        assert len(mat) == 1
        assert mat[0] == Matrix.from_str("1 1; -1 1")
        assert len(mat[0].columnspace()) == 2
        assert mat[0].nullspace() == []

    def test_question_1b(self):
        in_vec = Matrix.from_str("x; y")
        out_vec = Matrix.from_str("2**x; 0")
        with pytest.raises(ValueError):
            in_vec.standard_matrix(out_vec, matrices=1)

    def test_question_1c(self):
        in_vec = Matrix.from_str("x; y")
        out_vec = Matrix.from_str("x+y; 0; 0")
        mat = in_vec.standard_matrix(out_vec, matrices=1)
        assert len(mat) == 1
        assert mat[0].shape == (3, 2)
        assert len(mat[0].columnspace()) == 1
        assert len(mat[0].nullspace()) == 1

    def test_question_1d(self):
        in_vec = Matrix.from_str("x; y; z")
        out_vec = Matrix.from_str("1; y-x; y-z")
        with pytest.raises(ValueError):
            in_vec.standard_matrix(out_vec, matrices=1)

    def test_question_1e(self):
        in_vec = Matrix.create_unk_matrix(5, 1, "x", is_real=True)
        out_vec = Matrix.from_str("x3 + 2*x4 - x5", col_sep="?", is_real=True)
        result = in_vec.standard_matrix(out_vec, matrices=1)
        assert len(result) == 1
        assert result[0].shape == (1, 5)
        assert result[0] == Matrix.from_str("0 0 1 2 -1")

    def test_question_1f(self):
        x = Matrix.create_unk_matrix(3, 1)
        with pytest.raises(ValueError):
            x.standard_matrix(Matrix([x.dot(x)]), matrices=1)

    def test_question_2a(self):
        F_in = Matrix.from_str("x1; x2; x3")
        F_out = Matrix.from_str("x1 - 2*x2; x1 + x2 - 3*x3; 5*x2 - x3", col_sep="?")
        A_F = F_in.standard_matrix(F_out, matrices=1)[0]
        assert A_F == Matrix.from_str("1 -2 0; 1 1 -3; 0 5 -1")

        G_in = Matrix.from_str("x1; x2; x3")
        G_out = Matrix.from_str("x3 - x1; x2 + 5*x1; x1 + x2 + x3", col_sep="?")
        B_G = G_in.standard_matrix(G_out, matrices=1)[0]
        assert B_G == Matrix.from_str("-1 0 1; 5 1 0; 1 1 1")

    def test_question_2c(self):
        A_F = Matrix.from_str("1 -2 0; 1 1 -3; 0 5 -1")
        B_G = Matrix.from_str("-1 0 1; 5 1 0; 1 1 1")
        comp = A_F @ B_G
        expected = Matrix.from_str("-11 -2 1; 1 -2 -2; 24 4 -1")
        assert comp.shape == (3, 3)
        assert comp == expected

    def test_question_2d(self):
        B_G = Matrix.from_str("-1 0 1; 5 1 0; 1 1 1")
        H = B_G.inverse()
        assert (H @ B_G).applyfunc(sym.simplify) == Matrix.eye(3)

    def test_question_3a(self):
        T_in = Matrix.eye(3)
        T_out = Matrix.from_str("1 3 0 1; 2 2 -1 4; 0 4 1 6").T
        mat = T_in.standard_matrix(T_out, matrices=1)
        assert len(mat) == 1
        assert mat[0].shape == (4, 3)
        assert mat[0] == T_out

    def test_question_3b(self):
        # Overdetermined but consistent → still a unique solution
        T_in = Matrix.from_str("1 -1; 1 1; 2 0").T
        T_out = Matrix.from_str("2 0; 0 2; 2 2").T
        mat = T_in.standard_matrix(T_out, matrices=1)
        assert len(mat) == 1
        assert (mat[0] @ T_in).applyfunc(sym.simplify) == T_out

    def test_question_3c(self):
        T_in = Matrix.from_str("1 -1 0; 0 1 -1; -1 0 1").T
        T_out = Matrix.from_str("-1 1 0")
        mat = T_in.standard_matrix(T_out, matrices=1)
        assert len(mat[0].free_symbols) == 1
        free = mat[0].free_symbols.pop()
        expected = Matrix([[free, free + 1, free]])
        assert mat[0] == expected
