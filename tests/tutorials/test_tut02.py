import pytest
import sympy as sym

from ma1522 import Matrix


class TestTutorial02:
    def test_question_2a(self):
        A = Matrix.from_str("1 1 0 1; 0 1 1 0; 0 0 1 1", col_sep=" ", row_sep=";")
        assert A.shape == (3, 4)

        X = A.inverse("right", matrices=2)
        # Particular solution must satisfy A @ X_part == I_3
        assert (A @ X.part_sol).applyfunc(sym.simplify) == Matrix.eye(3)
        assert X.part_sol.shape == (4, 3)

    def test_question_2b(self):
        B = Matrix.from_str("1 0 1; 1 1 0; 0 1 1; 0 0 1", col_sep=" ", row_sep=";")
        assert B.shape == (4, 3)

        Y = B.inverse("left", matrices=2)
        # Particular solution must satisfy Y_part @ B == I_3
        assert (Y.part_sol @ B).applyfunc(sym.simplify) == Matrix.eye(3)
        assert Y.part_sol.shape == (3, 4)

    def test_question_3a(self):
        A = Matrix.from_str("5 -2 6 0; -2 1 3 1")
        result = (
            A.copy()
            .reduce_row(1, -2 / 5, 0)
            .scale_row(0, 1 / 5)
            .scale_row(1, 5)
            .reduce_row(0, -2 / 5, 1)
        )
        # Result must equal the RREF of A
        expected_rref = Matrix([[1, 0, 12, 2], [0, 1, 27, 5]])
        assert result.applyfunc(sym.simplify) == expected_rref

    def test_question_3b(self):
        A = Matrix.from_str("-1 3 -4; 2 4 1; -4 2 -9")
        result = (
            A.copy()
            .reduce_row(1, -2, 0)
            .reduce_row(2, 4, 0)
            .reduce_row(2, -1, 1)
            .scale_row(0, -1)
            .scale_row(1, sym.Rational(1, 10))
            .reduce_row(0, -3, 1)
        )
        expected_rref = A.rref(pivots=False)
        assert result.applyfunc(sym.simplify) == expected_rref

    def test_question_3c(self):
        A = Matrix.from_str("1 -1 0; 2 -2 1; 1 2 3")
        result = (
            A.copy()
            .reduce_row(1, 2, 0)
            .reduce_row(2, 1, 0)
            .swap_row(1, 2)
            .scale_row(1, sym.Rational(1, 3))
            .reduce_row(1, 1, 2)
            .reduce_row(0, -1, 1)
        )
        expected_rref = A.rref(pivots=False)
        assert result.applyfunc(sym.simplify) == expected_rref

    def test_question_4a(self):
        mat = Matrix.from_str("-1 3; 3 -2")
        inv = mat.inverse()
        assert mat.det() != 0
        assert (inv @ mat).applyfunc(sym.simplify) == Matrix.eye(2)

    def test_question_4b(self):
        mat = Matrix.from_str("-1 3 -4; 2 4 1; -4 2 -9")
        # Singular matrix → det == 0 → inverse raises ValueError
        assert mat.det() == 0
        with pytest.raises(ValueError):
            mat.inverse()

    def test_question_5(self):
        a, b, c = sym.symbols("a b c")
        mat = Matrix([[1, 1, 1], [a, b, c], [a**2, b**2, c**2]])

        det = mat.det()
        expected = (b - a) * (c - a) * (c - b)
        assert sym.simplify(det - expected) == 0
