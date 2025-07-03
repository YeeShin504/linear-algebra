import math

import pytest
import sympy as sym
import mpmath as mp


from symbolic import Matrix, PartGen, ScalarFactor


class TestMatrixManipulators:
    def test_copy(self):
        """Test the copy method preserves matrix data and aug_pos."""

        numeric_matrix = Matrix([[1, 2], [3, 4]], aug_pos=1)

        x, y = sym.symbols("x y")
        symbolic_matrix = Matrix([[x, y], [x + y, 2 * y]], aug_pos=set())

        special_matrix = Matrix([[2, 4, 6], [8, 10, 12]], aug_pos=set([0, 1]))

        matrices = [numeric_matrix, symbolic_matrix, special_matrix]

        for matrix in matrices:
            copied = matrix.copy()

            assert copied is not matrix
            assert copied.tolist() == matrix.tolist()
            assert copied._aug_pos == matrix._aug_pos

    def test_simplify_rational(self):
        """Test rational simplification in the simplify method."""
        matrix = Matrix([[1.0, 1 / 3], [2 / 3, 0.25]])
        matrix.simplify(rational=True)

        expected = Matrix(
            [[1, sym.Rational(1, 3)], [sym.Rational(2, 3), sym.Rational(1, 4)]]
        )
        assert matrix == expected

    def test_simplify_expand(self):
        """Test expansion in the simplify method."""
        x = sym.Symbol("x")
        matrix = Matrix(
            [
                [(x + 1) * (x + 2)],
            ]
        )
        matrix.simplify(expand=True)

        expected = Matrix(
            [
                [x**2 + 3 * x + 2],
            ]
        )
        assert matrix == expected

    def test_simplify_factor(self):
        """Test factoring in the simplify method."""
        x = sym.Symbol("x")
        matrix = Matrix(
            [
                [x**2 + 3 * x + 2],
            ]
        )
        matrix.simplify(expand=False)

        expected = Matrix(
            [
                [(x + 1) * (x + 2)],
            ]
        )
        print(matrix)
        print(expected)
        # raise InterruptedError
        assert matrix == expected

    def test_simplify_collect(self):
        """Test collecting terms in the simplify method."""
        x, y, z = sym.symbols("x y z")
        matrix = Matrix(
            [
                [x * y + x - 3 + 2 * x**2 - z * x**2 + x**3],
            ]
        )
        matrix.simplify(collect_sym=x)

        expected = Matrix(
            [
                [x**3 + x**2 * (2 - z) + x * (y + 1) - 3],
            ]
        )
        assert matrix == expected

    def test_identify(self):
        """Test the identify method."""
        matrix = Matrix([[math.sqrt(2), math.pi], [math.e, 1 / math.sqrt(2)]])
        result = matrix.identify(constants=["pi"])

        expected = Matrix([[sym.sqrt(2), sym.pi], [sym.E, sym.sqrt(2) / 2]])
        print(f"{result=}")
        print(f"{expected=}")
        assert (result - expected).norm().evalf() <= 1e-10

    def test_select_cols(self):
        """Test selecting columns from a matrix."""
        matrix = Matrix([[2, 4, 6], [8, 10, 12]])

        result = matrix.select_cols(0, 2)
        expected = Matrix([[2, 6], [8, 12]])
        assert result == expected

        result = matrix.select_cols(-2)
        expected = Matrix([[4], [10]])
        assert result == expected

        with pytest.raises(IndexError):
            matrix.select_cols(3)

    def test_select_rows(self):
        """Test selecting rows from a matrix."""
        matrix = Matrix([[2, 4, 6], [8, 10, 12]])

        result = matrix.select_rows(0)
        expected = Matrix([[2, 4, 6]])
        assert result == expected

        result = matrix.select_rows(0, 1)
        assert result == matrix

        with pytest.raises(IndexError):
            matrix.select_rows(2)

    def test_sep_part_gen(self):
        """Test separating a matrix into particular and general solutions."""
        x, y = sym.symbols("x y")
        matrix = Matrix([[1 + x, 2 - x], [3 + y, 4 * y]])
        result = matrix.sep_part_gen()

        assert isinstance(result, PartGen)
        assert result.part_sol == Matrix([[1, 2], [3, 0]])
        assert result.gen_sol == Matrix([[x, -x], [y, 4 * y]])

    def test_scalar_factor_column(self):
        """Test factorizing a matrix by columns."""
        matrix = Matrix([[1, 2], [3, 4]]).normalized(factor=False)
        result = matrix.scalar_factor(column=True)

        assert isinstance(result, ScalarFactor)
        assert result.order == "FD"
        assert result.full @ result.diag == matrix
        expected_F = Matrix([[1, 1], [3, 2]])
        expected_D = Matrix.diag(*[1 / sym.sqrt(10), 1 / sym.sqrt(5)])
        assert result.full == expected_F
        assert result.diag == expected_D

    def test_scalar_factor_row(self):
        """Test factorizing a matrix by rows."""
        x = sym.symbols("x")
        matrix = Matrix([[0, 0], [-x, 2 * x]])
        result = matrix.scalar_factor(column=False)

        assert isinstance(result, ScalarFactor)
        assert result.order == "DF"
        assert result.diag @ result.full == matrix

        # Check specific values (allowing for numeric approximation)
        expected_F = Matrix([[2, 3], [2, 3]])
        expected_D = Matrix.diag([3, 5])
