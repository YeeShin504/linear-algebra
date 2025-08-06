"""Include the following methods in the tests
- copy
- simplify
- identify
- select_cols
- select_rows
- sep_part_gen
- scalar_factor
"""

import pytest
import sympy as sym

from ma1522 import Matrix, PartGen, ScalarFactor


class TestBasicManipulators:
    def test_copy(self):
        mat = Matrix([[1, 2], [3, 4]], aug_pos={1})
        copied_mat = mat.copy()
        assert copied_mat == mat
        assert copied_mat is not mat

    def test_simplify(self):
        x = sym.symbols("x")
        mat = Matrix([[x**2 - 1, x + 1], [sym.sin(x) ** 2 + sym.cos(x) ** 2, 2]])  # type: ignore
        simplified_mat = mat.copy()
        simplified_mat.simplify()
        assert simplified_mat == Matrix([[x**2 - 1, x + 1], [1, 2]])

    def test_identify(self):
        mat = Matrix([[1.0000000000001, 2.0], [3.0, 4.0]])
        identified_mat = mat.identify(tolerance=1e-10)
        assert (identified_mat - Matrix([[1, 2], [3, 4]])).norm() < 1e-9

    def test_select_cols(self):
        """Test column selection"""
        matrix = Matrix([[2, 4, 6], [8, 10, 12]])

        # Select multiple columns
        result = matrix.select_cols(0, 2)
        expected = Matrix([[2, 6], [8, 12]])
        assert result == expected

        # Select single column
        result = matrix.select_cols(-2)
        expected = Matrix([[4], [10]])
        assert result == expected

        # Test invalid index
        with pytest.raises(IndexError):
            matrix.select_cols(3)

    def test_select_rows(self):
        """Test row selection"""
        matrix = Matrix([[2, 4, 6], [8, 10, 12]])

        # Select single row
        result = matrix.select_rows(0)
        expected = Matrix([[2, 4, 6]])
        assert result == expected

        # Select all rows
        result = matrix.select_rows(0, 1)
        assert result == matrix

        # Test invalid index
        with pytest.raises(IndexError):
            matrix.select_rows(2)

    def test_sep_part_gen(self):
        x, y = sym.symbols("x y")
        mat = Matrix([[x + 1, y], [2, x + y]])
        part_gen = mat.sep_part_gen()
        assert isinstance(part_gen, PartGen)
        assert part_gen.part_sol == Matrix([[1, 0], [2, 0]])
        assert part_gen.gen_sol == Matrix([[x, y], [0, x + y]])

    def test_scalar_factor(self):
        mat = Matrix([[2, 4], [6, 8]])
        scalar_factor = mat.scalar_factor()
        assert isinstance(scalar_factor, ScalarFactor)
        assert scalar_factor.full == Matrix([[1, 1], [3, 2]])
        assert scalar_factor.diag == Matrix([[2, 0], [0, 4]])
