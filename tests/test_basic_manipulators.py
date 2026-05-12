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
        with pytest.warns(RuntimeWarning):
            identified_mat = mat.identify(tol=1e-10)
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

    def test_sep_unk_splits_each_symbol(self):
        x, y = sym.symbols("x y")
        mat = Matrix([[x + y, 2 * x]])
        parts = mat.sep_unk()
        assert parts[x] == Matrix([[1, 2]])
        assert parts[y] == Matrix([[1, 0]])

    def test_scalar_factor(self):
        mat = Matrix([[2, 4], [6, 8]])
        scalar_factor = mat.scalar_factor()
        assert isinstance(scalar_factor, ScalarFactor)
        assert scalar_factor.full == Matrix([[1, 1], [3, 2]])
        assert scalar_factor.diag == Matrix([[2, 0], [0, 4]])

    def test_scalar_factor_symbolic_reconstructs(self):
        x = sym.symbols("x")
        mat = Matrix([[2 * x, 4], [6 * x, 8]])
        scalar_factor = mat.scalar_factor()
        assert scalar_factor.full @ scalar_factor.diag == mat

    def test_scalar_factor_by_rows(self):
        mat = Matrix([[2, 4], [6, 8]])
        scalar_factor = mat.scalar_factor(column=False)
        assert scalar_factor.order == "DF"
        assert scalar_factor.diag @ scalar_factor.full == mat

    def test_symbolic_reduce_row_updates_l_matrix(self):
        x = sym.symbols("x")
        mat = Matrix([[1, 0], [x, 1]])
        l_matrix = Matrix.eye(2)

        mat._symbolic_reduce_row(1, 0, 0, L_matrix=l_matrix)

        assert mat == Matrix([[1, 0], [0, 1]])
        assert l_matrix == Matrix([[1, 0], [x, 1]])

    def test_symbolic_reduce_row_warns_when_max_tries_exceeded(self):
        mat = Matrix([[0], [1]])
        with pytest.warns(RuntimeWarning, match="Max tries exceeded"):
            mat._symbolic_reduce_row(1, 0, 0, max_tries=0)

    def test_hermitian_transpose(self):
        """Verify the .H property for complex matrices."""
        A = Matrix([[1, 2 + sym.I], [3 - sym.I, 4]])
        expected = Matrix([[1, 3 + sym.I], [2 - sym.I, 4]])
        assert expected == A.H

    def test_is_mat_orthogonal_behavioral(self):
        """Verify is_mat_orthogonal calls .is_diagonal() as a method, not a property."""
        A = Matrix.eye(2)
        assert A.is_mat_orthogonal() is True
        B = Matrix([[1, 1], [0, 1]])
        assert B.is_mat_orthogonal() is False

    def test_is_mat_orthogonal_negative(self):
        """Verify is_mat_orthogonal returns False for non-orthogonal matrices."""
        # Non-square matrix (not orthogonal by definition)
        A = Matrix([[1, 0, 0], [0, 1, 0]])
        assert A.is_mat_orthogonal(verbosity=0) is False
        # Square but non-orthogonal
        B = Matrix([[1, 1], [0, 1]])
        assert B.is_mat_orthogonal(verbosity=0) is False

    def test_scale_row_warning(self):
        """Verify that scaling a row by zero triggers a UserWarning."""
        mat = Matrix([[1, 2], [3, 4]])
        with pytest.warns(UserWarning):
            mat.scale_row(0, 0)

    def test_normalized_zero_column(self):
        """Verify that normalized() handles zero columns without errors."""
        mat = Matrix([[0, 1], [0, 0]])
        result = mat.normalized()
        assert result[0, 0] == 0
        assert result[1, 0] == 0
