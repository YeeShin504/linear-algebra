"""Include the following methods in the tests
- aug_line
- rm_aug_line
- row_join
- col_join
- scale_row
- swap_row
- reduce_row
- get_pivot_row
- get_pivot_pos
- get_pivot_elements
- ref
- find_all_cases
- evaluate_cases
- rref
- solve
"""

import pytest
from sympy import S

from ma1522 import Matrix


class TestAugLine:
    """Tests for Matrix.aug_line() method

    Covers:
    - Adding augmentation lines at different positions
    - Multiple augmentation lines
    - Edge cases (invalid positions)
    """

    def test_single_aug_line(self):
        mat = Matrix([[1, 2], [3, 4]])
        mat2 = mat.aug_line(1)
        assert hasattr(mat2, "_aug_pos")
        assert mat2._aug_pos == set([1])

    def test_multiple_aug_lines(self):
        """Test adding multiple augmentation lines"""
        mat = Matrix([[1, 2], [3, 4]])
        mat = mat.aug_line(0).aug_line(1)
        assert mat._aug_pos == set([0, 1])

    def test_invalid_position(self):
        """Test invalid augmentation position raises error"""
        mat = Matrix([[1, 2], [3, 4]])
        with pytest.raises(IndexError):
            mat.aug_line(3)


class TestRmAugLine:
    """Tests for Matrix.rm_aug_line() method

    Covers:
    - Removing existing augmentation lines
    - Attempting to remove non-existent lines
    """

    def test_remove_existing_line(self):
        """Test removing an existing augmentation line"""
        mat = Matrix([[1, 2], [3, 4]], aug_pos={0, 1})
        mat = mat.rm_aug_line(0)
        assert mat._aug_pos == set([1])

    def test_remove_nonexistent_line(self):
        """Test removing a non-existent line has no effect"""
        mat = Matrix([[1, 2], [3, 4]], aug_pos={0})
        mat = mat.rm_aug_line(1)
        assert mat._aug_pos == set([0])


class TestJoinOperations:
    """Tests for basic join operations

    Covers:
    - row_join()
    - col_join()
    """

    # TODO


class TestRowOperations:
    """Tests for basic row operations

    Covers:
    - scale_row()
    - swap_row()
    - reduce_row()
    """

    def test_scale_row(self):
        """Test scaling a row by a scalar"""
        mat = Matrix([[1, 2], [3, 4]])
        mat.scale_row(0, S(2))
        assert mat == Matrix([[2, 4], [3, 4]])

    def test_swap_row(self):
        """Test swapping two rows"""
        mat = Matrix([[1, 2], [3, 4]])
        mat.swap_row(0, 1)
        assert mat == Matrix([[3, 4], [1, 2]])

    def test_reduce_row(self):
        """Test row reduction operation"""
        mat = Matrix([[1, 2], [3, 4]])
        mat.reduce_row(1, 3.0, 0)
        assert mat == Matrix([[1, 2], [0, -2]])


class TestPivotOperations:
    """Tests for pivot-related operations

    Covers:
    - get_pivot_row()
    - get_pivot_pos()
    - get_pivot_elements()
    """

    @pytest.mark.parametrize(
        "matrix, col, start_row, expected",
        [
            (Matrix([[1, 2], [0, 1]]), 0, 0, 0),  # First row is pivot
            (Matrix([[0, 2], [1, 1]]), 0, 0, 1),  # Second row is pivot
            (Matrix([[0, 2], [0, 1]]), 0, 0, None),  # No pivot
        ],
    )
    def test_get_pivot_row(self, matrix, col, start_row, expected):
        """Test finding pivot row in a column"""
        assert matrix.get_pivot_row(col, start_row) == expected

    def test_get_pivot_pos(self):
        """Test finding all pivot positions"""
        mat = Matrix([[1, 2, 3], [0, 1, 2], [0, 0, 1]])
        assert mat.get_pivot_pos() == [(0, 0), (1, 1), (2, 2)]

    def test_get_pivot_elements(self):
        """Test getting pivot elements"""
        mat = Matrix([[1, 2, 3], [0, 1, 2], [0, 0, 1]])
        assert mat.get_pivot_elements() == [1, 1, 1]


class TestReductions:
    """Tests for matrix reduction methods

    Covers:
    - ref()
    - rref()
    """

    def test_ref(self):
        """Test row echelon form"""
        mat = Matrix([[1, 2], [3, 4]])
        plu = mat.ref()
        assert plu.U[1, 0] == 0  # Verify lower triangular is zero
        assert plu.U == Matrix([[1, 2], [0, -2]])  # Verify exact REF result

    def test_rref(self):
        """Test reduced row echelon form"""
        mat = Matrix([[1, 2], [3, 4]])
        rref_mat, pivots = mat.rref()
        assert rref_mat[1, 0] == 0  # Below diagonal is zero
        assert rref_mat[0, 1] == 0  # Above pivot is zero
        assert rref_mat == Matrix([[1, 0], [0, 1]])  # Verify exact RREF result
        assert pivots == (0, 1)  # Verify pivot columns


class TestEvaluateCases:
    """Tests for rank of matrix with unknowns
    Covers:
    - find_all_cases
    - evaluate_cases
    - solve
    """

    # TODO
    def test_evaluate_cases(self):
        """Test evaluating special cases"""
        mat = Matrix([[1, 2], [3, 4]])
        rhs = Matrix([[1], [2]])
        # Just verify the method runs without error
        mat.evaluate_cases(rhs)

    def test_solve(self):
        A = Matrix([[1, 2], [3, 4]])
        b = Matrix([[1], [1]])
        x = A.solve(b)
        assert (A @ x - b) == Matrix([0, 0])
