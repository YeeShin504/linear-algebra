"""Include the following methods in the tests
- inverse
- elem
- adj
- column_constraints
"""

import pytest
import sympy as sym

from ma1522 import Matrix


class TestMatrixInverse:
    """Tests for Matrix inverse operations

    Covers:
    - Left/right inverses
    - Edge cases (singular matrices)
    - Special cases
    """

    def test_left_inverse(self):
        """Test left inverse for full column rank matrix"""
        mat = Matrix([[1, 0], [0, 1], [1, 1]])
        inv = mat.inverse(option="left", verbosity=2)
        assert inv is not None
        assert (inv @ mat) == Matrix.eye(2)

    def test_right_inverse(self):
        """Test right inverse for full row rank matrix"""
        mat = Matrix([[1, 0, 1], [0, 1, 1]])
        inv = mat.inverse(verbosity=2)
        assert inv is not None
        assert (mat @ inv) == Matrix.eye(2)

    def test_square_inverse(self):
        """Test inverse for square matrix"""
        mat = Matrix([[1, 2], [3, 4]])
        inv = mat.inverse(option="both", verbosity=2)
        assert (inv @ mat) == Matrix.eye(2)
        assert (mat @ inv) == Matrix.eye(2)

    def test_singular_matrix(self):
        """Test inverse fails for singular matrix"""
        mat = Matrix([[1, 2], [2, 4]])  # Determinant = 0
        with pytest.raises(ValueError):
            mat.inverse()

    def test_1x1_matrix(self):
        """Test special case for 1x1 matrix"""
        mat = Matrix([[5]])
        inv = mat.inverse()
        assert (inv @ mat) == Matrix.eye(1)
        assert (mat @ inv) == Matrix.eye(1)
        assert inv == Matrix([[sym.Rational(1, 5)]])


class TestAdjointOperations:
    """Tests for adjoint/adjugate operations"""

    def test_adjoint(self):
        """Test classical adjoint calculation"""
        mat = Matrix([[1, 2], [3, 4]])
        adj = mat.adj()
        assert adj == Matrix([[4, -2], [-3, 1]])

    def test_adjoint_property(self):
        """Test A * adj(A) = det(A)*I"""
        mat = Matrix([[1, 2], [3, 4]])
        assert (mat @ mat.adj()) == Matrix.diag(mat.det(), mat.det())


class TestElem:
    """Tests for Matrix.elem() method"""

    def test_elem(self):
        """Test elem() returns identity matrix of correct size"""
        mat = Matrix([[1, 2], [3, 4]])
        identity_matrix = mat.elem()
        assert identity_matrix == Matrix.eye(2)
        assert identity_matrix.shape == (2, 2)


class TestColumnConstraints:
    """Tests for column constraint calculations"""

    def test_column_constraints(self):
        """Test column constraints for non-invertible matrix"""
        x1, x2 = Matrix.create_unk_matrix(2, 1, "x")
        mat = Matrix([[1, 2], [2, 4]])
        constraints = mat.column_constraints()
        half_x2 = sym.Mul(x2, sym.Rational(1, 2))
        x1_minus_half = sym.Add(x1, sym.Mul(half_x2, -1))  # Equivalent to x1 - half_x2
        expected_rref = Matrix([[1, 2, half_x2], [0, 0, x1_minus_half]], aug_pos=2)
        assert constraints == expected_rref
