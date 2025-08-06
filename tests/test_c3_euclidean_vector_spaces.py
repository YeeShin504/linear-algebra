"""Include the following methods in the tests
- normalized
- is_linearly_independent
- simplify_basis
- extend_basis
- intersect_subspace
- is_same_subspace
- coords_relative
- transition_matrix
"""

import pytest

import sympy as sym

from ma1522 import Matrix


class TestEuclideanVectorSpaces:
    """Tests for euclidean vector space operations

    Covers:
    - Vector normalization
    - Linear independence
    - Same subspace
    - Relative coordinates
    - Transition matrix
    """

    def test_normalized(self):
        mat = Matrix([[1, 1], [0, 4]])
        normalized_mat = mat.normalized()
        assert normalized_mat == Matrix([[1 / sym.sqrt(2), 1 / sym.sqrt(2)], [0, 1]])  # type: ignore

    def test_is_linearly_independent(self):
        mat = Matrix([[1, 2], [3, 4]])
        assert mat.is_linearly_independent() is True
        mat = Matrix([[1, 2], [2, 4]])
        assert mat.is_linearly_independent() is False

    def test_is_same_subspace(self):
        mat1 = Matrix([[1, 0], [0, 1]])
        mat2 = Matrix([[1, 2], [3, 4]])
        assert mat1.is_same_subspace(mat2) is True

    def test_transition_matrix(self):
        mat1 = Matrix([[1, 0], [0, 1]])
        mat2 = Matrix([[2, 0], [0, 2]])
        transition = mat1.transition_matrix(to=mat2, verbosity=2)
        assert transition == Matrix([[2, 0], [0, 2]])


class TestSubspaceOperations:
    """Tests for subspace operations

    Covers:
    - Simplify column/row space operations
    - Basis extension
    - Subspace intersections
    """

    def test_simplify_basis(self):
        """Test simplify basis"""
        mat = Matrix([[1, 2, 2, 5], [3, 4, 6, 13], [0, 0, 0, 0]])
        col_basis = mat.simplify_basis(colspace=True)
        row_basis = mat.simplify_basis(colspace=False)
        assert col_basis == Matrix([[1, 0], [0, 1], [0, 0]])
        assert row_basis == Matrix([[1, 0, 2, 3], [0, 1, 0, 1]])

    def test_extend_basis(self):
        """Test extending to a basis"""
        mat = Matrix([[1, 1], [1, 1]])
        extended = mat.extend_basis()
        assert extended.cols == 2  # Original space dimension
        assert extended.rank() == 2  # Now full rank

    def test_subspace_intersection(self):
        """Test subspace intersection"""
        A = Matrix([[1, 0], [0, 1]])
        B = Matrix([[1, 1], [1, 1]])
        intersection = A.intersect_subspace(B)
        assert intersection.cols == 1  # One common direction


class TestCoordinatesBasis:
    def test_transition_matrix(self):
        """Test transition matrix between bases"""
        A = Matrix([[1, 0], [0, 1]])
        B = Matrix([[1, 1], [1, -1]])
        P = A.transition_matrix(B)
        assert (P @ A) == B

    def test_coordinate_transformation(self):
        """Test coordinate transformation"""
        A = Matrix([[1, 0], [0, 1]])
        v = Matrix([[1], [1]])
        coords = A.coords_relative(v)
        assert (A @ coords) == v
