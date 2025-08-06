"""Include the following methods in the tests
Basic Dunder Methods
    - __init__
    - __str__
    - __repr__
    - __eq__
"""

import pytest
import sympy as sym

from ma1522 import Matrix


class TestDunderMethods:
    def test_str_repr_eq(self):
        mat1 = Matrix([[1, 2], [3, 4]])
        mat2 = Matrix([[1, 2], [3, 4]])
        mat3 = Matrix([[5, 6], [7, 8]])
        assert str(mat1) == "Matrix([[1, 2], [3, 4]]), aug_pos: set()"
        assert repr(mat1) == "Matrix([\n[1, 2]\n[3, 4]\n])"
        assert mat1 == mat2
        assert mat1 != mat3

        mat_aug = Matrix([[1, 2], [3, 4]], aug_pos={1})
        assert str(mat_aug) == "Matrix([[1, 2], [3, 4]]), aug_pos: {1}"
        assert repr(mat_aug) == "Matrix([ \n[1 | 2]\n[3 | 4]\n])"
        assert mat1 != mat_aug
