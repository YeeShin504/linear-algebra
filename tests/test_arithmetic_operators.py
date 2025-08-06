"""Include the following methods in the tests
- __abs__
- __add__
- __mul__
- __rmul__
- __neg__
"""

import pytest
import sympy as sym

from ma1522 import Matrix


class TestArithmeticOperators:
    def test_add(self):
        mat1 = Matrix([[1, 2], [3, 4]])
        mat2 = Matrix([[5, 6], [7, 8]])
        assert mat1 + mat2 == Matrix([[6, 8], [10, 12]])
        mat1_aug = Matrix([[1, 2], [3, 4]], aug_pos={1})
        mat2_aug = Matrix([[5, 6], [7, 8]], aug_pos={2})
        assert mat1_aug + mat2_aug == Matrix([[6, 8], [10, 12]], aug_pos={1, 2})

    def test_subtraction(self):
        A = Matrix([[1, 2], [3, 4]])
        B = Matrix([[5, 6], [7, 8]])
        expected = Matrix([[-4, -4], [-4, -4]])
        assert A - B == expected

    def test_scalar_multiplication(self):
        A = Matrix([[1, 2], [3, 4]])
        expected = Matrix([[2, 4], [6, 8]])
        assert 2 * A == expected
        assert A * 2 == expected

    def test_matrix_multiplication(self):
        A = Matrix([[1, 2], [3, 4]])
        B = Matrix([[5, 6], [7, 8]])
        expected = Matrix([[19, 22], [43, 50]])
        assert A @ B == expected

    def test_neg(self):
        mat = Matrix([[1, 2], [3, 4]])
        assert -mat == Matrix([[-1, -2], [-3, -4]])
        mat_aug = Matrix([[-1, -2], [-3, -4]], aug_pos={1})
        assert -mat_aug == Matrix([[1, 2], [3, 4]], aug_pos={1})
