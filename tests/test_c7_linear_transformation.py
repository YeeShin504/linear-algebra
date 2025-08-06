"""Include the following methods in the tests
- standard_matrix
"""

import pytest
import sympy as sym

from ma1522 import Matrix


class TestLinearTransformations:
    def test_standard_matrix(self):
        """Test standard matrix representation"""
        A = Matrix([[1, 0], [0, 1]])
        b = Matrix([[1], [1]])
        X = A.standard_matrix(b)
        assert (X @ A) == b
