"""Include the following methods in the tests
- nullspace
"""

import pytest
import sympy as sym

from ma1522 import Matrix


class TestSubspacesAssocMatrix:
    def test_nullspace(self):
        mat = Matrix([[1, 2], [3, 6]])
        nullspace = Matrix.from_list(mat.nullspace())
        assert nullspace.cols == 1  # One free variable
        assert (mat @ nullspace).norm() == 0  # Av = 0
