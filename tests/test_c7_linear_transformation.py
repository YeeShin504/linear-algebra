"""Include the following methods in the tests
- standard_matrix
"""

import pytest
import sympy as sym

from ma1522 import Matrix


class TestLinearTransformations:
    def test_standard_matrix(self):
        """Test standard matrix representation"""
        standard_matrix = Matrix.create_rand_matrix(3, 3)
        input_vectors = Matrix.create_rand_matrix(3, 3)
        output_vectors = standard_matrix @ input_vectors
        sol = Matrix.standard_matrix(input_vectors, output_vectors)[0]
        assert sol == standard_matrix
