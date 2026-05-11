"""Include the following methods in the tests
- standard_matrix
"""

import pytest
import sympy as sym

from ma1522 import Matrix, PartGen


class TestLinearTransformations:
    def test_standard_matrix(self):
        """Test standard matrix representation"""
        standard_matrix = Matrix.create_rand_matrix(3, 3)
        input_vectors = Matrix.create_rand_matrix(3, 3)
        output_vectors = standard_matrix @ input_vectors
        sol = Matrix.standard_matrix(input_vectors, output_vectors)[0]
        assert sol == standard_matrix

    def test_standard_matrix_with_symbolic_input_vectors(self):
        x = sym.symbols("x")
        input_vectors = Matrix([[x, 0], [0, 1]])
        output_vectors = Matrix([[2 * x, 0], [0, 3]])

        sol = input_vectors.standard_matrix(output_vectors)[0]

        assert sol == Matrix([[2, 0], [0, 3]])

    def test_standard_matrix_returns_part_gen(self):
        input_vectors = Matrix([[1], [0]])
        output_vectors = Matrix([[2], [0]])

        sol = input_vectors.standard_matrix(output_vectors, matrices=2)[0]

        assert isinstance(sol, PartGen)
        assert sol.part_sol @ input_vectors == output_vectors

    def test_standard_matrix_raises_when_no_solution_exists(self):
        with pytest.raises(ValueError, match="No solution found"):
            Matrix.zeros(2, 1).standard_matrix(Matrix.ones(2, 1))

    def test_standard_matrix_rejects_invalid_return_mode(self):
        with pytest.raises(ValueError, match="Invalid value for matrices"):
            Matrix.eye(2).standard_matrix(Matrix.eye(2), matrices=3)
