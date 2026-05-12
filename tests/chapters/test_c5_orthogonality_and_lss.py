"""Include the following methods in the tests
- orthogonal_complement
- is_vec_orthogonal
- is_mat_orthogonal
- orthogonal_decomposition
- proj_comp
- norm_comp
- gram_schmidt
- QRdecomposition
- solve_least_squares
- create_vander
- apply_vander
"""

import sympy as sym

from ma1522 import Matrix, PartGen, VecDecomp


class TestChapter5:
    def test_orthogonal_complement(self):
        mat = Matrix([[1, 0], [0, 1], [0, 0]])
        ortho_comp = mat.orthogonal_complement()
        assert ortho_comp.cols == 1
        assert (mat.T @ ortho_comp).is_zero_matrix

        mat = Matrix([[1, 0], [0, 1]])
        orth_comp = mat.orthogonal_complement()
        assert orth_comp == Matrix([])

    def test_is_vec_orthogonal(self):
        mat = Matrix([[1, 0], [0, 1]])
        assert mat.is_vec_orthogonal() is True
        mat = Matrix([[2, 0], [0, 3]])
        assert mat.is_vec_orthogonal() is True
        mat = Matrix([[1, 1], [1, 1]])
        assert mat.is_vec_orthogonal() is False

    def test_is_mat_orthogonal(self):
        mat = Matrix([[1, 0], [0, 1]])
        assert mat.is_mat_orthogonal() is True
        mat = Matrix([[2, 0], [0, 3]])
        assert mat.is_mat_orthogonal() is False
        mat = Matrix([[1, 1], [1, 0]])
        assert mat.is_mat_orthogonal() is False

    def test_orthogonal_decomposition(self):
        """Test vector decomposition into orthogonal components"""
        A = Matrix([[1, 0], [1, 1]])
        b = Matrix([[-1], [2]])
        decomp = b.orthogonal_decomposition(A)
        assert isinstance(decomp, VecDecomp)
        assert (decomp.proj + decomp.norm) == b

    def test_proj_comp(self):
        to = Matrix([[1, 0], [0, 1]])
        vec = Matrix([[1], [1]])
        proj = vec.proj_comp(to)
        assert proj == vec

    def test_norm_comp(self):
        to = Matrix([[1, 0], [0, 1]])
        vec = Matrix([[1], [1]])
        norm = vec.norm_comp(to)
        assert norm == Matrix([[0], [0]])

    def test_gram_schmidt(self):
        mat = Matrix([[1, 1], [1, 0]])
        ortho_mat = mat.gram_schmidt(factor=True)
        assert ortho_mat.full.is_vec_orthogonal() is True  # type: ignore
        assert ortho_mat.full.select_cols(0).dot(ortho_mat.full.select_cols(1)) == 0  # type: ignore

    def test_gram_schmidt_complex_uses_hermitian_projection(self):
        mat = Matrix([[1, sym.I], [sym.I, 0]])
        ortho_mat = mat.gram_schmidt(factor=False, verbosity=0)
        assert isinstance(ortho_mat, Matrix)
        assert (
            ortho_mat.select_cols(0).dot(ortho_mat.select_cols(1), hermitian=True) == 0
        )

    def test_QRdecomposition(self):
        mat = Matrix([[1, 1], [1, 0]])
        q, r = mat.QRdecomposition()
        # Verify Q is orthogonal
        assert (q.T @ q) == Matrix.eye(2)
        # Verify R is upper triangular
        assert r[1, 0] == 0

    def test_QRdecomposition_verbose(self):
        mat = Matrix([[1, 1], [1, 0]])
        q, r = mat.QRdecomposition(verbosity=1)
        assert (q @ r).equals(mat)

    def test_full_QRdecomposition_rectangular(self):
        mat = Matrix([[1, 0], [0, 1], [1, 1]])
        q, r = mat.QRdecomposition(full=True)
        assert q.shape == (3, 3)
        assert r.shape == (3, 2)
        assert (q @ r).equals(mat)

    def test_solve_least_squares(self):
        """Test least squares solution"""
        A = Matrix([[0, 1], [1, 1], [2, 1]])
        b = Matrix([[6], [0], [0]])
        x = A.solve_least_squares(b)
        assert (A @ x - b).norm() == sym.sqrt(6)  # Minimized error

    def test_solve_least_squares_fallback_part_gen(self):
        A = Matrix([[1, 0], [0, 0]])
        b = Matrix([[1], [2]])
        result = A.solve_least_squares(b, verbosity=0, matrices=2)
        assert isinstance(result, PartGen)
        assert A @ result.part_sol == Matrix([[1], [0]])
