import pytest
import sympy as sym
from ma1522 import Matrix, SVD, PDP

pytest.skip("Skipping this entire file because it's under construction", allow_module_level=True)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _reconstruction_norm(svd: SVD, A: Matrix) -> float:
    """Numerical Frobenius norm of U S V^T - A."""
    U = svd.U.evalf()
    S = svd.S.evalf()
    V = svd.V.evalf()
    diff = U @ S @ V.T - A.evalf()
    return float(sym.re(diff.norm()))

def _diag_reconstruction_norm_no_inv(pdp: PDP, A: Matrix) -> float:
    """Verify A = P D P^{-1} without inverting P: check ||A P - P D|| numerically."""
    P = pdp.P.evalf()
    D = pdp.D.evalf()
    A_num = A.evalf()
    diff = A_num @ P - P @ D
    return float(sym.re(diff.norm()))

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestSVDDecomposition:
    """Regression and edge-case tests for Singular Value Decomposition."""
    
    def test_4x3_irrational_reconstruction(self):
        """Test SVD on matrix with irrational eigenvalues (Regression)."""
        A = Matrix([[1, -2, -1], [2, 0, 1], [2, -4, 2], [4, 0, 0]])
        # Note: verify=False because the symbolic norm check is extremely expensive for this matrix
        svd = A.singular_value_decomposition(verbosity=0, verify=False)
        assert _reconstruction_norm(svd, A) < 1e-8

    def test_rank_deficient_reconstruction(self):
        """Test SVD on rank-deficient matrix."""
        A = Matrix([[1, 2], [2, 4]])
        svd = A.singular_value_decomposition(verbosity=0, verify=False)
        assert _reconstruction_norm(svd, A) < 1e-8

    def test_fast_svd_tol_none(self):
        """Verify fast_svd handles tol=None correctly (Regression)."""
        mat = Matrix([[1.1, 1], [-0.1, 0]])
        svd = mat.fast_svd(option="sym", identify=True, tol=None)
        assert isinstance(svd, SVD)

class TestDiagonalization:
    """Regression tests for Matrix diagonalization."""

    def test_ata_irrational_diagonalization(self):
        """Test diagonalization of A^T A with irrational eigenvalues (Regression)."""
        A = Matrix([[1, -2, -1], [2, 0, 1], [2, -4, 2], [4, 0, 0]])
        ATA = A.T @ A
        pdp = ATA.diagonalize(verbosity=0)
        assert _diag_reconstruction_norm_no_inv(pdp, ATA) < 1e-8

class TestVectorSpaces:
    """Tests for Vector Space operations (Gram-Schmidt, Transition Matrix)."""

    def test_transition_matrix_indexing(self):
        """Verify transition_matrix uses correct slicing."""
        B = Matrix([[1, 0], [0, 1]])
        C = Matrix([[1, 1], [1, -1]])
        T = B.transition_matrix(to=C, verbosity=0)
        assert T == C.inv() @ B

    def test_gram_schmidt_orthonormal(self):
        """Verify Gram-Schmidt returns orthonormal vectors."""
        v1 = Matrix([1, 1, 0])
        v2 = Matrix([1, 0, 0])
        res = Matrix.from_list([v1, v2]).gram_schmidt(factor=False, verbosity=0)
        assert res.col(0) == v1.normalized()
        assert res.col(1) == (v2 - (v2.dot(v1)/v1.dot(v1))*v1).normalized()

class TestNegativeDecompositions:
    """Negative tests for decompositions (invalid inputs)."""
    
    def test_diagonalize_non_diagonalizable(self):
        """Verify diagonalize raises error for non-diagonalizable matrices."""
        # Shear matrix [[1, 1], [0, 1]] is not diagonalizable
        A = Matrix([[1, 1], [0, 1]])
        with pytest.raises(Exception): 
            A.diagonalize(verbosity=0)
