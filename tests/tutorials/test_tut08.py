import gen
import pytest
import sympy as sym

from ma1522 import Matrix
from ma1522.custom_types import ScalarFactor


class TestTutorial08:
    def test_question_1a(self):
        U = Matrix.from_str("1 1 1 1; 1 -1 1 0; 1 1 -1 -1; 1 2 0 1").T
        Q_SF = U.gram_schmidt(factor=True, verbosity=1)
        assert isinstance(Q_SF, ScalarFactor)
        Q = U.gram_schmidt(factor=False, verbosity=1)
        assert isinstance(Q, Matrix)
        assert (Q - Q_SF.eval()).applyfunc(sym.simplify) == Matrix.zeros(*Q.shape)

        assert (Q.T @ Q).applyfunc(lambda x: x.simplify()) == Matrix.eye(4)
        assert Q.cols == 4

    def test_question_1b(self):
        U = Matrix.from_str("1 2 2 1; 1 2 1 0; 1 0 1 0; 1 0 2 1").T
        with pytest.warns(UserWarning):  # linearly dependent set
            Q = U.gram_schmidt(factor=False, verbosity=1)
            assert not Q.is_mat_orthogonal()

    def test_question_2a(self):
        A = Matrix.from_str("0 1 1 0; 1 -1 1 -1; 1 0 1 0; 1 1 1 1")
        b = Matrix.from_str("6; 3; -1; 1")
        # System is inconsistent: b is NOT in col(A)
        assert not b.is_subspace_of(A)

    def test_question_2b(self):
        A = Matrix.from_str("0 1 1 0; 1 -1 1 -1; 1 0 1 0; 1 1 1 1")
        b = Matrix.from_str("6; 3; -1; 1")
        sol = A.solve_least_squares(b)
        part, gen = sol.sep_part_gen()
        assert part == Matrix([-6, -1, 7, 0])
        assert len(gen.free_symbols) == 1
        free = gen.free_symbols.pop()
        assert gen == free * Matrix([-1, -1, 1, 1])

    def test_question_2c(self):
        A = Matrix.from_str("0 1 1 0; 1 -1 1 -1; 1 0 1 0; 1 1 1 1")
        b = Matrix.from_str("6; 3; -1; 1")
        sol = A.solve_least_squares(b)
        proj = A @ sol

        assert proj == Matrix([6, 2, 1, 0])

    def test_question_3(self):
        A = Matrix.from_str("1 0.01; 1 0.012; 1 0.015; 1 0.02")
        b = Matrix.from_str("2.75E-4; 3.31E-4; 3.92E-4; 4.95E-4")

        A.simplify(tolerance=1e-10)
        b.simplify(tolerance=1e-10)

        assert not b.is_subspace_of(A, verbosity=2)

        sol = A.solve_least_squares(b)
        assert len(sol.free_symbols) == 0
        assert A.T @ A @ sol == A.T @ b

    def test_question_4(self):
        x = Matrix.from_str("4 4.5 5 5.5 6 6.5 7 8 8.5").T
        y = Matrix.from_str(
            "0.8651 0.4828 2.590 -4.389 -7.858 3.103 7.456 0.0965 4.326"
        ).T
        x.simplify(tolerance=1e-10)
        y.simplify(tolerance=1e-10)

        N = Matrix.create_vander(num_rows=x.rows, num_cols=5)
        N = N.apply_vander(x)
        sol = N.solve_least_squares(y, verbosity=2)
        assert len(sol.free_symbols) == 0
        assert N.T @ N @ sol == N.T @ y

    def test_question_5a(self):
        A = Matrix.from_str("1 1 0; 1 1 0; 1 1 1; 0 1 1")
        QR = A.QRdecomposition()
        Q, R = QR.Q, QR.R

        gram = (Q.T @ Q).applyfunc(lambda x: x.simplify())
        assert gram == Matrix.eye(Q.cols)
        assert R.is_upper
        assert R.shape == (A.cols, A.cols)
        assert (Q @ R).applyfunc(sym.simplify) == A

    def test_question_5b(self):
        A = Matrix.from_str("1 1 0; 1 1 0; 1 1 1; 0 1 1")
        b = Matrix.from_str("1; 1; 0; 0")
        QR = A.QRdecomposition()
        Q, R = QR.Q, QR.R

        sol_qr = (R.inverse() @ Q.T @ b).applyfunc(sym.simplify)
        sol_lsq = A.solve_least_squares(b).sep_part_gen()[0].applyfunc(sym.simplify)
        assert sol_qr == sol_lsq
        assert (A.T @ A @ sol_qr).applyfunc(sym.simplify) == (A.T @ b).applyfunc(
            sym.simplify
        )
