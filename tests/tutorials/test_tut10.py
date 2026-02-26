import sympy as sym

from ma1522 import Matrix, SVD


class TestTutorial10:
    def test_question_1a(self):
        A = Matrix.from_str("0.4 0.1 0.5; 0.2 0.6 0.2; 0.4 0.3 0.3")
        A.simplify(rational=True)
        assert A.is_stochastic()

    def test_question_1b(self):
        A = Matrix.from_str("0.4 0.1 0.5; 0.2 0.6 0.2; 0.4 0.3 0.3")
        A.simplify(rational=True)
        assert A.is_stochastic(verbosity=1)

        x_0 = Matrix.from_str("100; 0; 0")
        PDP = A.diagonalize(reals_only=True)

        x_3 = PDP.P @ PDP.D**3 @ PDP.P.inv() @ x_0

        print(x_3)
        assert x_3.applyfunc(sym.simplify) == Matrix(
            [35, sym.Rational(156, 5), sym.Rational(169, 5)]
        )

    def test_question_1d(self):
        A = Matrix.from_str("0.4 0.1 0.5; 0.2 0.6 0.2; 0.4 0.3 0.3")
        A.simplify(rational=True)
        x_0 = Matrix.from_str("100; 0; 0")
        PDP = A.diagonalize(reals_only=True)

        D_inf = Matrix.diag(0, 0, 1)
        diff = PDP.D**100 - D_inf
        assert sum(abs(diff.evalf(10))) < 1e-10

        steady = (PDP.P @ D_inf @ PDP.P.inv() @ x_0).applyfunc(sym.simplify)
        state = sym.Rational(100, 3)
        assert steady == Matrix([state, state, state])

    def test_question_1e(self):
        A = Matrix.from_str("0.4 0.1 0.5; 0.2 0.6 0.2; 0.4 0.3 0.3")
        A.simplify(rational=True)
        PDP = A.diagonalize(reals_only=True)
        D_inf = Matrix.diag(
            *[
                1 if abs(d - 1) < 1e-9 else 0
                for d in [float(v) for v in PDP.D.diagonal()]
            ]
        )
        alpha, beta, gamma = sym.symbols("alpha beta gamma")
        x = Matrix([alpha, beta, gamma])
        result = (PDP.P @ D_inf @ PDP.P.inv() @ x).applyfunc(sym.simplify)
        steady = (alpha + beta + gamma) / 3
        assert (result - Matrix([steady, steady, steady])).applyfunc(
            sym.simplify
        ) == Matrix.zeros(3, 1)

    def test_question_2(self):
        A = Matrix.from_str("1 0 3; 0 4 0; 0 0 4")
        PDP = A.diagonalize(reals_only=True)
        B = PDP.P @ sym.sqrt(PDP.D) @ PDP.P.inv()
        assert (B**2).applyfunc(sym.simplify) == A

    def test_question_3a(self):
        A = Matrix.from_str("3 1; 1 3")
        PDPT = A.orthogonally_diagonalize()
        P, D = PDPT.P, PDPT.D
        assert (P.T @ P).applyfunc(sym.simplify) == Matrix.eye(2)
        assert (P.T @ A @ P).applyfunc(sym.simplify) == D
        assert set(D.diagonal()) == {2, 4}

    def test_question_3b(self):
        A = Matrix.from_str("2 2 -2; 2 -1 4; -2 4 -1")
        PDPT = A.orthogonally_diagonalize()
        P, D = PDPT.P, PDPT.D
        assert (P.T @ P).applyfunc(sym.simplify) == Matrix.eye(3)
        assert (P.T @ A @ P).applyfunc(sym.simplify) == D

    def test_question_4(self):
        A = Matrix.from_str("1 -2 0 0; -2 1 0 0; 0 0 1 -2; 0 0 -2 1")
        PDPT = A.orthogonally_diagonalize(verbosity=1)
        P, D = PDPT.P, PDPT.D
        assert (P.T @ P).applyfunc(sym.simplify) == Matrix.eye(4)
        assert (P.T @ A @ P).applyfunc(sym.simplify) == D

    def test_question_5a(self):
        A = Matrix.from_str("3 2; 2 3; 2 -2")
        svd = A.singular_value_decomposition(verbosity=1)
        # Reconstruct A from SVD
        reconst = svd.eval()
        reconst.simplify()
        assert reconst == A
        # U has orthonormal columns
        assert (svd.U.T @ svd.U).applyfunc(sym.simplify) == Matrix.eye(svd.U.cols)
        # V is orthogonal
        assert (svd.V.T @ svd.V).applyfunc(sym.simplify) == Matrix.eye(svd.V.cols)

    def test_question_5b(self):
        A = Matrix.from_str("3 2; 2 3; 2 -2")
        svd = A.singular_value_decomposition()
        svdT = SVD(U=svd.V, S=svd.S.T, V=svd.U)
        AT_reconst = (svdT.eval()).applyfunc(sym.simplify)
        assert AT_reconst == A.T

    def test_question_5c(self):
        A = Matrix.from_str("1 0 1; 0 1 1; 1 1 2")
        # A is symmetric and rank-deficient (rank 2)
        assert A.rank() == 2
        svd = A.singular_value_decomposition()
        reconst = (svd.eval()).applyfunc(sym.simplify)
        assert reconst == A

    def test_question_6(self):
        A = Matrix.from_str("-18 13 -4 4; 2 19 -4 12; -14 11 -12 8; -2 21 4 8")
        svd = A.singular_value_decomposition(verbosity=1)
        reconst = (svd.eval()).applyfunc(sym.simplify)
        assert reconst == A
