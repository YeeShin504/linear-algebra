import sympy as sym

from ma1522 import Matrix


class TestTutorial03:
    def test_question_1(self):
        E_1 = Matrix.eye(4).scale_row(1, sym.Rational(1, 2))
        E_2 = Matrix.eye(4).reduce_row(0, 1, 1)
        E_3 = Matrix.eye(4).swap_row(1, 3)
        E_4 = Matrix.eye(4).reduce_row(2, -3, 0)

        A = E_4 @ E_3 @ E_2 @ E_1
        A_inv = E_1**-1 @ E_2**-1 @ E_3**-1 @ E_4**-1
        assert (A @ A_inv).applyfunc(sym.simplify) == Matrix.eye(4)
        assert (A_inv @ A).applyfunc(sym.simplify) == Matrix.eye(4)
        for E in [E_1, E_2, E_3, E_4]:
            assert E.det() != 0

    def test_question_2a(self):
        A = Matrix.from_str("2 -1 2; -6 0 -2; 8 -1 5")
        b = Matrix.from_str("1; 0; 4")

        plu = A.ref()
        assert plu.U.is_echelon
        assert plu.P == Matrix.eye(3)
        # Explicitly assert L and U from notebook
        assert plu.L == Matrix([[1, 0, 0], [-3, 1, 0], [4, -1, 1]])
        assert plu.U == Matrix([[2, -1, 2], [0, -3, 4], [0, 0, 1]])
        assert (plu.L @ plu.U).applyfunc(sym.simplify) == A

        # Solve the unique system A @ x = b
        x = A.solve(b)[0]
        assert (A @ x).applyfunc(sym.simplify) == b

    def test_question_2b(self):
        A = Matrix.from_str("2 -4 4 -2; 6 -9 7 -3; -1 -4 8 0")
        b = Matrix.from_str("0; 0; 5")

        plu = A.ref()
        assert plu.U.is_echelon
        # Explicitly assert L and U from notebook
        assert plu.L == Matrix([[1, 0, 0], [3, 1, 0], [sym.Rational(-1, 2), -2, 1]])
        assert plu.U == Matrix([[2, -4, 4, -2], [0, 3, -5, 3], [0, 0, 0, 5]])
        assert (plu.P.T @ plu.L @ plu.U).applyfunc(sym.simplify) == A

        # Solve underdetermined system A @ x = b (multiple solutions)
        sol = A.solve(b)[0]
        assert len(sol.free_symbols) > 0
        part = sol.subs({v: 0 for v in sol.free_symbols})
        assert (A @ part).applyfunc(sym.simplify) == b

    def test_question_3a(self):
        A = Matrix.from_str("2 -6 6; -4 5 -7; 3 5 -1; -6 4 -8; 8 -3 9")
        plu = A.ref()

        assert plu.U.is_echelon
        assert (plu.P.T @ plu.L @ plu.U).applyfunc(sym.simplify) == A

    def test_question_4(self):
        x = sym.Symbol("x")
        A = Matrix([[-x, 1, 0], [0, -x, 1], [2, -5, 4 - x]])

        det = A.det().factor()
        expected = -((x - 1) ** 2) * (x - 2)
        assert sym.simplify(det - expected) == 0

    def test_question_5(self):
        a, b, c, p, q, r, u, v, w, x = sym.symbols("a b c p q r u v w x")

        lhs_mat = Matrix(
            [
                [a + p * x, b + q * x, c + r * x],
                [p + u * x, q + v * x, r + w * x],
                [u + a * x, v + b * x, w + c * x],
            ]
        )
        rhs_mat = Matrix([[a, b, c], [p, q, r], [u, v, w]])

        ratio = sym.simplify(lhs_mat.det() / rhs_mat.det())
        assert sym.simplify(ratio - (1 + x**3)) == 0

    def test_question_7(self):
        mat = Matrix.from_str("1 5 3; 0 2 -2; 0 1 3")
        rhs = Matrix.from_str("1; 2; 0")

        sol = mat.cramer_solve(rhs)
        # Verify solution satisfies the system
        assert (mat @ sol).applyfunc(sym.simplify) == rhs
        assert sol == Matrix(
            [sym.Rational(-2, 1), sym.Rational(3, 4), sym.Rational(-1, 4)]
        )

    def test_question_8(self):
        A = Matrix.from_str("1 -1 2; 0 2 1; 3 0 6")

        adj = A.adj()
        det = A.det()

        assert det != 0

        A_inv = adj / det
        assert (A @ A_inv).applyfunc(sym.simplify) == Matrix.eye(3)
        assert (A_inv @ A).applyfunc(sym.simplify) == Matrix.eye(3)
        assert (adj / det).applyfunc(sym.simplify) == A.inverse()
