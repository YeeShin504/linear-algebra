import pytest
import sympy as sym

from ma1522 import Matrix


class TestTutorial06:
    def test_question_1a(self):
        S = Matrix.from_str("1 2 -1; 0 2 1; 0 -1 3").T
        basis = S.get_linearly_independent_vectors()
        assert basis.shape == (3, 3)
        assert S.rank() == 3

    def test_question_1b(self):
        S = Matrix.from_str("1 2 -1; 0 2 1; 0 -1 3").T
        w = Matrix.from_str("1; 1; 1")
        coords = w.coords_relative(S)
        # Check: S @ coords == w
        assert (S @ coords).applyfunc(sym.simplify) == w
        assert coords.applyfunc(sym.simplify) == Matrix(
            [
                sym.Integer(1),
                sym.Rational(-1, 7),
                sym.Rational(5, 7),
            ]
        )

    def test_question_1c(self):
        S = Matrix.from_str("1 2 -1; 0 2 1; 0 -1 3").T
        T = Matrix.from_str("1 5 4; -1 3 7; 2 2 4").T
        P_T_to_S = T.transition_matrix(S)
        assert (S @ P_T_to_S).applyfunc(sym.simplify) == T

    def test_question_1d(self):
        S = Matrix.from_str("1 2 -1; 0 2 1; 0 -1 3").T
        T = Matrix.from_str("1 5 4; -1 3 7; 2 2 4").T
        P_T_to_S = T.transition_matrix(S)
        P_S_to_T = S.transition_matrix(T)
        assert P_S_to_T.applyfunc(sym.simplify) == P_T_to_S.inverse()

    def test_question_1e(self):
        T = Matrix.from_str("1 5 4; -1 3 7; 2 2 4").T
        w = Matrix.from_str("1; 1; 1")
        w_T = w.coords_relative(T)
        assert (T @ w_T).applyfunc(sym.simplify) == w

    def test_question_3a(self):
        A = Matrix.from_str("1 -1 1; 1 1 -1; -1 -1 1")
        b = Matrix.from_str("2; 1; 0")
        # b is NOT in the column space of A
        with pytest.raises(ValueError):
            b.coords_relative(A)

    def test_question_3b(self):
        A = Matrix.from_str("1 9 1; -1 3 1; 1 1 1")
        b = Matrix.from_str("5 1 -1")
        scalars = b.T.coords_relative(A.T)
        assert scalars == Matrix([1, -3, 1])
        assert (A.T @ scalars).applyfunc(sym.simplify) == b.T

    def test_question_3c(self):
        A = Matrix.from_str("1 2 0 1; 0 1 2 1; 1 2 1 3; 0 1 2 2")
        assert A.rank() == 4
        assert A.is_same_subspace()
        assert A.T.is_same_subspace()

    def test_question_4a(self):
        A = Matrix.from_str("1 2 5 3; 1 -4 -1 -9; -1 0 -3 1; 2 1 7 0; 0 1 1 2")
        r = A.rank()
        n = A.nullity()
        assert r + n == A.cols
        rowspace = A.get_linearly_independent_vectors(colspace=False)
        assert rowspace.rows == r
        colspace = A.get_linearly_independent_vectors(colspace=True)
        assert colspace.cols == r
        nullvecs = A.nullspace()
        assert len(nullvecs) == n

    def test_question_4b(self):
        A = Matrix.from_str("1 3 7; 2 1 8; 3 -5 -1; 2 -2 2; 1 1 5")
        rowspace = A.simplify_basis(colspace=False)
        assert rowspace.shape == (3, 3)
        assert rowspace.rank() == 3
        colspace = A.simplify_basis(colspace=True)
        assert colspace.shape == (5, 3)
        assert colspace.rank() == 3
        with pytest.warns(UserWarning):  # trivial null space
            nullspace = A.nullspace(verbosity=2)
            assert len(nullspace) == 0

    def test_question_5a(self):
        W = Matrix.from_str("1 -2 0 0 3; 2 -5 -3 -2 6; 0 5 15 10 0; 2 1 15 8 6").T
        basis = W.get_linearly_independent_vectors(colspace=True)
        assert basis.shape == (5, 3)
        assert basis.rank() == 3

    def test_question_5c(self):
        W = Matrix.from_str("1 -2 0 0 3; 2 -5 -3 -2 6; 0 5 15 10 0; 2 1 15 8 6").T
        extended = W.extend_basis()
        assert extended.shape == (5, 5)
        assert extended.rank() == 5

    def test_question_6(self):
        V = Matrix.from_str("1 0 1 3; 2 -1 0 1; -1 3 5 12; 0 1 2 5; 3 -1 1 4").T
        S_prime = V.get_linearly_independent_vectors(colspace=True)
        assert S_prime.is_same_subspace(V)
        assert S_prime.shape == (4, 2)
        assert S_prime.rank() == 2
