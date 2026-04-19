import sympy as sym

from ma1522 import Matrix


class TestTutorial05:
    def test_question_1a(self):
        S = Matrix.from_str("2 -1 0; 0 3 2; 2 4 3; 3 6 6").T
        assert not S.is_linearly_independent()
        assert S.rank() == 3

    def test_question_1b(self):
        S = Matrix.from_str("1 1 0; 3 4 2").T
        assert S.is_linearly_independent()
        assert S.rank() == 2

    def test_question_1c(self):
        # Contains the zero vector → always linearly dependent
        S = Matrix.from_str("1 1 0; 3 4 2; 0 0 0").T
        assert not S.is_linearly_independent()

    def test_question_1d(self):
        S = Matrix.from_str("1 0 0; 0 1 1; 1 2 -1").T
        assert S.is_linearly_independent()
        assert S.rank() == 3

    def test_question_2d(self):
        a, b, c = sym.symbols("a b c")
        u, v, w = sym.symbols("u v w")
        v1, v2, v3 = u, u + v, u + v + w
        sols = sym.solve(a * v1 + b * v2 + c * v3, (a, b, c))
        assert sols == {a: 0, b: 0, c: 0}

    def test_question_2e(self):
        a, b, c, d = sym.symbols("a b c d")
        u, v, w = sym.symbols("u v w")
        v1, v2, v3, v4 = u + v, v + w, u + w, u + v + w
        sols = sym.solve(a * v1 + b * v2 + c * v3 + d * v4, (a, b, c, d))
        assert len(sols) == 1
        half = sym.Rational(1, 2)
        assert sols[0] == (-half * d, -half * d, -half * d, d)

    def test_question_3a(self):
        a, b, c, d = sym.symbols("a b c d")
        V = Matrix([a + b, a + c, c + d, b + d])
        vecs = V.sep_unk()

        expected_a = Matrix([1, 1, 0, 0])
        expected_b = Matrix([1, 0, 0, 1])
        expected_c = Matrix([0, 1, 1, 0])
        expected_d = Matrix([0, 0, 1, 1])

        assert vecs[a] == expected_a
        assert vecs[b] == expected_b
        assert vecs[c] == expected_c
        assert vecs[d] == expected_d

        V_mat = Matrix.from_list([*vecs.values()])
        basis = V_mat.get_linearly_independent_vectors()
        assert basis.shape == (4, 3)
        assert basis.rank() == 3

    def test_question_3b(self):
        V = Matrix.from_str("1 0 -1; -1 2 3; 0 3 0; 1 -1 1").T
        basis = V.get_linearly_independent_vectors()
        assert basis.shape == (3, 3)
        assert basis.rank() == 3

    def test_question_3c(self):
        mat = Matrix.from_str("1 0 1 1 -1; 0 1 1 2 1; 1 1 2 1 -2")
        null = mat.nullspace()
        assert len(null) == mat.cols - mat.rank()
        for v in null:
            assert (mat @ v).applyfunc(sym.simplify) == Matrix.zeros(mat.rows, 1)

    def test_question_4(self):
        a = sym.Symbol("a", real=True)
        U = Matrix([[a, 1, -1], [-1, a, 1], [1, -1, a]]).T
        singular_vals = sym.solve(U.det(), a)
        assert len(singular_vals) == 1
        assert singular_vals[0] == 0

    def test_question_5b(self):
        U = Matrix.from_str("1 1 1 1; 1 2 2 1").T
        V = Matrix.from_str("1 0 1 0; 1 0 2 -1").T
        uv_basis = U.row_join(V).get_linearly_independent_vectors()
        assert uv_basis.shape == (4, 3)
        assert uv_basis.rank() == 3

    def test_question_5e(self):
        U = Matrix.from_str("1 1 1 1; 1 2 2 1").T
        V = Matrix.from_str("1 0 1 0; 1 0 2 -1").T
        intersection = U.intersect_subspace(V)
        assert intersection.shape == (4, 1)
        for j in range(intersection.cols):
            col = intersection.select_cols(j)
            aug_U = U.row_join(col)
            aug_V = V.row_join(col)
            assert aug_U.rank() == U.rank()
            assert aug_V.rank() == V.rank()
