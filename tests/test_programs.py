"""Run the TI-BASIC programs in the simulator and check results mathematically.

    python3 -m unittest discover -s tests -v
"""

import random
import sys
import unittest
from pathlib import Path

TI = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(TI / "tools"))
from sim import Calc  # noqa: E402

SRC = TI / "src"
TOL = 1e-7


def mul(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(len(b))) for j in range(len(b[0]))] for i in range(len(a))]


def T(a):
    return [list(r) for r in zip(*a)]


def eye(n):
    return [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]


def close(a, b, tol=TOL):
    if len(a) != len(b) or len(a[0]) != len(b[0]):
        return False
    return all(abs(x - y) < tol for r, s in zip(a, b) for x, y in zip(r, s))


def rank(a):
    m = [list(map(float, r)) for r in a]
    r = 0
    for c in range(len(m[0])):
        p = max(range(r, len(m)), key=lambda i: abs(m[i][c]), default=None)
        if p is None or abs(m[p][c]) < 1e-9:
            continue
        m[r], m[p] = m[p], m[r]
        for i in range(len(m)):
            if i != r:
                f = m[i][c] / m[r][c]
                m[i] = [x - f * y for x, y in zip(m[i], m[r])]
        r += 1
        if r == len(m):
            break
    return r


def calc(**mats):
    c = Calc(SRC)
    for k, v in mats.items():
        c.set_matrix(k, v)
    return c


class TestRREF(unittest.TestCase):
    def test_rref_singular_with_steps(self):
        A = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
        for steps in (0, 1):
            c = calc(A=A).run("RREF", [steps, 0])
            self.assertTrue(close(c.matrix("C"), [[1, 0, -1], [0, 1, 2], [0, 0, 0]]))
            self.assertEqual(c.reals["R"], 2)
            self.assertEqual(c.list("PIV"), [1, 2])

    def test_ref_only(self):
        A = [[0, 2, 4], [1, 1, 1], [2, 4, 6]]
        c = calc(A=A).run("RREF", [1, 1])
        U = c.matrix("C")
        for i in range(3):
            for j in range(i):
                self.assertAlmostEqual(U[i][j], 0)
        self.assertEqual(c.reals["R"], 2)

    def test_zero_matrix(self):
        c = calc(A=[[0, 0], [0, 0]]).run("RREF", [0, 0])
        self.assertEqual(c.reals["R"], 0)


class TestSolve(unittest.TestCase):
    def test_general_solution(self):
        A = [[1, 1, 1, 1], [1, 2, 3, 4], [2, 3, 4, 5]]
        b = [[6], [14], [20]]
        c = calc(A=A, B=b).run("SOLVE", [1])
        xp, N = c.matrix("C"), c.matrix("D")
        self.assertEqual(c.reals["F"], 2)
        self.assertTrue(close(mul(A, xp), b))
        self.assertTrue(close(mul(A, N), [[0, 0]] * 3))
        self.assertEqual(rank(N), 2)

    def test_unique(self):
        A = [[2, 1], [1, 3]]
        c = calc(A=A, B=[[3], [5]]).run("SOLVE", [0])
        self.assertTrue(close(c.matrix("C"), [[0.8], [1.4]]))
        self.assertEqual(c.reals["F"], 0)

    def test_inconsistent(self):
        c = calc(A=[[1, 2], [2, 4]], B=[[1], [3]]).run("SOLVE", [0])
        self.assertIn("INCONSISTENT", c.screen)


class TestLU(unittest.TestCase):
    def check(self, A, swaps):
        c = calc(A=A).run("LU", [1])
        L, U, P = c.matrix("C"), c.matrix("D"), c.matrix("E")
        self.assertTrue(close(mul(P, mul(L, U)), A))
        n = len(A)
        for i in range(n):
            self.assertAlmostEqual(L[i][i], 1)
            for j in range(i + 1, n):
                self.assertAlmostEqual(L[i][j], 0)
        self.assertEqual(close(P, eye(n)), not swaps)

    def test_no_swap(self):
        self.check([[2, 4, 1], [4, 9, 3], [-2, -1, 7]], False)

    def test_swap(self):
        self.check([[0, 1, 2], [1, 2, 3], [2, 7, 1]], True)

    def test_rectangular(self):
        self.check([[1, 2, 3, 4], [2, 4, 7, 9], [3, 6, 10, 15]], False)


class TestInverse(unittest.TestCase):
    def test_invertible(self):
        A = [[1, 2, 3], [0, 1, 4], [5, 6, 0]]
        c = calc(A=A).run("INVERSE", [1])
        self.assertAlmostEqual(c.reals["O"], 1)
        self.assertTrue(close(mul(A, c.matrix("C")), eye(3)))
        self.assertTrue(close(c.matrix("D"), c.matrix("C")))  # det = 1
        self.assertTrue(close(T(c.matrix("E")), c.matrix("D")))

    def test_singular_adjoint(self):
        A = [[1, 2], [2, 4]]
        c = calc(A=A).run("INVERSE", [0])
        self.assertEqual(c.reals["O"], 0)
        self.assertTrue(close(c.matrix("D"), [[4, -2], [-2, 1]]))

    def test_cramer(self):
        A = [[2, 1, -1], [-3, -1, 2], [-2, 1, 2]]
        c = calc(A=A, B=[[8], [-11], [-3]]).run("CRAMER")
        self.assertTrue(close(c.matrix("C"), [[2], [3], [-1]]))


class TestSubspaces(unittest.TestCase):
    def test_spaces(self):
        A = [[1, 2, 0, 1], [2, 4, 1, 4], [3, 6, 2, 7]]
        c = calc(A=A).run("SPACES")
        self.assertEqual(c.reals["R"], 2)
        self.assertTrue(close(c.matrix("C"), [[1, 0], [2, 1], [3, 2]]))
        self.assertTrue(close(mul(A, c.matrix("E")), [[0, 0]] * 3))
        self.assertEqual(len(c.matrix("D")), 2)

    def test_colcon(self):
        A = [[1, 2], [2, 4], [3, 6]]
        c = calc(A=A).run("COLCON")
        C = c.matrix("C")
        self.assertEqual(len(C), 2)
        self.assertTrue(close(mul(C, A), [[0, 0]] * 2))
        self.assertEqual(rank(C), 2)

    def test_subset_and_equal(self):
        A = [[1, 0], [1, 1], [0, 1]]
        B = [[1, 1], [2, 1], [1, 0]]  # same plane
        c = calc(A=A, B=B).run("SUBSP", [2])
        self.assertIn("=> SAME SUBSPACE", c.screen)
        c = calc(A=[[1], [0], [0]], B=B).run("SUBSP", [1])
        self.assertIn("SPAN(A) NOT IN SPAN(B)", c.screen)

    def test_intersection(self):
        A = [[1, 0], [0, 1], [0, 0]]  # xy-plane
        B = [[0, 1], [1, 0], [0, 1]]  # contains e2
        c = calc(A=A, B=B).run("SUBSP", [3])
        C = c.matrix("C")
        self.assertEqual(c.reals["R"], 1)
        self.assertAlmostEqual(C[0][0], 0)
        self.assertAlmostEqual(C[2][0], 0)

    def test_extend_and_simplify(self):
        A = [[1], [1], [0]]
        c = calc(A=A).run("SUBSP", [4])
        C = c.matrix("C")
        self.assertEqual(rank(C), 3)
        self.assertEqual([r[0] for r in C], [1, 1, 0])
        c = calc(A=[[1, 2, 3], [1, 2, 4], [1, 2, 5]]).run("SUBSP", [5])
        self.assertEqual(c.reals["R"], 2)

    def test_coordinates(self):
        S = [[1, 1], [1, -1]]
        c = calc(A=S, B=[[3], [1]]).run("COORD", [1])
        self.assertTrue(close(c.matrix("C"), [[2], [1]]))
        Tm = [[1, 0], [1, 1]]
        c = calc(A=S, B=Tm).run("COORD", [2])
        P = c.matrix("C")
        self.assertTrue(close(mul(Tm, P), S))


class TestOrthogonality(unittest.TestCase):
    def test_check(self):
        c = calc(A=[[1, 1], [1, -1]]).run("ORTHO", [1])
        self.assertIn("COLS ARE ORTHOGONAL", c.screen)
        self.assertIn("NOT ORTHONORMAL", c.screen)

    def test_complement(self):
        A = [[1, 0], [1, 1], [1, 2]]
        c = calc(A=A).run("ORTHO", [2])
        self.assertTrue(close(mul(T(A), c.matrix("C")), [[0], [0]]))

    def test_gram_schmidt_qr(self):
        A = [[1, 1, 0], [1, 0, 1], [0, 1, 1], [1, 1, 1]]
        c = calc(A=A).run("GRAMSCH")
        Q, R = c.matrix("D"), c.matrix("E")
        self.assertTrue(close(mul(T(Q), Q), eye(3)))
        self.assertTrue(close(mul(Q, R), A))
        V = c.matrix("C")
        self.assertEqual([r[0] for r in V], [1, 1, 0, 1])
        self.assertAlmostEqual(c.list("NRM")[0], 3)
        # exact form: r_ij = F(i,j) / sqrt(|v_i|^2)
        F, nrm = c.matrix("F"), c.list("NRM")
        exact = [[F[i][j] / nrm[i] ** 0.5 for j in range(3)] for i in range(3)]
        self.assertTrue(close(exact, R))

    def test_gram_schmidt_dependent(self):
        A = [[1, 2, 0], [1, 2, 1], [0, 0, 1]]
        c = calc(A=A).run("GRAMSCH")
        self.assertEqual(c.reals["K"], 2)

    def test_least_squares(self):
        A = [[1, 0], [1, 1], [1, 2]]
        b = [[6], [0], [0]]
        c = calc(A=A, B=b).run("LEASTSQ")
        self.assertTrue(close(c.matrix("C"), [[5], [-3]]))
        r = c.matrix("F")
        self.assertTrue(close(mul(T(A), r), [[0], [0]]))

    def test_least_squares_not_unique(self):
        A = [[1, 1], [1, 1], [0, 0]]
        c = calc(A=A, B=[[1], [3], [5]]).run("LEASTSQ")
        self.assertEqual(c.reals["F"], 1)
        self.assertTrue(close(c.matrix("E"), [[2], [2], [0]]))

    def test_polyfit(self):
        c = calc()
        c.lists["L1"] = [0.0, 1.0, 2.0, 3.0]
        c.lists["L2"] = [1.0, 3.0, 7.0, 13.0]  # 1 + x + x^2
        c.run("POLYFIT", [2])
        self.assertTrue(close(c.matrix("C"), [[1], [1], [1]]))


def check_diag(tc, A, c):
    P, D = c.matrix("C"), c.matrix("D")
    tc.assertTrue(close(mul(A, P), mul(P, D), 1e-6), (P, D))
    tc.assertEqual(rank(P), len(A))


class TestEigen(unittest.TestCase):
    def test_distinct(self):
        A = [[4, 1], [2, 3]]
        c = calc(A=A).run("EIGEN")
        self.assertEqual(c.list("CP"), [1, -7, 10])
        self.assertEqual(c.list("EV"), [5, 2])
        check_diag(self, A, c)
        self.assertTrue(close(mul(c.matrix("C"), c.matrix("E")), eye(2)))  # [E] = P^-1
        self.assertIn("P^-1 = [E]", c.screen)

    def test_repeated_diagonalizable(self):
        A = [[2, 0, 0], [0, 3, 4], [0, 4, 9]]  # eigenvalues 11, 2, 1
        c = calc(A=A).run("EIGEN")
        check_diag(self, A, c)
        A = [[1, 1, 1], [1, 1, 1], [1, 1, 1]]  # 3, 0, 0
        c = calc(A=A).run("EIGEN")
        self.assertEqual(c.list("EV"), [3, 0])
        self.assertEqual(c.list("AM"), [1, 2])
        self.assertEqual(c.list("GM"), [1, 2])
        check_diag(self, A, c)

    def test_triple_root(self):
        A = [[2, 0, 0, 0], [0, 2, 0, 0], [0, 0, 2, 0], [1, 2, 3, 5]]
        c = calc(A=A).run("EIGEN")
        self.assertEqual(c.list("EV"), [5, 2])
        self.assertEqual(c.list("AM"), [1, 3])
        check_diag(self, A, c)

    def test_defective(self):
        A = [[1, 1], [0, 1]]
        c = calc(A=A).run("EIGEN")
        self.assertEqual(c.list("GM"), [1])
        self.assertIn("NOT DIAGONALIZABLE", c.screen)

    def test_complex(self):
        c = calc(A=[[0, -1], [1, 0]]).run("EIGEN")
        self.assertEqual(c.reals["Z"], 1)
        self.assertEqual(sorted(e.imag for e in c.lists["EC"]), [-1, 1])

    def test_complex_lists_are_complex_typed(self):
        # The calculator gives ERR:DATA TYPE when a nonreal number is stored into an
        # element of a real list, so ZEIG must create ʟZR and ʟEC as complex lists.
        c = calc(A=[[4, 1], [2, 3]]).run("EIGEN")
        self.assertIn("ZR", c.clists)
        self.assertIn("EC", c.clists)

    def test_irrational(self):
        A = [[1, 2], [2, 2]]
        c = calc(A=A).run("EIGEN")
        check_diag(self, A, c)

    def test_fractions(self):
        A = [[0.5, 0.25], [0.5, 0.75]]  # 1 and 1/4
        c = calc(A=A).run("EIGEN")
        self.assertEqual(c.list("EV"), [1, 0.25])
        check_diag(self, A, c)

    def test_random_integer(self):
        rnd = random.Random(1522)
        for _ in range(25):
            n = rnd.choice([2, 3, 4])
            # A = S D S^-1 with integer eigenvalues (some repeated); det S = +-1 keeps A integer
            while True:
                S = [[rnd.randint(-2, 2) for _ in range(n)] for _ in range(n)]
                if abs(abs(Calc(SRC)._det(S)) - 1) < 1e-9:
                    break
            d = [rnd.choice([-2, -1, 0, 1, 3]) for _ in range(n)]
            D = [[d[i] if i == j else 0 for j in range(n)] for i in range(n)]
            A = mul(mul(S, D), Calc(SRC)._inverse(S))
            A = [[float(round(x)) for x in r] for r in A]
            c = calc(A=A).run("EIGEN")
            check_diag(self, A, c)


class TestOrthDiagSVD(unittest.TestCase):
    def test_orthdiag(self):
        for A in ([[2, 1, 1], [1, 2, 1], [1, 1, 2]], [[1, 2], [2, 1]], [[5, 2, 0], [2, 2, 0], [0, 0, 3]]):
            c = calc(A=A).run("ORTDIAG")
            P, D = c.matrix("C"), c.matrix("D")
            self.assertTrue(close(mul(T(P), P), eye(len(A)), 1e-6))
            self.assertTrue(close(mul(mul(P, D), T(P)), A, 1e-6))

    def test_not_symmetric(self):
        c = calc(A=[[1, 2], [3, 4]]).run("ORTDIAG")
        self.assertIn("A IS NOT SYMMETRIC", c.screen)

    def test_svd(self):
        for A in (
            [[3, 2, 2], [2, 3, -2]],
            [[1, 1], [0, 1], [1, 0]],
            [[1, 2], [2, 4]],
            [[1, 0, 1, 0], [0, 1, 0, 1], [1, 1, 1, 1]],
            [[2]],
        ):
            c = calc(A=A).run("SVD")
            U, S, V = c.matrix("C"), c.matrix("D"), c.matrix("E")
            self.assertTrue(close(mul(T(U), U), eye(len(A)), 1e-6), U)
            self.assertTrue(close(mul(T(V), V), eye(len(A[0])), 1e-6), V)
            self.assertTrue(close(mul(mul(U, S), T(V)), A, 1e-6), (U, S, V))
        c = calc(A=[[3, 2, 2], [2, 3, -2]]).run("SVD")
        self.assertTrue(all(abs(x - y) < 1e-9 for x, y in zip(c.list("SV"), [5, 3, 0])))
        self.assertIn("V VECTORS UNNORMALISED", c.screen)
        self.assertIn("A*V = SIGMA*U, SO", c.screen)

    def test_svd_extra_u(self):
        A = [[1, 1], [1, 1], [0, 0]]  # rank 1, so U needs two extra vectors
        c = calc(A=A).run("SVD")
        U, S, V = c.matrix("C"), c.matrix("D"), c.matrix("E")
        self.assertTrue(close(mul(mul(U, S), T(V)), A, 1e-6))
        self.assertIn("EXTRA U VECTORS:", c.screen)


class TestVectors(unittest.TestCase):
    def test_square(self):
        A = [[1, 1], [1, 1]]
        B = [[1, 1, 1], [1, -1, 0]]  # eigenvector (2), null vector (eigenvalue 0), neither
        c = calc(A=A, B=B).run("TESTVEC")
        s = "\n".join(c.screen)
        self.assertIn("EIGENVECTOR, LAMBDA=", s)
        self.assertIn("A*V=0: IN NULL(A)", s)
        self.assertIn("EIGENVECTOR, LAMBDA=0", s)
        self.assertIn("NOT AN EIGENVECTOR", s)
        self.assertEqual(c.screen.count("NOT IN COL(A)"), 2)
        self.assertEqual(c.screen.count("IN COL(A)"), 1)
        i = c.screen.index("EIGENVECTOR, LAMBDA=")
        self.assertEqual(c.screen[i + 1].strip(), "2")  # lambda of the first vector
        self.assertEqual(c.reals["B"], 1)  # rank A

    def test_rectangular(self):
        A = [[1, 2, 3], [2, 4, 6]]
        c = calc(A=A, B=[[1, 0], [2, 1]]).run("TESTVEC")  # vectors in R^2: column space test only
        self.assertEqual(c.screen.count("NOT IN COL(A)"), 1)
        self.assertEqual(c.screen.count("IN COL(A)"), 1)
        self.assertNotIn("EIGENVECTOR", "\n".join(c.screen))
        c = calc(A=A, B=[[1], [1], [-1]]).run("TESTVEC")  # vector in R^3: null space test only
        self.assertIn("A*V=0: IN NULL(A)", c.screen)
        self.assertNotIn("COL(A)", "\n".join(c.screen))

    def test_bad_size(self):
        c = calc(A=[[1, 2], [3, 4]], B=[[1], [2], [3]]).run("TESTVEC")
        self.assertIn("ROWS OF [B] MUST EQUAL", c.screen)


class TestMisc(unittest.TestCase):
    def test_markov(self):
        A = [[0.7, 0.2], [0.3, 0.8]]
        c = calc(A=A, B=[[1], [0]]).run("MARKOV", [2])
        self.assertIn("A IS STOCHASTIC", c.screen)
        self.assertTrue(close(c.matrix("C"), [[0.4], [0.6]]))
        self.assertTrue(close(c.matrix("D"), mul(mul(A, A), [[1], [0]])))

    def test_lintran(self):
        V = [[1, 0, 1], [2, -1, 0], [0, 3, 1]]
        W = [[4, 2, 3], [5, -1, 0], [1, 4, 2]]
        c = calc(A=V, B=W).run("LINTRAN")
        self.assertTrue(close(c.matrix("C"), [[2, 1, 1], [-0.6, 2.8, 0.6], [0.6, 0.2, 1.4]]))

    def test_menu(self):
        c = calc(A=[[1, 2], [3, 4]])
        c.run("MA1522", [1, 0, 0])
        self.assertTrue(close(c.matrix("C"), eye(2)))
        c = calc(A=[[4, 1], [2, 3]])
        c.run("MA1522", [6, 6, 5])
        self.assertEqual(c.list("EV"), [5, 2])
        c.run("MA1522", [7])
        self.assertIn("  -> QUIT", c.screen)
        c = calc(A=[[1, 1], [1, 1]], B=[[1], [1]])
        c.run("MA1522", [6, 6, 6, 4])
        self.assertIn("EIGENVECTOR, LAMBDA=", c.screen)


if __name__ == "__main__":
    unittest.main()
