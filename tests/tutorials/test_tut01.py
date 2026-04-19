import pytest

import sympy as sym

from ma1522 import Matrix


class TestTutorial01:
    def test_question_3a(self):
        mat = Matrix([[3, 2, -4], [2, 3, 3], [5, -3, 1]])

        aug = Matrix([3, 15, 14])

        aug_mat = mat.row_join(aug)

        # Intermediate: REF should produce an upper-triangular echelon form
        plu = aug_mat.ref()
        assert plu.U.is_echelon

        # Intermediate: full RREF should be the identity augmented with the solution
        assert plu.U.rref(pivots=False) == Matrix(  # type: ignore
            [[1, 0, 0, 3], [0, 1, 0, 1], [0, 0, 1, 2]], aug_pos=2
        )
        # Solution column only (using the RREF namedtuple's .rref attribute)
        rref = plu.U.rref()
        assert rref.rref[:, -1] == Matrix([3, 1, 2])  # type: ignore

    def test_question_3b(self):
        mat = Matrix([[1, 1, -1, -2], [2, 1, -1, 1], [-1, 1, -3, 1]])
        aug = Matrix([0, -2, 4])
        aug_mat = mat.row_join(aug)

        # Intermediate: RREF shows the parametric structure — x4 is the free variable
        rref = aug_mat.rref(pivots=False)
        assert rref == Matrix(
            [
                [1, 0, 0, 3, -2],
                [0, 1, 0, sym.Rational(-19, 2), 2],
                [0, 0, 1, sym.Rational(-9, 2), 0],
            ],
            aug_pos=3,
        )

        # Solution expressed in terms of the free parameter
        sol = mat.solve(aug)[0]
        unk = sol[-1]
        assert sol == Matrix([-3 * unk - 2, 19 * unk / 2 + 2, 9 * unk / 2, unk])  # type: ignore

    def test_question_3c(self):
        aug_mat = Matrix([[1, -4, 2, -2], [1, 2, -2, -3], [1, -1, 0, 4]], aug_pos=2)

        # Intermediate: RREF has a pivot in the augmented column → inconsistent system
        rref = aug_mat.rref(pivots=False)
        # Last row must be all-zero on the coefficient side with a non-zero RHS entry
        assert rref[-1, :-1] == Matrix([[0, 0, 0]])
        assert rref[-1, -1] != 0

        mat = aug_mat.select_cols(0, 1, 2)
        aug = aug_mat.select_cols(3)

        with pytest.raises(ValueError):
            mat.solve(aug)[0]

    def test_question_4(self):
        a, b = sym.symbols("a b")
        mat = Matrix([[a, 0, b, 2], [a, a, 4, 4], [0, a, 2, b]], aug_pos=2)

        # Explicitly test evaluate_cases as shown in notebook
        A = mat.select_cols(0, 1, 2)
        rhs = mat.select_cols(3)
        cases = A.evaluate_cases(rhs=rhs, verbosity=0)

        assert len(cases) == 4
        # Case 1: b != 2, a = 0
        assert cases[0].conditions == {a: 0}
        assert cases[0].is_consistent is False
        # Case 2: b != 2, a != 0
        assert cases[1].conditions == {}
        assert cases[1].is_consistent is True
        # Case 3: b = 2, a != 0
        assert cases[2].conditions == {b: 2}
        assert cases[2].is_consistent is True
        # Case 4: b = 2, a = 0
        assert cases[3].conditions == {a: 0, b: 2}
        assert cases[3].is_consistent is True

        # Intermediate: REF always produces an upper-triangular form
        plu = mat.ref()
        assert plu.U.is_echelon

        # Case 1: b != 2, a = 0 → no solution (pivot in augmented column)
        assert plu.U.subs({a: 0}).rref(pivots=False) == Matrix(
            [[0, 0, 1, 0], [0, 0, 0, 1], [0, 0, 0, 0]], aug_pos=2
        )
        # Case 2: b != 2, a != 0 → unique solution
        assert plu.U.rref(pivots=False) == Matrix(
            [[1, 0, 0, (2 - b) / a], [0, 1, 0, (b - 2) / a], [0, 0, 1, 1]], aug_pos=2
        )
        # Case 3: b = 2, a != 0 → infinitely many solutions, 1 free variable (x3)
        assert plu.U.subs({b: 2}).rref(pivots=False) == Matrix(
            [[1, 0, 2 / a, 2 / a], [0, 1, 2 / a, 2 / a], [0, 0, 0, 0]], aug_pos=2
        )
        # Case 4: b = 2, a = 0 → infinitely many solutions, 2 free variables (x1, x2)
        assert plu.U.subs({a: 0, b: 2}).rref(pivots=False) == Matrix(
            [[0, 0, 1, 1], [0, 0, 0, 0], [0, 0, 0, 0]], aug_pos=2
        )

    def test_question_6(self):
        x, y, z = sym.symbols("x y z")

        # Intermediate: substituting u=x², v=y², w=z² linearises the system;
        # RREF gives a unique solution u=4, v=0, w=1
        mat = Matrix([[1, -1, 2, 6], [2, 2, -5, 3], [2, 5, 1, 9]], aug_pos=2)
        assert mat.rref(pivots=False) == Matrix(
            [[1, 0, 0, 4], [0, 1, 0, 0], [0, 0, 1, 1]], aug_pos=2
        )

        # Four real (x,y,z) solutions come from x=±2, y=0, z=±1
        vec = Matrix([x**2, y**2, z**2])
        coef_mat = Matrix([[1, -1, 2], [2, 2, -5], [2, 5, 1]])
        aug = Matrix([6, 3, 9])
        sols = sym.solve(coef_mat @ vec - aug)
        assert len(sols) == 4

        # Each solution must satisfy x²=4, y²=0, z²=1
        for sol in sols:
            assert sol[x] ** 2 == 4
            assert sol[y] ** 2 == 0
            assert sol[z] ** 2 == 1

    def test_question_7(self):
        aug_mat = Matrix(
            [
                [1, 0, 1, 0, 0, 0, 0, 800],
                [1, -1, 0, 1, 0, 0, 0, 200],
                [0, 1, 0, 0, -1, 0, 0, 500],
                [0, 0, 1, 0, 0, 1, 0, 750],
                [0, 0, 0, -1, 0, -1, 1, -600],
                [0, 0, 0, 0, 1, 0, -1, -50],
            ],
            aug_pos=6,
        )

        # Explicitly assert the matrix from notebook
        assert aug_mat.shape == (6, 8)

        # Intermediate (7a): RREF has 5 pivots among 7 coefficient columns
        # → 2 free variables, so the system is underdetermined
        rref = aug_mat.rref(pivots=False)
        assert rref.shape == (6, 8)

        # Assert specific RREF values from notebook
        assert rref[0, -1] == 50
        assert rref[1, -1] == 450
        assert rref[5, -1] == 0

        # 5 coefficient columns become identity columns; x6 and x7 are free
        pivot_cols = [
            c
            for c in range(7)
            if any(
                rref[r, c] == 1 and all(rref[rr, c] == 0 for rr in range(6) if rr != r)
                for r in range(6)
            )
        ]
        assert len(pivot_cols) == 5

        mat = aug_mat.select_cols(*range(7))
        aug = aug_mat.select_cols(7)

        sol = mat.solve(aug)[0]

        # Intermediate: the general solution has exactly 2 free symbols (x6 and x7)
        assert len(sol.free_symbols) == 2

        y = sol[5]
        z = sol[6]

        # (7b) Substitute x6=50, x7=100
        assert sol.subs({y: 50, z: 100}) == Matrix([100, 550, 700, 650, 50, 50, 100])  # type: ignore

        # (7c) Set x6=-50 so that x1=0; x7 remains free
        assert sol.subs({y: -50}) == Matrix([0, z + 450, 800, z + 650, z - 50, -50, z])  # type: ignore

        # Coverage: Run once with verbosity to hit print lines
        mat.solve(aug, verbosity=1)
