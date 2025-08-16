import pytest

import sympy as sym

from ma1522 import Matrix


class TestTutorial01:
    def test_question_3a(self):
        mat = Matrix([[3, 2, -4], [2, 3, 3], [5, -3, 1]])

        aug = Matrix([3, 15, 14])

        aug_mat = mat.row_join(aug)
        plu = aug_mat.ref()
        assert plu.U.is_echelon

        rref = plu.U.rref()
        assert rref.rref[:, -1] == Matrix([3, 1, 2])  # type: ignore

    def test_question_3b(self):
        mat = Matrix([[1, 1, -1, -2], [2, 1, -1, 1], [-1, 1, -3, 1]])
        aug = Matrix([0, -2, 4])
        sol = mat.solve(aug)[0]
        unk = sol[-1]
        assert sol == Matrix([-3 * unk - 2, 19 * unk / 2 + 2, 9 * unk / 2, unk])  # type: ignore

    def test_question_3c(self):
        aug_mat = Matrix([[1, -4, 2, -2], [1, 2, -2, -3], [1, -1, 0, 4]], aug_pos=2)
        mat = aug_mat.select_cols(0, 1, 2)  # Remove the augmented column
        aug = aug_mat.select_cols(3)  # Select the augmented column

        with pytest.raises(ValueError):
            mat.solve(aug)[0]

    def test_question_4(self):
        # TODO: Implement better evaluate all cases
        a, b = sym.symbols("a b")  # Define symbolic variables
        mat = Matrix([[a, 0, b, 2], [a, a, 4, 4], [0, a, 2, b]], aug_pos=2)
        plu = mat.ref()
        assert plu.U.is_echelon
        # Case 1: b != 2, a = 0
        assert plu.U.subs({a: 0}).rref(pivots=False) == Matrix(
            [[0, 0, 1, 0], [0, 0, 0, 1], [0, 0, 0, 0]], aug_pos=2
        )
        # Case 2: b != 2, a != 0
        assert plu.U.rref(pivots=False) == Matrix(
            [[1, 0, 0, (2 - b) / a], [0, 1, 0, (b - 2) / a], [0, 0, 1, 1]], aug_pos=2
        )
        # Case 3: b = 2, a != 0
        assert plu.U.subs({b: 2}).rref(pivots=False) == Matrix(
            [[1, 0, 2 / a, 2 / a], [0, 1, 2 / a, 2 / a], [0, 0, 0, 0]], aug_pos=2
        )
        # Case 4: b = 2, a = 0
        assert plu.U.subs({a: 0, b: 2}).rref(pivots=False) == Matrix(
            [[0, 0, 1, 1], [0, 0, 0, 0], [0, 0, 0, 0]], aug_pos=2
        )

    def test_question_6(self):
        x, y, z = sym.symbols("x y z")
        vec = Matrix([x**2, y**2, z**2])
        mat = Matrix([[1, -1, 2], [2, 2, -5], [2, 5, 1]])
        aug = Matrix([6, 3, 9])
        assert len(sym.solve(mat @ vec - aug)) == 4

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

        mat = aug_mat.select_cols(*range(7))
        aug = aug_mat.select_cols(7)

        sol = mat.solve(aug)[0]
        y = sol[5]
        z = sol[6]
        assert sol.subs({y: 50, z: 100}) == Matrix([100, 550, 700, 650, 50, 50, 100])  # type: ignore
        assert sol.subs({y: -50}) == Matrix([0, z + 450, 800, z + 650, z - 50, -50, z])  # type: ignore
