import pytest

from ma1522 import Matrix


class TestTutorial04:
    def test_question_2a(self):
        U = Matrix.from_str("2 1 0 3; 3 -1 5 2; -1 0 2 1").T
        actual_constraints = U.column_constraints()
        vec = Matrix.create_unk_matrix(4, 1, "x")
        # Explicitly assert the constraint from notebook
        expected = vec[0, 0] + 7 * vec[1, 0] + 2 * vec[2, 0] - 3 * vec[3, 0]
        assert actual_constraints[3, 3] == expected

        # Check vectors from (a)
        rhs = Matrix.from_str("2 3 -7 3; 0 0 0 0; 1 1 1 1; -4 6 -13 4").T
        # (i) 2 + 7(3) + 2(-7) - 3(3) = 2 + 21 - 14 - 9 = 0 -> Consistent
        assert rhs.select_cols(0).is_subspace_of(U)
        # (ii) 0 -> Consistent
        assert rhs.select_cols(1).is_subspace_of(U)
        # (iii) 1 + 7 + 2 - 3 = 7 != 0 -> Inconsistent
        assert not rhs.select_cols(2).is_subspace_of(U)
        # (iv) -4 + 7(6) + 2(-13) - 3(4) = -4 + 42 - 26 - 12 = 0 -> Consistent
        assert rhs.select_cols(3).is_subspace_of(U)

    def test_question_2b(self):
        U = Matrix.from_str("2 1 0 3; 3 -1 5 2; -1 0 2 1").T
        E = U.extend_basis()
        assert E.shape == (4, 4)
        assert E.rank() == 4

    def test_question_3a(self):
        V = Matrix([[1, -1, -1]]).nullspace()
        S = Matrix.from_str("1 1 0; 5 2 3").T
        assert S.is_same_subspace(Matrix.from_list(V))

    def test_question_3b(self):
        S = Matrix.from_str("1 1 0; 5 2 3").T
        T = S.row_join(Matrix.from_str("0 0 1").T)
        T.rm_aug_line()
        assert T.is_same_subspace()

    def test_question_4(self):
        S_i = Matrix.from_str("1 0 0 1; 0 1 0 0; 1 1 1 1; 1 1 1 0").T
        assert S_i.is_same_subspace()

        S_ii = Matrix.from_str("1 2 1 0; 1 1 -1 0; 0 0 0 1").T
        assert not S_ii.is_same_subspace()

        S_iii = Matrix.from_str("6 4 -2 4; 2 0 0 1; 3 2 -1 2; 5 6 -3 2; 0 4 -2 -1").T
        assert not S_iii.is_same_subspace()

        S_iv = Matrix.from_str("1 1 0 0; 1 2 -1 1; 0 0 1 1; 2 1 2 1; 1 2 3 4").T
        assert S_iv.is_same_subspace()

    def test_question_5a(self):
        U = Matrix.from_str("2 -2 0; -1 1 -1; 0 0 9").T
        V = Matrix.from_str("1 -1 -5; 0 1 1").T
        assert not U.is_same_subspace(V)

    def test_question_5b(self):
        U = Matrix.from_str("1 6 4; 2 4 -1; -1 2 5").T
        V = Matrix.from_str("1 -2 -5; 0 8 9").T
        assert U.is_same_subspace(V)
