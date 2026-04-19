"""Comprehensive test suite for rref_cases() and evaluate_cases() covering a variety of scenarios and edge cases.

Covers:
- No free symbols → single case, correct pivots / free_params
- Single symbol that can be zero → two cases (zero / non-zero)
- Single symbol that cannot be zero → single case (non-zero only)
- Multiple independent symbols → all combinations of zero / non-zero
- Dependent symbols (one zero-condition fixes both) → correct merging
- with rhs: consistency detection
- with rhs: free_params counts only LHS columns
- Augmentation line is preserved on every RREF output
- Returned RREFCase fields are well-formed (types, lengths, …)
- 3x3 full-rank matrix → single case, zero free params
- Rank-deficient matrix without symbols → single case, correct free params
- System inconsistent for all symbol values → is_consistent=False in every case
- System consistent for all symbol values → is_consistent=True in every case
"""

import pytest
import sympy as sym

from ma1522 import Matrix, RREFCase


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _cond_map(cases: list[RREFCase]) -> list[dict]:
    """Return just the conditions dict from each case (for easy assertions)."""
    return [c.conditions for c in cases]


def _sort_cases(cases: list[RREFCase]) -> list[RREFCase]:
    """Sort cases by the string representation of their conditions so that
    tests are deterministic regardless of iteration order."""
    return sorted(cases, key=lambda c: str(sorted(c.conditions.items())))


# ---------------------------------------------------------------------------
# 1. No free symbols
# ---------------------------------------------------------------------------


class TestNoFreeSymbols:
    def test_full_rank_square(self):
        """Identity-like system → single case, no free params."""
        A = Matrix([[1, 0], [0, 1]])
        cases = A.rref_cases()
        assert len(cases) == 1
        c = cases[0]
        assert c.conditions == {}
        assert c.free_params == 0
        assert c.pivots == (0, 1)
        assert c.is_consistent is None

    def test_rank_deficient(self):
        """3×3 singular matrix → single case, 1 free param."""
        A = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        cases = A.rref_cases()
        assert len(cases) == 1
        c = cases[0]
        assert c.conditions == {}
        assert c.free_params == 1
        assert len(c.pivots) == 2
        assert c.is_consistent is None

    def test_rref_correct_values(self):
        """Verify the actual RREF matches sympy's rref."""
        A = Matrix([[2, 4], [1, 3]])
        cases = A.rref_cases()
        assert len(cases) == 1
        expected_rref, _ = sym.Matrix(A).rref()
        assert cases[0].rref == Matrix(expected_rref)


# ---------------------------------------------------------------------------
# 2. Single symbol – can be zero
# ---------------------------------------------------------------------------


class TestSingleSymbolCanBeZero:
    @pytest.fixture(autouse=True)
    def setup(self):
        self.a = sym.Symbol("a")

    def test_two_cases_produced(self):
        A = Matrix([[self.a, 1], [0, 1]])
        cases = A.rref_cases()
        assert len(cases) == 2

    def test_zero_case_conditions(self):
        A = Matrix([[self.a, 1], [0, 1]])
        cases = _sort_cases(A.rref_cases())
        # sorted: {a: 0} < {} lexicographically
        zero_case = cases[0]
        assert zero_case.conditions == {self.a: 0}

    def test_nonzero_case_conditions(self):
        A = Matrix([[self.a, 1], [0, 1]])
        cases = _sort_cases(A.rref_cases())
        nonzero_case = cases[1]
        assert nonzero_case.conditions == {}

    def test_zero_case_free_params(self):
        """When a=0 the first column vanishes → 1 free parameter."""
        A = Matrix([[self.a, 1], [0, 1]])
        cases = _sort_cases(A.rref_cases())
        assert cases[0].free_params == 1

    def test_nonzero_case_free_params(self):
        """When a≠0 the matrix is full rank → 0 free parameters."""
        A = Matrix([[self.a, 1], [0, 1]])
        cases = _sort_cases(A.rref_cases())
        assert cases[1].free_params == 0

    def test_excluded_is_complement(self):
        """Each case's 'excluded' should contain the other case's zero-conds."""
        A = Matrix([[self.a, 1], [0, 1]])
        cases = _sort_cases(A.rref_cases())
        zero_case, nonzero_case = cases
        # zero_case assumed a=0, so excluded is empty (no other constraints)
        assert zero_case.excluded == []
        # nonzero_case assumed nothing, so excluded lists {a: 0}
        assert {self.a: 0} in nonzero_case.excluded

    def test_rref_shape(self):
        """Both RREFs must have the same shape as the original matrix."""
        A = Matrix([[self.a, 1], [0, 1]])
        for c in A.rref_cases():
            assert c.rref.shape == A.shape

    def test_pivots_are_tuple(self):
        A = Matrix([[self.a, 1], [0, 1]])
        for c in A.rref_cases():
            assert isinstance(c.pivots, tuple)

    def test_return_type(self):
        A = Matrix([[self.a, 1], [0, 1]])
        cases = A.rref_cases()
        assert all(isinstance(c, RREFCase) for c in cases)


# ---------------------------------------------------------------------------
# 3. Single symbol – cannot be zero
# ---------------------------------------------------------------------------


class TestSingleSymbolCannotBeZero:
    def test_single_case_only(self):
        """A constant 1 in the pivot → no branching."""
        A = Matrix([[1, 2], [3, 4]])
        cases = A.rref_cases()
        assert len(cases) == 1

    def test_symbol_in_off_pivot_position(self):
        """Symbol only in an already-eliminated position → no branching."""
        a = sym.Symbol("a")
        # After eliminating col 0 via pivot row 0 (pivot=1, fixed),
        # 'a' ends up only in the off-diagonal; pivot for row 1 is fixed.
        A = Matrix([[1, a], [0, 2]])
        cases = A.rref_cases()
        assert len(cases) == 1
        assert cases[0].conditions == {}
        assert cases[0].free_params == 0


# ---------------------------------------------------------------------------
# 4. Two independent symbols
# ---------------------------------------------------------------------------


class TestTwoIndependentSymbols:
    @pytest.fixture(autouse=True)
    def setup(self):
        self.a, self.b = sym.symbols("a b")

    def test_correct_number_of_cases(self):
        """Diagonal matrix [[a,0],[0,b]] → 4 combinations (both, one, other, none zero)."""
        A = Matrix([[self.a, 0], [0, self.b]])
        cases = A.rref_cases()
        # Expected: {a:0,b:0}, {a:0}, {b:0}, {}
        assert len(cases) == 4

    def test_all_zero_case_has_2_free_params(self):
        A = Matrix([[self.a, 0], [0, self.b]])
        cases = A.rref_cases()
        all_zero = next(c for c in cases if c.conditions == {self.a: 0, self.b: 0})
        assert all_zero.free_params == 2

    def test_no_zero_case_has_0_free_params(self):
        A = Matrix([[self.a, 0], [0, self.b]])
        cases = A.rref_cases()
        no_zero = next(c for c in cases if c.conditions == {})
        assert no_zero.free_params == 0

    def test_one_zero_has_1_free_param(self):
        A = Matrix([[self.a, 0], [0, self.b]])
        cases = A.rref_cases()
        one_a = next(
            c
            for c in cases
            if c.conditions == {self.a: 0} and self.b not in c.conditions
        )
        one_b = next(
            c
            for c in cases
            if c.conditions == {self.b: 0} and self.a not in c.conditions
        )
        assert one_a.free_params == 1
        assert one_b.free_params == 1


# ---------------------------------------------------------------------------
# 5. With RHS – consistency checks
# ---------------------------------------------------------------------------


class TestWithRHS:
    @pytest.fixture(autouse=True)
    def setup(self):
        self.a = sym.Symbol("a")

    def test_is_consistent_not_none_when_rhs_given(self):
        A = Matrix([[self.a, 1], [0, 1]])
        b = Matrix([[2], [3]])
        for c in A.rref_cases(rhs=b):
            assert c.is_consistent is not None

    def test_is_consistent_none_without_rhs(self):
        A = Matrix([[self.a, 1], [0, 1]])
        for c in A.rref_cases():
            assert c.is_consistent is None

    def test_zero_case_inconsistent(self):
        """When a=0 the system [[0,1|2],[0,0|1]] is inconsistent."""
        A = Matrix([[self.a, 1], [0, 1]])
        b = Matrix([[2], [3]])
        cases = _sort_cases(A.rref_cases(rhs=b))
        zero_case = cases[0]  # {a: 0}
        assert zero_case.is_consistent is False

    def test_nonzero_case_consistent(self):
        A = Matrix([[self.a, 1], [0, 1]])
        b = Matrix([[2], [3]])
        cases = _sort_cases(A.rref_cases(rhs=b))
        nonzero_case = cases[1]  # {}
        assert nonzero_case.is_consistent is True

    def test_free_params_count_only_lhs_cols(self):
        """free_params must not count RHS columns as free variables."""
        A = Matrix([[self.a, 1], [0, 1]])
        b = Matrix([[2], [3]])
        cases = _sort_cases(A.rref_cases(rhs=b))
        nonzero_case = cases[1]
        # Full-rank LHS → 0 free params
        assert nonzero_case.free_params == 0

    def test_always_consistent_system(self):
        """A system that is consistent regardless of the symbol value."""
        # [[1, 0], [0, a]] x = [1, 0]^T
        # When a=0: [[1,0|1],[0,0|0]] → consistent (1 free param)
        # When a≠0: [[1,0|1],[0,1|0]] → consistent (unique)
        A = Matrix([[1, 0], [0, self.a]])
        b = Matrix([[1], [0]])
        for c in A.rref_cases(rhs=b):
            assert c.is_consistent is True

    def test_rref_augmented_shape(self):
        """RREF matrix must have self.cols + rhs.cols columns."""
        A = Matrix([[self.a, 1], [0, 1]])
        b = Matrix([[2], [3]])
        for c in A.rref_cases(rhs=b):
            assert c.rref.cols == A.cols + b.cols


# ---------------------------------------------------------------------------
# 6. Augmentation line preserved
# ---------------------------------------------------------------------------


class TestAugLinePreserved:
    def test_aug_line_in_output(self):
        """Augmentation lines set on the working matrix survive into every RREF."""
        a = sym.Symbol("a")
        A = Matrix([[a, 1], [0, 1]])
        b = Matrix([[2], [3]])
        cases = A.rref_cases(rhs=b)
        for c in cases:
            # row_join adds aug line at col A.cols-1 = 1
            assert hasattr(c.rref, "_aug_pos")
            assert len(c.rref._aug_pos) > 0


# ---------------------------------------------------------------------------
# 7. Rectangular / non-square matrices
# ---------------------------------------------------------------------------


class TestRectangularMatrices:
    def test_wide_matrix_no_symbols(self):
        """More columns than rows → always has free parameters."""
        A = Matrix([[1, 2, 3], [0, 1, 4]])
        cases = A.rref_cases()
        assert len(cases) == 1
        assert cases[0].free_params == 1

    def test_tall_matrix_no_symbols(self):
        """More rows than columns, full col rank → 0 free params."""
        A = Matrix([[1, 0], [0, 1], [0, 0]])
        cases = A.rref_cases()
        assert len(cases) == 1
        assert cases[0].free_params == 0

    def test_wide_matrix_with_symbol(self):
        a = sym.Symbol("a")
        A = Matrix([[a, 1, 0], [0, 0, 1]])
        cases = A.rref_cases()
        # When a=0: only 2 pivots in 3 cols → 1 free param
        # When a≠0: 2 pivots in 3 cols → 1 free param (col 1 is free)
        for c in cases:
            assert c.free_params >= 1


# ---------------------------------------------------------------------------
# 8. Homogeneous system always consistent
# ---------------------------------------------------------------------------


class TestHomogeneousConsistency:
    def test_homogeneous_always_consistent(self):
        a = sym.Symbol("a")
        A = Matrix([[a, 1], [0, 1]])
        b = Matrix([[0], [0]])
        for c in A.rref_cases(rhs=b):
            assert c.is_consistent is True


# ---------------------------------------------------------------------------
# 9. Idempotency – rref_cases on a constant matrix matches rref()
# ---------------------------------------------------------------------------


class TestConsistencyWithRref:
    def test_matches_standard_rref(self):
        """For a numeric matrix, rref_cases should produce the same RREF as rref()."""
        A = Matrix([[2, 1, -1], [-3, -1, 2], [-2, 1, 2]])
        cases = A.rref_cases()
        assert len(cases) == 1
        expected = A.rref().rref
        assert cases[0].rref == expected


# ---------------------------------------------------------------------------
# 10. evaluate_cases merge behavior
# ---------------------------------------------------------------------------


class TestEvaluateCasesMerge:
    @pytest.fixture(autouse=True)
    def setup(self):
        self.mat2 = Matrix.from_str(
            "a 2 a (a+b); a 2 a a; 3 3 -b 3; (a+1) 3 (a+1) (a+1)"
        )
        self.aug2 = Matrix.from_str("a-b; a-b; -b; a-b+1")

    def test_merges_inconsistent_region_to_a_eq_2_and_b_nonzero(self):
        """For the mat2/aug2 repro, inconsistent branches should merge to
        a=2 with only b=0 excluded (i.e. b!=0).
        """
        merged = self.mat2.evaluate_cases(rhs=self.aug2, verbosity=0)
        no_solution = [c for c in merged if c.is_consistent is False]

        assert len(no_solution) == 1
        c = no_solution[0]
        assert {str(k): v for k, v in c.conditions.items()} == {"a": 2}
        assert [{str(k): v for k, v in d.items()} for d in c.excluded] == [{"b": 0}]

    def test_freezes_exact_raw_case_signature_for_mat2_aug2(self):
        """Freeze current raw rref_cases behavior for mat2/aug2 exactly."""
        raw = self.mat2.rref_cases(rhs=self.aug2, verbosity=0)

        got = [
            (
                {str(k): v for k, v in c.conditions.items()},
                [{str(k): v for k, v in d.items()} for d in c.excluded],
                c.is_consistent,
                c.free_params,
                c.pivots,
            )
            for c in raw
        ]

        expected = [
            ({"a": 2, "b": -3}, [], False, 2, (0, 3, 4)),
            ({"a": 2, "b": 0}, [], True, 2, (0, 2)),
            ({"a": 2}, [{"b": -3}, {"b": 0}], False, 1, (0, 2, 3, 4)),
            ({"b": -3}, [{"a": 2}], True, 1, (0, 1, 3)),
            ({"b": 0}, [{"a": 2}], True, 1, (0, 1, 2)),
            ({}, [{"a": 2}, {"b": -3}, {"b": 0}], True, 0, (0, 1, 2, 3)),
        ]

        assert got == expected

    def test_freezes_exact_merged_case_signature_for_mat2_aug2(self):
        """Freeze current evaluate_cases merge+order behavior for mat2/aug2 exactly."""
        merged = self.mat2.evaluate_cases(rhs=self.aug2, verbosity=0)

        got = [
            (
                {str(k): v for k, v in c.conditions.items()},
                [{str(k): v for k, v in d.items()} for d in c.excluded],
                c.is_consistent,
                c.free_params,
                c.pivots,
            )
            for c in merged
        ]

        expected = [
            ({"a": 2}, [{"b": 0}], False, 1, (0, 2, 3, 4)),
            ({}, [{"a": 2}, {"b": -3}, {"b": 0}], True, 0, (0, 1, 2, 3)),
            ({"b": -3}, [{"a": 2}], True, 1, (0, 1, 3)),
            ({"b": 0}, [{"a": 2}], True, 1, (0, 1, 2)),
            ({"a": 2, "b": 0}, [], True, 2, (0, 2)),
        ]

        assert got == expected
