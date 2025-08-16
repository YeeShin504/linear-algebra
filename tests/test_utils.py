import pytest

import sympy as sym

from ma1522.utils import _powerset, _is_zero


class TestUtils:
    def test_powerset(self):
        assert list(_powerset([1, 2])) == [(), (1,), (2,), (1, 2)]
        assert list(_powerset([])) == [()]

    def test_is_zero(self):
        a, b = sym.symbols("a b", real=True)
        assert _is_zero(-a * b / 2 - a + b) is True
        assert _is_zero(a - b) is True
        assert _is_zero(0) is True
        assert _is_zero(a * b) is True
        assert _is_zero(b) is True
        assert _is_zero(1) is False

        n, d = sym.fraction(1)
        assert _is_zero(d) is False
