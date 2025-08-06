import pytest

import sympy as sym

from ma1522.utils import _powerset, _is_zero


class TestUtils:
    def test_powerset(self):
        assert list(_powerset([1, 2])) == [(), (1,), (2,), (1, 2)]
        assert list(_powerset([])) == [()]

    def test_is_zero(self):
        x = sym.symbols("x")
        assert _is_zero(0) is True
        assert _is_zero(sym.Integer(0)) is True
        assert _is_zero(x - x) is True
        assert _is_zero(x) is False
        assert _is_zero(1) is False
