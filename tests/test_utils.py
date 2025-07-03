import pytest
import sympy as sym
from utils import _powerset, _is_zero


class TestPowerSet:
    def test_powerset_empty(self) -> None:
        assert _powerset([]) == [()]

    def test_powerset_3(self) -> None:
        out = _powerset([1, 2, 3])
        expected = [(), (1,), (2,), (3,), (1, 2), (1, 3), (2, 3), (1, 2, 3)]
        assert out == expected


class TestIsZero:
    def test_is_zero(self) -> None:
        a, b = sym.symbols("a, b")
        assert _is_zero(-a * b / 2 - a + b)
        assert _is_zero(a - b)
        assert _is_zero(0)
        assert _is_zero(a * b)
        assert _is_zero(b)
        assert not _is_zero(1)

        n, d = sym.fraction(1)
        assert not _is_zero(d)
