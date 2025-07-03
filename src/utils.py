from itertools import chain, combinations
import sympy as sym
from sympy.printing.latex import LatexPrinter
from typing import Iterable, List, NamedTuple

PRINTER = LatexPrinter()


def _powerset(args: Iterable) -> List:
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(args)
    return list(chain.from_iterable(combinations(s, r) for r in range(len(s) + 1)))


def _is_zero(expr) -> bool:
    """
    returns True if expr can be 0 for real values of unknowns
    use is_complex rather than is_real as symbol.is_complex and
    symbol.is_real returns None
    """

    if (not isinstance(expr, sym.Expr)) or isinstance(expr, sym.Number):
        return expr == 0

    # set symbols assumption to true
    real_symbols = sym.symbols(f"x:{len(expr.free_symbols)}")
    for symbol, real_symbol in zip(expr.free_symbols, real_symbols):
        expr = expr.subs({symbol: real_symbol})

    sol = sym.solve(sym.Eq(expr, 0), expr.free_symbols)
    return len(sol) > 0


def _gen_latex_repr(obj: NamedTuple, printer=None) -> str:
    def text(txt: str) -> str:
        return "\\text{" + txt + "}"

    list_repr = []
    for k, v in obj._asdict().items():
        k_repr = text(k)
        if hasattr(v, "_latex"):
            v_repr = v._latex(printer)
        elif hasattr(v, "_repr_latex_") and _unwrap_latex(v.__repr__()) != v.__repr__():
            v_repr = _unwrap_latex(v.__repr__())
        else:
            v_repr = text(v.__repr__())
        list_repr.append(k_repr + " = " + v_repr)

    merged = ", \\quad".join(list_repr)
    expr = text(type(obj).__name__) + " = \\left(" + merged + "\\right)"
    return expr


def _wrap_latex(expr: str | None) -> str:
    return f"${expr}$"


def _unwrap_latex(expr: str | None) -> str:
    if expr is None:
        return ""
    return expr.replace("$", "").rstrip()
