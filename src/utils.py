from __future__ import annotations
import dataclasses
from itertools import chain, combinations
from typing import Iterable, List, Literal, Optional, TYPE_CHECKING

import sympy as sym

if TYPE_CHECKING:
    from custom_types import Printable
    from sympy.printing.latex import LatexPrinter


def _powerset(args: Iterable) -> List:
    """Generates the powerset of an iterable.

    Args:
        args: The iterable to generate the powerset from.

    Returns:
        A list of tuples representing the powerset.

    Examples:
        >>> _powerset([1, 2, 3])
        [(), (1,), (2,), (3,), (1, 2), (1, 3), (2, 3), (1, 2, 3)]
    """
    s = list(args)
    return list(chain.from_iterable(combinations(s, r) for r in range(len(s) + 1)))


def _is_zero(expr) -> bool:
    """Checks if a symbolic expression can be zero.

    This function checks if a given symbolic expression can be equal to zero for
    some real values of its variables.

    Args:
        expr: The symbolic expression to check.

    Returns:
        True if the expression can be zero, False otherwise.
    """

    if (not isinstance(expr, sym.Expr)) or isinstance(expr, sym.Number):
        return expr == 0

    # set symbols assumption to true
    real_symbols = sym.symbols(f"x:{len(expr.free_symbols)}")
    for symbol, real_symbol in zip(expr.free_symbols, real_symbols):
        expr = expr.subs({symbol: real_symbol})

    sol = sym.solve(sym.Eq(expr, 0), expr.free_symbols)
    return len(sol) > 0


def _gen_latex_repr(obj: Printable, printer: LatexPrinter | None = None) -> str:
    """Generates a LaTeX representation of a printable object.

    Args:
        obj: The object to represent.
        printer: The LaTeX printer to use.

    Returns:
        A LaTeX string representation of the object.
    """
    def text(txt: str) -> str:
        return "\\text{" + txt + "}"

    list_repr = []
    for k, v in dataclasses.asdict(obj).items():
        k_repr = text(k)
        if hasattr(v, "_latex"):
            # used for overriding printing behaviour in sympy objects
            v_repr = v._latex(printer)
        elif hasattr(v, "_repr_latex_") and _unwrap_latex(v.__repr__()) != v.__repr__():
            # used for objects that support IPython printing in latex
            v_repr = _unwrap_latex(v.__repr__())
        else:
            v_repr = text(v.__repr__())
        list_repr.append(k_repr + " = " + v_repr)

    merged = ", \\quad".join(list_repr)
    expr = text(type(obj).__name__) + " = \\left(" + merged + "\\right)"
    return expr


def _wrap_latex(expr: str | None) -> str:
    """Wraps a string in LaTeX math delimiters.

    Args:
        expr: The string to wrap.

    Returns:
        The wrapped string.
    """
    return f"${expr}$"


def _unwrap_latex(expr: str | None) -> str:
    """Unwraps a string from LaTeX math delimiters.

    Args:
        expr: The string to unwrap.

    Returns:
        The unwrapped string.
    """
    if expr is None:
        return ""
    # return expr.replace("$", "").rstrip()
    return (
        expr.strip()
        .removeprefix("$")
        .removeprefix("$")  # repeated for $$
        .removesuffix("$")
        .removesuffix("$")  # repeated for $$
        .strip()
    )


def is_IPython() -> bool:
    """Checks if the code is running in an IPython environment.
    Used to determine the printing options for the objects.

    Returns:
        True if running in IPython, False otherwise.
    """
    # Adapted from https://stackoverflow.com/questions/15411967/how-can-i-check-if-code-is-executed-in-the-ipython-notebook
    try:
        from IPython.core.getipython import get_ipython

        shell = get_ipython().__class__.__name__
        if shell == "ZMQInteractiveShell" or shell == "TerminalInteractiveShell":
            return True  # Jupyter notebook, qtconsole or terminal running IPython
        else:
            return False  # Other type
    except NameError:
        return False  # Probably standard Python interpreter
    except ImportError:
        return False  # IPython module does not exist


def display(*args, opt: Optional[Literal["math"]] = None, **kwargs) -> None:
    """Displays objects in a rich format, depending on the environment.

    This function displays objects using IPython's display mechanism if available,
    otherwise it falls back to `sympy.pprint`.

    Args:
        *args: The objects to display.
        opt: If "math", displays the object as a math expression.
        **kwargs: Additional keyword arguments to pass to the display function.
    """
    if is_IPython():
        import IPython.display

        if opt == "math":
            IPython.display.display(IPython.display.Math(*args, **kwargs))
        else:
            IPython.display.display(*args, **kwargs)
    else:
        sym.pprint(*args, **kwargs)