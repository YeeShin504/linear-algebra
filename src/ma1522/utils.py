from __future__ import annotations
import dataclasses
from itertools import chain, combinations
from typing import TYPE_CHECKING

import sympy as sym

if TYPE_CHECKING:
    from typing import Iterable, Literal
    from sympy.printing.latex import LatexPrinter
    from ma1522.custom_types import Printable


def _powerset(args: Iterable) -> list:
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

    # def text(txt: str) -> str:
    #     return "\\text{" + txt + "}"

    # list_repr = []
    # for k, v in dataclasses.asdict(obj).items():
    #     k_repr = text(k)
    #     if hasattr(v, "_latex"):
    #         # used for overriding printing behaviour in sympy objects
    #         v_repr = v._latex(printer)
    #     elif hasattr(v, "_repr_latex_") and _unwrap_latex(v.__repr__()) != v.__repr__():
    #         # used for objects that support IPython printing in latex
    #         v_repr = _unwrap_latex(v.__repr__())
    #     else:
    #         v_repr = text(v.__repr__())
    #     list_repr.append(k_repr + " = " + v_repr)

    # merged = ", \\quad".join(list_repr)
    return _textify(type(obj).__name__) + _gen_latex_repr_dict(
        dataclasses.asdict(obj), printer=printer
    )


def _gen_latex_repr_dict(obj: dict, printer: LatexPrinter | None = None) -> str:
    """Generates a LaTeX representation of a dictionary.

    Args:
        obj (dict): The dictionary to represent.
        printer (LatexPrinter): The LaTeX printer to use.

    Returns:
        (str): A LaTeX string representation of the dictionary.
    """

    list_repr = []
    for k, v in obj.items():
        k_repr = _textify(k)
        if hasattr(v, "_latex"):
            # used for overriding printing behaviour in sympy objects
            v_repr = v._latex(printer)
        elif hasattr(v, "_repr_latex_") and _unwrap_latex(v.__repr__()) != v.__repr__():
            # used for objects that support IPython printing in latex
            v_repr = _unwrap_latex(v.__repr__())
        else:
            # either it does not have _repr_latex_ or its __repr__ has wrapped latex
            # representation, so we unwrap it
            # v_repr = _wrap_latex(v.__repr__())
            v_repr = _unwrap_latex(v.__repr__())
            # v_repr = textify(v.__repr__())
        list_repr.append(f"{k_repr} = {v_repr}")

    merged = ", \\quad".join(list_repr)
    return r"\left\{" + merged + r"\right\}"


def _textify(txt: str) -> str:
    """Converts a string to a LaTeX text representation."""
    txt = txt.replace("_", r"\_")  # escape underscores
    return r"\text{" + txt + "}"


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


def _is_IPython() -> bool:
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


def display(*args, opt: Literal["math", "dict"] | None = None, **kwargs) -> None:
    """Displays objects in a rich format, depending on the environment.

    This function displays objects using IPython's [`display`][IPython.display] mechanism if available,
    otherwise it falls back to [`sympy.pprint`][sympy.printing.pretty.pretty.PrettyPrinter].

    Args:
        *args: The objects to display.
        opt:

            - If "math", displays the object as a math expression.
            - If "dict", generates a LaTeX representation of the dictionary for display.
            - If none, assumes the object can be passed into IPython's [`display`][] function directly.
        **kwargs: Additional keyword arguments to pass to the display function.

    See Also:
        - [`IPython.display.display`][IPython.display]: The class used to display
            objects in IPython environments.
        - [`sympy.pprint`][sympy.printing.pretty.pretty.PrettyPrinter]: The class used to pretty-print
            objects in non-IPython environments.
    """
    if _is_IPython():
        import IPython.display

        if opt == "math":
            # Display the object as a math expression
            IPython.display.display(IPython.display.Math(*args, **kwargs))
        elif opt == "dict":
            # Generate a LaTeX representation of the dictionary
            from sympy.printing.latex import LatexPrinter

            printer = kwargs.pop("printer", LatexPrinter())
            for arg in args:
                if not isinstance(arg, dict):
                    continue
                IPython.display.display(
                    IPython.display.Math(_gen_latex_repr_dict(arg, printer=printer))
                )
        else:
            IPython.display.display(*args, **kwargs)
    else:
        sym.pprint(*args, **kwargs)
