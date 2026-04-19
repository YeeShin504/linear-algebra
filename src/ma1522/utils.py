from __future__ import annotations
import dataclasses
from itertools import chain, combinations
import re
from typing import TYPE_CHECKING

import sympy as sym

if TYPE_CHECKING:
    from typing import Iterable, Literal
    from sympy.printing.latex import LatexPrinter
    from ma1522.custom_types import Printable


def sympy_commands():
    commands = """
    # Note: zero-indexing, see https://docs.sympy.org/latest/modules/matrices/matrices.html
    import sympy as sym

    # Variables
    >>> a, b = sym.symbols("a b") # symbols
    >>> I_3 = eye(3) # identity matrix
    >>> zeros = sym.zeros(2, cols=2) # zero matrix
    >>> A = Matrix([[...]] ) # user-defined matrix
    >>> B = Matrix.vstack(A, I_3, ...) # repeated application of col_join
    >>> C = sym.nsimplify(A, tolerance=0.001, rational=True) # convert to fraction
    >>> D = C.evalf() # convert to decimal

    # Matrix Operations
    >>> A @ B # matrix multiplication
    >>> A + B # element wise addition
    >>> A - B # element wise subtraction
    >>> A.col_del(col) # delete column
    >>> A.col_insert(pos, B) # insert matrix B at pos in A, column wise
    >>> A.col_join(B) # insert matrix B below matrix A
    >>> A.dot(B) # dot product of A and B
    >>> A.exp() # exponential
    >>> A.flat() # flatten to row vector
    >>> A.pow(exp) # power
    >>> A.reshape(rows, cols) # reshape
    >>> A.rot90(k=1) # rotate 90deg, k times
    >>> A.row_del(row) # delete row
    >>> A.row_insert(pos, B) # insert matrix B at pos in A, row wise
    >>> A.row_join(B) # insert matrix B on RHS of matrix A
    >>> A.vec() # stack to column vector
    >>> A.xreplace(dict) # replace sym (key) with value

    # Symbolic Methods
    >>> A.T # transpose
    >>> A.inv() # inverse
    >>> A.adj() # adjoint (as per MA1522 definition)
    >>> A.cofactor(i, j) # (i, j) cofactor
    >>> A.cofactor_matrix() # cofactor matrix
    >>> A.columnspace(simplify=False) # list of vectors that span column space of A
    >>> A.conjugate() # conjugate
    >>> A.copy() # copy matrix
    >>> A.det(method='bareiss', iszerofunc=None) # determinant, use domain-ge/laplace as method
    >>> A.diag() # diagonal
    >>> A.echelon_form() # REF
    >>> A.eigenvals(rational=True) # eigenvalues
    >>> A.eigenvects() # eigenvectors
    >>> A.elementary_row_op(op='n->kn', row=None, k=None, row1=None, row2=None) # ERO, "n->kn"/"n<->m"/"n->n+km"
    >>> A.is_nilpotent() # check if nilpotent
    >>> A.is_symmetric() # check if symmetric
    >>> A.minor(i, j) # (i, j) minor (WRONG DEFINITION, uses determinant)
    >>> A.nullspace() # nullspace
    >>> A.rank(iszerofunc=<function _iszero>, simplify=False) # rank
    >>> A.rowspace(simplify=False) # list of vectors that span row space of A
    >>> A.rref() # returns [rref, list_of_pivot_cols]
    >>> A.LUdecomposition(iszerofunc=<function _iszero>, simpfunc=None, rankcheck=False) # LU decomposition
    >>> A.lower_triangular_solve(rhs) # solve Ax = rhs, A lower triangular
    >>> A.upper_triangular_solve(rhs) # solve Ax = rhs, A upper triangular

    # Custom Commands (verbosity >= 1 returns ERO (idx + 1), >= 2 returns matrix at each step)
    >>> is_zero(expr, symbolic: bool = True)
    >>> Matrix.from_latex(expr, row_join=True, norm=False, aug_pos=None) # Parse LaTeX matrix/vector to Matrix
    >>> Matrix.from_str(matrix_str, row_sep=';', col_sep=' ', aug_pos=None, is_real=True) # Parse string to Matrix
    >>> Matrix.from_list(vectors, row_join=True, aug_pos=None) # Create Matrix from list of vectors
    >>> Matrix.create_unk_matrix(num_rows: int, num_cols: int, symbol: str, is_real: bool) # Matrix with symbolic entries
    >>> Matrix.create_rand_matrix(num_rows: int, num_cols: int) # Matrix with random entries
    >>> A.simplify(rational=True, tolerance=1e-4, simplify=True, expand=True, collect_sym=None) # Simplify entries
    >>> A.identify(tolerance: float) # Identify symbolic/numeric constants
    >>> A.elem() # Identity matrix with same number of rows as A
    >>> A.select_rows(*idx) # New matrix with selected rows
    >>> A.select_cols(*idx) # New matrix with selected columns
    >>> A.scale_row(idx: int, scalar: float, verbosity=0) # Scale row
    >>> A.swap_row(idx_1: int, idx_2: int, verbosity=0) # Swap rows
    >>> A.reduce_row(idx_1: int, scalar: float, idx_2: int, verbosity=0) # Row reduction
    >>> A.get_pivot_row(col_idx: int, row_start_idx: int, follow_GE=False) # Find pivot row
    >>> A.ref(verbosity=2, max_tries=2, follow_GE=False, matrices=2) # Row echelon form (REF)
    >>> A.evaluate_cases(rhs: Matrix, verbosity=0) # Display solution cases for symbolic systems
    >>> A.column_constraints(use_id=False, use_ref=False) # RREF of [A | b], constraints for b
    >>> A.extend_basis(span_subspace=A.elem()) # Augment basis to span subspace
    >>> A.transition_matrix(to: Matrix) # Transition matrix between bases
    >>> A.intersect_subspace(other: Matrix, verbosity=1) # Basis for intersection of subspaces
    >>> A.is_same_subspace(other: Matrix, verbosity=1) # Check if subspaces are equal
    >>> A.inverse(option: str, verbosity=0) # Inverse (left/right/both)
    >>> A.orthogonal_complement(verbosity=0) # Null(A^T)
    >>> A.is_vec_orthogonal(verbosity=1) # Check if columns are orthogonal
    >>> A.normalized(factor=False) # Normalize columns
    >>> A.scalar_factor(column=True) # Factor out common divisors
    >>> A.gram_schmidt(factor=True, verbosity=1) # Orthonormalize columns
    >>> A.QRdecomposition(full=False) # QR decomposition
    >>> A.solve_least_squares(rhs: Matrix, verbosity=1) # Least squares solution
    >>> A.cpoly(force_factor=True) # Characteristic polynomial
    >>> A.is_diagonalizable(reals_only=True, verbosity=1) # Diagonalizability
    >>> A.diagonalize(reals_only=True, verbosity=0) # Diagonalization
    >>> A.is_orthogonally_diagonalizable # Check if symmetric
    >>> A.orthogonally_diagonalize(reals_only=True, factor=True, verbosity=1) # Orthogonal diagonalization
    >>> A.equilibrium_vectors() # Probability vectors Ax = x
    >>> A.fast_svd(option='np', identify=True, tolerance=None) # Fast SVD (numeric)
    >>> A.singular_value_decomposition(verbosity=0) # Full SVD (A = U @ S @ V.T)
    """
    print(commands)


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
        if shell in ["ZMQInteractiveShell", "TerminalInteractiveShell", "Interpreter"]:
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
            - If none, assumes the object can be passed into IPython's [`display`][IPython.display.display] function directly.
        **kwargs: Additional keyword arguments to pass to the display function.

    See Also:
        - [`IPython.display`][IPython.display]: The classes used to display objects in IPython environments.
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


def _standardise_symbol(
    symbols: set[sym.Symbol], is_real: bool | None = None
) -> list[sym.Symbol]:
    """Standardizes the subscripts of symbols by converting them to a consistent format.

    Args:
        symbols (set[sym.Symbol]): A set of SymPy's symbols to standardize.
        is_real (bool, optional): Whether the symbols are real numbers.
    Returns:
        list[sym.Symbol]: A list of standardized sympy symbols.
    """
    pattern = r"([a-zA-Z]+)(\d+)"
    replacement = r"\1_\2"

    res = []
    for symbol in symbols:
        new_symbol = re.sub(pattern, replacement, str(symbol))
        res.append(sym.symbols(new_symbol, real=is_real))
    return res
