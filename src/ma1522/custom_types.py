from __future__ import annotations
from abc import abstractmethod
import dataclasses
from enum import Enum
from typing import TYPE_CHECKING, NamedTuple

from sympy.printing.latex import LatexPrinter

from .utils import _gen_latex_repr
# from ma1522 import utils

if TYPE_CHECKING:
    from typing import Literal
    import numpy as np
    from symbolic import Matrix


class Shape(Enum):
    """Enumeration for different matrix shapes and structural properties.

    This enum defines various matrix shapes that can be used to specify
    the structure of matrices in mathematical operations and optimizations.
    Each shape represents a specific pattern of zero and non-zero elements.

    Examples:
        >>> shape = Shape.SYMMETRIC
        >>> print(shape.value)
        'SYMMETRIC'
    """

    DIAGONAL = "DIAGONAL"
    r"""Diagonal matrix.
    
    A matrix where all off-diagonal entries are zero. Only elements on the
    main diagonal ($i,j$ where $i = j$) can be non-zero. The diagonal entries
    can have different values (unlike SCALAR matrices).
    
    Example:
        $$
        \begin{pmatrix}
        a & 0 & 0 \\
        0 & b & 0 \\
        0 & 0 & c
        \end{pmatrix}
        $$
    """

    SCALAR = "SCALAR"
    r"""Scalar matrix (diagonal matrix with equal diagonal entries).
    
    A diagonal matrix where all diagonal entries are equal to the same scalar value,
    and all off-diagonal entries are zero. This is also known as a scalar matrix
    or scalar multiple of the identity matrix.
    
    Example:
        $$
        \begin{pmatrix}
        c & 0 & 0 \\
        0 & c & 0 \\
        0 & 0 & c
        \end{pmatrix}
        $$
    """

    STRICT_UPPER = "STRICT_UPPER"
    r"""Strictly upper triangular matrix.
    
    A matrix where all elements on and below the main diagonal are zero.
    Only elements above the main diagonal ($i,j$ where $i < j$) can be non-zero.
    
    Example:
        $$
        \begin{pmatrix}
        0 & a & b \\
        0 & 0 & c \\
        0 & 0 & 0
        \end{pmatrix}
        $$
    """

    STRICT_LOWER = "STRICT_LOWER"
    r"""Strictly lower triangular matrix.
    
    A matrix where all elements on and above the main diagonal are zero.
    Only elements below the main diagonal ($i,j$ where $i > j$) can be non-zero.
    
    Example:
        $$
        \begin{pmatrix}
        0 & 0 & 0 \\
        a & 0 & 0 \\
        b & c & 0
        \end{pmatrix}
        $$
    """

    UPPER = "UPPER"
    r"""Upper triangular matrix.
    
    A matrix where all elements below the main diagonal are zero.
    Elements on and above the main diagonal ($i,j$ where $i <= j$) can be non-zero.
    
    Example:
        $$
        \begin{pmatrix}
        a & b & c \\
        0 & d & e \\
        0 & 0 & f
        \end{pmatrix}
        $$
    """

    LOWER = "LOWER"
    r"""Lower triangular matrix.
    
    A matrix where all elements above the main diagonal are zero.
    Elements on and below the main diagonal ($i,j$ where $i >= j$) can be non-zero.
    
    Example:
        $$
        \begin{pmatrix}
        a & 0 & 0 \\
        b & c & 0 \\
        d & e & f
        \end{pmatrix}
        $$
    """

    SYMMETRIC = "SYMMETRIC"
    r"""Symmetric matrix.
    
    A square matrix where elements are symmetric about the main diagonal,
    meaning $A_{i, j} = A_{j,i}$ for all valid indices $i$ and $j$.
    
    Example:
        $$
        \begin{pmatrix}
        a & b & c \\
        b & d & e \\
        c & e & f
        \end{pmatrix}
        $$
    """


#####################
# PRINTABLE OBJECTS #
#####################

PRINTER = LatexPrinter()


# Base class that all LaTeX objects should inherit
@dataclasses.dataclass
class Printable:
    """Base class for objects that can be printed as LaTeX."""

    def _latex(self, printer=None) -> str:
        """Generates a LaTeX representation of the object."""
        return _gen_latex_repr(self, printer)

    def _repr_latex_(self) -> str:
        """Returns the LaTeX representation for IPython."""
        return f"${self._latex(PRINTER)}$"

    def __iter__(self):
        return iter(
            tuple(getattr(self, field.name) for field in dataclasses.fields(self))
        )

    def __getitem__(self, idx: int):
        fields = dataclasses.fields(self)
        return getattr(self, fields[idx].name)

    def __setitem__(self, idx: int, value) -> None:
        fields = dataclasses.fields(self)
        setattr(self, fields[idx].name, value)

    @abstractmethod
    def eval(self) -> Matrix:
        """Evaluates the object to a matrix."""
        ...

    def evalf(self, *args, **kwargs):
        """Evaluates the object to a matrix of floats."""
        return (self.eval()).evalf(*args, **kwargs)


@dataclasses.dataclass
class PartGen(Printable):
    """Represents a matrix as a sum of a particular and general solution."""

    part_sol: Matrix
    gen_sol: Matrix

    def _latex(self, printer=None) -> str:
        return (
            "\\left("
            + self.part_sol._latex(printer)
            + " + "
            + self.gen_sol._latex(printer)
            + "\\right)"
        )

    def eval(self) -> Matrix:
        return self.part_sol + self.gen_sol


@dataclasses.dataclass
class ScalarFactor(Printable):
    """Represents a matrix factored into a diagonal and a full matrix."""

    diag: Matrix
    full: Matrix
    order: Literal["FD", "DF"]

    def _latex(self, printer=None) -> str:
        if self.order == "FD":
            return self.full._latex(printer) + self.diag._latex(printer)
        else:
            return self.diag._latex(printer) + self.full._latex(printer)

    def eval(self) -> Matrix:
        if self.order == "FD":
            return self.full @ self.diag
        else:
            return self.diag @ self.full


@dataclasses.dataclass
class PLU(Printable):
    """Represents a PLU decomposition of a matrix."""

    P: Matrix
    L: Matrix
    U: Matrix

    def _latex(self, printer=None) -> str:
        return self.P._latex(printer) + self.L._latex(printer) + self.U._latex(printer)

    def eval(self) -> Matrix:
        return self.P @ self.L @ self.U


@dataclasses.dataclass
class RREF(Printable):
    """Represents the reduced row echelon form of a matrix."""

    rref: Matrix
    pivots: tuple[int, ...]

    def eval(self) -> Matrix:
        return self.rref


# class Inverse(NamedTuple):
#     left: Optional[Matrix]
#     right: Optional[Matrix]

#     def _latex(self, printer=None) -> str:
#         return _gen_latex_repr(self, printer)

#     def _ipython_display_(self) -> None:
#         from IPython.display import display, Math

#         display(Math(self._latex(PRINTER)))
#         # IPython.display.display(IPython.display.Latex("$" + self._latex(PRINTER) + "$"))


# @dataclasses.dataclass
# class Inverse(Printable):
#     """Represents the inverse of a matrix."""
#     both: Matrix | PartGen


# @dataclasses.dataclass
# class LeftInverse(Printable):
#     """Represents the left inverse of a matrix."""
#     left: Matrix | PartGen


# @dataclasses.dataclass
# class RightInverse(Printable):
#     """Represents the right inverse of a matrix."""
#     right: Matrix | PartGen


@dataclasses.dataclass
class VecDecomp(Printable):
    """Represents a vector decomposition into projection and normal components."""

    proj: Matrix
    norm: Matrix

    def eval(self) -> Matrix:
        return self.proj + self.norm


@dataclasses.dataclass
class QR(Printable):
    """Represents a QR decomposition of a matrix."""

    Q: Matrix
    R: Matrix

    def _latex(self, printer=None) -> str:
        return self.Q._latex(printer) + self.R._latex(printer)

    def eval(self) -> Matrix:
        return self.Q @ self.R


@dataclasses.dataclass
class PDP(Printable):
    """Represents a PDP diagonalization of a matrix."""

    P: Matrix
    D: Matrix

    def _latex(self, printer=None) -> str:
        P_inv = self.P.inverse(matrices=1)  # inv exists and is unique
        return self.P._latex(printer) + self.D._latex(printer) + P_inv._latex(printer)  # type: ignore

    def eval(self) -> Matrix:
        return self.P @ self.D @ self.P.inverse(matrices=1)


@dataclasses.dataclass
class SVD(Printable):
    """Represents a Singular Value Decomposition of a matrix."""

    U: Matrix
    S: Matrix
    V: Matrix

    def _latex(self, printer=None) -> str:
        return (
            self.U._latex(printer) + self.S._latex(printer) + self.V.T._latex(printer)
        )

    def eval(self) -> Matrix:
        return self.U @ self.S @ self.V.T


class NumSVD(NamedTuple):
    """Represents a numerical Singular Value Decomposition of a matrix."""

    U: np.typing.NDArray
    S: np.typing.NDArray
    V: np.typing.NDArray

    def __repr__(self) -> str:
        return f"""NumSVD(
        U = \n{self.U.__repr__()}, \n
        S = \n{self.S.__repr__()}, \n
        V = \n{self.V.__repr__()})"""
