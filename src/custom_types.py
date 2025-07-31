from __future__ import annotations
from abc import abstractmethod
import dataclasses
from enum import Enum
from typing import Literal, NamedTuple, TYPE_CHECKING

from sympy.printing.latex import LatexPrinter

from utils import _gen_latex_repr

if TYPE_CHECKING:
    from symbolic import Matrix
    import numpy as np


class Shape(Enum):
    """Enum for different matrix shapes."""
    SCALAR = "SCALAR"
    STRICT_UPPER = "STRICT_UPPER"
    STRICT_LOWER = "STRICT_LOWER"
    UPPER = "UPPER"
    LOWER = "LOWER"
    SYMMETRIC = "SYMMETRIC"


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
        """Returns an iterator over the fields of the dataclass."""
        return iter(
            tuple(getattr(self, field.name) for field in dataclasses.fields(self))
        )

    def __getitem__(self, idx: int):
        """Gets the value of a field by index."""
        fields = dataclasses.fields(self)
        return getattr(self, fields[idx].name)

    def __setitem__(self, idx: int, value) -> None:
        """Sets the value of a field by index."""
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
        """Generates a LaTeX representation of the particular and general solution."""
        return (
            "\\left("
            + self.part_sol._latex(printer)
            + " + "
            + self.gen_sol._latex(printer)
            + "\\right)"
        )

    def eval(self) -> Matrix:
        """Evaluates the sum of the particular and general solution."""
        return self.part_sol + self.gen_sol


@dataclasses.dataclass
class ScalarFactor(Printable):
    """Represents a matrix factored into a diagonal and a full matrix."""
    diag: Matrix
    full: Matrix
    order: Literal["FD", "DF"]

    def _latex(self, printer=None) -> str:
        """Generates a LaTeX representation of the factored matrix."""
        if self.order == "FD":
            return self.full._latex(printer) + self.diag._latex(printer)
        else:
            return self.diag._latex(printer) + self.full._latex(printer)

    def eval(self) -> Matrix:
        """Evaluates the product of the factored matrix."""
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
        """Generates a LaTeX representation of the PLU decomposition."""
        return self.P._latex(printer) + self.L._latex(printer) + self.U._latex(printer)

    def eval(self) -> Matrix:
        """Evaluates the product of the PLU decomposition."""
        return self.P @ self.L @ self.U


@dataclasses.dataclass
class RREF(Printable):
    """Represents the reduced row echelon form of a matrix."""
    rref: Matrix
    pivots: tuple[int, ...]

    def eval(self) -> Matrix:
        """Returns the reduced row echelon form matrix."""
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


@dataclasses.dataclass
class QR(Printable):
    """Represents a QR decomposition of a matrix."""
    Q: Matrix
    R: Matrix

    def _latex(self, printer=None) -> str:
        """Generates a LaTeX representation of the QR decomposition."""
        return self.Q._latex(printer) + self.R._latex(printer)

    def eval(self) -> Matrix:
        """Evaluates the product of the QR decomposition."""
        return self.Q @ self.R


@dataclasses.dataclass
class PDP(Printable):
    """Represents a PDP diagonalization of a matrix."""
    P: Matrix
    D: Matrix

    def _latex(self, printer=None) -> str:
        """Generates a LaTeX representation of the PDP diagonalization."""
        P_inv = self.P.inverse(matrices=1)  # inv exists and is unique
        return self.P._latex(printer) + self.D._latex(printer) + P_inv._latex(printer)  # type: ignore

    def eval(self) -> Matrix:
        """Evaluates the product of the PDP diagonalization."""
        return self.P @ self.D @ self.P.inverse(matrices=1)


@dataclasses.dataclass
class SVD(Printable):
    """Represents a Singular Value Decomposition of a matrix."""
    U: Matrix
    S: Matrix
    V: Matrix

    def _latex(self, printer=None) -> str:
        """Generates a LaTeX representation of the SVD."""
        return (
            self.U._latex(printer) + self.S._latex(printer) + self.V.T._latex(printer)
        )

    def eval(self) -> Matrix:
        """Evaluates the product of the SVD."""
        return self.U @ self.S @ self.V.T


class NumSVD(NamedTuple):
    """Represents a numerical Singular Value Decomposition of a matrix."""
    U: np.typing.NDArray
    S: np.typing.NDArray
    V: np.typing.NDArray

    def __repr__(self) -> str:
        """Returns a string representation of the numerical SVD."""
        return f"""NumSVD(
        U = \n{self.U.__repr__()}, \n
        S = \n{self.S.__repr__()}, \n
        V = \n{self.V.__repr__()})"""