"""NUS MA1522 Linear Algebra for Computing library.

This package provides step-by-step workings for linear algebra algorithms with a focus on
symbolic computation and custom matrix operations. It extends SymPy's Matrix class to support
educational visualization and computation of linear algebra decompositions and transformations.

Main components:
- Matrix: Extended matrix class with symbolic computation capabilities
- Decomposition types: PLU, RREF, QR, PDP, SVD, and specialized variants
- Utilities: Helper functions for display and symbolic computation
"""

from .custom_types import (
    PDP,
    PLU,
    QR,
    RREF,
    SVD,
    NumSVD,
    PartGen,
    RREFCase,
    ScalarFactor,
    Shape,
    VecDecomp,
)
from .symbolic import Matrix
from .utils import display, sympy_commands

__all__ = [
    "PDP",
    "PLU",
    "QR",
    "RREF",
    "SVD",
    "Matrix",
    "NumSVD",
    "PartGen",
    "RREFCase",
    "ScalarFactor",
    "Shape",
    "VecDecomp",
    "display",
    "sympy_commands",
]
