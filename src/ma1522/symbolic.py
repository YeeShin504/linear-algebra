from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable
from typing import TYPE_CHECKING
from warnings import warn

import sympy as sym
from sympy.parsing.sympy_parser import parse_expr

import numpy as np
import mpmath as mp
from latex2sympy2 import latex2sympy
import re

# import IPython.display

from .utils import _powerset, _is_zero, _textify, display

from .custom_types import (
    Shape,
    PartGen,
    ScalarFactor,
    PLU,
    RREF,
    VecDecomp,
    QR,
    PDP,
    SVD,
    NumSVD,
)

if TYPE_CHECKING:
    from typing import Callable, Literal, DefaultDict

from sympy.core.expr import Expr
from sympy.core.symbol import Symbol
from sympy.core.mul import Mul

########
# MISC #
########


def sympy_commands():
    commands = """
    # Note : zero-indexing, see https://docs.sympy.org/latest/modules/matrices/matrices.html
    import sympy as sym

    # Variables
    >>> a, b = sym.symbols("a b") # symbols
    >>> I_3 = eye(3) # identity matrix
    >>> zeros = sym.zeros(2, cols=2) # zero matrix
    >>> A = sym.Matrix(
        [
            [...]
        ]
    ) # user-defined matrix
    >>> B = sym.Matrix.vstack(A, I_3, ...) # repeated application of col_join
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
    >>> A.adj() # defined as adjoint as per MA1522 definition
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
    >>> A.LUdecomposition(iszerofunc=<function _iszero>, simpfunc=None, rankcheck=False) # LU decomposition, returns (L, U, perm)
    >>> A.lower_triangular_solve(rhs) # solve the matrix equation Ax = rhs, where A is an lower triangular matrix
    >>> A.upper_triangular_solve(rhs) # solve the matrix equation Ax = rhs, where A is an upper triangular matrix

    # Custom Commands (verbosity >= 1 returns ERO (idx + 1), >= 2 returns matrix at each step)
    >>> is_zero(expr, symbolic: bool = True)
    >>> Matrix.from_latex(expr, row_join: bool = True, norm: bool = False) # creates a Matrix object from the LaTeX expression
    >>> Matrix.from_list(vectors: List, row_join: bool = True) # creates a Matrix object from a list of column vectors
    >>> Matrix.create_unk_matrix(num_rows: int, num_cols: int, symbol: str, is_real: bool) # creates a Matrix object with symbolic entries
    >>> Matrix.create_rand_matrix(num_rows: int, num_cols: int) # creates a Matrix object with random entries
    >>> A.simplify(rational: bool = True, tolerance: float, simplify: bool = True, expand: bool = True, collect_sym: sym.Symbol = None)
    >>> A.identify(tolerance: float) # returns a matrix simplified using mp.identify 
    >>> A.elem() # returns the identity matrix with same number of rows as A
    >>> A.select_rows(*idx) # returns a new matrix with row vectors of *idx 
    >>> A.select_cols(*idx) # returns a new matrix with column vectors of *idx 
    >>> A.scale_row(idx: int, scalar: float, verbosity: int = 0) 
    >>> A.swap_row(idx_1: int, idx_2: int, verbosity: int = 0)
    >>> A.reduce_row(idx_1: int, scalar: float, idx_2: int, verbosity: int = 0)
    >>> A.get_pivot_row(col_idx: int, row_start_idx: int, follow_GE: bool = False)
    >>> A.ref(verbosity: int = 2, max_tries: int = 2, follow_GE: bool = False, matrices: int = 2) # matrices = 1 (U), 2 (LU), 3 (PLU)
    >>> A.evaluate_cases(rhs: Matrix) # display a list of matrices for unknowns with critical values
    >>> A.column_constraints(use_id: bool = False, use_ref: bool = False) # returns the rref of [A | b], where b can be any vector. Use it to find constraints for b if A is not invertible.
    >>> A.extend_basis(span_subspace: Matrix = A.elem()) # returns an augmented matrix with additional column vectors to span the subspace of the argument.
    >>> A.transition_matrix(to: Matrix) # returns a transition matrix from A to the other matrix
    >>> A.intersect_subspace(other: Matrix, verbosity: int = 1) # returns a basis for the subspace that is in the intersection of A and other columnspace.
    >>> A.is_same_subspace(other: Matrix, verbosity: int = 1) # returns True if columnspace of A and other are equal
    >>> A.inverse(option: str, verbosity: int = 0) # returns an inverse of A. If A is non-square return either a left inverse or right inverse.
    >>> A.orthogonal_complement(verbosity: int = 0) # returns a matrix whose column vectors are perpendicular (independent) to the column vectors of A (i.e. Null(A^T))
    >>> A.is_vec_orthogonal(verbosity: int = 1) # checks if the column vectors of A are orthogonal
    >>> A.normalized(factor: bool = False) # returns a matrix whose columnn vectors are normalized
    >>> A.scalar_factor(column: bool = True) # returns 2 matrices, one of which is a diagonal matrix, used to simplify matrix expressions
    >>> A.gram_schmidt(factor: bool = True, verbosity: int = 1)  # returns a matrix with orthonormal columns (and displays gram-schmidt process)
    >>> A.QRdecomposition(full: bool = False) # QR decomposition, returns (Q, R)
    >>> A.solve_least_squares(rhs: Matrix, verbosity: int = 1) # solve the least square solution min(|Ax - rhs|)
    >>> A.cpoly(force_factor: bool = True) # returns a factorised characterstic polynomial for a square matrix
    >>> A.is_diagonalizable(reals_only: bool = True, verbosity: int = 1) # modified diagonalizability criteria to align with MA1522 (reals_only)
    >>> A.diagonalize(reals_only: bool = True, verbosity:  int = 0) # returns tuple of invertible P and diagonal D
    >>> A.is_orthogonally_diagonalizable # checks if matrix is symmetric
    >>> A.orthogonally_diagonalize(reals_only: bool = True, factor: bool = True, verbosity = 1) # returns tuple of orthogonal matrix P and diagonal D
    >>> A.equilibrium_vectors() # returns a matrix whose column vectors are probability vectors such that Ax = x
    >>> A.fast_svd(option: str = 'np', identify: bool = True, tolerance: float = None) # returns a svd via numerical methods. Attempt to convert back to symbolic via identify not guaranteed.
    >>> A.singular_value_decomposition(verbosity: int = 0) # returns the full SVD required by MA1522 (A = U @ S @ V.T)
    """
    print(commands)


sym.init_printing(use_unicode=True)
np.set_printoptions(formatter={"float": lambda x: f"{x:10.7g}"})


class Matrix(sym.MutableDenseMatrix):
    r"""A symbolic matrix class extending [`MutableDenseMatrix`][sympy.matrices.dense.MutableDenseMatrix] with enhanced linear algebra operations.

    The inherited methods from [`MutableDenseMatrix`][sympy.matrices.dense.MutableDenseMatrix]
    can be found in the [SymPy Matrices Documentation](https://docs.sympy.org/latest/modules/matrices/matrices.html). A summary of the
    inherited attributes and methods is also available on the [Inherited Methods Summary](inherited.md) page.

    This class provides comprehensive linear algebra functionality with support for:
        - Matrix creation from various sources (lists, $\rm\LaTeX$, random values)
        - Matrix decompositions (REF, RREF, LU, QR, SVD, diagonalization)
        - Vector space operations (orthogonalization, projections, basis manipulation)
        - Eigenvalue/eigenvector computations
        - Custom printing and $\rm\LaTeX$ formatting with augmented matrix support

    Key Features:
        - Maintains symbolic expressions throughout operations
        - Follows MA1522 syllabus conventions for linear algebra
        - Provides detailed step-by-step output for learning purposes
        - Supports both exact symbolic and numerical computations

    Attributes:
        _aug_pos (set[int]): Set of column indices where augmentation lines should be drawn
            for displaying augmented matrices.

    Examples:
        Basic matrix operations:
        >>> A = Matrix([[1, 2], [3, 4]])
        >>> A.rref()
        RREF(rref=Matrix([
            [1, 0]
            [0, 1]
        ]), pivots=(0, 1))

        Creating from LaTeX:
        >>> B = Matrix.from_latex(r'\begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}')
        ⎡1  2⎤
        ⎢    ⎥
        ⎣3  4⎦

        Eigenvalue decomposition:
        >>> A.diagonalize()
        >>> P, D = A.diagonalize()
        >>> P, D
        (Matrix([
        [-sqrt(33)/6 - 1/2, -1/2 + sqrt(33)/6],
        [                1,                 1]]), Matrix([
        [5/2 - sqrt(33)/2,                0],
        [               0, 5/2 + sqrt(33)/2]]))
    """

    def __init__(
        self,
        *args,
        aug_pos: Iterable[int] | int | None = None,
        **kwargs,
    ) -> None:
        if aug_pos is None:
            self._aug_pos = set()
        elif isinstance(aug_pos, int):
            self._aug_pos = set([aug_pos])
        elif isinstance(aug_pos, Iterable) and all(isinstance(i, int) for i in aug_pos):
            self._aug_pos = set(aug_pos)
        else:
            raise TypeError(
                f"Invalid type for aug_pos: {type(aug_pos)}. Expected Iterable[int]."
            )

    def __str__(self) -> str:
        res = super().__str__()
        if not hasattr(self, "_aug_pos"):
            # Matrices produced by parent methods may not have _aug_pos
            return res
        aug_pos = f"aug_pos: {self._aug_pos}"
        return f"{res}, {aug_pos}"

    def __repr__(self) -> str:
        if not hasattr(self, "_aug_pos"):
            # Matrices produced by parent methods may not have _aug_pos
            return super().__repr__()

        def rep_row(row: str, pos_set: set[int]) -> str:
            repr = ""
            elems = row.removesuffix(",").split(",")
            for idx, elem in enumerate(elems):
                if idx in pos_set:
                    repr += elem + " |"
                else:
                    repr += elem + ","
            return repr.removesuffix("|").removesuffix(",")

        res = super().__repr__().removeprefix("Matrix([").removesuffix("])")
        res_row_list = (rep_row(row, self._aug_pos) for row in res.split("\n"))
        return "Matrix([" + "\n".join(res_row_list) + "\n])"

    def __eq__(self, other) -> bool:
        if not hasattr(self, "_aug_pos"):
            self._aug_pos = set()
        if not hasattr(other, "_aug_pos"):
            other._aug_pos = set()

        return super().__eq__(other) and (self._aug_pos == other._aug_pos)

    # Override
    def _latex(self, printer=None) -> str:
        raw = printer._print(sym.Matrix(self))  # type: ignore
        if not hasattr(self, "_aug_pos"):
            # Matrices produced by parent methods may not have _aug_pos
            return raw

        # get latex representation of matrix with "array" format
        raw = sym.latex(sym.Matrix(self), mat_str="array")
        array_c = "\\begin{array}{" + "c" * self.cols + "}"
        for mat_str in ["smallmatrix", "matrix"]:
            old_beg = "\\begin{" + mat_str + "}"
            old_end = "\\end{" + mat_str + "}"
            raw = raw.replace(old_beg, array_c).replace(old_end, "\\end{array}")

        ls = [
            pos for pos in self._aug_pos if 0 <= pos < self.cols - 1
        ]  # remove trailing pos
        if len(ls) == 0:
            # no valid _aug_pos found
            return raw
        ls.sort()

        # create formatting string s to insert augment line visually
        delta = [ls[0]]
        delta.extend([ls[i] - ls[i - 1] for i in range(1, len(ls))])
        remainder = self.cols - sum(delta) - 1
        delta.append(remainder)
        s = "\\begin{array}{c" + "|".join(("c" * i for i in delta)) + "}"
        array_c = "\\begin{array}{" + "c" * self.cols + "}"

        return raw.replace(array_c, s)

    ###################
    # FACTORY METHODS #
    ###################
    @staticmethod
    def from_latex(
        expr: str, verbosity: int = 1, row_join: bool = True, norm: bool = False
    ) -> Matrix:
        r"""Converts a $\rm\LaTeX$ matrix/vector expression into a Matrix object.

        Parses $\rm\LaTeX$ matrix environments (pmatrix, array) and vector lists into a Matrix.
        Handles matrix multiplication expressions and normalizes vectors when requested.

        Args:
            expr (str): $\rm\LaTeX$ string containing:

                - Matrix environments:
                    `\begin{pmatrix} ... \end{pmatrix}`, `\begin{array}{ccc} ... \end{array}`
                - Vector lists:
                    `\{ \begin{pmatrix} ... \end{pmatrix}, \begin{pmatrix} ... \end{pmatrix}\}`
                - Matrix products:
                    `\begin{pmatrix}...\end{pmatrix}\begin{pmatrix} ... \end{pmatrix}`

            verbosity (int): Controls output detail level:
                - 0: No output
                - 1: Display parsed matrix

            row_join (bool): If True, vector lists are treated as columns.
                If False, vectors are treated as rows.

            norm (bool): If True, normalizes vectors to unit length (default False).

        Returns:
            (Matrix): The parsed matrix with optional normalization.

        Raises:
            Exception: If the $\rm\LaTeX$ expression is empty, invalid or cannot be parsed.

        Examples:
            >>> Matrix.from_latex(r'\begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}', verbosity=0)
            Matrix([
            [1, 2]
            [3, 4]
            ])

            >>> Matrix.from_latex(r'\begin{pmatrix} 1 \\ 2 \\ 3 \end{pmatrix}', verbosity=0, norm=True)
            Matrix([
            [  sqrt(14)/14]
            [   sqrt(14)/7]
            [3*sqrt(14)/14]
            ])
        """

        # Step 1: Modify the LaTeX string to ensure compatibility with the parser.
        # Convert array-like LaTeX to pmatrix for proper matrix formatting
        # Replace \begin{array}{ccc*} with \begin{pmatrix}
        modified_latex = re.sub(
            r"\\begin\{array\}(?:\{[^}]*\})?", r"\\begin{pmatrix}", expr
        )
        # Replace \end{array} with \end{pmatrix}
        modified_latex = re.sub(r"\\end\{array\}", r"\\end{pmatrix}", modified_latex)
        # Remove LaTeX semicolon for cleaner parsing
        modified_latex = re.sub(r"\\;", "", modified_latex)

        # Step 2: Use latex2sympy to parse the modified LaTeX expression into SymPy Matrix
        res = latex2sympy(modified_latex)
        if verbosity >= 1:
            display(res)

        # Step 3: Handle the parsed result based on its type (list, MatMul, or Matrix)
        mat = Matrix()
        if isinstance(res, list):
            vector_list = []
            for vector in res:
                vector = vector.expand()
                vector_list.append(vector)
            mat = Matrix.from_list(vector_list, row_join)
        # elif isinstance(res, sym.matrices.expressions.matmul.MatMul):
        elif isinstance(res, sym.MatMul):
            # If the matrix is a product of matrices, evaluate the product directly
            mat = Matrix(res.doit())
        elif isinstance(res, sym.Matrix):
            # Directly converts the SymPy Matrix into the custom Matrix object to inherit the custom methods
            mat = Matrix(res)
        else:
            # If the result is neither a list nor a matrix expression, return the raw result
            mat = res

        # Step 4: Normalisation
        if norm and isinstance(mat, Matrix):
            return mat.normalized(factor=False)  # type: ignore
        return mat

    @staticmethod
    def from_str(matrix_str: str, row_sep: str = ";", col_sep: str = " ") -> Matrix:
        """Parses a string representation of a matrix into a Matrix.

        This method enables quick creation of a Matrix object from a string format similar to
        the one used in MATLAB. It supports both row and column separators and
        uses SymPy's [`parse_expr`][sympy.parsing.sympy_parser.parse_expr] to convert
        the entries of the matrix from a string format into a Matrix object.

        Args:
            matrix_str (str): The string representation of the matrix.
            row_sep (str, optional): The separator for rows in the string.
            col_sep (str): The separator for columns in the string.

        Returns:
            (Matrix): A Matrix object representing the parsed matrix.

        Raises:
            SyntaxError: If the string cannot be parsed into a matrix.

        Examples:
            >>> Matrix.from_str("[1 2; 3 4]")
            Matrix([
            [1, 2]
            [3, 4]
            ])
        """

        matrix_str = (
            matrix_str.strip().removeprefix("[").removesuffix("]")
        )  # remove surrounding brackets
        rows = matrix_str.strip().split(row_sep)
        matrix = []
        for row in rows:
            cols = row.strip().split(col_sep)
            matrix.append([parse_expr(col.strip()) for col in cols])
        return Matrix(matrix)

    @staticmethod
    def from_list(vectors: list[Matrix], row_join: bool = True) -> Matrix:
        """Creates a Matrix object from a list of vectors.

        This method takes a list of vectors (each represented as a Matrix object)
        and combines them into a single matrix.

        Args:
            vectors (list[Matrix]): A list of Matrix objects, where each Matrix
                represents a row or column vector.

            row_join (bool, optional): If True, the vectors are joined horizontally.
                If False, the vectors are stacked vertically.

        Returns:
            (Matrix): A matrix constructed from the list of vectors.

        Raises:
            sympy.matrices.exceptions.ShapeError: If the vectors do not have compatible dimensions for joining.

        Examples:
            >>> vec1 = Matrix([[1], [2]])
            >>> vec2 = Matrix([[3], [4]])
            >>> Matrix.from_list([vec1, vec2])
            Matrix([
            [1, 3]
            [2, 4]
            ])

            >>> Matrix.from_list([vec1, vec2], row_join=False)
            Matrix([
            [1]
            [2]
            [3]
            [4]
            ])
        """
        if not vectors:
            return Matrix([])
        res = vectors.pop(0)
        for vec in vectors:
            if row_join:
                res = res.row_join(vec, aug_line=False)
            else:
                res = res.col_join(vec)
        return res

        # res = Matrix([[elem for elem in vector] for vector in vectors])
        # if row_join:
        #     # If `row_join=True`, transpose the matrix (vectors treated as rows)
        #     return res.T
        # else:
        #     # If `row_join=False`, return the matrix as-is (vectors treated as columns)
        #     return res

    def _shape(self, shape: Shape) -> Matrix:
        # if self.rows != self.cols:
        #     raise sym.NonSquareMatrixError()
        match shape:
            case Shape.DIAGONAL:
                res = Matrix.diag(*self.diagonal())
                if self.rows > self.cols:
                    res = res.col_join(Matrix.zeros(self.rows - self.cols, self.cols))
                elif self.rows < self.cols:
                    res = res.row_join(
                        Matrix.zeros(self.rows, self.cols - self.rows), aug_line=False
                    )
                return res
            case Shape.SCALAR:
                if self.rows != self.cols:
                    raise sym.NonSquareMatrixError(
                        "Scalar shape is only defined for square matrices."
                    )
                return self.diagonal()[0] * self.elem()
            case Shape.UPPER:
                return self.upper_triangular()
            case Shape.LOWER:
                return self.lower_triangular()
            case Shape.STRICT_UPPER:
                return self._shape(Shape.UPPER) - self._shape(Shape.DIAGONAL)
            case Shape.STRICT_LOWER:
                return self._shape(Shape.LOWER) - self._shape(Shape.DIAGONAL)
            case Shape.SYMMETRIC:
                if self.rows != self.cols:
                    raise sym.NonSquareMatrixError(
                        "Symmetric shape is only defined for square matrices."
                    )
                return self._shape(Shape.UPPER) + self._shape(Shape.STRICT_UPPER).T

    @staticmethod
    def create_unk_matrix(
        r: int = 1,
        c: int = 1,
        symbol: str | None = None,
        is_real: bool = True,
        shape: Shape | None = None,
    ) -> Matrix:
        r"""Creates a symbolic matrix with unknown entries.

        This method generates a matrix of size $r \times c$ with symbolic
        entries. The entries are named based on the provided `symbol` parameter and
        indexed by their row and column positions. The `is_real` flag determines whether
        the symbols are real-valued.

        Note:
            - For a column vector without a specified symbol, the entries will be named
                following conventions, i.e., $\begin{pmatrix} x \end{pmatrix}$,
                $\begin{pmatrix} x \\ y \end{pmatrix}$,
                $\begin{pmatrix} x \\ y \\ z \end{pmatrix}$, for 1, 2 and 3 rows respectively.

        Args:
            r (int, optional): The number of rows in the matrix.
            c (int, optional): The number of columns in the matrix.
            symbol (str, optional): The base name for the symbols used in the matrix entries.
            is_real (bool, optional): If True, the symbols are real-valued;
                otherwise, they are complex.
            shape (Shape, optional): If provided, the matrix will be reshaped to this
                specific shape. Supported shapes include:

                - [DIAGONAL][(p).Shape.DIAGONAL]: Returns a diagonal matrix.
                - [SCALAR][(p).Shape.SCALAR]: Returns a scalar matrix.
                - [UPPER][(p).Shape.UPPER]: Returns an upper triangular matrix.
                - [LOWER][(p).Shape.LOWER]: Returns a lower triangular matrix.
                - [STRICT_UPPER][(p).Shape.STRICT_UPPER]: Returns an upper triangular matrix without the diagonal.
                - [STRICT_LOWER][(p).Shape.STRICT_LOWER]: Returns a lower triangular matrix without the diagonal.
                - [SYMMETRIC][(p).Shape.SYMMETRIC]: Returns a symmetric matrix.

        Returns:
            (Matrix): A matrix with symbolic entries of the specified size.

        Raises:
            sympy.matrices.exceptions.NonSquareMatrixError: If `shape` is ill-defined on a non-square matrix.

        Examples:
            >>> Matrix.create_unk_matrix(2, 2, symbol='a')
            Matrix([
            [a_1,1, a_1,2],
            [a_2,1, a_2,2]])

            >>> Matrix.create_unk_matrix(3, 1, symbol='y')
            Matrix([
            [y_1]
            [y_2]
            [y_3]
            ])

        See Also:
            - [`create_rand_matrix`][..]: Creates a matrix with random entries.
        """

        # Create a vector of size rows with entries (x, y) or (..., y, z)
        if r <= 26 and c == 1 and symbol is None:
            ls = []
            match r:
                case 1:
                    ls = list("x")
                case 2:
                    ls = list("xy")
                case 3:
                    ls = list("xyz")
                case _:
                    ascii_lowercase = list("abcdefghijklmnopqrstuvwxyz")
                    ls = ascii_lowercase[26 - r :]
            entries = sym.symbols(ls, real=is_real)
            return Matrix(entries)

        if symbol and c == 1:
            entries = sym.symbols(f"{symbol}_(1:{r + 1})", real=is_real)
            return Matrix(entries)

        # Creates a matrix of size rows * cols with entries symbol_i,j
        symbol = symbol or "x"  # default is "x"
        entries = sym.symbols(f"{symbol}_(1:{r + 1})\\,(1:{c + 1})", real=is_real)

        res = Matrix(entries).reshape(r, c)
        if shape:
            return res._shape(shape)
        else:
            return res

    @staticmethod
    def create_rand_matrix(
        r: int = 1,
        c: int = 1,
        shape: Shape | None = None,
        *args,
        **kwargs,
    ) -> Matrix:
        r"""Creates a matrix with random entries.

        This method generates a matrix of size $r \times c$ where the
        entries are real integers. The values in the matrix are generated using SymPy's
        [`randMatrix`][sympy.matrices.dense.randMatrix] function.

        Note:
            - The entries in the matrix are generated randomly and will change each time
                the function is called. Setting a random seed using `seed` in `**kwargs`
                will ensure reproducibility of the random values.

        Args:
            r (int, optional): The number of rows in the matrix.
            c (int, optional): The number of columns in the matrix.
            shape (Shape, optional): If provided, the matrix will be reshaped to this
                specific shape. Supported shapes include:

                - [DIAGONAL][(p).Shape.DIAGONAL]: Returns a diagonal matrix.
                - [SCALAR][(p).Shape.SCALAR]: Returns a scalar matrix.
                - [UPPER][(p).Shape.UPPER]: Returns an upper triangular matrix.
                - [LOWER][(p).Shape.LOWER]: Returns a lower triangular matrix.
                - [STRICT_UPPER][(p).Shape.STRICT_UPPER]: Returns an upper triangular matrix without the diagonal.
                - [STRICT_LOWER][(p).Shape.STRICT_LOWER]: Returns a lower triangular matrix without the diagonal.
                - [SYMMETRIC][(p).Shape.SYMMETRIC]: Returns a symmetric matrix.
            *args: Additional arguments passed to the [`randMatrix`][sympy.matrices.dense.randMatrix]
                function.
            **kwargs: Additional arguments passed to the [`randMatrix`][sympy.matrices.dense.randMatrix]
                function.

        Returns:
            (Matrix): A Matrix with random entries of the specified size.

        Raises:
            sympy.matrices.exceptions.NonSquareMatrixError: If `shape` is ill-defined on a non-square matrix.

        Examples:
            >>> Matrix.create_rand_matrix(2, 3, seed=42)
            Matrix([
            [81, 14,  3]
            [94, 35, 31]
            ])

        See Also:
            - [`create_unk_matrix`][..]: Creates a matrix with symbolic entries.
            - [`randMatrix`][sympy.matrices.dense.randMatrix]: SymPy function to create a random matrix.
        """
        res = Matrix(sym.randMatrix(*args, r=r, c=c, **kwargs))
        if shape:
            return res._shape(shape)
        else:
            return res

    # Override
    @staticmethod
    def eye(*args, **kwargs) -> Matrix:
        aug_pos = kwargs.pop("aug_pos", None)
        return Matrix(sym.eye(*args, **kwargs), aug_pos=aug_pos)

    # Override
    @staticmethod
    def zeros(*args, **kwargs) -> Matrix:
        aug_pos = kwargs.pop("aug_pos", None)
        return Matrix(sym.zeros(*args, **kwargs), aug_pos=aug_pos)

    # Override
    @staticmethod
    def ones(*args, **kwargs) -> Matrix:
        aug_pos = kwargs.pop("aug_pos", None)
        return Matrix(sym.ones(*args, **kwargs), aug_pos=aug_pos)

    # Override
    @staticmethod
    def diag(*args, **kwargs) -> Matrix:
        aug_pos = kwargs.pop("aug_pos", None)
        return Matrix(sym.diag(*args, **kwargs), aug_pos=aug_pos)

    # Override
    @property
    def T(self) -> Matrix:
        return Matrix(super().T)

    ###################################################
    # OVERRIDE OVERLOADED PYTHON ARITHMETIC OPERATORS #
    ###################################################

    def __abs__(self) -> Matrix:
        aug = self._aug_pos.copy() if hasattr(self, "_aug_pos") else set()
        return Matrix(super().__abs__(), aug_pos=aug)

    def __add__(self, other: Matrix) -> Matrix:
        if not hasattr(self, "_aug_pos"):
            self._aug_pos = set()
        if not hasattr(other, "_aug_pos"):
            other._aug_pos = set()

        aug = self._aug_pos | other._aug_pos
        return Matrix(super().__add__(other), aug_pos=aug)

    def __mul__(self, other) -> Matrix:
        res = super().__mul__(other)
        if hasattr(other, "shape") and self.shape != res.shape:
            return Matrix(res)
        else:
            # probably scalar multiplication
            aug = self._aug_pos.copy() if hasattr(self, "_aug_pos") else set()
            return Matrix(res, aug_pos=aug)

    def __rmul__(self, other) -> Matrix:
        res = super().__rmul__(other)
        if hasattr(other, "shape") and self.shape != res.shape:
            return Matrix(res)
        else:
            # probably scalar multiplication
            aug = self._aug_pos.copy() if hasattr(self, "_aug_pos") else set()
            return Matrix(res, aug_pos=aug)

    def __neg__(self) -> Matrix:
        aug = self._aug_pos.copy() if hasattr(self, "_aug_pos") else set()
        return Matrix(super().__neg__(), aug_pos=aug)

    ######################
    # BASIC MANIPULATORS #
    ######################

    # Override
    def copy(self) -> Matrix:
        """
        Creates a copy of the matrix, preserving augmentation lines.

        This method returns a new [`Matrix`][...] object that is a deep copy of the current matrix,
        including any augmentation line positions (used for displaying augmented matrices).

        Returns:
            (Matrix): A new matrix object with the same entries and augmentation lines as the original.

        Examples:
            >>> mat = Matrix([[1, 2], [3, 4]], aug_pos={1})
            >>> mat_copy = mat.copy()
            >>> mat_copy == mat
            True
            >>> mat_copy is mat
            False
        """
        new_mat = super().copy()
        aug = self._aug_pos.copy() if hasattr(self, "_aug_pos") else set()
        return Matrix(new_mat, aug_pos=aug)

    # Override
    def subs(self, *args, **kwargs) -> Matrix:
        """Substitutes values in the matrix entries.

        This method overrides SymPy's [`subs`][sympy.matrices.matrixbase.MatrixBase.subs]
        method to ensure that the augmentation lines are preserved after substitution.

        Args:
            *args: Positional arguments for substitution.
            **kwargs: Keyword arguments for substitution.

        Returns:
            (Matrix): A new matrix with substituted values, preserving augmentation lines.
        """
        new_mat = super().subs(*args, **kwargs)
        aug = self._aug_pos.copy() if hasattr(self, "_aug_pos") else set()
        return Matrix(new_mat, aug_pos=aug)

    # Override
    def simplify(
        self,
        rational: bool = True,
        suppress_warnings: bool = False,
        tolerance: float = 1e-4,
        simplify: bool = True,
        expand: bool = True,
        collect_sym: Symbol | None = None,
        *args,
        **kwargs,
    ) -> None:
        """Simplifies the matrix by applying various simplification techniques.

        This method performs several operations on the matrix to simplify its entries:
            - Rational simplification.
            - General symbolic simplification.
            - Expansion or factoring of expressions.
            - Collecting terms involving a specific symbol (if provided).

        Note:
            - Rational simplification attempts to convert entries into rational numbers if possible.
                If there is a residue (e.g. attempting to convert a non-rational number into a rational),
                a warning is printed with the approximation error.
            - Expansion and factoring can be controlled by the `expand` parameter.
            - The matrix is modified in place.

        Args:
            rational (bool, optional): If True, applies rational simplification
                to the matrix entries using [`sym.nsimplify`][sympy.simplify.simplify.nsimplify].
            suppress_warnings (bool, optional): If True, suppresses warnings about non-zero residues
                after rational simplification.
            tolerance (float, optional): The tolerance for rational simplification.
            simplify (bool, optional): If True, applies general symbolic simplification using [`sym.simplify`][sympy.simplify.simplify.simplify].
            expand (bool, optional): If True, applies expansion to the matrix entries. If False, applies factoring instead.
            collect_sym (Symbol, optional): A symbol to collect terms with. If provided,
                [`sym.collect`][sympy.simplify.radsimp.collect] will be applied to all entries of the matrix with respect to this symbol.
            *args: Additional arguments passed to the [`sym.simplify`][sympy.simplify.simplify.simplify] function.
            **kwargs: Additional arguments passed to the [`sym.simplify`][sympy.simplify.simplify.simplify] function.

        Returns:
            (Matrix): A new simplified matrix with the applied operations.

        Examples:
            >>> mat = Matrix([[sym.symbols('x') + 1, sym.symbols('x') + 2], [sym.symbols('x') + 3, sym.symbols('x') + 4]])
            >>> mat.simplify(rational=False, expand=True)
            Matrix([[x + 1, x + 2], [x + 3, x + 4]])
        """

        temp = self.copy()
        if rational:
            temp = sym.nsimplify(temp, tolerance=tolerance, rational=True)
            residues = (temp - self).norm()
            if residues != 0 and not suppress_warnings:
                res = residues.evalf()
                warn(
                    f"""
                    Non-zero Approximation Error: {res}
                    Rational approximation might have failed. Try lower tolerance.""",
                    UserWarning,
                )
        if simplify:
            temp = sym.simplify(temp, *args, **kwargs)
        if expand:
            temp = sym.expand(temp)
        else:
            temp = temp.applyfunc(lambda x: sym.factor(x))
        if collect_sym is not None:
            temp = temp.applyfunc(lambda x: sym.collect(x, collect_sym))

        # Create a new Matrix object from the simplified list and update the original object
        aug = self._aug_pos.copy() if hasattr(self, "_aug_pos") else set()
        temp = Matrix(temp, aug_pos=aug)  # prevent SymPy from losing _aug_pos
        self.__dict__.update(temp.__dict__)

    def identify(self, tol: float | None = None, *args, **kwargs) -> Matrix:
        r"""Identifies the matrix by applying a transformation function to each entry.

        This method applies a transformation to each element of the matrix using
        the [`identify`][mpmath.identify] function from the `mp` module. After identification,
        the method checks if there is any residue (i.e., if the matrix has been modified).

        Note:
            - If there is a residue (i.e., unable to identify an entry, such as $\pi$),
                a warning is printed with the approximation error. This can be resolved by
                lowering `tolerance` or supplying appropriate `constants` as `**kwargs` to
                the `identify` function.

        Args:
            tol (float, optional): A tolerance value that is passed to the
                [`identify`][mpmath.identify] function. If None, no tolerance is applied.
            *args: Additional positional arguments passed to the
                [`identify`][mpmath.identify] function.
            **kwargs: Additional keyword arguments passed to the
                [`identify`][mpmath.identify] function.

        Returns:
            (Matrix): A new matrix that results from applying the transformation to
                each element of the original matrix.

        Examples:
            >>> import math
            >>> mat = Matrix([[math.sqrt(4), math.e], [1/math.sqrt(2), 0.0]])
            >>> mat.identify()
            Matrix([
            [        1, E]
            [sqrt(2)/2, 0]
            ])

        See Also:
            - [`mpmath.identify`][mpmath.identify]: The function used to identify
              the entries of the matrix.
            - [`simplify`][..simplify]:
              For general simplification of the matrix entries.
        """

        temp = self.applyfunc(lambda x: mp.identify(x, tol=tol, *args, **kwargs))
        residues = (temp - self).norm()
        if residues != 0:
            res = residues.evalf()
            warn(f"Non-zero Identification Error: {res}", UserWarning)
        return temp

    def select_cols(self, *args: int) -> Matrix:
        """Selects columns from the matrix based on the provided column indices.

        This method returns a new matrix consisting of the columns specified by the
        provided indices. The columns are selected from the original matrix, and the
        result is returned as a new matrix.

        Args:
            *args (int): One or more column indices (0-based) to select from the matrix.

        Returns:
            (Matrix): A new matrix consisting of the selected columns.

        Examples:
            >>> mat = Matrix([[1, 2, 3], [4, 5, 6]])
            >>> mat.select_cols(0, 2)
            Matrix([
            [1, 3]
            [4, 6]
            ])
        """

        res = []
        for idx in args:
            res.append(list(self.col(idx)))
        return Matrix(res).T

    def select_rows(self, *args: int) -> Matrix:
        """Selects rows from the matrix based on the provided row indices.

        This method returns a new matrix consisting of the rows specified by the
        provided indices. The rows are selected from the original matrix, and the
        result is returned as a new matrix.

        Args:
            *args (int): One or more row indices (0-based) to select from the matrix.

        Returns:
            (Matrix): A new matrix consisting of the selected rows.

        Examples:
            >>> mat = Matrix([[1, 2, 3], [4, 5, 6]])
            >>> mat.select_rows(0)
            Matrix([[1, 2, 3]
            ])
        """

        res = []
        for idx in args:
            res.append(list(self.row(idx)))
        return Matrix(res)

    def sep_part_gen(self) -> PartGen:
        """
        Separates a matrix into its particular and general solution parts.

        This method separates the matrix into two components:
            - The **particular solution**, which is the solution to the system when
                all free variables are set to zero.
            - The **general solution**, which is the full solution including the
                homogeneous part.

        It assumes that the matrix is in symbolic form and contains free variables that can be set to zero.

        Returns:
            (PartGen): A dataclass containing two matrices:

                - `part_sol` ([`Matrix`][...]): The particular solution
                    (with free variables set to zero).
                - `gen_sol` ([`Matrix`][...]): The general solution (the original matrix
                    minus the particular solution).

        Examples:
            >>> from sympy import symbols
            >>> x = symbols('x')
            >>> mat = Matrix([[x + 2, 3], [3*x, x - 2]])
            >>> PG = mat.sep_part_gen()
            >>> PG.part_sol
            Matrix([
            [2,  3],
            [0, -2]]) # Particular solution
            >>> PG.gen_sol
            Matrix([
            [  x, 0]
            [3*x, x]
            ])  # General solution
        """

        set_0 = dict(((symbol, 0) for symbol in self.free_symbols))
        part_sol = self.subs(set_0)
        gen_sol = self - part_sol
        return PartGen(part_sol, gen_sol)

    def sep_unk(self) -> dict[Expr, Matrix]:
        """Separates the matrix into matrices with each free symbol set to 1.

        Returns:
            (dict[Expr, Matrix]): Returns a dictionary where the sum of the key*value pairs
                reconstructs the original matrix. Each key is a free symbol, and each value is a
                matrix with that symbol set to 1 and all other free symbols set to 0.
        """
        syms = self.free_symbols
        res: dict[Expr, Matrix] = defaultdict(Matrix)
        for s in syms:
            sub = dict(((symbol, 0) for symbol in syms if symbol != sym))
            sub[s] = 1  # type: ignore
            res[s] = self.subs(sub)
        return res

    def scalar_factor(self, column: bool = True) -> ScalarFactor:
        r"""Factorizes a matrix into the form $\mathbf{A} = \mathbf{FD}$, where $\mathbf{D}$ is a diagonal matrix
        and $\mathbf{F}$ contains the vectors with common divisors factored out (if `column=True`). If `column=False`,
        then returns $\mathbf{A} = \mathbf{DF}$ instead.

        Args:
            column (bool): If `True`, factorizes by columns. If `False`, factorizes by rows.

        Returns:
            (ScalarFactor): A dataclass of two matrices (F, D) and order (FD or DF)

                - `diag` ([`Matrix`][...]): The diagonal matrix containing the common divisors.
                - `full` ([`Matrix`][...]): The matrix with common divisors factored out.
                - `order` (str): The order of the factorization, either "FD" (for column factorization) or "DF" (for row factorization).

        Examples:
            >>> mat = Matrix([[6, 9], [12, 18]])
            >>> SF = mat.scalar_factor(column=True)
            >>> SF.full, SF.diag
            (Matrix([
            [1, 3]
            [2, 5]
            ]), Matrix([
            [6, 0]
            [0, 3]
            ]))
        """

        scalars = []
        F = self.copy()
        if column:
            for i in range(self.cols):
                g = sym.gcd(tuple(self.col(i)))
                F[:, i] /= g
                scalars.append(g)
            D = Matrix.diag(*scalars)
            assert self == F @ D
            return ScalarFactor(diag=D, full=F, order="FD")
        else:
            for i in range(self.rows):
                g = sym.gcd(tuple(self.row(i)))
                F[i, :] /= g
                scalars.append(g)
            D = Matrix.diag(*scalars)
            assert self == D @ F
            return ScalarFactor(diag=D, full=F, order="DF")

    #############################
    # CHAPTER 1: LINEAR SYSTEMS #
    #############################

    def aug_line(self, pos: int = -1) -> Matrix:
        """Inserts an augmented line at the specified position.

        This method adds an augmented line (i.e., a visual vertical line)
        to the matrix at the specified column position. If no position is provided (default: -1),
        the line is inserted at the last column.

        Note:
            - The method updates the `_aug_pos` attribute to track the position of the inserted line.
            - Negative `pos` will be converted before inserting it into `_aug_pos`.

        Args:
            pos (int, optional):
                The position (column index) where the augmented line will be inserted.
                Default -1 means the augmented line is added at the end of the matrix.

        Returns:
            (Matrix): The current matrix with the augmented line added at the specified position.

        Raises:
            IndexError: If the `pos` is out of range.

        Examples:
            >>> mat = Matrix([[1, 2], [3, 4]])
            >>> mat.aug_line(0)
            Matrix([
            [1 | 2]
            [3 | 4]
            ])
            >>> mat.aug_line().row_join(mat)
            Matrix([
            [1, 2 | 1, 2]
            [3, 4 | 3, 4]
            ])

        See Also:
            - [`rm_aug_line`][..]: Removes an augmentation line from the matrix.
        """

        new_pos = pos
        if new_pos < 0:
            new_pos += self.cols

        if not 0 <= new_pos < self.cols:
            raise IndexError(
                f"Position for augmented line ({pos}) out of range ({self.cols})."
            )

        if not hasattr(self, "_aug_pos"):
            self._aug_pos: set[int] = set()
        self._aug_pos.add(new_pos)
        return self

    def rm_aug_line(self, pos: int | None = None) -> Matrix:
        """Remove an augmentation line from the matrix.

        Removes the specified position from the matrix's augmentation line tracking.
        If the matrix has no augmentation lines tracked or the position is not
        currently marked as an augmentation line, the matrix remains unchanged.

        Note:
            This method modifies the matrix's internal `_aug_pos` attribute which
            tracks augmentation line positions. If no `_aug_pos` attribute exists,
            it will be initialized as an empty set.

        Args:
            pos (int, optional): The column position of the augmentation line to remove.
                If it is not set, all augmentation lines will be removed

        Returns:
            (Matrix): The matrix instance (supports method chaining).

        Examples:
            >>> mat = Matrix([[1, 2], [3, 4]], aug_pos=0)
            >>> mat.rm_aug_line(0)
            Matrix([
            [1, 2]
            [3, 4]
            ])
        """
        if not hasattr(self, "_aug_pos") or pos is None:
            self._aug_pos = set()
            return self
        if pos in self._aug_pos:
            self._aug_pos.remove(pos)
        return self

    # Override
    def row_join(self, other: Matrix, aug_line: bool = True) -> Matrix:
        """Joins two matrices horizontally (column-wise), preserving augmentation lines.

        This method concatenates the columns of `self` and `other` to form a new matrix. Any augmentation lines
        (vertical lines for augmented matrices) tracked in either matrix are preserved and adjusted for the new column positions.

        Args:
            other (Matrix): The matrix to join to the right of `self`.
            aug_line (bool, optional): If `True`, adds an augmentation line between the two matrices.

        Returns:
            (Matrix): A new matrix formed by joining `self` and `other` column-wise, with updated augmentation lines.

        Raises:
            sympy.matrices.exceptions.ShapeError: If the number of rows in `self` and `other` do not match.

        Examples:
            >>> A = Matrix([[1, 2], [3, 4]])
            >>> B = Matrix([[5], [6]])
            >>> A.row_join(B)
            Matrix([
            [1, 2 | 5]
            [3, 4 | 6]
            ])
        """
        if not hasattr(self, "_aug_pos"):
            self._aug_pos = set()
        if not hasattr(other, "_aug_pos"):
            other._aug_pos = set()

        offset = self.cols
        # Iterate over a copy of the set to prevent runtime errors
        for pos in list(other._aug_pos):
            self._aug_pos.add(pos + offset)
        if aug_line:
            self._aug_pos.add(offset - 1)
        return Matrix(super().row_join(other), aug_pos=self._aug_pos)

    # Override
    def col_join(self, other: Matrix) -> Matrix:
        """
        Joins two matrices vertically (row-wise), preserving augmentation lines.

        This method concatenates the rows of `self` and `other` to form a new matrix. Only augmentation lines
        (vertical lines for augmented matrices) that are present in both matrices at the same column positions are preserved.

        Note:
            - Both `self` and `other` matrices should have the same number of columns for the join to be valid.
            - The method updates the `_aug_pos` attribute to include only those positions that are common
              in both matrices, ensuring that the augmentation lines are correctly aligned after the join.

        Args:
            other (Matrix): The matrix to join below `self`.

        Returns:
            (Matrix): A new matrix formed by joining `self` and `other` row-wise, with preserved augmentation lines.

        Raises:
            sympy.matrices.exceptions.ShapeError: If the number of columns in `self` and `other` do not match.

        Examples:
            >>> A = Matrix([[1, 2]], aug_pos={0})
            >>> B = Matrix([[3, 4]], aug_pos={0})
            >>> A.col_join(B)
            Matrix([
            [1 | 2]
            [3 | 4]
            ])
        """
        if not hasattr(self, "_aug_pos"):
            self._aug_pos = set()
        if not hasattr(other, "_aug_pos"):
            other._aug_pos = set()

        aug = self._aug_pos & other._aug_pos
        return Matrix(super().col_join(other), aug_pos=aug)

    def scale_row(
        self, idx: int, scalar: Expr | float | int, verbosity: int = 2
    ) -> Matrix:
        """
        Scales a row of the matrix by a scalar and simplifies the result.

        This method scales a specified row of the matrix by multiplying it with a scalar
        and then simplifies the matrix. The result is stored back in the matrix. Optionally,
        the method can print information about the row scaling and display the matrix,
        depending on the verbosity level.

        Note:
            - The method modifies the matrix in-place and returns the updated matrix.
            - After scaling the row, the matrix is simplified using [`simplify`][..].

        Args:
            idx (int): The index of the row to scale (0-based).
            scalar (Expr, float, int): The scalar by which to multiply the row.
            verbosity (int, optional): The level of verbosity for output.

                - 0: No output.
                - 1: Print the row scaling operation.
                - 2: Print the row scaling operation and display the matrix.

        Returns:
            (Matrix): The modified matrix with the scaled row.

        Raises:
            IndexError: If the `idx` is out of range for the number of rows in the matrix.

        Examples:
            >>> mat = Matrix([[1, 2], [3, 4]])
            >>> mat.scale_row(0, 2, verbosity=0)
            Matrix([
            [2, 4]
            [3, 4]
            ])
        """

        if scalar == 0:
            warn("Matrix rows should not be scaled by 0", UserWarning)

        scalar = sym.sympify(scalar)
        self[idx, :] *= scalar  # type: ignore
        self.simplify(suppress_warnings=True)

        if verbosity >= 1:
            display(
                f"\\left({sym.latex(scalar)}\\right) R_{idx + 1} \\rightarrow R_{idx + 1}",
                opt="math",
            )
        if verbosity >= 2:
            display(self)
            print("\n")

        return self

    def swap_row(self, idx_1: int, idx_2: int, verbosity: int = 2) -> Matrix:
        """Swaps two rows of the matrix.

        This method swaps the contents of two rows in the matrix. The operation is performed
        in-place, and the modified matrix is returned. Optionally, the method can print
        information about the row swap and display the matrix, depending on the verbosity level.

        Note:
            - The method modifies the matrix in-place and returns the updated matrix.
            - After performing the row swaps, the matrix is simplified using [`simplify`][..].

        Args:
            idx_1 (int): The index of the first row to swap (0-based).
            idx_2 (int): The index of the second row to swap (0-based).
            verbosity (int, optional): The level of verbosity for output.

                - 0: No output.
                - 1: Print the row swap operation.
                - 2: Print the row swap operation and display the matrix.

        Returns:
            (Matrix): The modified matrix after the row swap.

        Raises:
            IndexError: If the `idx_1` or `idx_2` is out of range for the number of rows in the matrix.

        Examples:
            >>> mat = Matrix([[1, 2], [3, 4]])
            >>> mat.swap_row(0, 1, verbosity=0)
            Matrix([
            [3, 4]
            [1, 2]
            ])
        """

        self[idx_1, :], self[idx_2, :] = self[idx_2, :], self[idx_1, :]

        if verbosity >= 1:
            display(f"R_{idx_1 + 1} \\leftrightarrow R_{idx_2 + 1}", opt="math")
        if verbosity >= 2:
            display(self)
            print("\n")

        return self

    def reduce_row(
        self, idx_1: int, scalar: Expr | float | int, idx_2: int, verbosity: int = 2
    ) -> Matrix:
        """Reduces a row by subtracting a scalar multiple of another row.

        This method modifies a row by subtracting a specified scalar multiple of another row.
        The result is stored back in the matrix. Optionally, the method can print information
        about the row reduction and display the matrix, depending on the verbosity level.

        Note:
            - The method modifies the matrix in-place and returns the updated matrix.
            - After performing the row reduction, the matrix is simplified using [`simplify`][..].

        Args:
            idx_1 (int): The index of the row to reduce (0-based).
            scalar (Expr, float, int): The scalar by which to multiply the second row.
            idx_2 (int): The index of the row from which to subtract the scalar multiple (0-based).
            verbosity (int, optional): The level of verbosity for output.

                - 0: No output.
                - 1: Print the row reduction operation.
                - 2: Print the row reduction operation and display the matrix.

        Returns:
            (Matrix): The modified matrix after the row reduction.

        Raises:
            IndexError: If the `idx_1` or `idx_2` is out of range for the number of rows in the matrix.

        Examples:
            >>> mat = Matrix([[1, 2], [3, 4]])
            >>> mat.reduce_row(0, 2, 1, verbosity=0)
            Matrix([
            [-5, -6]
            [ 3,  4]
            ])
        """

        scalar = sym.sympify(scalar)
        self[idx_1, :] = self[idx_1, :] - scalar * self[idx_2, :]  # type: ignore
        self.simplify(suppress_warnings=True)

        if verbosity >= 1:
            display(
                f"R_{idx_1 + 1} - \\left({sym.latex(scalar)}\\right)R_{idx_2 + 1} \\rightarrow R_{idx_1 + 1}",
                opt="math",
            )
        if verbosity >= 2:
            display(self)
            print("\n")

        return self

    def get_pivot_row(
        self, col_idx: int, row_start_idx: int, follow_GE: bool = False
    ) -> int | None:
        """Finds the row index of the pivot element in a given column.

        This method attempts to find a row that contains a non-zero element in the
        specified column. If the `follow_GE` flag is `False`, it first looks for
        a non-zero constant that does not contain any symbolic expressions. If no
        such element is found, it will return the first non-zero element. If the
        entire column contains only zeros, the method returns -1.

        Args:
            col_idx (int): The index of the column to search for the pivot.
            row_start_idx (int): The row index to start searching from.
            follow_GE (bool, optional): Flag to control whether to follow Gaussian elimination strategy.

                - `True`: Always return the first non-zero element, even if it is symbolic.
                - `False`: First look for non-zero constants that are not symbolic expressions.

        Returns:
            (int): The index of the row containing the pivot element, or None if no pivot is found.

        Examples:
            >>> mat = Matrix([[1, 2, 3], [4, 5, 6], [0, 0, 0]])
            >>> mat.get_pivot_row(0, 0)
            0
        """

        # Step 1: Search for a non-zero constant that is not symbolic (if not following Gaussian elimination)
        # that it is easier to reduce other rows
        if not follow_GE:
            for row_idx in range(row_start_idx, self.rows):
                term = self[row_idx, col_idx]
                if term != 0:
                    # Check if it's not a symbolic expression
                    if not isinstance(term, Expr):
                        return row_idx
                    # Check if it's a non-symbolic constant
                    elif len(term.free_symbols) == 0:
                        return row_idx

        # Step 2: If no non-zero constant is found, return the first non-zero element (symbolic or not)
        for row_idx in range(row_start_idx, self.rows):
            term = self[row_idx, col_idx]
            if term != 0:
                return row_idx

        # Step 3: If no non-zero element is found, return None (indicating no pivot)
        return None

    def get_pivot_pos(self) -> list[tuple[int, int]]:
        """Finds the positions of the pivot elements in the matrix.

        This method checks the matrix to determine the positions of the pivots
        (the first non-zero entry in each row) by examining each column one-by-one.
        It assumes that the matrix is in Row Echelon Form (REF), as checked by the
        [`is_echelon`][sympy.matrices.matrixbase.MatrixBase.is_echelon] property.

        It uses [`get_pivot_row`][..] to find the pivot row for each column.
        For each pivot found, a tuple (row, column) is added to the result list.

        Returns:
            (list[tuple[int, int]]): A list of lists, where each sublist contains a
                tuple representing the position (row, column) of a pivot.

        Examples:
            >>> mat = Matrix([[1, 2, 3], [0, 0, 5], [0, 0, 0]])
            >>> mat.get_pivot_pos()
            [(0, 0), (1, 2)]
        """

        assert self.is_echelon  # check for REF

        pivot_pos: list[tuple[int, int]] = []
        cur_row_pos = 0
        for cur_col_pos in range(self.cols):
            pivot_row = self.get_pivot_row(cur_col_pos, cur_row_pos, follow_GE=False)

            if pivot_row is not None:
                pivot_pos.append((pivot_row, cur_col_pos))
                cur_row_pos += 1

        return pivot_pos

    def get_pivot_elements(self) -> list[Expr]:
        """Retrieves the pivot elements from the matrix.

        This method identifies the pivot positions (row, column) using the
        [`get_pivot_pos`][..] method and
        then extracts the elements at those positions in the matrix.

        Returns:
            (list[Expr]): A list of pivot elements corresponding
                to the positions identified by [`get_pivot_pos`][..].

        Examples:
            >>> mat = Matrix([[1, 2, 3], [0, 0, 5], [0, 0, 0]])
            >>> mat.get_pivot_elements()
            [1, 5]
        """

        pivot_elements: list[Expr] = []

        for i, j in self.get_pivot_pos():
            pivot_elements.append(self[i, j])  # type: ignore

        return pivot_elements

    def ref(
        self,
        verbosity: int = 2,
        max_tries: int = 2,
        follow_GE: bool = False,
    ) -> PLU:
        """Find the Row Echelon Form (REF) of the matrix.

        This method applies Gaussian elimination (or a similar approach) to bring
        the matrix to row echelon form.

        Note:
            - PLU decomposition is the generalisation of the LU decomposition. Unlike
                LU decomposition, PLU works for any matrix.
            - The REF is obtained from `PLU.U` and the matrix is LU factorisable if `PLU.P`
                is the identity matrix.
            -

        Args:
            verbosity (int, optional): Level of verbosity for the output.

                - 0: No output.
                - 1: Output basic information (e.g., row operations).
                - 2: Output detailed information (e.g., matrix states after each operation).
            max_tries (int, optional): Maximum number of tries to reduce a row in case of symbolic denominators.
            follow_GE (bool, optional): Whether to strictly follow Gaussian elimination rules.

                - `True`: Always return the first non-zero element, even if it is symbolic.
                - `False`: First look for non-zero constants that are not symbolic expressions.

        Returns:
            (PLU): A dataclass containing the following matrices:

                - `P` ([`Matrix`][...]): The permutation matrix used during the transformation.
                - `L` ([`Matrix`][...]): The lower triangular matrix representing the multipliers used in the elimination process.
                - `U` ([`Matrix`][...]): The upper triangular matrix in row echelon form.

        Examples:
            >>> mat = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
            >>> mat.ref(verbosity=0)
            PLU(P=Matrix([
            [1, 0, 0]
            [0, 1, 0]
            [0, 0, 1]
            ]), L=Matrix([
            [1, 0, 0]
            [4, 1, 0]
            [7, 2, 1]
            ]), U=Matrix([
            [1,  2,  3]
            [0, -3, -6]
            [0,  0,  0]
            ]))

        See Also:
            - [`Matrix.LUdecomposition`][sympy.matrices.matrixbase.MatrixBase.LUdecomposition]
        """

        U = self.copy()

        I = self.elem()
        L = self.elem()
        P = self.elem()

        # Loop over each column
        cur_row_pos = 0

        for cur_col_pos in range(self.cols):
            # Find the first non-zero row in the current column
            pivot_row = U.get_pivot_row(cur_col_pos, cur_row_pos, follow_GE)

            if pivot_row is None:
                # If no non-zero pivot is found, continue to the next column
                continue

            # Swap the current row with the pivot row if necessary
            if pivot_row != cur_row_pos:
                U.swap_row(cur_row_pos, pivot_row, verbosity=verbosity)
                P_elem = I.copy().swap_row(cur_row_pos, pivot_row, verbosity=0)
                P = P @ P_elem
                L = P_elem @ L @ P_elem

            # Eliminate the current column in rest of the rows below
            for row_idx in range(cur_row_pos + 1, self.rows):
                # reduce the row_idx iteratively via partial fractions to
                # prevent division by a possible 0 term
                tries = 0
                while U[row_idx, cur_col_pos] != 0:
                    tries += 1
                    if tries > max_tries:
                        warn(
                            f"ERROR: Max tries exceeded to reduce row {row_idx + 1} with row {cur_row_pos + 1}",
                            UserWarning,
                        )
                        break
                    try:
                        scalar = U[row_idx, cur_col_pos] / U[cur_row_pos, cur_col_pos]  # type: ignore
                        scalar = scalar.expand().simplify()

                        try:
                            decomp = sym.apart(scalar)  # partial fractions
                        except Exception as e:
                            decomp = scalar
                        if isinstance(decomp, sym.Add):
                            terms = decomp.args
                        else:
                            # there is only 1 term (could be integer or proper fraction)
                            terms = [decomp]

                        for term in terms:
                            _, d = sym.fraction(term)

                            # ensure denominator is non-zero so that reduction is valid
                            if not _is_zero(d):
                                U.reduce_row(
                                    row_idx,
                                    term,
                                    cur_row_pos,
                                    verbosity=verbosity,  # type: ignore
                                )
                                elem = I.copy().reduce_row(
                                    row_idx, -term, cur_row_pos, verbosity=0
                                )  # type: ignore
                                L = L @ elem

                        # Cases where pivot row contains symbols such that scalar is a
                        # fraction with symbolic denominator.
                        # To reduce further, can only scale row_idx accordingly
                        if U[row_idx, cur_col_pos] != 0:
                            scalar = (
                                U[cur_row_pos, cur_col_pos] / U[row_idx, cur_col_pos]
                            )  # type: ignore
                            tmp = scalar.simplify()
                            if tmp is not None:
                                scalar = tmp
                            n, d = sym.fraction(scalar)
                            # to scale by n, n cannot be 0 both numerically or symbolically
                            # to scale by 1/d, d cannot be 0, same argument as n
                            # if (n != 0) and (not _is_zero(d)):
                            if (not _is_zero(n)) and (not _is_zero(d)):
                                U.scale_row(row_idx, scalar, verbosity=verbosity)
                                elem = I.copy().scale_row(
                                    row_idx, 1 / scalar, verbosity=0
                                )
                                L = L @ elem

                    except Exception as error:
                        print(f"Exception encountered: {error}")
                        return PLU(P, L, U)
            cur_row_pos += 1

        # Return the appropriate number of matrices based on the `matrices` parameter
        return PLU(P, L, U)

    def find_all_cases(self) -> list:
        cases = []
        det = (self.T @ self).det()
        if len(det.free_symbols) == 0:
            # determinant of the matrix is fixed
            return []
        elif len(det.free_symbols) == 1:
            for sol in sym.solve(det):
                cases.extend([{det.free_symbols.pop(): sol}])
        else:
            for sol in sym.solve(det):
                cases.extend([sol])
        cases = [dict(case) for case in set(tuple(case.items()) for case in cases)]

        # if variable is not in dictionary, it can be treated as entire real
        # except for specific cases found in other combinations
        combinations = set()
        for subset in _powerset(cases):
            combined = dict()
            for d in subset:
                combined = {**combined, **d}
            combinations.add(tuple(sym.ordered(combined.items())))

        return list(sym.ordered(dict(combination) for combination in combinations))

    def evaluate_cases(self, rhs: Matrix) -> None:
        cases = self.find_all_cases()
        all_possible_values = set(
            possible_val for case in cases for possible_val in case.items()
        )

        for i, case in enumerate(cases, 1):
            print(
                f"Case {i}: {case}, not including {[dict([val]) for val in all_possible_values.symmetric_difference(set(case.items()))]}"
            )
            U = (
                self.row_join(rhs, aug_line=False)
                .subs(case)
                .ref(verbosity=0, follow_GE=False)
                .U
            )
            display(U)

    # Override
    def rref(self, *args, pivots: bool = True, **kwargs) -> RREF | Matrix:
        """Computes the Reduced Row Echelon Form (RREF) of the matrix.

        This method is a wrapper for [`rref`][sympy.matrices.matrixbase.MatrixBase.rref] method
        and returns the matrix in Reduced Row Echelon Form (RREF) along with the pivot positions.

        Args:
            pivots (bool, optional):
                If `True`, returns a tuple containing the RREF matrix and a list of pivot columns
            *args: Positional arguments passed to SymPy's [`rref`][sympy.matrices.matrixbase.MatrixBase.rref] method.
            **kwargs: Keyword arguments passed to SymPy's [`rref`][sympy.matrices.matrixbase.MatrixBase.rref] method.

        Returns:
            (RREF): A dataclass containing the following:

                - `rref` ([`Matrix`][...]): The matrix in Reduced Row Echelon Form.
                - `pivots` ([`list`][list][[`int`][int]]): A list of integers representing the indices of the pivot columns.
            (Matrix): If `pivots` is `False`, only the `rref` will be returned.

        Examples:
            >>> mat = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
            >>> mat.rref()
            RREF(rref=Matrix([
            [1, 0, -1]
            [0, 1,  2]
            [0, 0,  0]
            ]), pivots=(0, 1))
        """
        if pivots:
            rref_mat, pivot_pos = super().rref(*args, **kwargs)
        else:
            rref_mat = super().rref(*args, pivots=False, **kwargs)

        aug = self._aug_pos.copy() if hasattr(self, "_aug_pos") else set()
        rref_mat = Matrix(rref_mat, aug_pos=aug)

        if pivots:
            return RREF(rref_mat, pivot_pos)
        else:
            return rref_mat

    # Override
    def solve(self, rhs: Matrix) -> Matrix | list[Matrix]:
        """Solves the linear system `Ax = rhs` for `x`.

        This method uses SymPy's [`solve`][sympy.solvers.solvers.solve] method to find a solution vector `x` such that `self @ x = rhs`.
        If a unique solution exists, it is returned. If multiple solutions exist, a list of solution matrices is returned.
        If no solution exists (e.g., the system is inconsistent), a [ValueError][] is raised.

        Args:
            rhs (Matrix): The right-hand side matrix or vector in the equation `Ax = rhs`.

        Returns:
            (Matrix): The solution vector or matrix `x` that satisfies `Ax = rhs`.
            (list[Matrix]): If multiple solutions exist, a list of solution matrices is returned.

        Raises:
            ValueError: If no solution is found for the linear system, an exception is raised.

        Examples:
            >>> A = Matrix([[1, 2], [3, 4]])
            >>> b = Matrix([[5], [11]])
            >>> A.solve(b)
            Matrix([
            [1],
            [2]])

        See Also:
            - [`sympy.solve`][sympy.solvers.solvers.solve]: For solving equations in general.
            - [`rref`][..]: For finding the reduced row echelon form of the matrix.
            - [`solve_least_squares`][..]: For solving least squares problems.
        """
        # Use sympy's solve function directly
        x = Matrix.create_unk_matrix(r=self.cols, c=1)
        solution = sym.solve(self @ x - rhs, x.free_symbols)

        if isinstance(solution, dict):
            # If a single solution is found, substitute the values into the unknown matrix
            return x.subs(solution)
        elif isinstance(solution, list) and len(solution) > 0:
            # If multiple solutions are found, return all solutions found
            return [x.subs(sol) for sol in solution]
        else:
            # If no solution is found (e.g., inconsistent system or empty list), raise an error
            display(self.row_join(rhs, aug_line=True).rref())
            raise ValueError(
                "No solution found for the linear system. The system may be inconsistent."
            )

    #############################
    # CHAPTER 2: MATRIX ALGEBRA #
    #############################

    def inverse(
        self,
        option: Literal["left", "right", "both"] | None = None,
        matrices: int = 1,
        verbosity: int = 0,
    ) -> Matrix | PartGen | None:
        """Computes the left or right inverse of a matrix, depending on its rank and the specified option.

        The method checks whether the matrix has full row rank or full column rank and computes either:
            - The **left inverse** (if the matrix has full column rank).
            - The **right inverse** (if the matrix has full row rank).

        If neither option is provided, the method automatically determines which inverse to compute based on the matrix's rank.

        Args:
            option (str, optional): Specifies which inverse to compute:

                - `'left'` for the left inverse (requires the matrix to have full column rank).
                - `'right'` for the right inverse (requires the matrix to have full row rank).
                - `'both'` for the inverse of a square matrix (works on both sides).

            matrices (int, optional): Specifies the number of matrices to return:

                - 1: Returns only the inverse matrix.
                - 2: Returns the particular and general solutions of the inverse.

            verbosity (int, optional): Level of verbosity for displaying intermediate steps:

                - 0: No output.
                - 1: Display matrices before and after RREF.

        Returns:
            (Matrix): If `matrices = 1`, returns the inverse matrix.
            (PartGen): If `matrices = 2`, returns a dataclass containing the particular and general solutions of the inverse.

        Raises:
            ValueError: If no valid inverse (left or right or both) is found, an exception is raised.

        Examples:
            >>> mat = Matrix([[1, 2], [3, 4]])
            >>> mat.inverse()
            Matrix([
            [ -2,    1],
            [3/2, -1/2]])
        """

        if option is None:
            rank = self.rank()
            if rank == self.cols:
                if verbosity:
                    print("Left inverse found!")
                option = "left"
            if rank == self.rows:
                if verbosity:
                    print("Right inverse found!")
                option = "right"
                if self.rows == self.cols:
                    # square matrix inverse works on both sides
                    option = "both"
            else:
                raise ValueError(
                    f"No inverse found! Rank: {rank}, Rows: {self.rows}, Columns: {self.cols}. Try pseudo-inverse: .pinv()"
                )

        if option == "both" and self.rows != self.cols:
            raise ValueError(
                "Cannot compute both left and right inverse for non-square matrices!"
            )

        if (option is not None) and (verbosity >= 1):
            if option == "left":
                aug = self.T.copy().row_join(
                    Matrix.eye(self.cols, aug_pos=range(self.cols))
                )
                print("Before RREF: [self^T | eye]")
                display(aug)
                print("\nAfter RREF:")
                display(aug.rref())
            else:
                aug = self.copy().row_join(
                    Matrix.eye(self.rows, aug_pos=range(self.rows))
                )
                print("Before RREF: [self | eye]")
                display(aug)
                print("\nAfter RREF:")
                display(aug.rref())

        if option is not None:
            X = Matrix.create_unk_matrix(r=self.cols, c=self.rows, symbol="x")
            if option == "left":
                eqn = X @ self - sym.eye(self.cols)
            else:
                eqn = self @ X - sym.eye(self.rows)

            sol = sym.solve(eqn, X.free_symbols)
            if isinstance(sol, list) and len(sol) > 0:
                # Multiple sets of solutions found, picks the first 1
                X = X.subs(sol[0])
            elif isinstance(sol, dict):
                X = X.subs(sol)
            else:
                raise ValueError(
                    f"No {option} inverse found! Try pseudo-inverse: .pinv()"
                )

            if matrices == 1:
                return X
            elif matrices == 2:
                return X.sep_part_gen()

    def elem(self) -> Matrix:
        """Returns the identity matrix with the same number of rows as the current matrix.

        This method creates an identity matrix to be used for elementary row operations,
        i.e., `A = I A` where `I` is the identity matrix.

        Returns:
            (Matrix): An identity matrix with the same number of rows as the current matrix.

        Examples:
            >>> mat = Matrix([[1, 2], [3, 4], [5, 6]])
            >>> mat.elem()
            Matrix([
            [1, 0, 0]
            [0, 1, 0]
            [0, 0, 1]
            ])
        """
        return Matrix.eye(self.rows)

    # override
    def adjoint(self) -> Matrix:
        """Computes the adjugate (classical adjoint) of the matrix.

        This method calculates the classical adjoint (also known as the [adjugate](https://en.wikipedia.org/wiki/Adjugate_matrix)
        in literature) of the matrix. The adjoint of a matrix (as defined in MA1522 syllabus) is the transpose of its cofactor matrix.

        Note:
            If you wish to compute the conjugate transpose of the matrix (SymPy's definition for adjoint),
            use `self.H` directly or `super(symbolic.Matrix, self).adjoint()` to call the parent method.

        Returns:
            (Matrix): The classical adjoint (or adjugate) matrix of the current matrix.

        Examples:
            >>> mat = Matrix([[1, 2], [3, 4]])
            >>> mat.adjoint()
            Matrix([
            [ 4, -2],
            [-3,  1]])

        See Also:
            - SymPy's [`Matrix.adjugate`][sympy.matrices.matrixbase.MatrixBase.adjugate]
            - SymPy's [`Matrix.adjoint`][sympy.matrices.matrixbase.MatrixBase.adjoint] for conjugate transpose.
        """
        warn(
            """The classical adjoint of the matrix is computed rather than the conjugate transpose.
            Please use self.adj() instead to remove ambiguity.""",
            DeprecationWarning,
            stacklevel=2,
        )
        return self.adjugate()

    def adj(
        self,
        method: Literal["bareiss", "berkowitz", "bird", "laplace", "lu"] = "berkowitz",
    ) -> Matrix:
        """Alias for the [`adjoint`][..] method.

        It returns the classical adjoint (or [adjugate](https://en.wikipedia.org/wiki/Adjugate_matrix)) of the matrix.

        Args:
            method (str, optional): Method to use to find the cofactors, can be "bareiss", "berkowitz",
                "bird", "laplace" or "lu".
        Returns:
            (Matrix): The classical adjoint of the current matrix.

        Examples:
            >>> mat = Matrix([[1, 2], [3, 4]])
            >>> mat.adj()
            Matrix([
            [ 4, -2],
            [-3,  1]])

        See Also:
            - SymPy's [`adjugate`][sympy.matrices.matrixbase.MatrixBase.adjugate]
        """
        return self.adjugate(method=method)

    # override
    def cramer_solve(
        self, rhs: Matrix, det_method: str = "laplace", verbosity: int = 2
    ) -> Matrix:
        """Solves the linear system using Cramer's Rule.

        This method applies Cramer's Rule to solve the linear system represented by the matrix and the right-hand side vector.
        It computes the determinant of the matrix and uses it to find the solution vector.

        Args:
            rhs (Matrix): The right-hand side vector in the equation `Ax = rhs`.
            det_method (str, optional): The method to use for computing the determinant. Options include:

                - `'laplace'`: Uses the Laplace expansion method.
                - `'berkowitz'`: Uses the Berkowitz algorithm.
                - `'bird'`: Uses the Bird's algorithm.
                - `'bareiss'`: Uses the Bareiss algorithm.
                - `'lu'`: Uses LU decomposition.

            verbosity (int, optional): Level of verbosity for displaying intermediate steps:

                - 0: No output.
                - 1: Display basic information.
                - 2: Display detailed information.

        Returns:
            (Matrix): The solution vector `x` that satisfies `self @ x = rhs`.

        Raises:
            sympy.matrices.exceptions.NonSquareMatrixError: If the matrix is not square.
            sympy.matrices.exceptions.ShapeError: If the matrix and the right-hand side vector have incompatible dimensions.
            ValueError: If the determinant is zero, indicating that the system has no unique solution.

        Examples:
            >>> A = Matrix([[1, 2], [3, 4]])
            >>> b = Matrix([[5], [11]])
            >>> A.cramer_solve(b, verbosity=0)
            Matrix([
            [1],
            [2]])

        See Also:
            - [`solve`][..]: For solving linear systems using other methods.
            - SymPy's [`Matrix.det`][sympy.matrices.matrixbase.MatrixBase.det] for computing the determinant.
            - SymPy's [`Matrix.cramer_solve`][sympy.matrices.matrixbase.MatrixBase.cramer_solve]
        """
        if self.rows != rhs.rows:
            raise sym.ShapeError(
                "The right-hand side vector must have the same number of rows as the matrix."
            )
        if rhs.cols != 1:
            raise sym.ShapeError(
                "The right-hand side vector must be a column vector (1 column)."
            )
        if self.rows != self.cols:
            raise sym.NonSquareMatrixError(
                "Cramer's Rule can only be applied to square matrices."
            )
        det = self.det(method=det_method)
        if det == 0:
            raise ValueError("Determinant is zero, no unique solution exists.")

        entries = []
        for i in range(self.cols):
            # Create a copy of the matrix and replace the i-th column with the rhs vector
            modified_matrix = self.copy()
            modified_matrix[:, i] = rhs[:, 0]
            if verbosity >= 2:
                print(f"Modified matrix for column {i + 1}:")
                display(modified_matrix)
            det_i = modified_matrix.det(method=det_method) / det
            if verbosity >= 1:
                display(
                    "\\text{Determinant for column }"
                    + str(i + 1)
                    + ": "
                    + sym.latex(det_i),
                    opt="math",
                )
            entries.append(det_i)
        return Matrix(entries)

    def column_constraints(self, use_ref: bool = False, verbosity: int = 1) -> Matrix:
        r"""Computes the column constraints for the matrix by appending a symbolic vector.

        This method creates a matrix where a random column vector $\begin{pmatrix} x_1 \\ \vdots \\ x_m \end{pmatrix}$
        is added to the matrix as an additional column. It then constructs a larger augmented matrix
        and optionally computes its Row Echelon Form (REF) or Reduced Row Echelon Form (RREF).

        The method modifies the matrix to ensure that the unknown vector is not reduced in RREF,
        and the constraints for the matrix columns are calculated accordingly.

        Args:
            use_ref (bool, optional): Whether to use Row Echelon Form (REF) instead of Reduced Row Echelon Form (RREF).
                If `False`, RREF will be used.
            verbosity (int, optional): Verbosity level for displaying information.

                - 0: No output
                - 1: Display all information

        Returns:
            (Matrix): A new matrix containing the result after applying REF or RREF to the augmented matrix.

        Examples:
            >>> mat = Matrix([[1, 2], [2, 4]]) # linearly dependent columns
            >>> mat.column_constraints(verbosity=0)
            Matrix([
            [1, 2 |       x_2/2]
            [0, 0 | x_1 - x_2/2]
            ])

            >>> mat = Matrix([[1, 2], [3, 4]]) # linearly independent columns
            >>> mat.column_constraints(verbosity=0)
            Matrix([
            [1, 0 |    -2*x_1 + x_2]
            [0, 1 | 3*x_1/2 - x_2/2]
            ])
        """

        # write a random vector as x_1, ..., x_m, given m rows
        vector = Matrix.create_unk_matrix(self.rows, 1, "x")

        # insert hidden column vectors so that the unknown vector is not reduced in rref
        hidden = self.elem()

        M = self.copy().row_join(hidden).row_join(vector)
        if use_ref:
            res = M.ref().U
        else:
            res = M.rref(pivots=False)

        visible_cols = (*range(self.cols), -1)
        res_matrix = res.select_cols(*visible_cols).aug_line(-2)  # type: ignore

        if verbosity:
            print("Before RREF: [self | vec]")
            display(M.select_cols(*visible_cols).aug_line(-2))
            print("After RREF")
            display(res_matrix)
            print(
                "For the system to be consistent, the following constraints must be satisfied."
            )
            for i in range(res_matrix.rows):
                # check for zero row
                if res_matrix[i, : self.cols].norm() == 0:  # type: ignore
                    display(sym.Eq(res_matrix[i, -1], 0))
        return res_matrix

    ######################################
    # CHAPTER 3: EUCLIDEAN VECTOR SPACES #
    ######################################

    def normalized(
        self,
        iszerofunc: Callable[[Expr], int] | None = None,
        factor: bool = False,
    ) -> Matrix | ScalarFactor:
        """Normalizes the column vectors of the matrix (scaling each vector to have a unit norm).

        Args:
            iszerofunc (Callable[[Expr], int], optional): Function to determine if a value should be treated as zero
            factor (bool, optional): If `True`, returns the [`ScalarFactor`][(p).] dataclass.
                If `False`, only the normalized matrix is returned.

        Returns:
            (ScalarFactor): A dataclass containing the normalized matrix and the diagonal matrix of norms if `factor=True`.
            (Matrix): The normalized matrix if `factor=False`.

        Examples:
            >>> mat = Matrix([[0, 3], [0, 4]])
            >>> mat.normalized()
            Matrix([
            [0, 3/5]
            [0, 4/5]
            ])

            >>> mat = Matrix([[3], [4]])
            >>> mat.normalized(factor=True)
            ScalarFactor(diag=Matrix([
            [0,   0]
            [0, 1/5]
            ]), full=Matrix([
            [0, 3]
            [0, 4]
            ]), order='FD')
        """

        for i in range(self.cols):
            col = self.col(i)
            scalar = sym.sqrt(sum(x**2 for x in col))  # Manual norm calculation
            if iszerofunc is None:
                if scalar != 0:
                    self[:, i] /= scalar  # type: ignore
            else:
                if iszerofunc(scalar) != 0:
                    self[:, i] /= scalar  # type: ignore

        if factor:
            return self.scalar_factor(column=True)
        else:
            return self

    def is_linearly_independent(
        self,
        colspace: bool = True,
        verbosity: int = 0,
    ) -> bool:
        """Determines if the vectors in the matrix are linearly independent.

        This method checks whether the columns (or rows) of the matrix are linearly independent
        by computing its reduced row echelon form (RREF) and comparing the number of pivot columns
        to the number of columns (or rows).

        Args:
            colspace (bool, optional): If `True`, checks linear independence of columns (column space).
                If `False`, checks linear independence of rows (row space).

            verbosity (int, optional): Level of output during the RREF and check.

                - 0: No output (default).
                - 1: Print summary of the RREF and check.
                - 2: Print the matrix before and after RREF, and the check details.

        Returns:
            (bool): `True` if the matrix's columns (or rows) are linearly independent, `False` otherwise.

        Examples:
            >>> mat = Matrix([[1, 2], [3, 4]])
            >>> mat.is_linearly_independent()
            True

            >>> mat = Matrix([[1, 2], [2, 4]])
            >>> mat.is_linearly_independent()
            False
        """
        rref = self.rref(pivots=True)
        assert isinstance(rref, RREF), "RREF should return a RREF dataclass"
        rref_mat, pivots = rref.rref, rref.pivots

        if verbosity == 1:
            print("rref(self)")
        elif verbosity >= 2:
            print("Before RREF: self")
            display(self)
            print("\nAfter RREF:")
            display(rref_mat)

        if colspace:
            if verbosity >= 1:
                print(
                    f"Check if Number of columns ({self.cols}) == Number of pivot columns ({len(pivots)})"  # type: ignore
                )
            return self.cols == len(pivots)
        else:
            if verbosity >= 1:
                print(
                    f"Check if Number of rows ({self.rows}) == Number of pivot columns ({len(pivots)})"  # type: ignore
                )
            return self.rows == len(pivots)

    def get_linearly_independent_vectors(
        self, colspace: bool = True, verbosity: int = 2
    ) -> Matrix:
        """Returns a matrix containing the linearly independent vectors from the column space or row space.

        This method computes the reduced row echelon form (RREF) of the matrix and selects the non-zero rows
        as linearly independent vectors. The result is a matrix whose columns (or rows) are linearly independent.

        Args:
            colspace (bool, optional): If `True`, returns linearly independent vectors from the column space.
                If `False`, returns from the row space.

            verbosity (int, optional): Level of output verbosity.

                - 0: No output.
                - 1: Print a summary of the RREF and selection.

        Returns:
            (Matrix): A matrix whose columns (if colspace=True) or rows (if colspace=False) are linearly independent vectors.

        Examples:
            >>> mat = Matrix([[1, 2], [2, 4]])
            >>> mat.get_linearly_independent_vectors(colspace=True, verbosity=0)
            Matrix([
            [1]
            [2]
            ])
        """
        if colspace:
            rref = self.rref(pivots=True)
            assert isinstance(rref, RREF), "RREF should return a RREF dataclass"
            if verbosity >= 1:
                print("Before RREF: [self]")
                display(self)
                print("\nAfter RREF:")
                display(rref)
                print("Select columns of self corresponding to pivot positions.")
            return self.select_cols(*rref.pivots)
        else:
            rref = self.T.rref(pivots=True)
            assert isinstance(rref, RREF), "RREF should return a RREF dataclass"
            if verbosity >= 1:
                print("Before RREF: [self^T]")
                display(self.T)
                print("\nAfter RREF:")
                display(rref)
                print("Select rows of self corresponding to pivot positions.")
            return self.select_rows(*rref.pivots)

    def simplify_basis(self, colspace: bool = True, verbosity: int = 2) -> Matrix:
        """Returns a simplified basis for the column space or row space of the matrix.

        This method computes a basis for either the column space or the row space of the matrix
        by reducing the matrix (or its transpose) to reduced row echelon form (RREF) and selecting
        the nonzero rows as basis vectors. The result is a matrix whose columns (or rows) form a basis
        for the specified subspace.

        Args:
            colspace (bool, optional): If `True`, returns a basis for the column space.
                If `False`, returns a basis for the row space.

            verbosity (int, optional): Level of output verbosity.

                - 0: No output.
                - 1: Print a summary of the RREF and basis selection.
                - 2: Print the matrix before and after RREF, and show the selected basis vectors.

        Returns:
            (Matrix): A matrix whose columns (if colspace=True) or rows (if colspace=False) form a basis
                for the corresponding subspace.

        Examples:
            >>> mat = Matrix([[1, 2], [2, 4]])
            >>> mat.simplify_basis(colspace=True, verbosity=0)
            Matrix([
            [1]
            [2]
            ])

            >>> mat = Matrix([[1, 2], [3, 4]])
            >>> mat.simplify_basis(colspace=False, verbosity=0)
            Matrix([
            [1, 0]
            [0, 1]
            ])
        """
        if colspace:
            rref_mat = self.T.rref(pivots=False)
            assert isinstance(rref_mat, Matrix), "RREF should return a Matrix"
            if verbosity == 1:
                print("Select non-zero rows of rref(self.T) as basis vectors.")
            if verbosity >= 2:
                print("Before RREF: self^T")
                display(self.T)
                print("\nAfter RREF:")
                display(rref_mat)
        else:
            rref_mat = self.rref(pivots=False)
            assert isinstance(rref_mat, Matrix), "RREF should return a Matrix"
            if verbosity == 1:
                print("Select non-zero rows of rref(self) as basis vectors.")
            if verbosity >= 2:
                print("Before RREF: self")
                display(self)
                print("\nAfter RREF:")
                display(rref_mat)

        idxs = []
        for i in range(rref_mat.rows):
            if any(not _is_zero(x) for x in rref_mat[i, :]):  # type: ignore
                idxs.append(i)

        mat = rref_mat.select_rows(*idxs)
        if colspace:
            return mat.T
        else:
            return mat

    def extend_basis(
        self, span_subspace: Matrix | None = None, verbosity: int = 2
    ) -> Matrix:
        r"""Extends the matrix to form a basis for the span of the given subspace.

        This method extends the column space of the current matrix to include the columns of the provided
        `span_subspace`, computes the Reduced Row Echelon Form (RREF) of the augmented matrix,
        and then selects the pivot columns to return the extended basis.

        If no `span_subspace` is provided, the identity matrix (i.e. $\mathrm{span} \left(\mathbb{R}^n \right)$) is used as the default.
        The result is a matrix with the extended basis that spans the combined space of the `self`
        matrix and the `span_subspace`.

        Args:
            span_subspace (Matrix, optional): A matrix whose columns represent the subspace to
                be added to the current matrix. If `None`, the identity matrix is used.

            verbosity (int, optional): Verbosity level for displaying information.

                - 0: No output.
                - 1: Display steps.
                - 2: Display the matrix before and after RREF.

        Returns:
            (Matrix): A matrix whose column space represents the extended basis, consisting of the pivot columns
                    from the RREF of the augmented matrix.

        Examples:
            >>> mat = Matrix([[1, 2], [2, 4]])
            >>> mat.extend_basis()
            Matrix([
            [1, 1]
            [2, 0]
            ])
        """

        if span_subspace is None:
            span_subspace = self.elem()
        aug = self.copy().row_join(span_subspace)
        rref = aug.rref(pivots=True)
        assert isinstance(rref, RREF), "RREF should return a RREF dataclass"
        rref_mat, pivots = rref.rref, rref.pivots

        if verbosity == 1:
            print("rref([self | span_subspace])")
        elif verbosity >= 2:
            print("Before RREF: [self | span_subspace]")
            display(aug)
            print("\nAfter RREF:")
            display(rref_mat)

        return aug.select_cols(*pivots)

    def intersect_subspace(self, other: Matrix, verbosity: int = 2) -> Matrix:
        """Computes the intersection of two subspaces by finding the nullspace of their orthogonal complements.

        This method computes the intersection of the subspaces spanned by the columns of the current matrix
        (`self`) and the provided matrix (`other`). The intersection is computed by finding the union of the nullspace of
        the row space of the two matrices, and then finding its orthogonal complement.

        Args:
            other (Matrix): The second matrix representing the other subspace to intersect with the current matrix.
            verbosity (int, optional): Level of verbosity for displaying intermediate steps:

                - 0: No output.
                - 1: Display steps.
                - 2: Display the relevant matrices.
                Defaults to 2.

        Returns:
            (Matrix): A matrix whose columns form a basis for the intersection of the two subspaces.

        Examples:
            >>> mat1 = Matrix([[1, 0], [0, 1]])
            >>> mat2 = Matrix([[1, 1], [0, 0]])
            >>> mat1.intersect_subspace(mat2, verbosity=0)
            Matrix([
            [1],
            [0]])
        """

        # Construct 2 matrices A and B, whose solution space (ie nullspace) is
        # the subspace self and other respectively. Observe that the solution
        # space is orthogonal to the row space, so it is the orthogonal complement.

        A = self.orthogonal_complement().T
        B = other.orthogonal_complement().T

        # Now we obtain A and B which represent the linear system of 2 different
        # subspaces. When we solve these simultaneously, we will find the solution
        # space which contains vectors which are solutions to both linear systems.
        aug = A.col_join(B)
        if verbosity == 1:
            print("A = Null(self^T)^T")
            print("B = Null(other^T)^T")
            print("Null([A ; B])")

        if verbosity >= 2:
            print(
                "A linear system whose solution space is the subspace of self. Null(self^T)^T"
            )
            display(A)
            print(
                "\nA linear system whose solution space is the subspace of other. Null(other^T)^T"
            )
            display(B)
            print("\nBefore RREF: [self ; other]")
            display(aug)
            print("\nAfter RREF:")
            display(aug.rref())

        return Matrix.from_list(aug.nullspace())

    def is_same_subspace(self, other: Matrix | None = None, verbosity: int = 2) -> bool:
        r"""Checks if two subspaces are the same by verifying if each subspace is a subspace of the other.

        This method determines whether the subspaces spanned by the columns of the current matrix (`self`)
        and the provided matrix (`other`) are the same. It does so by checking if the row-reduced echelon forms
        (RREF) of the augmented matrices `[self | other]` and `[other | self]` reveal that both subspaces
        span each other.

        Args:
            other (Matrix, optional): The second matrix representing the other subspace to compare with the current matrix.
                If `None`, the identity matrix is used to check if it spans $\mathbb{R}^\text{self.rows}}$
            verbosity (int, optional): Level of verbosity for displaying intermediate steps:

                - 0: No output.
                - 1: Display the steps.
                - 2: Display the relevant matrices.
                Defaults to 2.

        Returns:
            (bool): `True` if the subspaces spanned by `self` and `other` are the same, `False` otherwise.

        Raises:
            sympy.matrices.exceptions.ShapeError: If the number of rows in the current matrix and the target matrix are different.

        Examples:
            >>> mat1 = Matrix([[1, 0], [0, 1]])
            >>> mat2 = Matrix([[1, 2], [3, 4]])
            >>> mat1.is_same_subspace(mat2, verbosity=0)
            True
        """
        if other is None:
            rref = self.rref(pivots=True)
            assert isinstance(rref, RREF), "RREF should return a RREF dataclass"
            if verbosity >= 1:
                print("Check rref(self) does not have zero rows")
            if verbosity >= 2:
                print("Before RREF: self")
                display(self)
                print("\nAfter RREF:")
                display(rref.rref)
            return len(rref.pivots) == self.rows  # no zero rows

        if self.rows != other.rows:
            raise sym.ShapeError(
                f"The matrices have incompatible number of rows ({self.rows}, {other.rows})"
            )

        sub = other.copy().row_join(self).rref(pivots=True)
        assert isinstance(sub, RREF), "RREF should return a RREF dataclass"
        sub_rref, sub_pivots = sub.rref, sub.pivots
        super = self.copy().row_join(other).rref(pivots=True)
        assert isinstance(super, RREF), "RREF should return a RREF dataclass"
        super_rref, super_pivots = super.rref, super.pivots

        if verbosity == 1:
            print("Check rref([other | self])")
            print("Check rref([self | other])")
        elif verbosity >= 2:
            print("Check if span(self) is subspace of span(other)")
            print("\nBefore RREF: [other | self]")
            display(other.copy().row_join(self))
            print("\nAfter RREF:")
            display(sub_rref)
            if max(sub_pivots) >= other.cols:
                print(
                    f"Span(self) is not a subspace of span(other). Pivot columns: {sub_pivots}"
                )
            else:
                print(
                    f"Span(self) is a subspace of span(other). Pivot columns: {sub_pivots}"
                )

            print("\nCheck if span(other) is subspace of span(self)")
            print("\nBefore RREF: [self | other]")
            display(self.copy().row_join(other))
            print("\nAfter RREF:")
            display(super_rref)
            if max(super_pivots) >= self.cols:
                print(
                    f"Span(other) is not a subspace of span(self). Pivot columns: {super_pivots}"
                )
            else:
                print(
                    f"Span(other) is a subspace of span(self). Pivot columns: {super_pivots}"
                )

        return (max(sub_pivots) < other.cols) and (max(super_pivots) < self.cols)

    def coords_relative(self, to: Matrix, verbosity: int = 2) -> Matrix:
        """Computes the coordinates of the current vector relative to a given basis.

        This method finds the coordinate vector `c` such that `self = to @ c`,
        where `to` is a matrix whose columns form a basis, and `self` is a column vector. The method
        achieves this by augmenting the target matrix with the current matrix, performing
        Reduced Row Echelon Form (RREF), and extracting the appropriate part of the resulting matrix.

        Args:
            to (Matrix): The matrix whose columns form the target basis.

            verbosity (int, optional): Level of output verbosity.

                - 0: No output.
                - 1: Print summary of the RREF and solution.
                - 2: Print the matrix before and after RREF, and show the solution details.

        Returns:
            (Matrix): The coordinate vector of `self` relative to the basis `to`.

        Raises:
            sympy.matrices.exceptions.ShapeError: If `self` is not a column vector
                or if the number of rows in `self` and `to` do not match.
            Exception: If the system is inconsistent and no solution exists.

        Examples:
            >>> v = Matrix([[3], [7]])
            >>> B = Matrix([[1, 2], [1, 3]])
            >>> v.coords_relative(B, verbosity=0)
            Matrix([
            [-5],
            [ 4]])
        """
        if self.cols != 1:
            raise sym.ShapeError(
                f"self should be a vector with 1 column. ({self.cols})"
            )
        if self.rows != to.rows:
            raise sym.ShapeError(
                f"The matrices have incompatible number of rows ({self.rows}, {to.rows})"
            )

        M = to.copy().row_join(self)
        rref = M.rref(pivots=True)
        assert isinstance(rref, RREF), "RREF should return a RREF dataclass"
        rref_mat, pivots = rref.rref, rref.pivots

        if verbosity == 1:
            print("Solve system via rref([to | self])")
        elif verbosity >= 2:
            print("Before RREF: [to | self]")
            display(M)
            print("\nAfter RREF:")
            display(rref_mat)

        if to.cols in pivots:
            raise Exception("No solution found due to inconsistent system.")

        sol = to.solve(self)
        assert isinstance(sol, Matrix), "Solution should be a Matrix (vector) object"
        return sol

    def transition_matrix(self, to: Matrix, verbosity: int = 2) -> Matrix:
        """Computes the transition matrix that transforms this matrix to another matrix.

        This method computes the transition matrix `P` such that `self = P @ to`,
        where `to` is the target basis, and `self` is the current basis. The method
        achieves this by augmenting the target matrix with the current matrix, performing
        Reduced Row Echelon Form (RREF), and extracting the appropriate part of the resulting matrix.

        Args:
            to (Matrix): The matrix to which the current matrix should be transformed.

            verbosity (int, optional): Verbosity level for displaying information.

                - 0: No output.
                - 1: Display the steps.
                - 2: Display the relevant matrices.

        Returns:
            (Matrix): The transition matrix `P` that satisfies `self = P @ to`.

        Raises:
            AssertionError: If the columns of the `self` matrix and `to` matrix do not span the same subspace.

        Examples:
            >>> mat1 = Matrix([[1, 0], [0, 1]])
            >>> mat2 = Matrix([[2, 0], [0, 2]])
            >>> mat1.transition_matrix(to=mat2, verbosity=0)
            Matrix([
            [1/2,   0],
            [  0, 1/2]])
        """
        assert self.is_same_subspace(to, verbosity=0), (
            "Column vectors of both matrices must span the same subspace."
        )

        M = to.copy().row_join(self)
        res = M.rref(pivots=False)
        assert isinstance(res, Matrix), "RREF should return a Matrix object"
        if verbosity == 1:
            print("rref([to | self])")
        elif verbosity >= 2:
            print("Before RREF: [to | self]")
            display(M)
            print("\nAfter RREF:")
            display(res)
        P = res[: self.cols, self.cols :]
        return P  # type: ignore

    ###############################################
    # CHAPTER 4: SUBSPACES ASSOCIATED TO A MATRIX #
    ###############################################

    # Override
    def nullspace(self, verbosity: int = 0, *args, **kwargs) -> list[Matrix]:
        """
        Computes the null space (kernel) of the matrix, i.e., the set of vectors that satisfy `self @ x = 0`.

        This method utilizes the rank-nullity theorem to determine if the null space exists. Fixes the
        issue with SymPy [implementation][sympy.matrices.matrixbase.MatrixBase.nullspace] of nullspace where
        it raises an exception if the nullspace is trivial (only contain the zero vector).
        If the matrix has full column rank (i.e., rank = number of columns), it has no non-trivial null space,
        and an empty list is returned.

        Args:
            verbosity (int, optional): Level of verbosity for displaying intermediate steps.
                - 0: No output.
                - 1: Display the matrix before and after row-reduction (RREF).
                Default is 0.
            *args: Additional positional arguments passed to SymPy's [`nullspace`][sympy.matrices.matrixbase.MatrixBase.nullspace] method.
            **kwargs: Additional keyword arguments passed to SymPy's [`nullspace`][sympy.matrices.matrixbase.MatrixBase.nullspace] method.

        Returns:
            list[Matrix]: A list of [`Matrix`][...] objects representing the null space vectors. Returns an empty list if the null space is trivial.

        Examples:
            >>> mat = Matrix([[1, 2], [3, 6]])
            >>> mat.nullspace()
            [Matrix([
            [-2],
            [ 1]])]

        See Also:
            - SymPy's [`Matrix.nullspace`][sympy.matrices.matrixbase.MatrixBase.nullspace]
            - [`orthogonal_complement`][..] for computing the orthogonal complement of the matrix.
        """

        # Issue with SymPy implementation of nullspace when there is None
        # Using rank nullity theorem to verify there are vectors spanning nullspace
        if verbosity >= 1:
            print("Before RREF: [self]")
            display(self)
            print("\nAfter RREF:")
            display(self.rref())

        if self.rank() == self.cols:
            if verbosity >= 1:
                warn("No nullspace detected!", UserWarning)
            return []
        else:
            return super().nullspace(*args, **kwargs)

    #######################################################
    # CHAPTER 5: ORTHOGONALITY AND LEAST SQUARES SOLUTION #
    #######################################################

    def orthogonal_complement(self, verbosity: int = 0) -> Matrix:
        """Computes the orthogonal complement of the matrix (the null space of its transpose).

        The orthogonal complement consists of all vectors that are orthogonal to the column space of the matrix.
        This method computes the null space of the transpose of the matrix, which gives a basis for the orthogonal complement.

        Note:
            The orthogonal complement is the set of all vectors `v` such that `self^T @ v = 0`.

        Args:
            verbosity (int, optional): Level of verbosity for debugging.

        Returns:
            (Matrix): A matrix whose columns form a basis for the orthogonal complement.

        Examples:
            >>> mat = Matrix([[1, 0], [0, 0]])
            >>> mat.orthogonal_complement()
            Matrix([
            [0],
            [1]])
        """

        return Matrix.from_list(self.transpose().nullspace(verbosity))

    def is_vec_orthogonal(self, verbosity: int = 1) -> bool:
        r"""Checks if the column vectors of the matrix are orthogonal to each other.

        This method computes `self^T @ self` and checks if the result is diagonal.
        If the result is diagonal, the vectors are orthogonal (i.e., $u_i \dot u_j = 0 \forall i != j$).

        Note:
            This method checks for orthogonality, not orthonormality. For orthonormality, use [`is_mat_orthogonal`][..].

        Args:
            verbosity (int, optional): Level of verbosity for displaying intermediate results.

                - 0: No output.
                - 1: Display the matrix product `self.T @ self`.

        Returns:
            (bool): `True` if the column vectors are orthogonal, `False` otherwise.

        Examples:
            >>> mat = Matrix([[1, 0], [0, 2]])
            >>> mat.is_vec_orthogonal(verbosity=0)
            True
        """

        res = self.T @ self
        if verbosity >= 1:
            print("self^T @ self")
            display(res)
        return res.is_diagonal()

    def is_mat_orthogonal(self, verbosity: int = 1) -> bool:
        """
        Checks if the matrix is orthogonal (i.e., its columns are orthonormal).

        A matrix is orthogonal if its columns are orthonormal, i.e., if `self.T @ self` is the identity matrix.
        This method computes `self.T @ self` and checks if the result is an identity matrix.

        Args:
            verbosity (int, optional): Level of verbosity for displaying intermediate results.

                - 0: No output.
                - 1: Display the matrix product `self.T @ self`.

        Returns:
            (bool): True if the matrix is orthogonal, False otherwise.

        Examples:
            >>> mat = Matrix([[1, 0], [0, 1]])
            >>> mat.is_mat_orthogonal(verbosity=0)
            True

            >>> mat = Matrix([[1, 0], [0, 2]])
            >>> mat.is_mat_orthogonal(verbosity=0)
            False
        """

        res = self.T @ self
        if verbosity >= 1:
            print("self^T @ self")
            display(res)
        return res.is_diagonal and all(entry == 1 for entry in res.diagonal())

    def orthogonal_decomposition(self, to: Matrix, verbosity: int = 0) -> VecDecomp:
        """Decomposes the current vector (or matrix) into its orthogonal projection onto a subspace and its orthogonal complement.

        This method computes the orthogonal decomposition of `self` relative to the subspace spanned by the columns of `to`.
        It finds the projection of `self` onto the subspace (`proj`) and the component orthogonal to the subspace (`norm`), such that:
            `self = proj + norm`

        The projection is computed using the least squares solution.

        Args:
            to (Matrix): The matrix whose columns form the subspace onto which to project `self`.
            verbosity (int, optional): Level of verbosity for displaying intermediate results.

                - 0: No output.
                - 1: Display the projected and normal components.
                - 2: Display detailed steps.

        Returns:
            (VecDecomp): A dataclass with fields:

                - proj ([`Matrix`][...]): The projection of `self` onto the subspace spanned by `to`.
                - norm ([`Matrix`][...]): The component of `self` orthogonal to the subspace spanned by `to`.

        Examples:
            >>> v = Matrix([[3], [4]])
            >>> B = Matrix([[1], [0]])
            >>> v.orthogonal_decomposition(B)
            VecDecomp(proj=Matrix([
            [3],
            [0]]), norm=Matrix([
            [0]
            [4]
            ]))
        """

        sol = to.solve_least_squares(self, verbosity=verbosity)
        proj = to @ sol
        norm = self - proj

        if verbosity >= 1:
            print("Projected component: Au")
            display(proj)
            print("Normal component: b - b_proj")
            display(norm)

        assert proj + norm == self
        return VecDecomp(proj, norm)

    def proj_comp(self, to: Matrix, verbosity: int = 0) -> Matrix:
        """Computes the orthogonal projection of the current vector (or matrix) onto the subspace spanned by the columns of another matrix.

        This method returns the component of `self` that lies in the subspace defined by the columns of `to`.
        It is equivalent to the projection of `self` onto the subspace, as computed by the orthogonal decomposition.

        Args:
            to (Matrix): The matrix whose columns form the subspace onto which to project `self`.
            verbosity (int, optional): Level of verbosity for displaying intermediate results.

                - 0: No output.
                - 1: Display the projected component.

        Returns:
            (Matrix): The projection of `self` onto the subspace spanned by `to`.

        Examples:
            >>> v = Matrix([[3], [4]])
            >>> B = Matrix([[1], [0]])
            >>> v.proj_comp(B)
            Matrix([
            [3],
            [0]])
        """

        return self.orthogonal_decomposition(to=to, verbosity=verbosity).proj

    def norm_comp(self, to: Matrix, verbosity: int = 0) -> Matrix:
        """Computes the component of the current vector (or matrix) orthogonal to the subspace spanned by the columns of another matrix.

                This method returns the part of `self` that is perpendicular to the subspace defined by the columns of `to`.
                It is equivalent to the normal component from the orthogonal decomposition.

                Args:
                    to (Matrix): The matrix whose columns form the subspace to which the orthogonal component is computed.
                    verbosity (int, optional): Level of verbosity for displaying intermediate results.

                        - 0: No output (default).
                        - 1: Display the normal component.

                Returns:
                    (Matrix): The component of `self` orthogonal to the subspace spanned by `to`.

                Examples:
                    >>> v = Matrix([[3], [4]])
                    >>> B = Matrix([[1], [0]])
                    >>> v.norm_comp(B)
                    Matrix([
        [0],
        [4]])
        """

        return self.orthogonal_decomposition(to=to, verbosity=verbosity).norm

    def gram_schmidt(
        self, factor: bool = True, verbosity: int = 1
    ) -> Matrix | ScalarFactor:
        """Performs Gram-Schmidt orthogonalization to convert a set of vectors (columns of the matrix) into
        an orthogonal set (that includes 0 vectors if any).

        Args:
            factor (bool): If `True`, the resulting orthogonal vectors will be scaled to have integer factors.
            verbosity (int): Level of verbosity (default is `1`):

                - 0: No output.
                - 1: Display intermediate results for each step of the process.

        Returns:
            (Matrix): A matrix whose columns are the orthogonalized vectors.

        Examples:
            >>> mat = Matrix([[1, 2], [3, 4]])
            >>> mat.gram_schmidt(factor=False, verbosity=0)
            Matrix([
            [  sqrt(10)/10, 3*sqrt(10)/10]
            [3*sqrt(10)/10,  -sqrt(10)/10]
            ])
        """

        if self.cols == 0:
            return Matrix()
        if verbosity >= 1:
            display(f"v_{1} = {sym.latex(self.select_cols(0))}", opt="math")

        orthogonal_set = [self.select_cols(0)]
        for i in range(1, self.cols):
            u = self.select_cols(i)
            latex_eq = f"v_{i + 1} = {sym.latex(u)}"
            for _, v in enumerate(orthogonal_set, start=1):
                if v.norm() != 0:
                    latex_eq += f"- ({sym.latex(v.dot(u))}/{sym.latex(v.dot(v))}) {sym.latex(v)}"
                    u -= (v.dot(u) / v.dot(v)) * v

            if verbosity >= 1:
                disp_u = u.copy()
                if factor:
                    scalar = sym.gcd(tuple(u))  # type: ignore
                    disp_u = sym.MatMul(scalar, u / scalar, evaluate=False)
                latex_eq += f" = {sym.latex(disp_u)}"
                display(latex_eq, opt="math")

            if u.norm() == 0 and (verbosity >= 1):
                warn(
                    "Vectors are linearly dependent. Note that there is no QR factorisation",
                    UserWarning,
                )
            orthogonal_set.append(u)

        return Matrix.from_list(orthogonal_set).normalized(factor=factor)

    # Override
    def QRdecomposition(self, full: bool = False) -> QR:
        """Computes the QR decomposition of the matrix. Optionally computes the full QR decomposition.

        A full QR decomposition returns an **orthogonal (square) matrix** `Q` and an upper triangular matrix `R`
        such that `self = Q @ R`. On the other hand, a reduced QR decomposition returns `Q` and `R` such that
        `self = Q @ R`, where `Q` has **orthonormal columns**.

        Args:
            full (bool): If `True`, computes the full QR decomposition.

        Returns:
            (QR): A dataclass containing:

                - `Q` ([`Matrix`][...]): A matrix with orthonormal columns (or orthogonal matrix if full QR decomposition).
                - `R` ([`Matrix`][...]): An upper triangular matrix.

        Examples:
            >>> mat = Matrix([[1, 2], [3, 4]])
            >>> mat.QRdecomposition()
            QR(Q=Matrix([
            [  sqrt(10)/10, 3*sqrt(10)/10]
            [3*sqrt(10)/10,  -sqrt(10)/10]
            ]), R=Matrix([
            [sqrt(10), 7*sqrt(10)/5]
            [       0,   sqrt(10)/5]
            ]))

        See Also:
            - SymPy's [`Matrix.QRdecomposition`][sympy.matrices.matrixbase.MatrixBase.QRdecomposition]
        """

        # Modified SymPy's implementation to compute full QR decomposition if required.
        Q, R = super().QRdecomposition()
        if full and Q.rows != Q.cols:
            Q = Matrix(Q)
            Q_aug = Q.row_join(Q.elem(), aug_line=False).QRdecomposition()[0]
            R_aug = Matrix(R.col_join(sym.zeros(Q_aug.cols - R.rows, R.cols)))
            assert Q_aug @ R_aug == self
            return QR(Q_aug, R_aug)
        return QR(Q, R)

    def solve_least_squares(
        self, rhs: Matrix, verbosity: int = 1, matrices: int = 1, *args, **kwargs
    ) -> Matrix | PartGen:
        r"""Solves the least squares problem $\min || \mathrm{self} \, \mathbf{x} - \mathrm{rhs}||^2$.

        Uses SymPy's built-in method for least squares when the rank condition is met, otherwise uses a custom
        solution approach using the normal equations: $\mathrm{self}^{\top} \mathrm{self} \, \mathbf{x} = \mathrm{self}^{\top} \mathrm{rhs}$

        Args:
            rhs (Matrix): The right-hand side matrix/vector `b` in `Ax = b`.
            verbosity (int, optional): Level of verbosity (default is `1`):

                - 0: No output.
                - 1: Display intermediate steps.
            matrices (int, optional):

                - 1: Returns the least squares solution matrix.
                - 2: Returns a [`PartGen`][(p).] with the part solution and general solution.
            *args: Additional positional arguments passed to SymPy's
                [`solve_least_squares`][sympy.matrices.matrixbase.MatrixBase.solve_least_squares] method.
            **kwargs: Additional arguments passed to to SymPy's
                [`solve_least_squares`][sympy.matrices.matrixbase.MatrixBase.solve_least_squares] method.

        Returns:
            (Matrix): If `matrices=1`, returns the least squares solution matrix.
            (PartGen): If `matrices=2`, returns a [`PartGen`][(p).] with the part solution and general solution.

        Examples:
            >>> A = Matrix([[1, 0], [0, 0]])
            >>> b = Matrix([1, 2])
            >>> A.solve_least_squares(b, verbsosity=0, matrices=1)
            Exception Encountered: Matrix must be non-singular.
            Attempting custom solve...
            Matrix([
            [1],
            [y]]) # in this case, least squares solution is not unique, so it returns a general solution.
        """

        if verbosity == 0:
            try:
                A, b = sym.Matrix(self), sym.Matrix(rhs)
                return A.solve_least_squares(rhs=b, *args, **kwargs)
            except Exception as e:
                print(f"Exception Encountered: {str(e)}")
                print("Attempting custom solve...")

        ATA, ATb = self.T @ self, self.T @ rhs
        sol = Matrix.create_unk_matrix(ATb.rows, 1)
        sol = sol.subs(sym.solve(ATA @ sol - ATb))

        if verbosity >= 1:
            print("Before RREF: [self.T self | self.T rhs]")
            aug_matrix = ATA.copy().row_join(ATb)
            display(aug_matrix)
            print("\nAfter RREF")
            display(aug_matrix.rref())

        if matrices == 1:
            return sol
        else:
            return sol.sep_part_gen()

    @staticmethod
    def create_vander(
        num_rows: int = 1, num_cols: int = 1, symbol: str = "x", is_real: bool = True
    ) -> Matrix:
        """Creates a Vandermonde matrix with symbolic entries.

        This method generates a Vandermonde matrix of size `num_rows` x `num_cols`
        where the entries are symbolic expressions. Each row in the matrix is formed
        by raising a symbolic variable (indexed by row) to increasing powers (from 0
        to `num_cols-1`). The `is_real` flag determines whether the symbols are real-valued.

        Args:
            num_rows (int, optional): The number of rows in the Vandermonde matrix.
            num_cols (int, optional): The number of columns in the Vandermonde matrix.
            symbol (str, optional): The base name for the symbols used in the matrix entries.
            is_real (bool, optional): If True (default), the symbols are real-valued;
                otherwise, they are complex.

        Returns:
            (Matrix): A Vandermonde matrix with symbolic entries.

        Examples:
            >>> Matrix.create_vander(2, 4, symbol='a')
            Matrix([
            [1, a_1, a_1**2, a_1**3]
            [1, a_2, a_2**2, a_2**3]
            ])

        See Also:
            - [`apply_vander`][..] for applying the Vandermonde transformation to a matrix.
        """

        entries = sym.symbols(f"{symbol}_(1:{num_rows + 1})", is_real=is_real)
        res = []
        for entry in entries:
            sub_res = []
            for col_idx in range(num_cols):
                # Raise the symbol to the power of the column index
                sub_res.append(sym.Pow(entry, col_idx))
            res.append(sub_res)
        return Matrix(res)

    def apply_vander(self, x: Matrix) -> Matrix:
        """
        Applies a Vandermonde transformation to the current matrix using the given vector.

        This method applies a Vandermonde transformation to the current matrix by
        substituting the free symbols in the last column with corresponding values
        from the provided vector `x`. The number of rows in `self` must match the
        number of elements in `x`, and `x` must be a column vector.

        Note:
            - The matrix `self` is expected to be created via [`Matrix.create_vander()`][..create_vander].
            - The `x` vector provides the values to substitute in place of these symbols.

        Args:
            x (Matrix): A column vector (Matrix object with a single column) containing
                the values to substitute into the last column of the matrix.

        Returns:
            (Matrix): A new Matrix object where the free symbols in the last column of
                the original matrix are substituted by the corresponding values from `x`.

        Raises:
            sympy.matrices.exceptions.ShapeError: If `x` is not a column vector or if
                the number of rows in `self` does not match the size of `x`.

        Examples:
            >>> mat = Matrix.create_vander(2, 2)
            >>> x = Matrix([1, 2])
            >>> mat.apply_vander(x)
            Matrix([
            [1, 1],
            [1, 2]])
        """
        # Validate the size of the vector x
        if x.cols != 1:
            raise sym.ShapeError(
                f"Input vector x must be a column vector. ({self.cols})"
            )
        if self.rows != x.rows:
            raise sym.ShapeError(
                f"Number of rows in matrix ({self.rows}) must match the size of the input vector ({x.rows})"
            )

        # Get the free symbols from the last column of the matrix
        ordered_syms = [entry.free_symbols.pop() for entry in self.select_cols(-1)]  # type: ignore

        # Create a substitution dictionary mapping symbols to values from vector x
        substitution = {var: val for var, val in zip(ordered_syms, x)}  # type: ignore
        return self.subs(substitution)

    #############################
    # CHAPTER 6: EIGEN-ANALYSIS #
    #############################

    def cpoly(self, force_factor: bool = True) -> Mul | tuple[Mul, Mul]:
        """Computes the characteristic polynomial of the matrix and attempts to factor it into real and complex parts.

        The characteristic polynomial is defined as `det(x * I - self)`, where `I` is the identity matrix of the same size as `self`.

        Args:
            force_factor (bool): If `True`, the polynomial is fully factored, even if it doesn't have real factors.
                If `False`, the polynomial is returned in its factored form if possible.

        Returns:
            (Mul): If the polynomial factors only into real terms, returns a single factored polynomial.
            (tuple[Mul, Mul]): If the polynomial has both real and complex factors, returns a tuple of two polynomials
                one with real factors and the other with complex factors.

        Examples:
            >>> mat = Matrix([[-1, 0], [0, 4]])
            >>> mat.cpoly()
            (x - 4)*(x + 1)
        """
        x = sym.symbols("x", real=True)
        poly = (x * self.elem() - self).det()
        if not force_factor:
            return poly.factor()
        # Attempt to factor poly into real factors
        try:
            roots = sym.roots(
                poly
            )  # TODO: FIX sym.roots NotImplementedError for multi variable
            real_fact = []
            for root, mult in roots.items():
                term = x - root
                if mult != 1:
                    term = sym.Pow(term, mult, evaluate=False)
                if root.is_real:
                    real_fact.append(term)
                    poly /= term

            linear_fact = Mul(*real_fact, evaluate=False)
            complex_fact = poly.expand().cancel().factor()

            if complex_fact == 1:
                return linear_fact  # type: ignore
            else:
                return linear_fact, complex_fact  # type: ignore
        except Exception as error:
            print(f"Encountered Error: {error}")
            return poly.factor()

    # Override
    def is_diagonalizable(
        self, reals_only: bool = True, verbosity: int = 1, *args, **kwargs
    ) -> bool:
        """Checks if the matrix is diagonalizable, with the option to focus only on real eigenvalues.

        A matrix is diagonalizable if it has enough linearly independent eigenvectors to form a basis for the space.

        Args:
            reals_only (bool, optional): If True, diagonalization will focus on real eigenvalues.
            verbosity (int, optional): Controls the level of output during the diagonalization process.
            *args: Additional positional arguments passed to SymPy's
                [`is_diagonalizable`][sympy.matrices.matrixbase.MatrixBase.is_diagonalizable] method.
            **kwargs: Additional arguments passed to SymPy's
                [`is_diagonalizable`][sympy.matrices.matrixbase.MatrixBase.is_diagonalizable] method.

        Returns:
            (bool): True if the matrix is diagonalizable, False otherwise.

        Examples:
            >>> mat = Matrix([[1, 2, 0], [0, 3, 0], [2, -4, 2]])
            >>> mat.is_diagonalizable(reals_only=True, verbosity=0)
            True

        See Also:
            - SymPy's [`Matrix.is_diagonalizable`][sympy.matrices.matrixbase.MatrixBase.is_diagonalizable]
            - Sympy's [`Matrix.eigenvects`][sympy.matrices.matrixbase.MatrixBase.eigenvects] for computing eigenvalues and their multiplicities.
            - [`diagonalize`][..] for diagonalizing the matrix.
        """

        # Changed default for reals_only to True to align with MA1522 syllabus
        if verbosity >= 1:
            print("Characteristic Polynomial is: ")
            display(self.cpoly())
            print("Eigenvectors are:")
            for val, mult, space in self.eigenvects():
                if (val.is_real and reals_only) or not reals_only:
                    res = {
                        "eigenvalue": val,
                        "algebraic multiplicity": mult,
                        "eigenspace": Matrix.from_list(space),
                    }
                    display(res, opt="dict")

        return super().is_diagonalizable(reals_only, *args, **kwargs)

    def eigenvects_associated(
        self, eigenvalue: Expr | int | float
    ) -> list[Matrix] | None:
        """Computes the eigenvectors associated with a given eigenvalue.

        This method finds all (nonzero) vectors `v` such that `(eigenvalue * I - self) * v = 0`
        where `I` is the identity matrix of the same size as `self`.

        Args:
            eigenvalue (Expr | int | float): The eigenvalue for which to compute the associated eigenvectors.

        Returns:
            (list[Matrix]): A list of eigenvectors (as Matrix objects) associated with the given eigenvalue,
            (None): If the eigenvalue does not correspond to any eigenvectors.

        Examples:
            >>> mat = Matrix([[2, 0], [0, 3]])
            >>> mat.eigenvects_associated(2)
            [Matrix([
            [1],
            [0]])]
        """
        return (eigenvalue * self.elem() - self).nullspace()

    # Override
    def diagonalize(
        self, reals_only: bool = True, verbosity: int = 0, *args, **kwargs
    ) -> PDP:
        """Diagonalizes the matrix if possible, focusing on real eigenvalues unless specified otherwise.

        Args:
            reals_only (bool, optional): If `True`, diagonalization will focus on real eigenvalues.
            verbosity (int, optional): Controls the level of output during the diagonalization process.

                - 0: No output.
                - 1: Displays the characteristic polynomial and eigenvectors.
            *args: Additional positional arguments passed to SymPy's
                [`diagonalize`][sympy.matrices.matrixbase.MatrixBase.diagonalize] method.
            **kwargs: Additional arguments passed to SymPy's
                [`diagonalize`][sympy.matrices.matrixbase.MatrixBase.diagonalize] method.

        Returns:
            (PDP): A dataclass containing:

                - `P` ([`Matrix`][...]): The matrix of eigenvectors.
                - `D` ([`Matrix`][...]): The diagonal matrix of eigenvalues.

        Examples:
            >>> mat = Matrix([[1, 2], [3, 4]])
            >>> mat.diagonalize()
            PDP(P=Matrix([
            [-sqrt(33)/6 - 1/2, -1/2 + sqrt(33)/6]
            [                1,                 1]
            ]), D=Matrix([
            [5/2 - sqrt(33)/2,                0]
            [               0, 5/2 + sqrt(33)/2]
            ]))
        """

        # Changed default for reals_only to True to align with MA1522 syllabus
        if verbosity >= 1:
            print("Characteristic Polynomial")
            poly = self.cpoly()
            display(poly)
            for root, _ in sym.roots(poly).items():
                if root.is_real:
                    display(
                        _textify("Before RREF: ")
                        + sym.latex(root)
                        + r"\mathbb{I} - \mathrm{self}",
                        opt="math",
                    )
                    expr = root * self.elem() - self
                    display(expr)

                    print("\nAfter RREF:")
                    display(expr.rref())

                    print("\nEigenvectors:")
                    display(expr.nullspace())
                    print("\n")

        P, D = super().diagonalize(reals_only, *args, **kwargs)
        P.rm_aug_line()  # Remove augmented line if exists
        return PDP(P, D)

    def is_orthogonally_diagonalizable(self, verbosity: int = 2) -> bool:
        """Determines whether the matrix is orthogonally diagonalizable.

        A matrix is orthogonally diagonalizable if and only if it is symmetric.
        This method checks the symmetry of the matrix and optionally displays
        diagnostic information based on the verbosity level.

        Args:
            verbosity (int, optional): Level of diagnostic output.

                - 0: No output.
                - 1: Displays the matrix.
                - 2: Displays the result of the symmetry check.

        Returns:
            (bool): True if the matrix is symmetric (orthogonally diagonalizable), False otherwise.

        Examples:
            >>> mat = Matrix([[1, 2], [2, 1]])
            >>> mat.is_orthogonally_diagonalizable(verbosity=0)
            True

        See Also:
            - SymPy's [`Matrix.is_symmetric`][sympy.matrices.matrixbase.MatrixBase.is_symmetric]

        """
        if verbosity:
            print("Check if matrix is symmetric:")
        if verbosity >= 2:
            print("\nCheck if self == self^T:")
            display(self == self.T)
        return self.is_symmetric()

    # Override
    def orthogonally_diagonalize(
        self, reals_only: bool = True, factor: bool = True, verbosity=1, *args, **kwargs
    ) -> PDP:
        """Orthogonally diagonalizes the matrix, ensuring that eigenvectors corresponding to different eigenvalues are orthogonal.

        Args:
            reals_only (bool): If True, only real eigenvalues are considered.
            factor (bool): If True, the eigenvectors are orthogonalized using the Gram-Schmidt process.
            verbosity (int): Controls the verbosity of output during the process.
            *args: Additional positional arguments passed to the [`diagonalize`][..] method.
            **kwargs: Additional arguments passed to the [`diagonalize`][..] method.

        Returns:
            (PDP): A dataclass containing:

                - `P` ([`Matrix`][...]): The orthogonal matrix of eigenvectors.
                - `D` ([`Matrix`][...]): The diagonal matrix of eigenvalues.

        Raises:
            AssertionError: If the matrix is not orthogonally diagonalizable (i.e., not symmetric).

        Examples:
            >>> mat = Matrix([[1, 2], [2, 1]])
            >>> mat.orthogonally_diagonalize(factor=False, verbosity=0)
            PDP(P=Matrix([
            [-sqrt(2)/2, sqrt(2)/2]
            [ sqrt(2)/2, sqrt(2)/2]
            ]), D=Matrix([
            [-1, 0]
            [ 0, 3]
            ]))

        See Also:
            - [`is_orthogonally_diagonalizable`][..] to check if the matrix is orthogonally diagonalizable.
            - [`diagonalize`][..] for diagonalizing the matrix.
        """

        # Changed default for reals_only to True to align with MA1522 syllabus
        # Note that you can just apply GSP on P directly, since eigenspace associated to different eigenvalues are orthogonal
        # However, we follow the steps given in MA1522 syllabus here
        assert self.is_orthogonally_diagonalizable(verbosity=verbosity)
        # P, D = super().diagonalize(reals_only, *args, **kwargs)
        P, D = self.diagonalize(
            reals_only=reals_only, verbosity=verbosity, *args, **kwargs
        )

        d: DefaultDict[Expr, list[Matrix]] = defaultdict(list)
        for vec, val in zip(P.columnspace(), D.diagonal()):
            d[val].append(vec)

        result = []
        for val, vecs in d.items():
            if len(vecs) > 1:
                # Require Gram Schmidt to ensure eigenvectors are orthogonal
                if verbosity >= 1:
                    print("Eigenvalue: ", val)
                    print("[Gram Schmidt Process]")
                if factor:
                    gram_result = Matrix.from_list(vecs).gram_schmidt(
                        factor=True, verbosity=verbosity
                    )
                    if isinstance(gram_result, ScalarFactor):
                        result.append(gram_result.full)
                    elif isinstance(gram_result, Matrix):
                        result.append(gram_result)
                    else:
                        raise TypeError(
                            f"Unexpected return type from gram_schmidt: {type(gram_result)}"
                        )
                else:
                    result.append(
                        Matrix.from_list(vecs).gram_schmidt(factor, verbosity)
                    )
            else:
                result.append(vecs[0].normalized())

        if len(result) == 0:
            ortho_P = P
        else:
            ortho_P = result[0]
            for m in result[1:]:
                ortho_P = ortho_P.row_join(m, aug_line=False)

        assert (ortho_P @ D @ ortho_P.T - self).norm() == 0
        return PDP(ortho_P, D)

    def is_stochastic(self, verbosity: int = 1) -> bool:
        """Checks if the matrix is stochastic.

        A matrix is stochastic if all its entries are non-negative and each column sums to 1.
        This property is commonly used to identify transition matrices in Markov chains.

        Args:
            verbosity (int, optional): Level of diagnostic output.

                - 0: No output.
                - 1: Displays the result of the checks.

        Returns:
            (bool): True if the matrix is stochastic, False otherwise.

        Examples:
            >>> mat = Matrix([[0.5, 0.1], [0.5, 0.9]])
            >>> mat.is_stochastic(verbosity=0)
            True

            >>> mat = Matrix([[1.1, 1], [-0.1, 0]])
            >>> mat.is_stochastic(verbosity=0)
            False
        """
        is_square = self.rows == self.cols
        is_non_negative = all(entry >= 0 for entry in self.flat())
        is_prob_vectors = all(sum(self[:, i]) == 1 for i in range(self.cols))  # type: ignore
        if verbosity >= 1:
            print(f"Check if matrix is square: {is_square}")
            print(f"Check if all entries are non-negative: {is_non_negative}")
            print(f"Check if each column sums to 1: {is_prob_vectors}")
        return is_square and is_non_negative and is_prob_vectors

    def equilibrium_vectors(self) -> Matrix:
        """Computes the equilibrium vectors of the matrix, i.e., the nullspace of (I - A).

        Note:
            - A matrix `P` has a unique equilibrium vector if it is stochastic and
                there exists some positive integer `k` such that `P^k` only has positive entries.

        Returns:
            (Matrix): A matrix containing equilibrium vectors normalized so that their
                column sums to 1.

        Examples:
            >>> mat = Matrix([[0.1, 0.9], [0.9, 0.1]])
            >>> mat.simplify()
            >>> mat.equilibrium_vectors()
            Matrix([
            [1/2],
            [1/2]])
        """

        P = Matrix.from_list((self.elem() - self).nullspace())
        for i in range(P.cols):
            if sum(P[:, i]) != 0:  # type: ignore
                P[:, i] /= sum(P[:, i])  # type: ignore
        return P

    def singular_value_decomposition(
        self, verbosity: int = 0, tol: float = 0.0, verify: bool = False
    ) -> SVD:
        """Performs Singular Value Decomposition (SVD) on the matrix, following the MA1522 syllabus.

        Note:
            - This function is known to take too much time and may kill Jupyter's kernel. Please use it with caution.
                A workaround is to set `verify=False` to skip the verification step, or use the faster numerical SVD
                method [`fast_svd`][..fast_svd] instead.

        Args:
            verbosity (int, optional): Controls the verbosity of the output.

                - 0: No output.
                - 1: Displays intermediate steps and results of the SVD process.
            tol (float, optional): Tolerance for verification of the SVD result.
            verify (bool): If `True`, verifies the result of the SVD by checking if `self = U @ S @ V.T`.
                If `False`, skips the verification step for performance reasons.

        Returns:
            (SVD): A dataclass containing:

                - `U` ([`Matrix`][...]): The left singular vectors.
                - `S` ([`Matrix`][...]): The diagonal matrix of singular values.
                - `V` ([`Matrix`][...]): The right singular vectors.

                Such that `self = U @ S @ V.T`.

        Examples:
            >>> mat = Matrix([[3, 2, 2], [2, 3, -2]])
            >>> U, S, V = mat.singular_value_decomposition(verbosity=0, verify=False)
            SVD(U=Matrix([
            [sqrt(2)/2, -sqrt(2)/2]
            [sqrt(2)/2,  sqrt(2)/2]
            ]), S=Matrix([
            [5, 0, 0]
            [0, 3, 0]
            ]), V=Matrix([
            [sqrt(2)/2,   -sqrt(2)/6,  2/3]
            [sqrt(2)/2,    sqrt(2)/6, -2/3]
            [        0, -2*sqrt(2)/3, -1/3]
            ]))

        See Also:
            - [`fast_svd`][..fast_svd] for a faster numerical SVD
            - SymPy's [`Matrix.singular_value_decomposition`][sympy.matrices.matrixbase.MatrixBase.singular_value_decomposition]
        """

        if verbosity >= 1:
            AT_A = self.T @ self
            print("A^T A")
            display(AT_A)
            P, D = AT_A.orthogonally_diagonalize(verbosity=verbosity)
            # Reverse index such that singular values are in decreasing order
            sigma = [sym.sqrt(val) for val in D.diagonal()][::-1]
            S = Matrix.diag(*[singular for singular in sigma if (singular != 0)])
            V = P.select_cols(*[i for i in range(P.cols)][::-1])

            u_list = []
            for idx, vec, val in zip(range(1, S.rows + 1), V.columnspace(), sigma):
                if val != 0:
                    u_i = self @ vec / val
                    u_list.append(u_i)
                    display(
                        f"u_{idx} = (1/{sym.latex(val)})A{sym.latex(vec)} = {sym.latex(u_i)}",
                        opt="math",
                    )

            U = Matrix.from_list(u_list)
            # Extend basis using orthogonal complement and gram-schmidt if insufficient vectors
            if U.cols < self.rows:
                print("\nExtending U with its orthogonal complement.")
                if U.cols == 0:
                    # Pad edge case with identity
                    orth = Matrix.eye(self.rows)
                else:
                    complement = U.orthogonal_complement(verbosity=verbosity)
                    gram_result = complement.gram_schmidt(
                        factor=True, verbosity=verbosity
                    )
                    if isinstance(gram_result, ScalarFactor):
                        orth = gram_result.full
                    else:
                        orth = gram_result

                    orth = orth.normalized(factor=False)
                    assert isinstance(orth, Matrix), (
                        f"Expected orth to be a Matrix, got {type(orth)}"
                    )
                U = U.row_join(orth, aug_line=False)

            # Add zero rows and columns to S so that matrix multiplication is defined
            m, n = self.shape
            r, c = S.shape
            S = S.row_join(sym.zeros(r, n - c), aug_line=False).col_join(
                sym.zeros(m - r, n)
            )

            assert (U @ S @ V.T - self).norm() == 0
            return SVD(U, S, V)

        m, n = self.shape
        U, S, V = super().singular_value_decomposition()
        # Reverse index such that singular values are in decreasing order
        new_S = Matrix.diag(*S.diagonal()[::-1])

        S_index = [i for i in range(S.cols)][::-1]
        new_U = Matrix(U).select_cols(*S_index)
        new_V = Matrix(V).select_cols(*S_index)

        # new_U = Matrix(U).
        # Add orthonormal columns to U and V so that they are square matrices
        new_U = new_U.QRdecomposition(full=True).Q
        new_V = new_V.QRdecomposition(full=True).Q

        # Add zero rows and columns to S so that matrix multiplication is defined
        r, c = new_S.shape
        new_S = new_S.row_join(sym.zeros(r, n - c), aug_line=False).col_join(
            sym.zeros(m - r, n)
        )

        if verify and (residues := (new_U @ new_S @ new_V.T - self).norm()) > tol:
            res = residues.evalf()
            warn(
                f"Verification failed: norm of residual is {res} > {tol}",
                UserWarning,
            )
        return SVD(new_U, new_S, new_V)

    def fast_svd(
        self,
        option: Literal["np", "sym"] = "np",
        identify: bool = True,
        tol: float | None = None,
    ) -> SVD | NumSVD:
        """A faster version of SVD that computes numerically using NumPy's SVD function.

        This method is designed to be efficient and suitable for large matrices, but it does not guarantee
        exact symbolic results like the [`singular_value_decomposition`][..singular_value_decomposition] method.
        It uses [`numpy.linalg.svd`][numpy.linalg.svd] function to compute the singular value decomposition and
        [`mpmath.identify`][mpmath.identify] function to identify rational numbers or surds if requested.

        Note:
            - This method might not return exact values, even if identification is enabled as it does
                not use SymPy's symbolic computation for SVD.
            - Use this method when performance is a concern and exact symbolic results are not required.

        Args:
            option (Literal["np", "sym"], optional): Whether to return numpy arrays or sympy matrices.
            identify (bool, optional): Whether to attempt identification of rational numbers or surds.
                If `True`, `option` must be `"sym"` to return symbolic matrices.
            tol (float, optional): Tolerance for [`mpmath.identify`][mpmath.identify] function as well as
                for verifying the SVD result.

        Returns:
            (SVD): A dataclass containing:

                - `U` ([`Matrix`][...]): The left singular vectors.
                - `S` ([`Matrix`][...]): The diagonal matrix of singular values.
                - `V` ([`Matrix`][...]): The right singular vectors.

            (NumSVD): A named tuple containing:

                - `U` ([`numpy.ndarray`][numpy.ndarray]): The left singular vectors as a NumPy array.
                - `S` ([`numpy.ndarray`][numpy.ndarray]): The diagonal matrix of singular values as a NumPy array.
                - `V` ([`numpy.ndarray`][numpy.ndarray]): The right singular vectors as a NumPy array.

        Examples:
            >>> mat = Matrix([[3, 2, 2], [2, 3, -2]])
            >>> U, S, V = mat.fast_svd(option="sym", identify=False)
            SVD(U=Matrix([
            [-0.707106781186548, -0.707106781186548]
            [-0.707106781186547,  0.707106781186548]
            ]), S=Matrix([
            [5.0, 0.0, 0.0]
            [0.0, 3.0, 0.0]
            ]), V=Matrix([
            [   -0.707106781186548, -0.235702260395516, -0.666666666666667]
            [   -0.707106781186548,  0.235702260395516,  0.666666666666667]
            [-6.47932334256779e-17, -0.942809041582063,  0.333333333333333]
            ]))

        See Also:
            - [`singular_value_decomposition`][..singular_value_decomposition] for the
                symbolic version of SVD.
            - NumPy's [`numpy.linalg.svd`][numpy.linalg.svd] for the underlying numerical
                SVD implementation.
        """

        m, n = self.shape
        U, S, Vh = np.linalg.svd(np.array(self, dtype=np.float64))
        # To align with MA1522 Syllabus, return V instead of V.T
        # Need not use conjugate transpose as MA1522 deals with real matrices
        V = Vh.T

        # Create sigma matrix from singular values
        S = np.diag(S)
        r, c = S.shape
        S = np.concat((S, np.zeros((r, n - c))), axis=1)
        S = np.concat((S, np.zeros((m - r, n))), axis=0)
        if option == "np":
            return NumSVD(U, S, V)
        elif option == "sym":
            U, S, V = Matrix(U), Matrix(S), Matrix(V)
            if identify:
                U = U.identify(tol=tol)
                S = S.identify(tol=tol)
                V = V.identify(tol=tol)
                residues = (self - U @ S @ V.T).norm()
                if residues > tol:
                    res = residues.evalf()
                    warn(
                        f"Non-zero Identification Error: {res}",
                        UserWarning,
                    )
                return SVD(U, S, V)
            else:
                return SVD(U, S, V)
        else:
            warn(
                f"Invalid option '{option}'. Expected 'np' or 'sym'. Returning NumSVD.",
                UserWarning,
            )
            return NumSVD(U, S, V)

    ####################################
    # CHAPTER 7: LINEAR TRANSFORMATION #
    ####################################

    def standard_matrix(self, out: Matrix, matrices: int = 1) -> Matrix | PartGen:
        """Returns the standard matrix for the transformation from self to out.

        The standard matrix is a matrix `T` such that `T @ self = out`, where `self` is the matrix
        whos columns represent the input vectors and `out` is the matrix whose columns represent the output vectors.

        Note:
            - The standard matrix may not be unique if the transformation is not injective.
            - If multiple solutions are found, the first solution is returned.

        Args:
            out (Matrix): The target matrix for the transformation.
            matrices (int): The number of matrices to return:

                - 1: Returns the standard matrix.
                - 2: Returns a [`PartGen`][(p).] with the part solution and general solution.

        Returns:
            (Matrix): If `matrices=1`, returns the standard matrix for the transformation.
            (PartGen): If `matrices=2`, returns a [`PartGen`][(p).] with the part solution and general solution.

        Examples:
            >>> input = Matrix([[1, 0, 1], [2, -1, 0], [0, 3, 1]])
            >>> output = Matrix([[4, 2, 3], [5, -1, 0], [1, 4, 2]])
            >>> input.standard_matrix(output)
            Matrix([
            [   2,    1,   1]
            [-3/5, 14/5, 3/5]
            [ 3/5,  1/5, 7/5]
            ])
        """

        X = Matrix.create_unk_matrix(r=out.rows, c=self.rows)
        sol = sym.solve(X @ self - out, X.free_symbols)
        if isinstance(sol, list) and len(sol) > 0:
            print("Multiple solutions found:")
            display(X.subs(s) for s in sol)
            print("Returning the first solution as the standard matrix.")
            X = X.subs(sol[0])
        elif isinstance(sol, dict):
            X = X.subs(sol)

        if matrices == 1:
            return X
        else:
            return X.sep_part_gen()
