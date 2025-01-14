import sympy as sym
import numpy as np
import mpmath as mp
from latex2sympy2 import latex2sympy
import re

import IPython.display
from itertools import chain, combinations
from collections import defaultdict
from typing import List, Literal, Union


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


def is_IPython():
    # Adapted from https://stackoverflow.com/questions/15411967/how-can-i-check-if-code-is-executed-in-the-ipython-notebook
    try:
        shell = get_ipython().__class__.__name__
        if shell == "ZMQInteractiveShell" or shell == "TerminalInteractiveShell":
            return True  # Jupyter notebook, qtconsole or terminal running IPython
        else:
            return False  # Other type
    except NameError:
        return False  # Probably standard Python interpreter

def aug_print(matrix: "Matrix") -> None:
    # get latex representation of matrix
    raw = sym.latex(matrix, mat_str="array")

    # create formatting string s to insert augment line visually
    ls = [pos for pos in matrix.aug_pos if 0 < pos < matrix.cols]
    ls.sort()
    delta = [ls[0]]
    delta.extend([ls[i] - ls[i-1] for i in range(1, len(ls))])
    remainder = matrix.cols - sum(delta)
    delta.append(remainder)
    s = '{' + '|'.join(['c' * i for i in delta]) + '}'
    default_s = '{' + 'c' * matrix.cols + '}'

    formatted = raw.replace(default_s, s)
    display(IPython.display.Math(formatted))

def display(input) -> None:
    if is_IPython():
        if hasattr(input, "aug_pos") and len(input.aug_pos) > 0:
            aug_print(input)
        else:
            IPython.display.display(input)
    else:
        sym.pprint(input)            

sym.init_printing(use_unicode=True)
np.set_printoptions(formatter={"float": lambda x: f"{x:10.7g}"})

####################
# CUSTOM FUNCTIONS #
####################

def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return list(chain.from_iterable(combinations(s, r) for r in range(len(s) + 1)))

def is_zero(expr) -> bool:
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

class Matrix(sym.MutableDenseMatrix):
    def __init__(self, matrix) -> None:
        super().__init__()
        self.aug_pos = set()

    def copy(self) -> "Matrix":
        new_mat = super().copy()
        new_mat.aug_pos = self.aug_pos if hasattr(self, "aug_pos") else set()
        return new_mat

    ###########################
    # MISCELLANEOUS FUNCTIONS #
    ###########################
    def from_latex(expr: str, row_join=True, norm=False) -> "Matrix":
        """
        Converts a LaTeX expression into a SymPy Matrix object.

        This function parses a LaTeX matrix or vector expression and returns
        a Matrix object. If the expression represents a list of vectors,
        the function will treat them as column vectors by default based on
        the `row_join` flag. Additionally, if `norm = True`, the vectors will be normalized.

        Parameters:
        - expr (str): The raw LaTeX string representing the matrix or vector.
        - row_join (bool, optional): If True (default), the matrix will be treated
          as a list of column vectors. If False, it will treat them as row vectors.
        - norm (bool, optional): If True, normalizes the column vectors to unit length.
          Default is False.

        Returns:
        - Matrix: A custom Matrix object (subclass of sympy.Matrix) representing the parsed LaTeX expression.

        Example:
        >>> Matrix.from_latex('\\begin{pmatrix} 1 & 2 \\\\ 3 & 4 \\end{pmatrix}')
        Matrix([[1, 2], [3, 4]])

        >>> Matrix.from_latex('\\begin{pmatrix} 1 \\\\ 2 \\\\ 3 \\end{pmatrix}', norm=True)
        Matrix([[1/√14], [2/√14], [3/√14]])
        """

        # Step 1: Modify the LaTeX string to ensure compatibility with the parser.
        # Convert array-like LaTeX to pmatrix for proper matrix formatting
        # Replace \begin{array}{ccc*} with \begin{pmatrix}
        modified_latex = re.sub(r"\\begin\{array\}\{[^}]*\}", r"\\begin{pmatrix}", expr)
        # Replace \end{array} with \end{pmatrix}
        modified_latex = re.sub(r"\\end\{array\}", r"\\end{pmatrix}", modified_latex)
        # Remove LaTeX semicolon for cleaner parsing
        modified_latex = re.sub(r"\\;", "", modified_latex)

        # Step 2: Use latex2sympy to parse the modified LaTeX expression into SymPy Matrix
        res = latex2sympy(modified_latex)
        display(res)

        # Step 3: Handle the parsed result based on its type (list, MatMul, or Matrix)
        if isinstance(res, list):
            vector_list = []
            for vector in res:
                vector = vector.expand()
                if norm:
                    vector /= vector.norm()
                vector_list.append(vector)
            return Matrix.from_list(vector_list, row_join)

        elif isinstance(res, sym.matrices.expressions.matmul.MatMul):
            # If the matrix is a product of matrices, evaluate the product directly
            return Matrix(res.doit())
        elif isinstance(res, sym.Matrix):
            # Directly converts the SymPy Matrix into the custom Matrix object to inherit the custom methods
            return Matrix(res)
        else:
            # If the result is neither a list nor a matrix expression, return the raw result
            return res

    def from_list(vectors: list["Matrix"], row_join: bool = True) -> "Matrix":
        """
        Creates a Matrix object from a list of vectors.

        This method takes a list of vectors (each represented as a Matrix object)
        and combines them into a single Matrix. If `row_join=True`, the vectors
        will be treated as rows and the resulting matrix will be transposed.
        If `row_join=False`, the vectors will be treated as columns in the matrix.

        Parameters:
        - vectors (list["Matrix"]): A list of Matrix objects, where each Matrix
          represents a row or column vector.
        - row_join (bool, optional): If True (default), the vectors are treated as columns,
          and the resulting matrix is transposed. If False, the vectors are treated as
          rows. Default is True.

        Returns:
        - Matrix: A custom Matrix object (subclass of sympy.Matrix) constructed from
          the list of vectors.

        Example:
        >>> vec1 = Matrix([[1], [2]])
        >>> vec2 = Matrix([[3], [4]])
        >>> Matrix.from_list([vec1, vec2])
        Matrix([[1, 3], [2, 4]])

        >>> Matrix.from_list([vec1, vec2], row_join=False)
        Matrix([[1, 2], [3, 4]])
        """

        res = Matrix([list(vector) for vector in vectors])
        if row_join:
            # If `row_join=True`, transpose the matrix (vectors treated as rows)
            return res.transpose()
        else:
            # If `row_join=False`, return the matrix as-is (vectors treated as columns)
            return res

    def create_unk_matrix(
        num_rows: int = 1, num_cols: int = 1, symbol: str = "x", is_real: bool = True
    ) -> "Matrix":
        """
        Creates a symbolic matrix with unknown entries.

        This method generates a matrix of size `num_rows` x `num_cols` with symbolic
        entries. The entries are named based on the provided `symbol` parameter and
        indexed by their row and column positions. The `is_real` flag determines whether
        the symbols are real-valued.

        Parameters:
        - num_rows (int, optional): The number of rows in the matrix. Default is 1.
        - num_cols (int, optional): The number of columns in the matrix. Default is 1.
        - symbol (str, optional): The base name for the symbols used in the matrix entries.
          Default is 'x'.
        - is_real (bool, optional): If True (default), the symbols are real-valued;
          otherwise, they are complex.

        Returns:
        - Matrix: A custom Matrix object (subclass of sympy.Matrix) with symbolic entries
          of the specified size.

        Example:
        >>> Matrix.create_unk_matrix(2, 2, symbol='a')
        Matrix([[a_(1, 1), a_(1, 2)], [a_(2, 1), a_(2, 2)]])

        >>> Matrix.create_unk_matrix(3, 1, symbol='y', is_real=False)
        Matrix([[y_(1, 1)], [y_(2, 1)], [y_(3, 1)]])
        """

        # Creates a matrix of size rows * cols with entries symbol_i,j
        entries = sym.symbols(
            f"{symbol}_(1:{num_rows+1})\\,(1:{num_cols+1})", is_real=is_real
        )
        return Matrix(entries).reshape(num_rows, num_cols)

    def create_rand_matrix(num_rows: int = 1, num_cols: int = 1) -> "Matrix":
        """
        Creates a matrix with random entries.

        This method generates a matrix of size `num_rows` x `num_cols` where the
        entries are randomly chosen. The values in the matrix are generated
        using SymPy's `randMatrix` function.

        Parameters:
        - num_rows (int, optional): The number of rows in the matrix. Default is 1.
        - num_cols (int, optional): The number of columns in the matrix. Default is 1.

        Returns:
        - Matrix: A custom Matrix object (subclass of sympy.Matrix) with random entries.

        Example:
        >>> Matrix.create_rand_matrix(2, 3)
        Matrix([[a random value, a random value, a random value],
                [a random value, a random value, a random value]])

        Notes:
        - The entries in the matrix are generated randomly and will change each time
          the function is called.
        """
        return Matrix(sym.randMatrix(num_rows, num_cols))

    def simplify(
        self,
        rational: bool = True,
        tolerance: float = 1e-4,
        simplify: bool = True,
        expand: bool = True,
        collect_sym: sym.Symbol = None,
        *args,
        **kwargs,
    ) -> "Matrix":
        """
        Simplifies the matrix by applying various simplification techniques.

        This method performs several operations on the matrix to simplify its entries:
        - Rational simplification using `sym.nsimplify`.
        - General symbolic simplification using `symsimplify`.
        - Expansion of expressions using `sym.expand` or factoring with `sym.factor`.
        - Collecting terms involving a specific symbol (if provided).

        Parameters:
        - rational (bool, optional): If True (default), applies rational simplification
          to the matrix entries using `sym.nsimplify`.
        - tolerance (float, optional): The tolerance for rational simplification (default is 1E-4).
        - simplify (bool, optional): If True (default), applies general symbolic simplification
          using `sym.simplify`.
        - expand (bool, optional): If True (default), applies expansion to the matrix entries.
          If False, applies factoring instead.
        - collect_sym (sym.Symbol, optional): A symbol to collect terms with. If provided,
          `sym.collect` will be applied to all entries of the matrix with respect to this symbol.
        - *args, **kwargs: Additional arguments passed to the `simplify` function.

        Returns:
        - Matrix: A new simplified matrix with the applied operations.

        Example:
        >>> mat = Matrix([[sym.symbols('x') + 1, sym.symbols('x') + 2], [sym.symbols('x') + 3, sym.symbols('x') + 4]])
        >>> mat.simplify(rational=False, expand=True)
        Matrix([[x + 1, x + 2], [x + 3, x + 4]])

        Notes:
        - Rational simplification attempts to convert entries into rational numbers if possible.
          If there is a residue (i.e., if the matrix was changed during the transformation),
          a warning is printed with the approximation error.
        - Expansion and factoring can be controlled by the `expand` parameter.
        - The matrix is modified in place and returned.
        """

        temp = self.copy()
        if rational:
            temp = sym.nsimplify(temp, tolerance=tolerance, rational=True)
            residues = (temp - self).norm()
            if residues != 0:
                print(f"Non-zero Approximation Error: {residues.evalf()}")
                print("Rational approximation might have failed. Try lower tolerance.")
        if simplify:
            temp = sym.simplify(temp, *args, **kwargs)
        if expand:
            temp = sym.expand(temp)
        else:
            temp = sym.factor(temp)
        temp = temp.tolist()
        if collect_sym is not None:
            temp = [[sym.collect(item, collect_sym) for item in row] for row in temp]

        # Create a new Matrix object from the simplified list and update the original object
        temp = Matrix(temp)
        temp.aug_pos = self.aug_pos if hasattr(self, "aug_pos") else set()
        self.__dict__.update(temp.__dict__)
        return self

    def identify(self, tolerance: float = None) -> "Matrix":
        """
        Identifies the matrix by applying a transformation function to each entry.

        This method applies a transformation to each element of the matrix using
        the `identify` function from the `mp` module. If a `tolerance` is provided,
        it will be passed to the transformation function. After applying the transformation,
        the method checks if there is any residue (i.e., if the matrix has been modified).

        Parameters:
        - tolerance (float, optional): A tolerance value that is passed to the
          `identify` function. If None (default), no tolerance is applied.

        Returns:
        - Matrix: A new matrix that results from applying the transformation to
          each element of the original matrix.

        Example:
        >>> mat = Matrix([[1, 2], [3, 4]])
        >>> mat.identify(tolerance=1E-5)
        Matrix([[1, 2], [3, 4]])  # Depending on the behavior of mp.identify

        Notes:
        - The `identify` function is assumed to perform some operation that "identifies" or
          transforms the entries in the matrix.
        - If there is a residue (i.e., if the matrix was changed during the transformation),
          a warning is printed with the approximation error.
        """

        temp = self.applyfunc(lambda x: mp.identify(x, tol=tolerance))
        residues = (temp - self).norm()
        if residues != 0:
            print(f"Non-zero Identification Error: {residues.evalf()}")
        return temp

    def select_cols(self, *args: int) -> "Matrix":
        """
        Selects columns from the matrix based on the provided column indices.

        This method returns a new matrix consisting of the columns specified by the
        provided indices. The columns are selected from the original matrix, and the
        result is returned as a new matrix.

        Parameters:
        - *args (int): One or more column indices (0-based) to select from the matrix.

        Returns:
        - Matrix: A new matrix consisting of the selected columns, transposed to
          match the original row-column structure.

        Example:
        >>> mat = Matrix([[1, 2, 3], [4, 5, 6]])
        >>> mat.select_cols(0, 2)
        Matrix([[1, 3], [4, 6]])

        Notes:
        - The method uses `self.col(idx)` to extract columns and then transposes the result
          to maintain the correct row-column structure.
        """

        res = []
        for idx in args:
            res.append(list(self.col(idx)))
        return Matrix(res).transpose()

    def select_rows(self, *args: int) -> "Matrix":
        """
        Selects rows from the matrix based on the provided row indices.

        This method returns a new matrix consisting of the rows specified by the
        provided indices. The rows are selected from the original matrix, and the
        result is returned as a new matrix.

        Parameters:
        - *args (int): One or more row indices (0-based) to select from the matrix.

        Returns:
        - Matrix: A new matrix consisting of the selected rows.

        Example:
        >>> mat = Matrix([[1, 2, 3], [4, 5, 6]])
        >>> mat.select_rows(0)
        Matrix([[1, 2, 3]])

        Notes:
        - The method uses `self.row(idx)` to extract rows directly and returns them in a new matrix.
        """

        res = []
        for idx in args:
            res.append(list(self.row(idx)))
        return Matrix(res)

    def sep_part_gen(self) -> tuple["Matrix", "Matrix"]:
        """
        Separates a matrix into its particular and general solution parts.

        This method separates the matrix into two components:
        - The **particular solution**, which is the solution to the system when all free variables are set to zero.
        - The **general solution**, which is the full solution including the homogeneous part.

        It assumes that the matrix is in symbolic form and contains free variables that can be set to zero.

        Returns:
            tuple: A tuple containing:
                - `part_sol` (Matrix): The particular solution (with free variables set to zero).
                - `gen_sol` (Matrix): The general solution (the original matrix minus the particular solution).
        
        Example:
            >>> mat = Matrix([[x, y], [x + y, 2*y]])
            >>> part_sol, gen_sol = mat.sep_part_gen()
            >>> part_sol
            Matrix([[0, 0], [0, 0]])  # Particular solution
            >>> gen_sol
            Matrix([[x, y], [x + y, 2*y]])  # General solution
        """

        set_0 = dict([(symbol, 0) for symbol in self.free_symbols])
        part_sol = self.subs(set_0)
        gen_sol = self - part_sol
        return part_sol, gen_sol

    def scalar_factor(self, column=True) -> tuple["Matrix", "Matrix"]:
        """
        Factorizes a matrix into the form A = BD, where D is a diagonal matrix and B contains the vectors
        with common divisors factored out (if `column=True`). If `column=False`, then returns A = DB instead.
        
        Args:
            column (bool): If `True`, factorizes by columns (default is `True`). If `False`, factorizes by rows.
            
        Returns:
            tuple["Matrix", "Matrix"]: A tuple of two matrices (B, D) such that A = BD, 
            where D is diagonal and B contains the scaled vectors.
        
        Example:
            >>> mat = Matrix([[6, 9], [12, 15]])
            >>> B, D = mat.scalar_factor(column=True)
            >>> B, D
            (Matrix([[1, 1.5], [2, 2.5]]), Matrix([6, 6]))
        """

        scalars = []
        B = self.copy()
        if column:
            for i in range(self.cols):
                g = sym.gcd(tuple(self.col(i)))
                B[:, i] /= g
                scalars.append(g)
            D = Matrix.diag(scalars)
            assert self == B @ D
            return B, D
        else:
            for i in range(self.rows):
                g = sym.gcd(tuple(self.row(i)))
                B[i, :] /= g
                scalars.append(g)
            D = Matrix.diag(scalars)
            assert self == D @ B
            return D, B

    #############################
    # CHAPTER 1: LINEAR SYSTEMS #
    #############################
    def aug_line(self, pos: int = -1) -> "Matrix":
        """
        Inserts an augmented line at the specified position.

        This method adds an augmented line (i.e., a line for augmented matrix operations)
        to the matrix at the specified column position. If no position is provided (default: -1),
        the line is inserted at the last column.

        Parameters:
        - pos (int, optional): The position (column index) where the augmented line
          will be inserted. Default is -1, which means the augmented line is added
          at the end of the matrix.

        Returns:
        - Matrix: The current matrix with the augmented line added at the specified position.

        Example:
        >>> mat = Matrix([[1, 2], [3, 4]])
        >>> mat.aug_line(1)
        Matrix([[1, 2], [3, 4], [0, 0]])

        Notes:
        - The method updates the `aug_pos` attribute to track the position of the inserted line.
        """

        new_pos = pos
        if new_pos < 0:
            new_pos += self.cols + 1

        if not 0 <= new_pos <= self.cols:
            raise IndexError(f'Position for augmented line ({pos}) out of range ({self.cols}).')

        self.aug_pos.add(new_pos)
        return self
    
    def row_join(self, other: "Matrix") -> "Matrix":
        offset = self.cols
        new_aug_pos = self.aug_pos if hasattr(self, "aug_pos") else {}
        other_aug_pos = other.aug_pos if hasattr(other, "aug_pos") else {}
        for pos in other_aug_pos:
            new_aug_pos.add(pos + offset)
        new_mat = Matrix(super().row_join(other))
        new_mat.aug_pos = new_aug_pos
        return new_mat

    def scale_row(self, idx: int, scalar: float, verbosity: int = 0) -> "Matrix":
        """
        Scales a row of the matrix by a scalar and simplifies the result.

        This method scales a specified row of the matrix by multiplying it with a scalar
        and then simplifies the matrix. The result is stored back in the matrix. Optionally,
        the method can print information about the row scaling and display the matrix,
        depending on the verbosity level.

        Parameters:
        - idx (int): The index of the row to scale (0-based).
        - scalar (float): The scalar by which to multiply the row.
        - verbosity (int, optional): The level of verbosity for output.
          - 0: No output.
          - 1: Print the row scaling operation.
          - 2: Print the row scaling operation and display the matrix.
          Default is 0 (no output).

        Returns:
        - Matrix: The modified matrix with the scaled row.

        Example:
        >>> mat = Matrix([[1, 2], [3, 4]])
        >>> mat.scale_row(0, 2)
        Matrix([[2, 4], [3, 4]])

        Notes:
        - The method modifies the matrix in-place and returns the updated matrix.
        - After scaling the row, the matrix is simplified using `self.simplify()`.
        """

        self[idx, :] = scalar * self[idx, :]
        self.simplify()

        if verbosity >= 1:
            print(f"R_{idx+1} <- ({scalar})R_{idx+1}")
        if verbosity >= 2:
            display(self)
            print("\n")

        return self

    def swap_row(self, idx_1: int, idx_2: int, verbosity: int = 0) -> "Matrix":
        """
        Swaps two rows of the matrix.

        This method swaps the contents of two rows in the matrix. The operation is performed
        in-place, and the modified matrix is returned. Optionally, the method can print
        information about the row swap and display the matrix, depending on the verbosity level.

        Parameters:
        - idx_1 (int): The index of the first row to swap (0-based).
        - idx_2 (int): The index of the second row to swap (0-based).
        - verbosity (int, optional): The level of verbosity for output.
          - 0: No output.
          - 1: Print the row swap operation.
          - 2: Print the row swap operation and display the matrix.
          Default is 0 (no output).

        Returns:
        - Matrix: The modified matrix after the row swap.

        Example:
        >>> mat = Matrix([[1, 2], [3, 4]])
        >>> mat.swap_row(0, 1)
        Matrix([[3, 4], [1, 2]])

        Notes:
        - The method modifies the matrix in-place and returns the updated matrix.
        - After performing the row swaps, the matrix is simplified using `self.simplify()`.
        """

        self[idx_1, :], self[idx_2, :] = self[idx_2, :], self[idx_1, :]

        if verbosity >= 1:
            print(f"R_{idx_1+1} <-> R_{idx_2+1}")
        if verbosity >= 2:
            display(self)
            print("\n")

        return self

    def reduce_row(
        self, idx_1: int, scalar: float, idx_2: int, verbosity: int = 0
    ) -> "Matrix":
        """
        Reduces a row by subtracting a scalar multiple of another row.

        This method modifies a row by subtracting a specified scalar multiple of another row.
        The result is stored back in the matrix. Optionally, the method can print information
        about the row reduction and display the matrix, depending on the verbosity level.

        Parameters:
        - idx_1 (int): The index of the row to reduce (0-based).
        - scalar (float): The scalar by which to multiply the second row.
        - idx_2 (int): The index of the row from which to subtract the scalar multiple (0-based).
        - verbosity (int, optional): The level of verbosity for output.
          - 0: No output.
          - 1: Print the row reduction operation.
          - 2: Print the row reduction operation and display the matrix.
          Default is 0 (no output).

        Returns:
        - Matrix: The modified matrix after the row reduction.

        Example:
        >>> mat = Matrix([[1, 2], [3, 4]])
        >>> mat.reduce_row(0, 2, 1)
        Matrix([[-5, -6], [3, 4]])

        Notes:
        - The method modifies the matrix in-place and returns the updated matrix.
        - After performing the row reduction, the matrix is simplified using `self.simplify()`.
        """

        self[idx_1, :] = self[idx_1, :] - scalar * self[idx_2, :]
        self.simplify()

        if verbosity >= 1:
            print(f"R_{idx_1+1} <- R_{idx_1+1} - ({scalar})R_{idx_2+1}")
        if verbosity >= 2:
            display(self)
            print("\n")

        return self

    def get_pivot_row(
        self, col_idx: int, row_start_idx: int, follow_GE: bool = False
    ) -> int:
        """
        Finds the row index of the pivot element in a given column.

        This method attempts to find a row that contains a non-zero element in the
        specified column. If the `follow_GE` flag is `False`, it first looks for
        a non-zero constant that does not contain any symbolic expressions. If no
        such element is found, it will return the first non-zero element. If the
        entire column contains only zeros, the method returns -1.

        Parameters:
        - col_idx (int): The index of the column to search for the pivot.
        - row_start_idx (int): The row index to start searching from.
        - follow_GE (bool, optional): Flag to control whether to follow Gaussian elimination strategy.
          - `False`: First look for non-zero constants that are not symbolic expressions.
          - `True`: Always return the first non-zero element, even if it is symbolic.
          Default is `False`.

        Returns:
        - int: The index of the row containing the pivot element, or -1 if no pivot is found.

        Example:
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
                    if not isinstance(term, sym.Expr):
                        return row_idx
                    # Check if it's a non-symbolic constant
                    elif len(term.free_symbols) == 0:
                        return row_idx

        # Step 2: If no non-zero constant is found, return the first non-zero element (symbolic or not)
        for row_idx in range(row_start_idx, self.rows):
            term = self[row_idx, col_idx]
            if term != 0:
                return row_idx

        # Step 3: If no non-zero element is found, return -1 (indicating no pivot)
        return -1

    def get_pivot_pos(self) -> list[tuple[int, int]]:
        """
        Finds the positions of the pivot elements in the matrix.

        This method checks the matrix to determine the positions of the pivots
        (the first non-zero entry in each row) by examining each column one-by-one.
        It assumes that the matrix is in Row Echelon Form (REF), as checked by the
        `is_echelon` property.

        It uses `get_pivot_row` to find the pivot row for each column. For each
        pivot found, a tuple (row, column) is added to the result list.

        Returns:
            list[tuple[int, int]]: A list of lists, where each sublist contains a
                                   tuple representing the position (row, column) of a pivot.

        Example:
        >>> mat = Matrix([[1, 2, 3], [0, 4, 5], [0, 0, 6]])
        >>> mat.get_pivot_pos()
        [(0, 0), (1, 1), (2, 2)]
        """

        assert self.is_echelon  # check for REF

        pivot_pos = list[tuple[int, int]] = []
        cur_row_pos = 0
        for cur_col_pos in range(self.cols):
            pivot_row = self.get_pivot_row(cur_col_pos, cur_row_pos, follow_GE=False)

            if pivot_row != -1:
                pivot_pos.append((pivot_row, cur_col_pos))
                cur_row_pos += 1

        return pivot_pos

    def get_pivot_elements(self) -> list[sym.Expr]:
        """
        Retrieves the pivot elements from the matrix.

        This method identifies the pivot positions (row, column) using the
        `get_pivot_pos` method and then extracts the elements at those positions
        in the matrix.

        Returns:
            list[sym.Expr]: A list of pivot elements corresponding to the positions
                            identified by `get_pivot_pos`.

        Example:
        >>> mat = Matrix([[1, 2, 3], [0, 4, 5], [0, 0, 6]])
        >>> mat.get_pivot_elements()
        [1, 4, 6]
        """

        pivot_elements: list[sym.Expr] = []

        for i, j in self.get_pivot_pos():
            pivot_elements.append(self[i, j])

        return pivot_elements

    def ref(
        self,
        verbosity: int = 2,
        max_tries: int = 2,
        follow_GE: bool = False,
        matrices: int = 2,
    ) -> Union["Matrix", tuple["Matrix", "Matrix"], tuple["Matrix", "Matrix", "Matrix"]]:
        """
        Performs the Row Echelon Form (REF) transformation on the matrix.

        This method applies Gaussian elimination (or a similar approach) to bring
        the matrix to row echelon form. Optionally, it can also return the LU decomposition
        and/or permutation matrix used during the process.

        Parameters:
        - verbosity (int, optional): Level of verbosity for the output.
          - 0: No output.
          - 1: Output basic information (e.g., row operations).
          - 2: Output detailed information (e.g., matrix states after each operation).
          Default is 2.
        - max_tries (int, optional): Maximum number of tries to reduce a row in case of symbolic denominators.
          Default is 2.
        - follow_GE (bool, optional): Whether to strictly follow Gaussian elimination rules.
          Default is False.
        - matrices (int, optional): Number of matrices to return.
          - 1: Returns only the upper triangular matrix (U).
          - 2: Returns the standard matrix (A) and the upper matrix (U) such that the matrix equals to AU.
               A is not guaranteed to be a lower triangular matrix since not all matrices are LU-factorisable.
          - 3: Returns the permutation matrix (P), lower triangular matrix (L), and upper matrix (U).
          Default is 2.

        Returns:
        - Matrix | tuple[Matrix]: The matrix or matrices depending on the `matrices` parameter.

        Example:
        >>> mat = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        >>> mat.ref(verbosity=1, matrices=2)
        (PermutationMatrix, UpperTriangularMatrix)
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

            if pivot_row == -1:
                # If no non-zero pivot is found, continue to the next column
                continue

            # Swap the current row with the pivot row if necessary
            if pivot_row != cur_row_pos:
                U.swap_row(cur_row_pos, pivot_row, verbosity=verbosity)
                P_elem = I.copy().swap_row(cur_row_pos, pivot_row)
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
                        print(
                            f"ERROR: Max tries exceeded to reduce row {row_idx+1} with row {cur_row_pos+1}"
                        )
                        break
                    try:
                        scalar = U[row_idx, cur_col_pos] / U[cur_row_pos, cur_col_pos]
                        scalar = scalar.expand().simplify()

                        try:
                            decomp = sym.apart(scalar)  # partial fractions
                        except:
                            decomp = scalar
                        # if isinstance(decomp, sym.core.add.Add) and (scalar != decomp):
                        if isinstance(decomp, sym.core.add.Add):
                            terms = decomp.args
                        else:
                            # there is only 1 term (could be integer or proper fraction)
                            terms = [decomp]

                        for term in terms:
                            _, d = sym.fraction(term)

                            # ensure denominator is non-zero so that reduction is valid
                            if not is_zero(d):
                                U.reduce_row(
                                    row_idx, term, cur_row_pos, verbosity=verbosity
                                )
                                elem = I.copy().reduce_row(row_idx, -term, cur_row_pos)
                                L = L @ elem

                        # Cases where pivot row contains symbols such that scalar is a
                        # fraction with symbolic denominator.
                        # To reduce further, can only scale row_idx accordingly
                        if U[row_idx, cur_col_pos] != 0:
                            scalar = (
                                U[cur_row_pos, cur_col_pos] / U[row_idx, cur_col_pos]
                            )
                            scalar.simplify()
                            n, d = sym.fraction(scalar)
                            # to scale by n, n cannot be 0 both numerically or symbolically
                            # to scale by 1/d, d cannot be 0, same argument as n
                            # if (n != 0) and (not is_zero(d)):
                            if (not is_zero(n)) and (not is_zero(d)):
                                U.scale_row(row_idx, scalar, verbosity=verbosity)
                                elem = I.copy().scale_row(row_idx, 1 / scalar)
                                L = L @ elem

                    except Exception as error:
                        print(f"Exception encountered: {error}")
                        if matrices == 1:
                            return U
                        elif matrices == 2:
                            return P @ L, U
                        elif matrices == 3:
                            return P, L, U
                        else:
                            return U

            cur_row_pos += 1

        # Return the appropriate number of matrices based on the `matrices` parameter
        if matrices == 1:
            return U
        elif matrices == 2:
            return P @ L, U
        elif matrices == 3:
            return P, L, U
        else:
            return U

    def find_all_cases(self) -> List:
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
        for subset in powerset(cases):
            combined = dict()
            for d in subset:
                combined = {**combined, **d}
            combinations.add(tuple(sym.ordered(combined.items())))

        return list(sym.ordered(dict(combination) for combination in combinations))

    def evaluate_cases(self, rhs: "Matrix") -> None:
        cases = self.find_all_cases()
        all_possible_values = set(
            possible_val for case in cases for possible_val in case.items()
        )

        for i, case in enumerate(cases, 1):
            print(
                f"Case {i}: {case}, not including {[dict([val]) for val in all_possible_values.symmetric_difference(set(case.items()))]}"
            )
            U = (
                self.row_join(rhs)
                .subs(case)
                .ref(verbosity=0, matrices=1, follow_GE=False)
            )
            display(U)

    def rref(self, *args, **kwargs) -> tuple["Matrix", tuple[int]]:
        """
        Computes the Reduced Row Echelon Form (RREF) of the matrix.

        This method wraps around the `rref` method of the superclass and returns
        the matrix in Reduced Row Echelon Form (RREF) along with the pivot positions.

        Args:
            *args: Positional arguments passed to the superclass's `rref` method.
            **kwargs: Keyword arguments passed to the superclass's `rref` method.

        Returns:
            tuple["Matrix", tuple[int]]: A tuple where the first element is a
                                        `Matrix` object in RREF, and the second element
                                        is a tuple containing the pivot positions (col indices).

        Example:
            >>> mat = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
            >>> rref_matrix, pivots = mat.rref()
            >>> print(rref_matrix)  # Reduced Row Echelon Form of the matrix
            >>> print(pivots)       # Indices of the pivot columns
        """

        rref_mat, pivots = super().rref(*args, **kwargs)
        rref_mat = Matrix(rref_mat)
        if hasattr(self, "aug_pos") and len(self.aug_pos) > 0:
            rref_mat.aug_pos = self.aug_pos
        return rref_mat, pivots

    def solve(self, rhs: "Matrix", *args, **kwargs) -> "Matrix":
        try:
            A, b = sym.Matrix(self), sym.Matrix(rhs)
            return A.solve(b, *args, **kwargs)
        except Exception as e:
            print(f"Exception Encountered: {str(e)}")
            print("Attempting custom solve...")
        x = Matrix.create_unk_matrix(num_rows=self.cols, num_cols=1)
        return x.subs(sym.solve(A@x - b, x.free_symbols))

    #############################
    # CHAPTER 2: MATRIX ALGEBRA #
    #############################

    def inverse(
        self, option: Union[Literal['left', 'right'], None], matrices: int = 1, verbosity: int = 0
    ) -> tuple["Matrix", "Matrix"]:
        """
        Computes the left or right inverse of a matrix, depending on its rank and the specified option.

        The method checks whether the matrix has full row rank or full column rank and computes either:
        - The **left inverse** (if the matrix has full column rank).
        - The **right inverse** (if the matrix has full row rank).

        If neither option is provided, the method automatically determines which inverse to compute based on the matrix's rank.

        Args:
            option (str, optional): Specifies which inverse to compute:
                - 'left' for the left inverse (requires the matrix to have full column rank).
                - 'right' for the right inverse (requires the matrix to have full row rank).
                Default is None, in which case the method attempts to determine the inverse type based on the matrix rank.
            matrices (int, optional): Specifies the number of matrices to return:
                - 1: Returns only the inverse matrix.
                - 2: Returns the particular and general solutions of the inverse.
            verbosity (int, optional): Level of verbosity for displaying intermediate steps:
                - 0: No output.
                - 1: Display matrices before and after RREF.
                Default is 0.

        Returns:
            tuple: A tuple containing:
                - The inverse matrix (if found).
                - If `matrices == 2`, it also returns the particular and general solutions of the inverse.

        Raises:
            Exception: If no valid inverse (left or right) is found, an exception is raised.

        Example:
            >>> mat = Matrix([[1, 2], [3, 4]])
            >>> mat.inverse(option="left")
            Matrix([[1, 0], [0, 1]])
        """

        if option is None:
            if self.rank() == self.cols:
                if verbosity:
                    print("Left inverse found!")
                option = "left"
            if self.rank() == self.rows:
                if verbosity:
                    print("Right inverse found!")
                option = "right"

        if (option is not None) and (verbosity >= 1):
            if option == "left":
                aug = self.T.aug_line().row_join(sym.eye(self.cols))
                print("Before RREF: [self^T | eye]")
                display(aug)
                print("\nAfter RREF:")
                display(aug.rref()[0])
            else:
                aug = self.aug_line().row_join(sym.eye(self.rows))
                print("Before RREF: [self | eye]")
                display(aug)
                print("\nAfter RREF:")
                display(aug.rref()[0])

        if option is not None:
            X = Matrix.create_unk_matrix(
                num_rows=self.cols, num_cols=self.rows, symbol="x"
            )
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
                print(f"No {option} inverse found! Try pseudo-inverse: .pinv()")
                return None

            if matrices == 1:
                return X
            else:
                return X.sep_part_gen()

    def elem(self) -> "Matrix":
        """
        Returns the identity matrix of the same size as the current matrix.

        `A = IA`

        This method creates and returns an identity matrix with the same number
        of rows as the current matrix. The identity matrix has ones on the diagonal 
        and zeros elsewhere.

        Returns:
        - Matrix: A custom Matrix object (subclass of sympy.Matrix) that represents
          the identity matrix of the same size as the current matrix.

        Example:
        >>> mat = Matrix([[1, 2], [3, 4]])
        >>> mat.elem()
        Matrix([[1, 0], [0, 1]])
        """
        elem_mat = Matrix(sym.eye(self.rows))
        elem_mat.aug_pos = set()
        return elem_mat

    def adjoint(self) -> "Matrix":
        """
        Computes the adjugate (classical adjoint) of the matrix.

        This method calculates the classical adjoint (also known as the adjugate in SymPy) of the matrix.
        The adjoint of a matrix is the transpose of its cofactor matrix. This is an alias
        for the classical adjoint, which aligns with the definition used in the MA1522 syllabus.

        Returns:
            Matrix: The classical adjoint (or adjugate) matrix of the current matrix.

        Example:
            >>> mat = Matrix([[1, 2], [3, 4]])
            >>> mat.adjoint()
            Matrix([[ -4,  2], [ 3, -1]])
        """
        print(
            "Warning: The classical adjoint of the matrix is computed rather than the conjugate transpose."
        )
        print("Please use self.adj() instead to remove ambiguity.")
        return self.adjugate()

    def adj(self) -> "Matrix":
        """
        Alias for the adjoint method.

        This method is an alias for the `adjoint` method. It returns the classical adjoint (or adjugate)
        of the matrix.

        Returns:
            Matrix: The classical adjoint (or adjugate) matrix of the current matrix.

        Example:
            >>> mat = Matrix([[1, 2], [3, 4]])
            >>> mat.adj()
            Matrix([[ -4,  2], [ 3, -1]])
        """
        return self.adjugate()

    def column_constraints(self, use_ref: bool = False, **kwargs) -> "Matrix":
        """
        Computes the column constraints for the matrix by appending a symbolic vector.

        This method creates a matrix where a random column vector (with variables `x_1, x_2, ..., x_m`)
        is added to the matrix as an additional column. It then constructs a larger augmented matrix
        and optionally computes its Row Echelon Form (REF) or Reduced Row Echelon Form (RREF).

        The method modifies the matrix to ensure that the unknown vector is not reduced in RREF,
        and the constraints for the matrix columns are calculated accordingly.

        Args:
            use_ref (bool, optional): Whether to use Row Echelon Form (REF) instead of Reduced Row Echelon Form (RREF).
                                      Defaults to False, in which case RREF will be used.
            **kwargs: Additional keyword arguments passed to the `ref` or `rref` method.

        Returns:
            Matrix: A new matrix containing the result after applying REF or RREF to the augmented matrix.

        Example:
            >>> mat = Matrix([[1, 2], [3, 4], [5, 6]])
            >>> mat.column_constraints(use_ref=True)
            Matrix([[1, 2, x_1], [0, 0, x_2], [0, 0, x_3]])
        """

        # write a random vector as x_1, ..., x_m, given m rows
        vector = Matrix.create_unk_matrix(self.rows, 1, "x")

        # insert hidden column vectors so that the unknown vector is not reduced in rref
        hidden = self.elem()

        M = self.copy().row_join(hidden).row_join(vector)
        if use_ref:
            res = M.ref(**kwargs, matrices=1)
        else:
            res, _ = M.rref()
        res_matrix = res[:, : self.cols].row_join(res[:, -1])
        # might try return res_matrix instead
        return Matrix(res_matrix)

    ######################################
    # CHAPTER 3: EUCLIDEAN VECTOR SPACES #
    ######################################

    def normalized(self, factor: bool = False) -> "Matrix":
        """
        Normalizes the column vectors of the matrix (scaling each vector to have a unit norm).

        Optionally, returns the matrix with the column scaling factors if `factor=True`.

        Args:
            factor (bool, optional): If `True`, returns the scaling factors used for normalization.
                Default is `False`, in which case the normalized matrix is returned.

        Returns:
            Matrix: The normalized matrix if `factor=False`, or a matrix of column scaling factors if `factor=True`.

        Example:
            >>> mat = Matrix([[3, 4], [1, 2]])
            >>> norm_mat = mat.normalized()
            >>> norm_mat
            Matrix([[0.6, 0.8], [0.2, 0.4]])

            >>> norm_factors = mat.normalized(factor=True)
            >>> norm_factors
            Matrix([sqrt(1/5), sqrt(1/5)])
        """

        for i in range(self.cols):
            scalar = self.col(i).norm()
            if scalar != 0:
                self[:, i] /= scalar

        if factor:
            return self.scalar_factor(column=True)
        else:
            return self

    def is_linearly_independent(
        # TODO
        self, colspace: bool = True, verbosity: int = 0
        ) -> bool:
        rref_mat, pivots = self.rref()

        if verbosity == 1:
            print("rref(self)")
        elif verbosity >= 2:
            print("Before RREF: self")
            display(self)
            print("\nAfter RREF:")
            display(rref_mat)

        if colspace:
            if verbosity >= 1:
                print(f"Check if Number of columns ({self.cols}) == Number of pivot columns ({len(pivots)})")
            return self.cols == len(pivots)
        else:
            if verbosity >= 1:
                print(f"Check if Number of rows ({self.rows}) == Number of pivot columns ({len(pivots)})")
            return self.rows == len(pivots)
    
    def simplify_basis(self, colspace: bool = True, verbosity: int = 2) -> "Matrix":
        # TODO
        if colspace:
            rref_mat, _ = self.T.rref()
            if verbosity == 1:
                print("Select non-zero rows of rref(self.T) as basis vectors.")
            if verbosity >= 2:
                print("Before RREF: self^T")
                display(self)
                print("\nAfter RREF:")
                display(rref_mat)
        else:
            rref_mat, _ = self.rref()
            if verbosity == 1:
                print("Select non-zero rows of rref(self) as basis vectors.")
            if verbosity >= 2:
                print("Before RREF: self")
                display(self)
                print("\nAfter RREF:")
                display(rref_mat)

        idxs = []
        for i in range(rref_mat.rows):
            if rref_mat[i, :].norm() > 0:
                idxs.append(i)
        return rref_mat.select_rows(*idxs).T
    
    def extend_basis(
        self, span_subspace: "Matrix" = None, verbosity: int = 2
    ) -> "Matrix":
        """
        Extends the matrix to form a basis for the span of the given subspace.

        This method extends the current matrix to include the columns of the provided
        `span_subspace`, computes the Reduced Row Echelon Form (RREF) of the augmented matrix,
        and then selects the pivot columns to return the extended basis.

        If no `span_subspace` is provided, the identity matrix is used as the default.
        The result is a matrix with the extended basis that spans the combined space of the original
        matrix and the `span_subspace`.

        Args:
            span_subspace (Matrix, optional): A matrix whose columns represent the subspace to
                                              be added to the current matrix. Defaults to `None`, in
                                              which case the identity matrix is used.
            verbosity (int, optional): Verbosity level for displaying information.
                - 0: No output.
                - 1: Display steps.
                - 2: Display matrix before and after RREF.
                Defaults to 2.

        Returns:
            Matrix: A matrix representing the extended basis, consisting of the pivot columns
                    from the RREF of the augmented matrix.

        Example:
            >>> mat = Matrix([[1, 2], [3, 4]])
            >>> subspace = Matrix([[1, 0], [0, 1]])
            >>> mat.extend_basis(span_subspace=subspace)
            Matrix([[1, 2, 1, 0], [3, 4, 0, 1]])
        """

        if span_subspace is None:
            span_subspace = self.elem()
        aug = self.aug_line().row_join(span_subspace)
        rref_mat, pivots = aug.rref()

        if verbosity == 1:
            print("rref([self | span_subspace])")
        elif verbosity >= 2:
            print("Before RREF: [self | span_subspace]")
            display(aug)
            print("\nAfter RREF:")
            display(rref_mat)

        return aug.select_cols(*pivots)

    def intersect_subspace(self, other: "Matrix", verbosity: int = 2) -> "Matrix":
        """
        Computes the intersection of two subspaces by finding the nullspace of their orthogonal complements.

        This method computes the intersection of the subspaces spanned by the columns of the current matrix
        (`self`) and the provided matrix (`other`). The intersection is computed by finding the nullspace of
        the row space of the two matrices, which is equivalent to the nullspace of their orthogonal complements.

        Args:
            other (Matrix): The second matrix representing the other subspace to intersect with the current matrix.
            verbosity (int, optional): Level of verbosity for displaying intermediate steps:
                - 0: No output.
                - 1: Display steps.
                - 2: Display the relevant matrices.
                Defaults to 2.

        Returns:
            Matrix: A matrix whose columns form a basis for the intersection of the two subspaces.

        Example:
            >>> mat1 = Matrix([[1, 0], [0, 1]])
            >>> mat2 = Matrix([[1, 1], [0, 0]])
            >>> mat1.intersect_subspace(mat2)
            Matrix([[0], [0]])
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
            display(aug.rref()[0])

        return Matrix.from_list(aug.nullspace())

    def is_same_subspace(self, other: "Matrix", verbosity: int = 2) -> bool:
        """
        Checks if two subspaces are the same by verifying if each subspace is a subspace of the other.

        This method determines whether the subspaces spanned by the columns of the current matrix (`self`)
        and the provided matrix (`other`) are the same. It does so by checking if the row-reduced echelon forms
        (RREF) of the augmented matrices `[self | other]` and `[other | self]` reveal that both subspaces
        span each other.

        Args:
            other (Matrix): The second matrix representing the subspace to compare with the current matrix.
            verbosity (int, optional): Level of verbosity for displaying intermediate steps:
                - 0: No output.
                - 1: Display the steps.
                - 2: Display the relevant matrices.
                Defaults to 2.

        Returns:
            bool: `True` if the subspaces spanned by `self` and `other` are the same, `False` otherwise.

        Raises:
            ShapeError: If the number of rows in the current matrix and the target matrix are different.
        
        Example:
            >>> mat1 = Matrix([[1, 0], [0, 1]])
            >>> mat2 = Matrix([[1, 0], [0, 1]])
            >>> mat1.is_same_subspace(mat2)
            True
        """

        if self.rows != other.rows:
            raise sym.ShapeError(f"The matrices have incompatible number of rows ({self.rows}, {other.rows})")

        sub_rref, sub_pivots = other.aug_line().row_join(self).rref()
        super_rref, super_pivots = self.aug_line().row_join(other).rref()

        if verbosity == 1:
            print('Check rref([other | self])')
            print('Check rref([self | other])')
        elif verbosity >= 2:
            print("Check if span(self) is subspace of span(other)")
            print("\nBefore RREF: [other | self]")
            display(other.aug_line().row_join(self))
            print("\nAfter RREF:")
            display(sub_rref)

            print("\nCheck if span(other) is subspace of span(self)")
            print("\nBefore RREF: [self | other]")
            display(self.aug_line().row_join(other))
            print("\nAfter RREF:")
            display(super_rref)

        return (max(sub_pivots) < other.cols) and (max(super_pivots) < self.cols)

    def coords_relative(self, to: "Matrix", verbosity: int = 2) -> "Matrix":
        # TODO

        if self.cols != 1:
            raise sym.ShapeError(f"self should be a vector with 1 column. ({self.cols})")
        if self.rows != to.rows:
            raise sym.ShapeError(f"The matrices have incompatible number of rows ({self.rows}, {to.rows})")

        M = to.aug_line().row_join(self)
        rref_mat, pivots = M.rref()
        if verbosity == 1:
            print("Solve system via rref([to | self])")
        elif verbosity >= 2:
            print("Before RREF: [to | self]")
            display(M)
            print("\nAfter RREF:")
            display(rref_mat)

        if to.cols in pivots:
            raise Exception("No solution found due to inconsistent system.")
        return to.solve(self)
    
    def transition_matrix(self, to: "Matrix", verbosity: int = 2) -> "Matrix":
        """
        Computes the transition matrix that transforms this matrix to another matrix.

        This method computes the transition matrix `P` such that:
        ```
        to = P * self
        ```
        where `to` is the target matrix, and `self` is the current matrix. The method
        achieves this by augmenting the target matrix with the current matrix, performing
        Reduced Row Echelon Form (RREF), and extracting the appropriate part of the resulting matrix.

        Args:
            to (Matrix): The matrix to which the current matrix should be transformed.
            verbosity (int, optional): Verbosity level for displaying information.
                - 0: No output.
                - 1: Display the steps.
                - 2: Display the relevant matrices.
                Defaults to 2.

        Returns:
            Matrix: The transition matrix `P` that satisfies `to = P * self`.

        Raises:
            AssertionError: If the columns of the current matrix and to matrix do not span the same subspace.

        Example:
            >>> mat1 = Matrix([[1, 0], [0, 1]])
            >>> mat2 = Matrix([[2, 0], [0, 2]])
            >>> mat1.transition_matrix(to=mat2)
            Matrix([[2, 0], [0, 2]])
        """
        assert self.is_same_subspace(
            to, verbosity=0
        ), "Column vectors of both matrices must span the same subspace."

        M = to.aug_line().row_join(self)
        res = M.rref()[0]
        if verbosity == 1:
            print("rref([to | self])")
        elif verbosity >= 2:
            print("Before RREF: [to | self]")
            display(M)
            print("\nAfter RREF:")
            display(res)
        P = res[: self.cols, self.cols :]
        return P

    ###############################################
    # CHAPTER 4: SUBSPACES ASSOCIATED TO A MATRIX #
    ###############################################

    def nullspace(self, verbosity: int = 0, *args, **kwargs) -> list["Matrix"]:
        """
        Computes the null space (kernel) of the matrix, i.e., the set of vectors that satisfy `Ax = 0`.
        
        This method utilizes the rank-nullity theorem to determine if the null space exists. Fixes the 
        issue with SymPy implementation of nullspace where it raises an exception if the nullspace does not exist. 
        If the matrix has full column rank (i.e., rank == number of columns), it has no non-trivial null space,
        and an empty list is returned.

        Args:
            verbosity (int, optional): Level of verbosity for displaying intermediate steps.
                - 0: No output.
                - 1: Display the matrix before and after row-reduction (RREF).
                Default is 0.
            *args: Additional arguments passed to the `nullspace` method of the superclass.
            **kwargs: Additional keyword arguments passed to the `nullspace` method of the superclass.

        Returns:
            list[Matrix]: A list of `Matrix` objects representing the null space vectors.

        Example:
            >>> mat = Matrix([[1, 2], [3, 6]])
            >>> nullspace = mat.nullspace()
            >>> nullspace
            [Matrix([[-2], [1]])]
        """

        # Issue with SymPy implementation of nullspace when there is None
        # Using rank nullity theorem to verify there are vectors spanning nullspace
        if verbosity >= 1:
            print("Before RREF: [self]")
            display(self)
            print("\nAfter RREF:")
            display(self.rref()[0])

        if self.rank() == self.cols:
            if verbosity >= 1:
                print("No nullspace detected!")
            return []
        else:
            return super().nullspace(*args, **kwargs)

    #######################################################
    # CHAPTER 5: ORTHOGONALITY AND LEAST SQUARES SOLUTION #
    #######################################################

    def orthogonal_complement(self, verbosity: int = 0) -> "Matrix":
        """
        Computes the orthogonal complement of the matrix. This is the null space of its transpose.

        The orthogonal complement consists of all vectors that are orthogonal to the column space of the matrix.

        Args:
            verbosity (int, optional): Level of verbosity for debugging (default is 0, no output).

        Returns:
            Matrix: A `Matrix` object representing the orthogonal complement.

        Example:
            >>> mat = Matrix([[1, 0], [0, 1]])
            >>> ortho_comp = mat.orthogonal_complement()
            >>> ortho_comp
        """

        return Matrix.from_list(self.transpose().nullspace(verbosity))

    def is_vec_orthogonal(self, verbosity: int = 1) -> bool:
        """
        Checks if the column vectors of the matrix are orthogonal to each other.

        This method computes `self^T @ self` and checks if the result is diagonal.
        - If the result is diagonal, the vectors are orthogonal (i.e., `u_i . u_j = 0` for all `i != j`).
        - Note: This method checks for orthogonality, not orthonormality. For orthonormality, `self^T @ self = I`.

        Args:
            verbosity (int, optional): Level of verbosity for displaying intermediate results.
                - 0: No output.
                - 1: Display the matrix product `self^T @ self`.
                Default is 0.

        Returns:
            bool: `True` if the column vectors are orthogonal, `False` otherwise.

        Example:
            >>> mat = Matrix([[1, 0], [0, 1]])
            >>> mat.is_vec_orthogonal()
            True
        """

        res = self.T @ self
        if verbosity >= 1:
            print("self^T @ self")
            display(res)
        return res.is_diagonal()
    
    def is_mat_orthogonal(self, verbosity: int = 1) -> bool:
        # TODO
        res = self.T @ self
        if verbosity >= 1:
            print("self^T @ self")
            display(res)
        return res.is_diagonal and all(entry == 1 for entry in res.diagonal())
    
    def orthogonal_decomposition(
        self, to: "Matrix", verbosity: int = 0
        ) -> tuple["Matrix", "Matrix"]:
        #  TODO
        sol = to.solve_least_squares(self, verbosity=verbosity)
        proj = to @ sol
        norm = self - proj

        if verbosity >= 1:
            print("Projected component: Au")
            display(proj)
            print("Normal component: b - b_proj")
            display(norm)

        return proj, norm
    
    def proj_comp(
        self, to: "Matrix", verbosity: int = 0
        ) -> "Matrix":
        # TODO
        return self.orthogonal_decomposition(to=to, verbosity=verbosity)[0]
    
    def norm_comp(
        self, to: "Matrix", verbosity: int = 0
        ) -> "Matrix":
        # TODO
        return self.orthogonal_decomposition(to=to, verbosity=verbosity)[1]

    def gram_schmidt(self, factor: bool = True, verbosity: int = 1) -> "Matrix":
        """
        Performs Gram-Schmidt orthogonalization to convert a set of vectors (columns of the matrix) into 
        an orthogonal set (that includes 0 vectors if any).
        
        Args:
            factor (bool): If `True`, the resulting orthogonal vectors will be scaled to have integer factors.
            verbosity (int): Level of verbosity (default is `1`):
                - 0: No output.
                - 1: Display intermediate results for each step of the process.
            
        Returns:
            Matrix: A matrix whose columns are the orthogonalized vectors.
        
        Example:
            >>> mat = Matrix([[1, 2], [2, 4]])
            >>> ortho_mat = mat.gram_schmidt()
            >>> ortho_mat
            Matrix([[1, 0], [2, 0]])
        """

        if self.cols == 0:
            return []
        if is_IPython() and (verbosity >= 1):
            display(IPython.display.Math(f"v_{1} = {sym.latex(self[:, 0])}"))

        orthogonal_set = [self[:, 0]]
        for i in range(1, self.cols):
            u = self[:, i]
            latex_eq = f"v_{i + 1} = {sym.latex(u)}"
            for _, v in enumerate(orthogonal_set, start=1):
                if v.norm() != 0:
                    latex_eq += f"- ({sym.latex(v.dot(u))}/{sym.latex(v.dot(v))}) {sym.latex(v)}"
                    u -= (v.dot(u) / v.dot(v)) * v

            if is_IPython() and (verbosity >= 1):
                disp_u = u.copy()
                if factor:
                    scalar = sym.gcd(tuple(u))
                    disp_u = sym.MatMul(scalar, u / scalar, evaluate=False)
                latex_eq += f" = {sym.latex(disp_u)}"
                display(IPython.display.Math(latex_eq))

            if u.norm() == 0 and (verbosity >= 1):
                print(
                    "Warning: vectors are linearly dependent. Note that there is no QR factorisation"
                )
            orthogonal_set.append(u)

        return Matrix.from_list(orthogonal_set).normalized(factor=factor)

    def QRdecomposition(
        self, full: bool = False, *args, **kwargs
    ) -> tuple["Matrix", "Matrix"]:
        """
        Computes the QR decomposition of the matrix. Optionally computes the full QR decomposition.
        
        Args:
            full (bool): If `True`, computes the full QR decomposition (default is `False`).
            *args, **kwargs: Additional arguments passed to the SymPy `QRdecomposition`.
            
        Returns:
            tuple["Matrix", "Matrix"]: A tuple (Q, R) where:
                - Q is an orthogonal matrix.
                - R is an upper triangular matrix.
        
        Example:
            >>> mat = Matrix([[1, 2], [3, 4]])
            >>> Q, R = mat.QRdecomposition()
            >>> Q, R
        """

        # Modified SymPy's implementation to compute full QR decomposition if required. 
        Q, R = super().QRdecomposition(*args, **kwargs)
        if full and Q.rows != Q.cols:
            Q = Matrix(Q)
            Q_aug = Q.row_join(Q.elem()).QRdecomposition()[0]
            R_aug = Matrix(R.col_join(sym.zeros(Q_aug.cols - R.rows, R.cols)))
            assert Q_aug @ R_aug == self
            return Q_aug, R_aug
        return Q, R

    def solve_least_squares(
        self, rhs: "Matrix", verbosity: int = 1, matrices: int = 1, *args, **kwargs
    ) -> Union["Matrix", tuple["Matrix"]]:
        """
        Solves the least squares problem `Ax = b`, where `A` is the matrix and `b` is the right-hand side.
        Uses SymPy's built-in method for least squares when the rank condition is met, otherwise uses a custom
        solution approach.
        
        Args:
            rhs (Matrix): The right-hand side matrix/vector `b` in `Ax = b`.
            verbosity (int): Level of verbosity (default is `1`):
                - 0: No output.
                - 1: Display intermediate steps.
            matrices (int): If `1`, returns the solution matrix; if `2`, returns a tuple with part and general solutions.
            *args, **kwargs: Additional arguments passed to the `solve_least_squares` method of SymPy.
            
        Returns:
            "Matrix" | tuple["Matrix"]: The least squares solution matrix or a tuple with part and general solutions.
        
        Example:
            >>> A = Matrix([[1, 1], [2, 2]])
            >>> b = Matrix([1, 2])
            >>> A.solve_least_squares(b)
            Matrix([1, 1])
        """

        if verbosity == 0:
            try:
                A, b = sym.Matrix(self), sym.Matrix(rhs)
                return A.solve_least_squares(b, *args, **kwargs)
            except Exception as e:
                print(f"Exception Encountered: {str(e)}")
                print("Attempting custom solve...")

        ATA, ATb = self.T @ self, self.T @ rhs
        sol = Matrix.create_unk_matrix(ATb.rows, 1)
        sol = sol.subs(sym.solve(ATA @ sol - ATb))

        if verbosity >= 1:
            print("[A.T A | A.T b]")
            aug_matrix = ATA.aug_line().row_join(ATb)
            display(aug_matrix)
            print("\nAfter RREF")
            display(aug_matrix.rref()[0])

        if matrices == 1:
            return sol
        else:
            return sol.sep_part_gen()
        
    def create_vander(
        num_rows: int = 1, num_cols: int = 1, symbol: str = "x", is_real: bool = True
    ) -> "Matrix":
        """
        Creates a Vandermonde matrix with symbolic entries.

        This method generates a Vandermonde matrix of size `num_rows` x `num_cols`
        where the entries are symbolic expressions. Each row in the matrix is formed
        by raising a symbolic variable (indexed by row) to increasing powers (from 0
        to `num_cols-1`). The `is_real` flag determines whether the symbols are real-valued.

        Parameters:
        - num_rows (int, optional): The number of rows in the Vandermonde matrix. Default is 1.
        - num_cols (int, optional): The number of columns in the Vandermonde matrix. Default is 1.
        - symbol (str, optional): The base name for the symbols used in the matrix entries.
          Default is 'x'.
        - is_real (bool, optional): If True (default), the symbols are real-valued;
          otherwise, they are complex.

        Returns:
        - Matrix: A custom Matrix object (subclass of sympy.Matrix) that represents
          the Vandermonde matrix with symbolic entries.

        Example:
        >>> Matrix.create_vander(3, 3, symbol='a')
        Matrix([[a_1^0, a_1^1, a_1^2], [a_2^0, a_2^1, a_2^2], [a_3^0, a_3^1, a_3^2]])

        >>> Matrix.create_vander(2, 4, symbol='b', is_real=False)
        Matrix([[b_1^0, b_1^1, b_1^2, b_1^3], [b_2^0, b_2^1, b_2^2, b_2^3]])
        """

        entries = sym.symbols(f"{symbol}_(1:{num_rows+1})", is_real=is_real)
        res = []
        for entry in entries:
            sub_res = []
            for col_idx in range(num_cols):
                # Raise the symbol to the power of the column index
                sub_res.append(sym.Pow(entry, col_idx))
            res.append(sub_res)
        return Matrix(res)

    def apply_vander(self, x: "Matrix") -> "Matrix":
        """
        Applies a Vandermonde transformation to the current matrix using the given vector.

        This method applies a Vandermonde transformation to the current matrix by
        substituting the free symbols in the last column with corresponding values
        from the provided vector `x`. The number of rows in `self` must match the
        number of elements in `x`, and `x` must be a column vector.

        Parameters:
        - x (Matrix): A column vector (Matrix object with a single column) containing
          the values to substitute into the last column of the matrix.

        Returns:
        - Matrix: A new Matrix object where the free symbols in the last column of
          the original matrix are substituted by the corresponding values from `x`.

        Example:
        >>> mat = Matrix([[x_1, x_2], [x_3, x_4]])
        >>> x = Matrix([[1], [2]])
        >>> mat.apply_vander(x)
        Matrix([[1, 2], [1, 4]])

        Notes:
        - The matrix `self` is expected to be created via Matrix.create_vander().
        - The `x` vector provides the values to substitute in place of these symbols.
        """
        # Step 1: Validate the size of the vector x
        if x.cols != 1:
            raise sym.ShapeError(f"Input vector x must be a column vector. ({self.cols})")
        if self.rows != x.rows:
            raise sym.ShapeError(f"Number of rows in matrix ({self.rows}) must match the size of the input vector ({x.rows})")

        # Step 2: Get the free symbols from the last column of the matrix
        ordered_syms = [entry.free_symbols.pop() for entry in self.select_cols(-1)]

        # Step 3: Create a substitution dictionary mapping symbols to values from vector x
        substitution = {var: val for var, val in zip(ordered_syms, x)}
        return self.subs(substitution)

    #############################
    # CHAPTER 6: EIGEN-ANALYSIS #
    #############################

    def cpoly(self, force_factor: bool = True) -> Union[sym.Mul, tuple[sym.Mul, sym.Mul]]:
        """
        Computes the characteristic polynomial of the matrix and attempts to factor it into real and complex parts.
        
        Args:
            force_factor (bool): If `True`, the polynomial is fully factored, even if it doesn't have real factors.
                                 If `False`, the polynomial is returned in its factored form if possible. Default is `True`.
        
        Returns:
            sym.Mul | tuple[sym.Mul, sym.Mul]: 
                - If the polynomial factors only into real terms, returns a single factored polynomial.
                - If the polynomial has both real and complex factors, returns a tuple of two polynomials: 
                one with real factors and the other with complex factors.
        
        Example:
            >>> mat = Matrix([[1, 2], [3, 4]])
            >>> mat.cpoly()
            (x - 5, x + 1)
        """
        x = sym.symbols("x", real=True)
        poly = (x * self.elem() - self).det()
        if not force_factor:
            return poly.factor()
        # Attempt to factor poly into real factors
        try:
            roots = sym.roots(poly) # TODO: FIX sym.roots NotImplementedError for multi variable
            real_fact = []
            for root, mult in roots.items():
                term = x - root
                if mult != 1:
                    term = sym.Pow(term, mult, evaluate=False)
                if root.is_real:
                    real_fact.append(term)
                    poly /= term

            linear_fact = sym.Mul(*real_fact, evaluate=False)
            complex_fact = poly.expand().cancel().factor()

            if complex_fact == 1:
                return linear_fact
            else:
                return linear_fact, complex_fact
        except Exception as error:
            print(f"Encountered Error: {error}")
            return poly.factor()

    def is_diagonalizable(
        self, reals_only: bool = True, verbosity: int = 1, *args, **kwargs
    ) -> bool:
        """
        Checks if the matrix is diagonalizable, with the option to focus only on real eigenvalues.

        Args:
            reals_only (bool): If True, only real eigenvalues are considered for diagonalizability.
            verbosity (int): The level of verbosity for debug output. 1 prints the characteristic polynomial 
                             and eigenvectors.
            *args, **kwargs: Additional arguments passed to the base method.

        Returns:
            bool: True if the matrix is diagonalizable, False otherwise.
        """

        # Changed default for reals_only to True to align with MA1522 syllabus
        if verbosity >= 1:
            print("Characteristic Polynomial is: ")
            display(self.cpoly())
            print("Eigenvectors are: (eigenval, algebraic multiplicity, eigenspace) ")
            for val, mult, space in self.eigenvects():
                if val.is_real:
                    display([val, mult, space])

        return super().is_diagonalizable(reals_only, *args, **kwargs)
    
    def eigenvects_associated(self, eigenvalue: sym.Expr) -> Union[list["Matrix"], None]:
        # TODO
        return (eigenvalue @ self.elem() - self).nullspace()
    
    def diagonalize(
        self, reals_only: bool = True, verbosity: int = 0, *args, **kwargs
    ) -> tuple["Matrix", "Matrix"]:
        """
        Diagonalizes the matrix if possible, focusing on real eigenvalues if specified.

        Args:
            reals_only (bool): If True, diagonalization will focus on real eigenvalues.
            verbosity (int): Controls the level of output during the diagonalization process.
            *args, **kwargs: Additional arguments passed to the base method.

        Returns:
            tuple[Matrix, Matrix]: A tuple (P, D), where P is the matrix of eigenvectors and D is the diagonal matrix of eigenvalues.
        """

        # Changed default for reals_only to True to align with MA1522 syllabus
        if verbosity >= 1:
            print("Characteristic Polynomial")
            poly = self.cpoly()
            display(poly)
            for root, _ in sym.roots(poly).items():
                if root.is_real:
                    print(f"Original: ({root})I - A")
                    expr = root * self.elem() - self
                    display(expr)

                    print("\nAfter RREF:")
                    display(expr.rref())

                    print("\nEigenvectors:")
                    display(expr.nullspace())
                    print("\n")

        return super().diagonalize(reals_only, *args, **kwargs)

    @property
    def is_orthogonally_diagonalizable(self) -> bool:
        """
        Checks if the matrix is orthogonally diagonalizable (i.e., symmetric).

        Returns:
            bool: True if the matrix is symmetric and hence orthogonally diagonalizable, False otherwise.
        """
        return self.is_symmetric()

    def orthogonally_diagonalize(
        self, reals_only: bool = True, factor: bool = True, verbosity=1, *args, **kwargs
    ) -> tuple["Matrix", "Matrix"]:
        """
        Orthogonally diagonalizes the matrix, ensuring that eigenvectors corresponding to different eigenvalues are orthogonal.

        Args:
            reals_only (bool): If True, only real eigenvalues are considered.
            factor (bool): If True, the eigenvectors are orthogonalized using the Gram-Schmidt process.
            verbosity (int): Controls the verbosity of output during the process.
            *args, **kwargs: Additional arguments passed to the base method.

        Returns:
            tuple[Matrix, Matrix]: A tuple (P, D) where P is the orthogonal matrix of eigenvectors and D is the diagonal matrix of eigenvalues.
        """

        # Changed default for reals_only to True to align with MA1522 syllabus
        # Note that you can just apply GSP on P directly, since eigenspace associated to different eigenvalues are orthogonal
        # However, we follow the steps given in MA1522 syllabus here
        assert self.is_orthogonally_diagonalizable
        # P, D = super().diagonalize(reals_only, *args, **kwargs)
        P, D = self.diagonalize(
            reals_only=reals_only, verbosity=verbosity, *args, **kwargs
        )
        d = defaultdict(list)
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
                    ortho, norm = Matrix.from_list(vecs).gram_schmidt(factor, verbosity)
                    result.append(ortho @ norm)
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
                ortho_P = ortho_P.row_join(m)

        assert (ortho_P @ D @ ortho_P.T - self).norm() == 0
        return ortho_P, D

    @property
    def is_stochastic(self) -> bool:
        # TODO
        is_non_negative = all(entry >= 0 for entry in self.flat())
        is_prob_vectors = all(sum(self[:, i]) == 1 for i in range(self.cols))
        return is_non_negative and is_prob_vectors
    
    def equilibrium_vectors(self) -> "Matrix":
        """
        Computes the equilibrium vectors of the matrix, i.e., the nullspace of (I - A).
        Ensures that no negative entries are present in the matrix before computing.

        Returns:
            Matrix: A matrix containing equilibrium vectors normalized so that their
            column sums are zero.
        """

        assert all(entry >= 0 for entry in self.flat())
        P = Matrix.from_list((self.elem() - self).nullspace())
        for i in range(P.cols):
            if sum(P[:, i]) != 0:
                P[:, i] /= sum(P[:, i])
        return P

    def singular_value_decomposition(
        self, verbosity: int = 0, *args, **kwargs
    ) -> tuple["Matrix", "Matrix", "Matrix"]:
        """
        Performs Singular Value Decomposition (SVD) on the matrix, following the MA1522 syllabus.
        Returns U, S, and V matrices for the decomposition.

        Args:
            verbosity (int): Controls the verbosity of the output.
            *args, **kwargs: Additional arguments passed to the base method.

        Returns:
            Tuple[Matrix, Matrix, Matrix]: The U, S, and V matrices from the decomposition.
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
                        IPython.display.Math(
                            f"u_{idx} = (1/{sym.latex(val)})A{sym.latex(vec)} = {sym.latex(u_i)}"
                        )
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
                    orth, _ = complement.gram_schmidt(verbosity=verbosity)
                U = U.row_join(orth.normalized())

            # Add zero rows and columns to S so that matrix multiplication is defined
            m, n = self.shape
            r, c = S.shape
            S = S.row_join(sym.zeros(r, n - c)).col_join(sym.zeros(m - r, n))

            print("U, S, V")
            display([U, S, V])
            assert (U @ S @ V.T - self).norm() == 0
            return U, S, V

        m, n = self.shape
        U, S, V = super().singular_value_decomposition(*args, **kwargs)
        # Reverse index such that singular values are in decreasing order
        new_S = Matrix.diag(*S.diagonal()[::-1])

        S_index = [i for i in range(S.cols)][::-1]
        new_U = Matrix(U).select_cols(*S_index)
        new_V = Matrix(V).select_cols(*S_index)

        # new_U = Matrix(U).
        # Add orthonormal columns to U and V so that they are square matrices
        new_U = new_U.QRdecomposition(full=True)[0]
        new_V = new_V.QRdecomposition(full=True)[0]

        # Add zero rows and columns to S so that matrix multiplication is defined
        r, c = new_S.shape
        new_S = new_S.row_join(sym.zeros(r, n - c)).col_join(sym.zeros(m - r, n))
        assert self == new_U @ new_S @ new_V.T

        return new_U, new_S, new_V

    def fast_svd(
        self, option: Literal["np", "sym"] = "np", identify: bool = True, tolerance: float = None
    ) -> Union[tuple[np.matrix, np.matrix, np.matrix], tuple["Matrix", "Matrix", "Matrix"]]:
        """
        A faster version of SVD that returns a numerical decomposition of the matrix. 
        Optionally identifies surds/rationals from the float.

        Args:
            option (Literal["np", "sym"]): Whether to return numpy arrays or sympy matrices.
            identify (bool): Whether to attempt identification of rational numbers or surds.
            tolerance (float): Tolerance for numerical identification.

        Returns:
            Tuple: A tuple (U, S, V) where U, S, and V are the SVD components.
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
            return U, S, V
        elif option == "sym":
            U, S, V = Matrix(U), Matrix(S), Matrix(V)
            if identify:
                U = U.identify(tolerance=tolerance)
                S = S.identify(tolerance=tolerance)
                V = V.identify(tolerance=tolerance)
                residues = (self - U @ S @ V.T).norm()
                if residues != 0:
                    print(f"Non-zero Identification Error: {residues.evalf()}")
                return U, S, V
            else:
                return U, S, V
        else:
            print("Error: Unknown option type!")

    ####################################
    # CHAPTER 7: LINEAR TRANSFORMATION #
    ####################################

    def standard_matrix(self, out: "Matrix", matrices: int = 1) -> "Matrix":
        """
        Returns the standard matrix for the transformation from self to out.

        Args:
            out (Matrix): The target matrix for the transformation.
            matrices (int): The number of matrices to return (default is 1).

        Returns:
            Matrix: The matrix representing the transformation from self to out.
        """
        
        X = Matrix.create_unk_matrix(num_rows=out.rows, num_cols=self.rows)
        sol = sym.solve(X @ self - out, X.free_symbols)
        if isinstance(sol, list) and len(sol) > 0:
            # Multiple sets of solutions found, picks the first 1
            X = X.subs(sol[0])
        elif isinstance(sol, dict):
            X = X.subs(sol)

        if matrices == 1:
            return X
        else:
            return X.sep_part_gen()
