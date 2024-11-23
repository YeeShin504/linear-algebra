import sympy as sym
import numpy as np
from latex2sympy2 import latex2sympy
import re

import IPython.display
from itertools import chain, combinations
from typing import List

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
    >>> A.adjoint() # redefined behaviour as per MA1522 definition
    >>> A.adjugate() # classical adjoint (defined as adjoint in MA1522)
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
    >>> A.QRdecomposition(full: bool = False) # QR decomposition, returns (Q, R)
    >>> A.solve_least_squares(rhs: sym.Matrix, verbosity: int = 1) # solve the least square solution min(|Ax - rhs|)
    >>> A.is_diagonalizable(reals_only: bool = True, verbosity: int = 1) # modified diagonalizability criteria to align with MA1522 (reals_only)
    >>> A.singular_value_decomposition() # returns the full SVD required by MA1522 (A = U @ S @ V.T)

    # Custom Commands (verbosity >= 1 returns ERO (idx + 1), >= 2 returns matrix at each step)
    >>> is_zero(expr, symbolic: bool = True)
    >>> Matrix.from_latex(expr, row_join: bool = True, norm: bool = False) # creates a Matrix object from the LaTeX expression
    >>> Matrix.from_list(vectors: List, row_join: bool = True) # creates a Matrix object from a list of column vectors
    >>> Matrix.create_unk_matrix(num_rows: int, num_cols: int, symbol: str, is_real: bool) # creates a Matrix object with symbolic entries
    >>> Matrix.create_rand_matrix(num_rows: int, num_cols: int) # creates a Matrix object with random entries
    >>> A.simplify(rational: bool = True, tolerance: float, simplify: bool = True, expand: bool = True, collect_sym: sym.Symbol = None)
    >>> A.id() # returns the identity matrix with same number of rows as A
    >>> A.select_rows(*idx) # returns a new matrix with row vectors of *idx 
    >>> A.select_cols(*idx) # returns a new matrix with column vectors of *idx 
    >>> A.scale_row(idx: int, scalar: float, verbosity: int = 0) 
    >>> A.swap_row(idx_1: int, idx_2: int, verbosity: int = 0)
    >>> A.reduce_row(idx_1: int, scalar: float, idx_2: int, verbosity: int = 0)
    >>> A.get_pivot_row(col_idx: int, row_start_idx: int, follow_GE: bool = False)
    >>> A.ref(verbosity: int = 2, max_tries: int = 2, follow_GE: bool = False, matrices: int = 2) # matrices = 1 (U), 2 (LU), 3 (PLU)
    >>> A.evaluate_cases() # display a list of matrices for unknowns with critical values
    >>> A.column_constraints(use_id: bool = False, use_ref: bool = False) # returns the rref of [A | b], where b can be any vector. Use it to find constraints for b if A is not invertible.
    >>> A.extend_basis(span_subspace: sym.Matrix = A.id()) # returns an augmented matrix with additional column vectors to span the subspace of the argument.
    >>> A.transition_matrix(to: sym.Matrix) # returns a transition matrix from A to the other matrix
    >>> A.intersect_subspace(other: sym.Matrix, verbosity: int = 1) # returns a basis for the subspace that is in the intersection of A and other columnspace.
    >>> A.is_same_subspace(other: sym.Matrix, verbosity: int = 1) # returns True if columnspace of A and other are equal
    >>> A.inverse(option: str, verbosity: int = 0) # returns an inverse of A. If A is non-square return either a left inverse or right inverse.
    >>> A.orthogonal_complement(verbosity: int = 0) # returns a matrix whose column vectors are perpendicular (independent) to the column vectors of A (i.e. Null(A^T))
    >>> A.is_vec_orthogonal(verbosity: int = 1) # checks if the column vectors of A are orthogonal
    >>> A.normalized(factor: bool = False) # returns a matrix whose columnn vectors are normalized
    >>> A.scalar_factor(column: bool = True) # returns 2 matrices, one of which is a diagonal matrix, used to simplify matrix expressions
    >>> A.gram_schmidt(factor: bool = True, verbosity: int = 1)  # returns a matrix with orthonormal columns (and displays gram-schmidt process)
    >>> A.cpoly() # returns a factorised characterstic polynomial for a square matrix
    >>> A.is_diagonalizable
    >>> A.equilibrium_vectors() # returns a matrix whose column vectors are probability vectors such that Ax = x
    >>> A.fast_svd(option: str = 'np', identify: bool = True, tolerance: float = None) # returns a svd via numerical methods. Attempt to convert back to symbolic via identify not guaranteed.
    """
    print(commands)

def is_IPython():
    # Adapted from https://stackoverflow.com/questions/15411967/how-can-i-check-if-code-is-executed-in-the-ipython-notebook
    try:
        shell = get_ipython().__class__.__name__
        if shell == 'ZMQInteractiveShell' or shell == 'TerminalInteractiveShell':
            return True   # Jupyter notebook, qtconsole or terminal running IPython
        else:
            return False  # Other type
    except NameError:
        return False      # Probably standard Python interpreter


def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return list(chain.from_iterable(combinations(s, r) for r in range(len(s)+1)))


display = IPython.display.display if is_IPython() else sym.pprint
sym.init_printing(use_unicode=True)
np.set_printoptions(formatter={'float': lambda x: f"{x:10.7g}"})

####################
# CUSTOM FUNCTIONS #
####################

def is_zero(expr) -> bool:
    # returns True if expr can be 0 for real values of unknowns
    # use is_complex rather than is_real as symbol.is_complex and 
    # symbol.is_real returns None

    def is_real_dict(case) -> bool:
        for _, value in case.items():
            # zero is considered complex
            if value.is_complex and value != 0:
                return False
        return True

    if (not isinstance(expr, sym.Expr)) or isinstance(expr, sym.Number):
        return expr == 0
    
    sol = sym.solve(sym.Eq(expr, 0), expr.free_symbols)
    # display(f"{sol=}")
    for case in sol:
        if isinstance(case, dict):
            if is_real_dict(case):
                return True
            else:
                continue
        elif case.is_symbol or case.is_real:
            return True
        elif not case.is_complex:
            return True
    return False


class Matrix(sym.MutableDenseMatrix):
    def __init__(self, matrix) -> None:
        super().__init__()

    def from_latex(expr: str,
                   row_join = True,
                   norm = False) -> sym.Matrix:
        # returns a Matrix object by parsing the LaTeX expression
        # row_join = True will treat the list of vectors as column vectors
        # norm = True will normalise the column vectors if they exist

        # Replace \begin{array}{ccc*} with \begin{pmatrix}
        modified_latex = re.sub(r'\\begin\{array\}\{[^}]*\}', r'\\begin{pmatrix}', expr)
        # Replace \end{array} with \end{pmatrix}
        modified_latex = re.sub(r'\\end\{array\}', r'\\end{pmatrix}', modified_latex)

        modified_latex = re.sub(r'\\;', "", modified_latex)

        res = latex2sympy(modified_latex)
        display(res)
        if isinstance(res, list):
            vector_list = []
            for vector in res:
                vector = vector.expand()
                if norm:
                    vector /= vector.norm()
                vector_list.append(vector)
            return Matrix.from_list(vector_list, row_join)

        elif isinstance(res, sym.matrices.expressions.matmul.MatMul):
            # If the matrix is a product of matrices, evaluate directly
            return Matrix(res.doit())
        elif isinstance(res, sym.Matrix):
            return Matrix(res)
        else:
            return res
        
    def from_list(vectors: List[sym.Matrix],
                  row_join: bool = True) -> sym.Matrix:
        res = Matrix([list(vector) for vector in vectors])
        if row_join:
            return res.transpose()
        else:
            return res
    
    def id(self) -> sym.Matrix:
        return Matrix(sym.eye(self.rows))
    
    def create_unk_matrix(num_rows: int = 1,
                          num_cols: int = 1,
                          symbol: str = 'x',
                          is_real: bool = True) -> sym.Matrix:
        # Creates a matrix of size rows * cols with entries symbol_ij
        entries = sym.symbols(f'{symbol}(1:{num_rows+1})(1:{num_cols+1})', is_real=is_real)
        return Matrix(entries).reshape(num_rows, num_cols)
    
    def create_rand_matrix(num_rows: int = 1,
                           num_cols: int = 1) -> sym.Matrix:
        return Matrix(sym.randMatrix(num_rows, num_cols))
    
    def simplify(self,
                 rational: bool = True,
                 tolerance: float = 1E-4,
                 simplify: bool = True,
                 expand: bool = True,
                 collect_sym: sym.Symbol = None,
                 *args,
                 **kwargs) -> sym.Matrix:
        temp = self.copy()
        if rational:
            temp = sym.nsimplify(temp, tolerance = tolerance, rational = True)
            assert (temp - self).norm() == 0, 'Rational approximation failed. Try lower tolerance.'
        if simplify:
            temp = sym.simplify(temp, *args, **kwargs)
        if expand:
            temp = sym.expand(temp)
        else:
            temp = sym.factor(temp)
        temp = temp.tolist()
        if collect_sym is not None:
            temp = [[sym.collect(item, collect_sym) for item in row] for row in temp]
        temp = Matrix(temp)
        self.__dict__.update(temp.__dict__)
        return self

    def select_cols(self, *args: int) -> sym.Matrix:
        res = []
        for idx in args:
            res.append(list(self.col(idx)))
        return Matrix(res).transpose()
    
    def select_rows(self, *args: int) -> sym.Matrix:
        res = []
        for idx in args:
            res.append(list(self.row(idx)))
        return Matrix(res)

    def scale_row(self,
                  idx: int,
                  scalar: float,
                  verbosity: int = 0):
        self[idx, :] = scalar * self[idx, :]
        self.simplify()

        if verbosity >= 1:
            print(f"R_{idx+1} <- ({scalar})R_{idx+1}")
        if verbosity >= 2:
            display(self)
            print('\n')
        
        return self

    def swap_row(self,
                 idx_1: int,
                 idx_2: int,
                 verbosity: int = 0):
        self[idx_1, :], self[idx_2, :] = self[idx_2, :], self[idx_1, :]

        if verbosity >= 1:
            print(f"R_{idx_1+1} <-> R_{idx_2+1}")
        if verbosity >= 2:
            display(self)
            print('\n')
        
        return self

    def reduce_row(self,
                   idx_1: int,
                   scalar: float,
                   idx_2: int,
                   verbosity: int = 0):
        self[idx_1, :] = self[idx_1, :] - scalar * self[idx_2, :]
        self.simplify()

        if verbosity >= 1:
            print(f"R_{idx_1+1} <- R_{idx_1+1} - ({scalar})R_{idx_2+1}")
        if verbosity >= 2:
            display(self)
            print('\n')

        return self

    def get_pivot_row(self,
                      col_idx: int,
                      row_start_idx: int,
                      follow_GE: bool = False) -> int:

        # Attempt to pick a pivot column that is a non-zero constant that do
        # not contain any symbols so that it is easier to reduce other rows
        if not follow_GE:
            for row_idx in range(row_start_idx, self.rows):
                term = self[row_idx, col_idx]
                if term != 0:
                    if not isinstance(term, sym.Expr):
                        return row_idx
                    elif len(term.free_symbols) == 0:
                        return row_idx

        # Attempt to pick the first non-zero row if all rows contain symbols
        for row_idx in range(row_start_idx, self.rows):
            term = self[row_idx, col_idx]
            if term != 0:
                return row_idx

        # if entire col is 0, ie no pivot_rows found, return -1
        return -1

    def ref(self, 
            verbosity: int = 2,
            max_tries: int = 2,
            follow_GE: bool = False,
            matrices: int = 2) -> sym.Matrix:
        # follow_GE can be set to True to follow Gaussian Elimination strictly
        # LU decomposition
        # matrices is the number of matrices return.
        # 1. Upper
        # 2. Perm @ Lower Upper
        # 3. Perm Lower Upper
        U = self.copy()

        I = self.id()
        L = self.id()
        P = self.id()

        # Loop over each column
        cur_row_pos = 0
        for cur_col_pos in range(min(self.rows, self.cols)):
            # Find the first non-zero row in the current column
            pivot_row = U.get_pivot_row(cur_col_pos, cur_row_pos, follow_GE)
            # print(f"{pivot_row=}")

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
            for row_idx in range(cur_row_pos+1, self.rows):
                # reduce the row_idx iteratively via partial fractions to
                # prevent division by a possible 0 term
                tries = 0
                while U[row_idx, cur_col_pos] != 0:
                    tries += 1
                    if tries > max_tries:
                        print(f"ERROR: Max tries exceeded to reduce row {row_idx+1} with row {cur_row_pos+1}")
                        break
                    try:
                        scalar = U[row_idx, cur_col_pos] / U[cur_row_pos, cur_col_pos] 
                        scalar = scalar.expand().simplify()
                        # print(f'{scalar=}')
                        try:
                            decomp = sym.apart(scalar) # partial fractions
                        except:
                            decomp = scalar
                        # if isinstance(decomp, sym.core.add.Add) and (scalar != decomp):
                        if isinstance(decomp, sym.core.add.Add):
                            terms = decomp.args
                        else:
                            # there is only 1 term (could be integer or proper fraction)
                            terms = [decomp]
                        # print(f'{terms=}')
                        for term in terms:
                            _, d = sym.fraction(term)
                            # print(f'{d=}')
                            # ensure denominator is non-zero so that reduction is valid
                            if not is_zero(d):
                                U.reduce_row(row_idx, term, cur_row_pos, verbosity=verbosity)
                                elem = I.copy().reduce_row(row_idx, -term, cur_row_pos)
                                L = L @ elem

                        # Cases where pivot row contains symbols such that scalar is a
                        # fraction with symbolic denominator.
                        # To reduce further, can only scale row_idx accordingly
                        if U[row_idx, cur_col_pos] != 0:
                            scalar = U[cur_row_pos, cur_col_pos] / U[row_idx, cur_col_pos]
                            scalar.simplify()
                            n, d = sym.fraction(scalar)
                            # to scale by n, n cannot be 0 both numerically or symbolically
                            # to scale by 1/d, d cannot be 0, same argument as n
                            # if (n != 0) and (not is_zero(d)):
                            if (not is_zero(n)) and (not is_zero(d)):
                                U.scale_row(row_idx, scalar, verbosity=verbosity)
                                elem = I.copy().scale_row(row_idx, 1/scalar)
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

        if matrices == 1:
            return U
        elif matrices == 2:
            return P @ L, U
        elif matrices == 3:
            return P, L, U
        else:
            return U

    def get_pivot_pos(self) -> List[List]:
        assert self.is_echelon # check for REF
        
        pivot_pos = []
        cur_row_pos = 0
        for cur_col_pos in range(self.cols):
            pivot_row = self.get_pivot_row(cur_col_pos, cur_row_pos, follow_GE = False)
            
            if pivot_row != -1:
                pivot_pos.append((pivot_row, cur_col_pos))
                cur_row_pos += 1
            
        return pivot_pos

    def get_pivot_elements(self) -> List:
        pivot_elements = []
        
        for i, j in self.get_pivot_pos():
            pivot_elements.append(self[i, j])
            
        return pivot_elements
        
    def find_all_cases(self) -> List:
        cases = []
        for pivot in self.get_pivot_elements():
            if len(pivot.free_symbols) > 0: # case found
                cases.extend(sym.solve(pivot, dict = True))
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

    def evaluate_cases(self) -> None:
        cases = self.find_all_cases()
        all_possible_values = set(possible_val for case in cases for possible_val in case.items())
            
        for i, case in enumerate(cases, 1):
            print(f"Case {i}: {case}, not including {[dict([val]) for val in all_possible_values.symmetric_difference(set(case.items()))]}")
            display(self.copy().xreplace(case).ref(verbosity = 0, matrices = 1, follow_GE = False))    
    
    def rref(self, *args, **kwargs) -> tuple[sym.Matrix, tuple[int]]:
        my_rref, pivots = super().rref(*args, **kwargs)
        return Matrix(my_rref), pivots
    
    def adjoint(self, *args, **kwargs) -> sym.Matrix:
        # point adjoint to classical adjoint (which is defined as adjugate in SymPy) to align with MA1522 Syllabus
        return Matrix(super().adjugate(*args, **kwargs))
    
    def column_constraints(self, 
                           use_id: bool = False,
                           use_ref: bool = False,
                           **kwargs) -> sym.Matrix:
        # write a random vector as x_1, ..., x_m, given m rows
        vector = Matrix.create_unk_matrix(self.rows, 1, 'x')

        # insert hidden column vectors so that the unknown vector is not reduced in rref
        if use_id:
            hidden = self.id()
        else:
            orth_vecs = self.orthogonal_complement()
            if len(orth_vecs) == 0:
                # No orth_vecs found (should be full rank) and default to id to be safe
                hidden = self.id()
            else:
                hidden = orth_vecs

        # display(hidden)
        M = self.copy().row_join(hidden).row_join(vector)
        if use_ref:
            res = M.ref(**kwargs, matrices=1)
        else:
            res, _ = M.rref()
        res_matrix = res[:, :self.cols].row_join(res[:, -1])
        return Matrix(res_matrix)
    
    def extend_basis(self,
                     span_subspace: sym.Matrix = None,
                     verbosity: int = 1) -> sym.Matrix:
        if span_subspace is None:
            span_subspace = self.id()
        aug = self.row_join(span_subspace)
        my_rref, pivots = aug.rref()

        if verbosity >= 1:
            display(my_rref)

        return aug.select_cols(*pivots)

    def transition_matrix(self, 
                          to: sym.Matrix,
                          verbosity: int = 1) -> sym.Matrix:
        assert self.rows == to.rows

        M = to.row_join(self)
        res = M.rref()[0]
        if verbosity >= 1:
            display(res)
        P = res[:self.cols, self.cols:]
        return P
    
    def intersect_subspace(self,
                           other: sym.Matrix,
                           verbosity: int = 1) -> sym.Matrix:
        # This method uses the orthogonal complement to find a matrix whose 
        # columns form a basis for the subspace (self n other)

        # Construct 2 matrices A and B, whose solution space (ie nullspace) is
        # the subspace self and other respectively. Observe that the solution 
        # space is orthogonal to the row space, so it is the orthogonal complement.

        A = self.orthogonal_complement().T
        B = other.orthogonal_complement().T
        # Now we obtain A and B which represent the linear system of 2 different 
        # subspaces. When we solve these simultaneously, we will find the solution 
        # space which contains vectors which are solutions to both linear systems.
        if verbosity >= 1:
            print('A linear system whose solution space is the subspace of self. Null(self^T)^T')
            display(A)
            print('\nA linear system whose solution space is the subspace of other. Null(other^T)^T')
            display(B)
            display(A.col_join(B).rref())

        return Matrix.from_list(A.col_join(B).nullspace())

    def is_same_subspace(self,
                         other: sym.Matrix,
                         verbosity: int = 1) -> bool:
        assert self.rows == other.rows

        sub_rref, sub_pivots = other.row_join(self).rref()
        super_rref, super_pivots = self.row_join(other).rref()
        
        if verbosity >= 1:
            print('Check if span(self) is subspace of span(other)')
            display(sub_rref)
            print('\nCheck if span(other) is subspace of span(self)')
            display(super_rref)
        
        return (max(sub_pivots) < self.cols) and (max(super_pivots) < other.cols)
        
    def inverse(self, 
                option: str = None, 
                verbosity: int = 0) -> tuple[sym.Matrix, sym.Matrix]:
        if option is None:
            if self.rank() == self.cols:
                if verbosity:
                    print("Left inverse found!")
                option = 'left'
            if self.rank() == self.rows:
                if verbosity:
                    print("Right inverse found!")
                option = 'right'
        
        if option is not None:
            X = Matrix.create_unk_matrix(self.cols, self.rows, 'x')
            if option == 'left':
                eqn = X @ self - sym.eye(self.cols)
            else:
                eqn = self @ X - sym.eye(self.rows)
            
            sol = sym.solve(eqn)
            if isinstance(sol, list) and len(sol) > 0:
                # Multiple sets of solutions found, picks the first 1
                X = X.subs(sol[0])
            elif isinstance(sol, dict):
                X = X.subs(sol)
            else:
                print(f'No {option} inverse found!')
                return None
            set_0 = dict([(symbol, 0) for symbol in X.free_symbols])
            part_sol = X.subs(set_0)
            gen_sol = X - part_sol
            return part_sol, gen_sol

        if verbosity:
            print("No inverse found! Try pseudo-inverse: .pinv()")
    
    def nullspace(self, 
                  verbosity: int = 1,
                  *args,
                  **kwargs) -> List[sym.Matrix]:
        # Issue with SymPy implementation of nullspace when there is None
        # Using rank nullity theorem to verify there are vectors spanning nullspace
        if self.rank() == self.cols:
            if verbosity >= 1:
                print("No nullspace detected!")
            return []
        else:
            return super().nullspace(*args, **kwargs)

    def orthogonal_complement(self, 
                              verbosity: int = 0) -> sym.Matrix:
        return Matrix.from_list(self.transpose().nullspace(verbosity))
    
    def is_vec_orthogonal(self,
                          verbosity: int = 1) -> bool:
        # Returns true if COLUMN vectors are orthoGONAL (ie u_i \dot u_j = 0)   
        # Note that this is different from the definition of orthoGONAL MATRIX 
        # where its column vectors are orthoNORMAL, such that M^T M = I
        res = self.T @ self
        if verbosity:
            display(res)
        return res.is_diagonal()
    
    def normalized(self, factor: bool = False) -> sym.Matrix:
        # normalize the column vectors of the matrix
        for i in range(self.cols):
            scalar = self.col(i).norm()
            if scalar != 0:
                self[:, i] /= scalar
            
        if factor:
            return self.scalar_factor(column=True)
        else:
            return self
    
    def scalar_factor(self,
                      column = True) -> tuple[sym.Matrix, sym.Matrix]:
        # Factorises a matrix into A = BD, where D is a diagonal matrix (for column vectors)
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

    def gram_schmidt(self,
                     factor: bool = True,
                     verbosity: int = 1) -> sym.Matrix:
        if self.cols == 0:
            return []
        if is_IPython():
            display(IPython.display.Math(f"v_{1} = {sym.latex(self[:, 0])}"))
        
        orthogonal_set = [self[:, 0]]
        for i in range(1, self.cols):
            u = self[:, i] # slicing gives a copy of what is sliced, so modifications of one object do not affect the other
            latex_eq = f"v_{i + 1} = {sym.latex(u)}"
            for _, v in enumerate(orthogonal_set, start = 1):
                if v.norm() != 0:
                    latex_eq += f"- ({sym.latex(v.dot(u))}/{sym.latex(v.dot(v))}) {sym.latex(v)}"
                    u -= (v.dot(u) / v.dot(v)) * v

            if is_IPython() and (verbosity >= 1):
                disp_u = u.copy()
                if factor:
                    scalar = sym.gcd(tuple(u))
                    disp_u = sym.MatMul(scalar, u/scalar, evaluate=False)
                latex_eq += f" = {sym.latex(disp_u)}"
                display(IPython.display.Math(latex_eq))

            if u.norm() == 0 and (verbosity >= 1):
                print("Warning: vectors are linearly dependent. Note that there is no QR factorisation")
            orthogonal_set.append(u)
        
        return Matrix.from_list(orthogonal_set).normalized(factor=factor)
    
    def QRdecomposition(self, 
                        *args, 
                        full: bool = False,
                        **kwargs) -> tuple[sym.Matrix, sym.Matrix]:
        # full is set to False as prefer not to override SymPy QRdecomposition as 
        # their default behaviour might be used for other algorithms like SVD
        Q, R = super().QRdecomposition(*args, **kwargs)
        if full and Q.rows != Q.cols:
            Q = Matrix(Q)
            Q_aug = Q.row_join(Q.id()).QRdecomposition()[0]
            R_aug = Matrix(R.col_join(sym.zeros(Q_aug.cols - R.rows, R.cols)))
            assert Q_aug @ R_aug == self
            return Q_aug, R_aug
        return Q, R
    
    def solve_least_squares(self, 
                            rhs: sym.Matrix, 
                            *args, 
                            verbosity: int = 1, **kwargs) -> sym.Matrix:
        # SymPy solve_least_squares requires unique solution, ie rank(self) > self.cols
        if verbosity == 0:
            try:
                A, b = sym.Matrix(self), sym.Matrix(rhs)
                return A.solve_least_squares(b, *args, **kwargs)
            # return super().solve_least_squares(rhs=rhs, method=method)
            except Exception as e:
                print(f'Exception Encountered: {str(e)}')
                print('Attempting custom solve...')

        ATA, ATb = self.T @ self, self.T @ rhs
        sol = Matrix.create_unk_matrix(ATb.rows)
        sol = sol.subs(sym.solve(ATA @ sol - ATb))

        if verbosity >= 1:
            print('[A.T A | A.T b]')
            aug_matrix = ATA.row_join(ATb)
            display(aug_matrix)
            print('After RREF')
            display(aug_matrix.rref()[0])

        set_0 = dict([(symbol, 0) for symbol in sol.free_symbols])
        part_sol = sol.subs(set_0)
        gen_sol = sol - part_sol
        return part_sol, gen_sol
    
    def cpoly(self,
              force_factor: bool = True):
        x = sym.symbols('x')
        poly = (x*self.id() - self).det()
        if not force_factor:
            return poly.factor()
        # Attempt to factor poly into real factors
        roots = sym.roots(poly)
        linear_fact = []
        complex_fact = []
        for root, mult in roots.items():
            term = (x - root)
            if mult != 1:
                term = sym.Pow(term, mult, evaluate=False)
            if root.is_real:
                linear_fact.append(term)
            else:
                complex_fact.append(term)
        linear_fact = sym.Mul(*linear_fact, evaluate=False)

        if len(complex_fact) == 0:
            return linear_fact
        
        complex_fact = sym.Mul(*complex_fact).expand().simplify().factor()
        return linear_fact, complex_fact
    
    def is_diagonalizable(self, 
                          reals_only: bool = True, 
                          verbosity: int = 1,
                          *args, **kwargs) -> bool:
        # Changed default for reals_only to True to align with MA1522 syllabus 
        if verbosity >= 1:
            print('Characteristic Polynomial is: ')
            display(self.cpoly())
            print('Eigenvectors are: (eigenval, multiplicity, eigenspace) ')
            for val, mult, space in self.eigenvects():
                if val.is_real:
                    display([val, mult, space])
        
        return super().is_diagonalizable(reals_only, *args, **kwargs)
    
    def equilibrium_vectors(self) -> sym.Matrix:
        neg_entries = [entry for entry in self.flat() if entry < 0]
        assert len(neg_entries) == 0
        P = Matrix.from_list((self.id() - self).nullspace())
        for i in range(P.cols):
            if sum(P[:, i]) != 0:
                P[:, i] /= sum(P[:, i])
        return P

    def singular_value_decomposition(self) -> tuple[sym.Matrix, sym.Matrix, sym.Matrix]:
        # return the full SVD that is required by MA1522 syllabus
        m, n = self.shape
        U, S, V = super().singular_value_decomposition()
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
        new_S = new_S.row_join(sym.zeros(r,n-c)).col_join(sym.zeros(m-r, n))
        assert self == new_U @ new_S @ new_V.T

        return new_U, new_S, new_V
    
    def fast_svd(self, 
                 option: str = 'np', 
                 identify: bool = True, 
                 tolerance: float = None):
        # return a faster version of SVD, but not symbolically
        # if identify, attempt to identify surds and rationals from the float
        m, n = self.shape 
        U, S, Vh = np.linalg.svd(np.array(self, dtype=float))

        # Create sigma matrix from singular values
        S = np.diag(S)
        r, c = S.shape
        S = np.concat((S, np.zeros((r, n - c))), axis=1)
        S = np.concat((S, np.zeros((m - r, n))), axis=0)
        if option == 'np':
            return U, S, Vh
        elif option == 'sym':
            U, S, V = Matrix(U), Matrix(S), Matrix(Vh)
            if identify:
                U = Matrix(sym.nsimplify(U, tolerance=tolerance))
                S = Matrix(sym.nsimplify(S, tolerance=tolerance))
                V = Matrix(sym.nsimplify(V, tolerance=tolerance))
                assert self == U @ S @ V.T, 'Identification failed!'
                return U, S, V
            else:
                return U, S, V
        else:
            print('Error: Unknown option type!')