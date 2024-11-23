import sympy as sym
from latex2sympy2 import latex2sympy

import IPython.display
from itertools import chain, combinations
from typing import List, Dict

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
    >>> A.adjoint() # adjoint
    >>> A.adjugate() # adjugate
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
    >>> A.simplify() # simplify
    >>> A.LUdecomposition(iszerofunc=<function _iszero>, simpfunc=None, rankcheck=False) # LU decomposition, returns (L, U, perm)
    >>> A.lower_triangular_solve(rhs) # solve the matrix equation Ax = rhs, where A is an lower triangular matrix
    >>> A.upper_triangular_solve(rhs) # solve the matrix equation Ax = rhs, where A is an upper triangular matrix

    # Custom Commands (verbosity >= 1 returns ERO (idx + 1), >= 2 returns matrix at each step)
    >>> Matrix.from_latex(expr, row_join: bool = True, norm: bool = False) # creates a Matrix object from the LaTeX expression
    >>> Matrix.create_unk_matrix(num_rows: int, num_cols: int, symbol: str, is_real: bool) # creates a Matrix object with unknown entries
    >>> A.id() # returns the identity matrix with same number of rows as A
    >>> A.select_rows(*idx) # returns a new matrix with row vectors of *idx 
    >>> A.select_cols(*idx) # returns a new matrix with column vectors of *idx 
    >>> A.scale_row(idx: int, scalar: float, verbosity: int = 0) 
    >>> A.swap_row(idx_1: int, idx_2: int, verbosity: int = 0)
    >>> A.reduce_row(idx_1: int, scalar: float, idx_2: int, verbosity: int = 0)
    >>> is_zero(expr, symbolic: bool = True)
    >>> A.get_pivot_row(col_idx: int, row_start_idx: int, follow_GE: bool = False)
    >>> A.ref(verbosity: int = 2, max_tries: int = 2, follow_GE: bool = False, matrices: int = 2) # matrices = 1 (U), 2 (LU), 3 (PLU)
    >>> A.evaluate_cases() # display a list of matrices for unknowns with critical values
    >>> A.column_constraints(use_id: bool = False) # returns the rref of [A | b], where b can be any vector. Use it to find constraints for b if A is not invertible.
    >>> A.extend_basis(span_subspace: sym.Matrix = A.id()) # returns an augmented matrix with additional column vectors to span the subspace of the argument.
    >>> A.transition_matrix(to: sym.Matrix) # returns a transition matrix from A to the other matrix
    >>> A.inverse(option: str) # returns an inverse of A. If A is non-square return either a left inverse or right inverse.
    >>> A.orthogonal_complement() # returns a list of column vectors that are perpendicular (independent) to the column vectors of A
    >>> A.is_orthogonal() # checks if the column vectors of A are orthogonal
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

####################
# CUSTOM FUNCTIONS #
####################

def is_zero(expr) -> bool:
    if not isinstance(expr, sym.Expr):
        return expr == 0
    sol = sym.solve(sym.Eq(expr, 0), expr.free_symbols)
    display(sol)
    real_sol = [v for v in sol if v.is_real]
    return len(real_sol) != 0

class Matrix(sym.MutableDenseMatrix):
    def __init__(self, matrix) -> None:
        super().__init__()

    def from_latex(expr: str,
                   row_join = True,
                   norm = False) -> sym.Matrix:
        # returns a Matrix object by parsing the LaTeX expression
        # row_join = True will treat the list of vectors as column vectors
        # norm = True will normalise the column vectors if they exist 
        res = latex2sympy(expr)
        display(res)
        if isinstance(res, list):
            vector_list = []
            for vector in res:
                vector = vector.expand()
                if norm:
                    vector = vector.normalized()
                vector_list.append(list(vector))
            
            if row_join:
                return Matrix(vector_list).transpose()
            else:
                return Matrix(vector_list)
            
        elif isinstance(res, sym.Matrix):
            return Matrix(res)
        else:
            return res
    
    def id(self) -> sym.Matrix:
        return Matrix(sym.eye(self.rows))
    
    def create_unk_matrix(num_rows: int = 1,
                          num_cols: int = 1,
                          symbol: str = 'x',
                          is_real: bool = True) -> sym.Matrix:
        # Creates a matrix of size rows * cols with entries symbol_ij
        entries = sym.symbols(f'{symbol}(1:{num_rows+1})(1:{num_cols+1})', is_real=True)
        return Matrix(entries).reshape(num_rows, num_cols)
    
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
                        display(U[row_idx, cur_col_pos])
                        print(f"ERROR: Max tries exceeded to reduce row {row_idx+1} with row {cur_row_pos+1}")
                        break
                    try:
                        scalar = U[row_idx, cur_col_pos] / U[cur_row_pos, cur_col_pos] 
                        # simplify scalar is needed to check if apart has effect
                        scalar.simplify()
                        try:
                            decomp = sym.apart(scalar) # partial fractions
                        except:
                            decomp = scalar
                        if isinstance(decomp, sym.core.add.Add) and (scalar != decomp):
                            terms = decomp.args
                        else:
                            # there is only 1 term (could be integer or proper fraction)
                            terms = [decomp]
                        for term in terms:
                            display(term)
                            _, d = sym.fraction(term)
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
                            # elif (not is_zero(n)):
                            #     U.scale_row(row_idx, n, verbosity=verbosity)
                            #     elem = I.copy().scale_row(row_idx, 1/n)
                            #     L = L @ elem
                    except Exception as e:
                        print(str(e))
                        return P, L, U

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
    
    def column_constraints(self, 
                           use_id = False,
                           use_ref = False,
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
            elif len(orth_vecs) == 1:
                hidden = orth_vecs[0]
            else:
                # Matrix.row_join requires at least 2 vectors
                hidden = Matrix.row_join(*self.orthogonal_complement())

        display(hidden)
        M = self.copy().row_join(hidden).row_join(vector)
        if use_ref:
            res = M.ref(*kwargs, matrices=1)
        else:
            res, _ = M.rref()
        res_matrix = res[:, :self.cols].row_join(res[:, -1])
        return Matrix(res_matrix)
    
    def extend_basis(self,
                     span_subspace: sym.Matrix = None,
                     verbosity: int = 1) -> List[sym.Matrix]:
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

        M = to.copy().row_join(self)
        res = M.rref()[0]
        if verbosity >= 1:
            display(res)
        P = res[:self.cols, self.cols:]
        # non_zero_rows = [row for row in P.tolist() if any(row)]
        return P
    
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
            X = X.subs(sym.solve(eqn))
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
                              verbosity: int = 0) -> List[sym.Matrix]:
        return self.transpose().nullspace(verbosity)
    
    def is_vec_orthogonal(self,
                          verbosity: int = 0) -> bool:
        # Returns true if COLUMN vectors are orthoGONAL (ie u_i \dot u_j = 0)   
        # Note that this is different from the definition of orthoGONAL MATRIX 
        # where its column vectors are orthoNORMAL, such that M M^T = I
        for i in range(self.cols):
            u = self.col(i)
            for j in range(i+1, self.cols):
                v = self.col(j)
                if u.dot(v) != 0:
                    if verbosity >= 1:
                        display([u, v])
                    return False
        return True
    