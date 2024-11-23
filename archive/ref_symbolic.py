import sympy as sym
import IPython.display
from itertools import chain, combinations

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
    >>> A.rref() # RREF
    >>> A.simplify() # simplify
    >>> A.LUdecomposition(iszerofunc=<function _iszero>, simpfunc=None, rankcheck=False) # LU decomposition, returns (L, U, perm)

    # Custom Commands (verbosity >= 1 returns ERO (idx + 1), >= 2 returns matrix at each step)
    scale_row(A: sym.Matrix, idx: int, scalar: float, verbosity: int = 0)
    swap_row(A: sym.Matrix, idx_1: int, idx_2: int, verbosity: int = 0)
    reduce_row(A: sym.Matrix, idx_1: int, scalar: float, idx_2: int, verbosity: int = 0)
    is_zero(expr, symbolic: bool = True)
    get_pivot_row(A: sym.Matrix, col_idx: int, row_start_idx: int, follow_GE: bool = False)
    def ref(A: sym.Matrix, verbosity: int = 2, max_tries: int = 2, follow_GE: bool = False, matrices: int = 2) # matrices = 1 (U), 2 (LU), 3 (PLU)
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
def scale_row(A: sym.Matrix,
              idx: int,
              scalar: float,
              verbosity: int = 0) -> sym.Matrix:
    A = sym.Matrix(A)
    A[idx, :] = scalar * A[idx, :]
    A = sym.simplify(A)
    if verbosity >= 1:
        print(f"R_{idx+1} <- ({scalar})R_{idx+1}")
    if verbosity >= 2:
        display(A)
        print('\n')
    return A

def swap_row(A: sym.Matrix,
             idx_1: int,
             idx_2: int,
             verbosity: int = 0) -> sym.Matrix:
    A = sym.Matrix(A)
    A[idx_1, :], A[idx_2, :] = A[idx_2, :], A[idx_1, :]

    if verbosity >= 1:
        print(f"R_{idx_1+1} <-> R_{idx_2+1}")
    if verbosity >= 2:
        display(A)
        print('\n')
    return A

def reduce_row(A: sym.Matrix,
               idx_1: int,
               scalar: float,
               idx_2: int,
               verbosity: int = 0) -> sym.Matrix:
    A = sym.Matrix(A)
    A[idx_1, :] = A[idx_1, :] - scalar * A[idx_2, :]
    A = sym.simplify(A)

    if verbosity >= 1:
        print(f"R_{idx_1+1} <- R_{idx_1+1} - ({scalar})R_{idx_2+1}")
    if verbosity >= 2:
        display(A)
        print('\n')
    return A

def is_zero(expr) -> bool:
    if not isinstance(expr, sym.Expr):
        return expr == 0
    sol = sym.solve(sym.Eq(expr, 0), expr.free_symbols)
    real_sol = [v for v in sol if v.is_real]
    return len(real_sol) != 0

def get_pivot_row(A: sym.Matrix,
                  col_idx: int,
                  row_start_idx: int,
                  follow_GE: bool = False) -> int:
    m = A.rows

    # Attempt to pick a pivot column that is a non-zero constant that do
    # not contain any symbols so that it is easier to reduce other rows
    if not follow_GE:
        for row_idx in range(row_start_idx, m):
            term = A[row_idx, col_idx]
            if term != 0:
                if not isinstance(term, sym.Expr):
                    return row_idx
                elif len(term.free_symbols) == 0:
                    return row_idx

    # Attempt to pick the first non-zero row if all rows contain symbols
    for row_idx in range(row_start_idx, m):
        term = A[row_idx, col_idx]
        if term != 0:
            return row_idx

    # if entire col is 0, ie no pivot_rows found, return -1
    return -1

def ref(A: sym.Matrix,
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
    U = A.copy()
    m, n = sym.shape(U)

    I = sym.eye(m)
    L = sym.eye(m)
    P = sym.eye(m)

    # Loop over each column
    cur_row_pos = 0
    for cur_col_pos in range(min(m, n)):
        # Find the first non-zero row in the current column
        pivot_row = get_pivot_row(U, cur_col_pos, cur_row_pos, follow_GE)
        # print(f"{pivot_row=}")

        if pivot_row == -1:
            # If no non-zero pivot is found, continue to the next column
            continue

        # Swap the current row with the pivot row if necessary
        if pivot_row != cur_row_pos:
            U = swap_row(U, cur_row_pos, pivot_row, verbosity=verbosity)
            P_elem = swap_row(I.copy(), cur_row_pos, pivot_row)
            P = P @ P_elem
            L = P_elem @ L @ P_elem

        # Eliminate the current column in rest of the rows below
        for row_idx in range(cur_row_pos+1, m):
            # reduce the row_idx iteratively via partial fractions to
            # prevent division by a possible 0 term
            tries = 0
            while U[row_idx, cur_col_pos] != 0:
                tries += 1
                if tries > max_tries:
                    print(f"ERROR: Max tries exceeded to reduce row {row_idx} with row {cur_row_pos}")
                    break
                scalar = sym.simplify(U[row_idx, cur_col_pos]/U[cur_row_pos, cur_col_pos])
                decomp = sym.apart(scalar) # partial fractions
                # simplify scalar is needed to check if apart has effect
                if scalar != decomp:
                    terms = decomp.args
                else:
                    # there is only 1 term (could be integer or proper fraction)
                    terms = [decomp]
                for term in terms:
                    _, d = sym.fraction(term)
                    # ensure denominator is non-zero so that reduction is valid
                    if not is_zero(d):
                        U = reduce_row(U, row_idx, term, cur_row_pos, verbosity=verbosity)
                        # L = sym.Matrix(L)
                        # L[row_idx, col_idx] = -term
                        # L = reduce_row(L, row_idx, -term, col_idx)
                        elem = reduce_row(I.copy(), row_idx, -term, cur_row_pos)
                        L = L @ elem

                # Cases where pivot row contains symbols such that scalar is a
                # fraction with symbolic denominator.
                # To reduce further, can only scale row_idx accordingly
                if U[row_idx, cur_col_pos] != 0:
                    scalar = sym.simplify(U[cur_row_pos, cur_col_pos]/U[row_idx, cur_col_pos])
                    n, d = sym.fraction(scalar)
                    # to scale by n, n can be symbols but not 0
                    # to scale by 1/d, d cannot be 0 both numerically or symbolically
                    if (n != 0) and (not is_zero(d)):
                        U = scale_row(U, row_idx, scalar, verbosity=verbosity)
                        elem = scale_row(I.copy(), row_idx, 1/scalar)
                        L = L @ elem

        cur_row_pos += 1

    if matrices == 1:
        return U
    elif matrices == 2:
        return P @ L, U
    elif matrices == 3:
        return P, L, U
    else:
        return U

def get_pivot_pos(A: sym.Matrix) -> sym.Matrix:
    assert A.is_echelon # check for REF
    
    pivot_pos = []
    cur_row_pos = 0
    for cur_col_pos in range(A.cols):
        pivot_row = get_pivot_row(A, cur_col_pos, cur_row_pos, follow_GE = False)
        
        if pivot_row != -1:
            pivot_pos.append((pivot_row, cur_col_pos))
            cur_row_pos += 1
        
    return pivot_pos

def get_pivot_elements(A: sym.Matrix) -> sym.Matrix:
    pivot_elements = []
    
    for i, j in get_pivot_pos(A):
        pivot_elements.append(A[i, j])
        
    return pivot_elements
    
def find_all_cases(A: sym.Matrix) -> sym.Matrix:
    cases = []
    for pivot in get_pivot_elements(A):
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

def evaluate_cases(A: sym.Matrix) -> sym.Matrix:
    cases = find_all_cases(A)
    all_possible_values = set(possible_val for case in cases for possible_val in case.items())
        
    for i, case in enumerate(cases, 1):
        print(f"Case {i}: {case}, not including {[dict([val]) for val in all_possible_values.symmetric_difference(set(case.items()))]}")
        display(ref(A.xreplace(case), verbosity = 0, matrices = 1, follow_GE = False))