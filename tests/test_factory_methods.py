"""Include the following methods in the tests
- from_latex
- from_list
- _shape
- create_unk_matrix
- create_rand_matrix
- eye
- zeros
- ones
- diag
- T
"""

import pytest
import sympy as sym

from ma1522 import Matrix
from ma1522.custom_types import Shape


class TestFromLatex:
    """Tests for Matrix.from_latex() factory method

    Covers:
    - Parsing matrices from different LaTeX environments (pmatrix, array)
    - Handling vector lists with row_join=True/False
    - Matrix multiplication expressions
    - Vector normalization
    - Edge cases (empty matrices, invalid input)
    - Special characters and symbols in matrices
    """

    def test_from_matrix(self):
        result = Matrix.from_latex(r"""
            \begin{pmatrix}
            1 & 2 \\
            3 & 4
            \end{pmatrix}    
        """)
        expected = Matrix([[1, 2], [3, 4]])
        assert result == expected

    @pytest.mark.parametrize(
        "latex_input",
        [r"\begin{pmatrix} \end{pmatrix}", r"\begin{array}{} \end{array}", r"{}", ""],
    )
    def test_empty_matrix(self, latex_input):
        """Test parsing empty matrix/vector inputs

        Verifies that various forms of empty input raise exceptions
        """
        with pytest.raises(Exception):
            Matrix.from_latex(latex_input)

    def test_single_element_matrix(self):
        """Test parsing a 1x1 matrix"""
        result = Matrix.from_latex(r"\begin{pmatrix} 5 \end{pmatrix}")
        expected = Matrix([[5]])
        assert result == expected

    def test_array_conversion(self):
        """Test conversion from array environment to pmatrix

        Verifies that array environments without column specifiers
        are properly converted to matrix format.
        """
        result = Matrix.from_latex(r"""
            \begin{array} 
            1 & 2 \\
            3 & 4 
            \end{array}
        """)
        expected = Matrix([[1, 2], [3, 4]])
        assert result == expected

    def test_array_conversion2(self):
        """Test conversion from array environment to pmatrix"""
        result = Matrix.from_latex(r"""
            \begin{array}{cc}
            1 & 2 \\
            3 & 4 
            \end{array}
        """)
        expected = Matrix([[1, 2], [3, 4]])
        assert result == expected

    def test_array_conversion3(self):
        """Test conversion from array environment to pmatrix"""
        result = Matrix.from_latex(r"""
            \begin{array}{} 
            1 & 2 \\
            3 & 4 
            \end{array}
        """)
        expected = Matrix([[1, 2], [3, 4]])
        assert result == expected

    def test_matmul_expression(self):
        """Test parsing a matrix multiplication expression"""
        result = Matrix.from_latex(r"""
            \begin{pmatrix} 
            1 & 2 \\
            3 & 4 
            \end{pmatrix}
            \begin{pmatrix} 
            5 & 6 \\
            7 & 8 
            \end{pmatrix}
        """)
        expected = Matrix([[1, 2], [3, 4]]) * Matrix([[5, 6], [7, 8]])
        assert result == expected

    @pytest.mark.parametrize(
        "latex_input,row_join,expected",
        [
            # row_join=True means vectors are treated as columns
            (
                r"\{ \begin{pmatrix} 1 \\ 3 \end{pmatrix}, \begin{pmatrix} 2 \\ 4 \end{pmatrix} \}",
                True,
                Matrix([[1, 2], [3, 4]]),
            ),
            # row_join=False means vectors are treated as rows
            (
                r"\{ \begin{pmatrix} 1 \\ 3 \end{pmatrix}, \begin{pmatrix} 2 \\ 4 \end{pmatrix} \}",
                False,
                Matrix([[1], [3], [2], [4]]),
            ),
        ],
    )
    def test_vector_list(self, latex_input, row_join, expected):
        """Test parsing a list of vectors"""
        result = Matrix.from_latex(expr=latex_input, row_join=row_join)
        assert result == expected

    def test_invalid_latex(self):
        """Test handling of invalid LaTeX input"""
        with pytest.raises(Exception):
            Matrix.from_latex(r"""
                \begin{pmatrix}
                1 & 2 \\
                3
                \end{pmatrix}
            """)

    @pytest.mark.parametrize(
        "latex_input,expected",
        [
            # Zero vector remains zero
            (r"\begin{pmatrix} 0 \\ 0 \end{pmatrix}", Matrix([[0], [0]])),
            # Already normalized vector
            (
                r"\begin{pmatrix} \sqrt{1/2} \\ \sqrt{1/2} \end{pmatrix}",
                Matrix(
                    [[sym.sqrt(sym.Rational(1, 2))], [sym.sqrt(sym.Rational(1, 2))]]
                ),
            ),
            # Matrix with mixed norms
            (
                r"\begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}",
                Matrix(
                    [
                        [sym.sqrt(sym.Rational(1, 10)), sym.sqrt(sym.Rational(1, 5))],
                        [
                            sym.Mul(3, sym.sqrt(sym.Rational(1, 10))),
                            sym.Mul(2, sym.sqrt(sym.Rational(1, 5))),
                        ],
                    ]
                ),
            ),
        ],
    )
    def test_vector_normalization(self, latex_input, expected):
        """Test vector normalization with various inputs

        Parameters:
            latex_input: LaTeX string to parse
            expected: Expected normalized matrix output
        """
        result = Matrix.from_latex(latex_input, norm=True)
        result.simplify()
        assert result == expected


class TestFromList:
    """Tests for Matrix.from_list() factory method

    Covers:
    - Creating matrices from column vectors (row_join=True)
    - Creating matrices from row vectors (row_join=False)
    - Edge cases (empty list, single vector)
    """

    @pytest.mark.parametrize(
        "vectors,row_join,expected",
        [
            # Column vectors
            ([Matrix([[1], [2]]), Matrix([[3], [4]])], True, Matrix([[1, 3], [2, 4]])),
            # Row vectors with row_join=False
            ([Matrix([[1, 2]]), Matrix([[3, 4]])], False, Matrix([[1, 2], [3, 4]])),
            # Single vector
            ([Matrix([[1], [2], [3]])], True, Matrix([[1], [2], [3]])),
        ],
    )
    def test_valid_vectors(self, vectors, row_join, expected):
        """Test creating matrix from valid vector lists"""
        result = Matrix.from_list(vectors, row_join)
        assert result == expected

    def test_from_list_non_mutation(self):
        """Verify from_list is independent of the input list and its contents (Regression)."""
        v1 = Matrix([1, 2])
        v2 = Matrix([3, 4])
        vecs = [v1, v2]
        res = Matrix.from_list(vecs)
        
        # Verify the list itself was not mutated (no .pop() occurred)
        assert len(vecs) == 2
        
        # Verify defensive copying: changing v1 should NOT change 'res'
        v1[0, 0] = 99
        assert res[0, 0] == 1, "Matrix should be independent of future mutations to the input vectors"

    @pytest.mark.parametrize(
        "vectors",
        [
            # Vectors with different lengths
            [Matrix([[1], [2]]), Matrix([[3], [4], [5]])],
            # Mixed dimensions
            [Matrix([[1, 2]]), Matrix([[3], [4], [5]])],
        ],
    )
    def test_invalid_vectors(self, vectors):
        """Test invalid vector inputs raise appropriate errors"""
        with pytest.raises((ValueError, IndexError)):
            Matrix.from_list(vectors)

    def test_empty_vectors_returns_empty_matrix(self):
        assert Matrix.from_list([]) == Matrix([])

    def test_augmented_position_is_preserved(self):
        result = Matrix.from_list([Matrix([1, 2]), Matrix([3, 4])], aug_pos=0)
        assert result == Matrix([[1, 3], [2, 4]], aug_pos=0)
        assert "|" in repr(result)


class TestShape:
    def test_diagonal_shape_pads_rectangular_matrix(self):
        wide = Matrix([[1, 2, 3], [4, 5, 6]])
        tall = Matrix([[1, 2], [3, 4], [5, 6]])

        assert wide._shape(Shape.DIAGONAL) == Matrix([[1, 0, 0], [0, 5, 0]])
        assert tall._shape(Shape.DIAGONAL) == Matrix([[1, 0], [0, 4], [0, 0]])

    def test_scalar_shape_rejects_non_square_matrix(self):
        with pytest.raises(sym.NonSquareMatrixError):
            Matrix([[1, 2, 3], [4, 5, 6]])._shape(Shape.SCALAR)

    def test_strict_and_symmetric_shapes(self):
        mat = Matrix([[1, 2], [3, 4]])

        assert mat._shape(Shape.STRICT_UPPER) == Matrix([[0, 2], [0, 0]])
        assert mat._shape(Shape.STRICT_LOWER) == Matrix([[0, 0], [3, 0]])
        assert mat._shape(Shape.SYMMETRIC) == Matrix([[1, 2], [2, 4]])

    def test_symmetric_shape_rejects_non_square_matrix(self):
        with pytest.raises(sym.NonSquareMatrixError):
            Matrix([[1, 2, 3], [4, 5, 6]])._shape(Shape.SYMMETRIC)


class TestCreateUnkMatrix:
    """Tests for Matrix.create_unk_matrix() factory method

    Covers:
    - Creating matrices with default parameters
    - Custom dimensions and symbols
    - Real vs complex entries
    - Shape constraints (diagonal, triangular, etc.)
    """

    @pytest.mark.parametrize(
        "rows,cols,symbol,is_real,shape",
        [
            # Default parameters
            (1, 1, "x", True, None),
            # Custom dimensions
            (3, 2, "a", True, None),
            # Complex entries
            (2, 2, "z", False, None),
            # Diagonal shape
            (3, 3, "d", True, Shape.SCALAR),
            # Upper triangular
            (3, 3, "u", True, Shape.UPPER),
            # Lower triangular
            (3, 3, "l", True, Shape.LOWER),
        ],
    )
    def test_create_unk_matrix(self, rows, cols, symbol, is_real, shape):
        """Test creating unknown matrices with various parameters"""
        result = Matrix.create_unk_matrix(
            r=rows, c=cols, symbol=symbol, is_real=is_real, shape=shape
        )

        # Verify dimensions
        assert result.shape == (rows, cols)

        # Verify symbol naming pattern
        assert all(
            str(entry).startswith(f"{symbol}_") for entry in result.flat() if entry != 0
        )

        # Verify real/complex type
        if is_real:
            assert all(entry.is_real for entry in result.flat())

        # Verify shape constraints if specified
        if shape:
            expected = result._shape(shape)
            assert result == expected

    def test_default_parameters(self):
        """Test that default parameters create expected matrix"""
        result = Matrix.create_unk_matrix()
        assert result.shape == (1, 1)
        assert str(result[0, 0]) == "x"
        assert sym.re(result[0, 0]) == result[0, 0]  # Verify real number


class TestCreateRandMatrix:
    def test_create_rand_matrix(self):
        mat = Matrix.create_rand_matrix(2, 2, seed=42)
        assert mat.shape == (2, 2)
        assert mat == Matrix([[81, 14], [3, 94]])

    def test_create_rand_matrix_with_shape(self):
        mat = Matrix.create_rand_matrix(2, 2, shape=Shape.STRICT_UPPER, seed=42)
        assert mat == Matrix([[0, 14], [0, 0]])

class TestApplyVander:
    """Regression tests for Vandermonde matrix applications."""
    def test_create_vander(self):
        result = Matrix.create_vander(2, 4)
        assert result.shape == (2, 4)
        assert [[str(entry) for entry in row] for row in result.tolist()] == [
            ["1", "x_1", "x_1**2", "x_1**3"],
            ["1", "x_2", "x_2**2", "x_2**3"],
        ]

    def test_basic_substitution(self):
        V = Matrix.create_vander(3, 3)
        x_vec = Matrix([[2], [3], [5]])
        result = V.apply_vander(x_vec)
        assert result == Matrix([[1, 2, 4], [1, 3, 9], [1, 5, 25]])
        assert result.free_symbols == set()

    def test_free_symbols_not_mutated(self):
        V = Matrix.create_vander(3, 3)
        syms_before = frozenset(V.free_symbols)
        x_vec = Matrix([[2], [3], [5]])
        V.apply_vander(x_vec)
        assert frozenset(V.free_symbols) == syms_before

    def test_apply_vander_rejects_non_column_vector(self):
        with pytest.raises(sym.ShapeError):
            Matrix.create_vander(2, 2).apply_vander(Matrix([[1, 2]]))

    def test_apply_vander_rejects_row_mismatch(self):
        with pytest.raises(sym.ShapeError):
            Matrix.create_vander(2, 2).apply_vander(Matrix([[1], [2], [3]]))


class TestOverriddenFactoryMethods:
    def test_eye(self):
        mat = Matrix.eye(3)
        assert mat == sym.eye(3)
        assert isinstance(mat, Matrix)

    def test_zeros(self):
        mat = Matrix.zeros(2, 3)
        assert mat == sym.zeros(2, 3)
        assert isinstance(mat, Matrix)

    def test_ones(self):
        mat = Matrix.ones(2, 2)
        assert mat == sym.ones(2, 2)
        assert isinstance(mat, Matrix)

    def test_diag(self):
        mat = Matrix.diag(1, 2, 3)
        assert mat == sym.diag(1, 2, 3)
        assert isinstance(mat, Matrix)

    def test_T_property(self):
        mat = Matrix([[1, 2], [3, 4]])
        assert mat.T == Matrix([[1, 3], [2, 4]])
        assert isinstance(mat.T, Matrix)

    def test_H_property(self):
        mat = Matrix([[1, 2*sym.I], [3+4*sym.I, 4]])
        assert mat.H == Matrix([[1, 3-4*sym.I], [-2*sym.I, 4]])
        assert isinstance(mat.H, Matrix)
