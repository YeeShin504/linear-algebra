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
        result = Matrix.from_latex(latex_input, row_join)
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

    @pytest.mark.parametrize(
        "vectors",
        [
            # Vectors with different lengths
            [Matrix([[1], [2]]), Matrix([[3], [4], [5]])],
            # Empty list
            [],
            # Mixed dimensions
            [Matrix([[1, 2]]), Matrix([[3], [4], [5]])],
        ],
    )
    def test_invalid_vectors(self, vectors):
        """Test invalid vector inputs raise appropriate errors"""
        with pytest.raises((ValueError, IndexError)):
            Matrix.from_list(vectors)


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
            num_rows=rows, num_cols=cols, symbol=symbol, is_real=is_real, shape=shape
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
        assert str(result[0, 0]) == "x_1,1"
        assert sym.re(result[0, 0]) == result[0, 0]  # Verify real number


class TestCreateRandMatrix:
    def test_create_rand_matrix(self):
        mat = Matrix.create_rand_matrix(2, 2)
        assert mat.shape == (2, 2)
        # Check if elements are numbers (not symbolic by default)
        assert all(isinstance(elem, (int, float)) for elem in mat.flat())


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
