# Tutorial: Getting Started

This guide is designed for undergraduate students who are new to Python and want to use this library for linear algebra computations. The `symbolic.Matrix` class is built on top of SymPy, a powerful Python library for symbolic mathematics, and is tailored for the NUS MA1522 course (AY 24/25 Sem 1).

## 1. Installation and Setup

First, you need to set up your Python environment. We recommend using a Jupyter Notebook for an interactive experience.

!!! tip "Prerequisites"
    - Python 3.10+
    - Jupyter Notebook or JupyterLab

Follow these steps to get everything set up:

1.  **Create and activate a virtual environment:**
    - On Windows:
      ```bash
      python -m venv venv
      venv\Scripts\activate
      ```
    - On macOS/Linux:
      ```bash
      source venv/bin/activate
      ```

2.  **Install dependencies:**
    ```bash
    pip install ma1522-linear-algebra notebook
    ```

3.  **Start Jupyter:**
    ```bash
    jupyter notebook
    ```

Now, create a new notebook and you're ready to go!

## 2. Creating Matrices

Let's start by importing the necessary functions and creating our first matrices.

```python
from symbolic import Matrix
from utils import display
```

### From a List of Lists

The most straightforward way to create a matrix is from a list of lists, where each inner list represents a row.

```python
A = Matrix([[1, 2, 3],
            [4, 5, 6],
            [7, 8, 9]])
display(A)
```

### From LaTeX

A key feature of this library is the ability to create a matrix directly from a LaTeX string. This is incredibly useful for copying matrices from online resources or textbooks.

```python
latex_expr = r'\begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}'
B = Matrix.from_latex(latex_expr)
display(B)
```

### Special Matrices

You can also create special matrices easily.

- **Identity Matrix:**
  ```python
  I = Matrix.eye(3)
  display(I)
  ```

- **Zero Matrix:**
  ```python
  Z = Matrix.zeros(2, 3)
  display(Z)
  ```

- **Symbolic Matrix:** Create a matrix with symbolic entries.
  ```python
  S = Matrix.create_unk_matrix(2, 2, 's')
  display(S)
  ```

## 3. Basic Operations

The `Matrix` class supports standard matrix operations.

```python
A = Matrix([[1, 2], [3, 4]])
B = Matrix([[5, 6], [7, 8]])

# Addition
display(A + B)

# Scalar Multiplication
display(2 * A)

# Matrix Multiplication
display(A @ B)

# Transpose
display(A.T)
```

## 4. Solving Linear Systems

Let's solve the matrix equation \(Ax = b\).

```python
A = Matrix([[1, 2, 3],
            [4, 5, 5],
            [7, 8, 9]])

b = Matrix([[1],
            [2],
            [3]])

# Solve directly
x = A.solve(rhs=b)
display(x)
```

### Using Row Reduction

Alternatively, you can use row reduction on an augmented matrix.

1.  **Create an augmented matrix:**
    ```python
    augmented_matrix = A.aug_line().row_join(b)
    display(augmented_matrix)
    ```
    The `aug_line()` method adds a visual separator.

2.  **Compute the Reduced Row Echelon Form (RREF):**
    ```python
    rref_matrix, pivots = augmented_matrix.rref()
    display(rref_matrix)
    ```

3.  **Step-by-step Row Echelon Form (REF):**
    For a detailed, step-by-step reduction, use the `ref()` method with `verbosity`.

    ```python
    # Create a new augmented matrix for demonstration
    augmented_matrix_2 = A.aug_line().row_join(b)

    # verbosity=1 shows the operations
    # verbosity=2 shows the matrix at each step
    plu = augmented_matrix_2.ref(verbosity=2)
    ```

## 5. Advanced Topics

The library provides functions for various advanced linear algebra concepts.

### Eigenvalues and Eigenvectors

```python
A = Matrix([[4, -1, 6],
            [2,  1, 6],
            [2, -1, 8]])

# Eigenvalues
eigenvals = A.eigenvals()
print(eigenvals)

# Eigenvectors
eigenvects = A.eigenvects()
display(eigenvects)
```

### Diagonalization

You can check if a matrix is diagonalizable and perform the diagonalization.

```python
# Perform diagonalization (P, D)
P, D = A.diagonalize()
display(P)
display(D)

# Verify A = PDP^-1
display(P @ D @ P.inv())
```

### Orthogonal Diagonalization

For symmetric matrices, you can perform orthogonal diagonalization.

```python
S = Matrix([[5, -2, 0],
            [-2, 8, 0],
            [0, 0, 4]])

if S.is_orthogonally_diagonalizable:
    P, D = S.orthogonally_diagonalize()
    display(P)
display(D)
    # Verify S = PDP^T
display(P @ D @ P.T)
```

### Singular Value Decomposition (SVD)

```python
A = Matrix([[1, 2], [2, 2], [2, 1]])
U, S, V = A.singular_value_decomposition()

print("U:")
display(U)
print("S:")
display(S)
print("V:")
display(V)

# Verify A = U * S * V.T
display(U @ S @ V.T)
```

## 6. Subspace Analysis

Explore fundamental subspaces of a matrix.

```python
A = Matrix([[1, 2, 3, 4],
            [5, 6, 7, 8],
            [9, 10, 11, 12]])

# Column Space
col_space = A.columnspace()
display(Matrix.from_list(col_space))

# Null Space
null_space = A.nullspace()
display(Matrix.from_list(null_space))

# Row Space
row_space = A.rowspace()
display(Matrix.from_list(row_space, row_join=False))
```

---

This tutorial covers the core functionalities of the `symbolic.Matrix` class. For more details on specific functions, you can refer to the code references.