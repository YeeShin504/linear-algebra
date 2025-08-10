# linear-algebra

## About

This project builds on SymPy's Matrix class and is designed for students taking NUS MA1522 Linear Algebra for Computing. It has implementations of most of the algorithms taught during the course (as of Sem 1 AY24/25). Documentation of the code is still a work in progress. 

### Key Features

1. Import matrices from $\rm\LaTeX$ directly. $\rm\LaTeX$ expression can be obtained from Canvas by right-clicking on the equations.
2. Step-by-step workings for most algorithms (including LU Factorisation and SVD)

## Installation and Usage

### Installation

#### Prerequisites

This project is best supported in a Jupyter Notebook environment with Python 3.10+. You can download Python from [here](https://www.python.org/downloads/).


#### Install Dependencies

It is recommended to use a virtual environment for managing dependencies.

1. Create a virtual environment:
    ```bash
    python -m venv venv
    ```

2. Activate the virtual environment:
    - On Windows:
      ```bash
      venv\Scripts\activate
      ```
    - On macOS/Linux:
      ```bash
      source venv/bin/activate
      ```

3. Install the library:
    ```bash
    pip install ma1522-linear-algebra
    ```
   It is recommended to use a Jupyter Notebook environment to run the code.
   ```bash
   pip install notebook
   ```

### Usage

Create a Jupyter Notebook `test.ipynb`. Within the notebook, run the following code.
```python
from ma1522 import *

# Create Matrix objects
A = Matrix([[1, 2, 3],
            [4, 5, 5],
            [7, 8, 9]])

b = Matrix([[1], 
            [2], 
            [3]])

# Join matrices along the columns via `row_join`. 
augmented_matrix = A.aug_line().row_join(b)

# `aug_line` adds a visual line that can be seen using `display`
display(augmented_matrix)

# Solution to the matrix equation Ax = b can be found using `solve`.
A.solve(rhs=b)

# Alternatively, the full steps with LU Factorisation can be found using `ref` with the appropriate options.
augmented_matrix.ref(verbosity=2)
```

Documentation of more functions can be found [here](https://yeeshin504.github.io/linear-algebra/).

More usage examples can be found under the examples folder. Proper documentation would be coming soon for all functions found in `symbolic.py`.

### FAQ

## Development

### Work in Progress

- [ ] Include better documentations for how to use the functions.
- [ ] Better `ref` algorithm to determine boundary conditions for a given matrix.

### Credits

I would like to thank [@DenseLance](https://github.com/DenseLance) for his contributions.

### License

This project is licensed under the MIT License - see the [LICENSE](LICENSE.txt) file for details.