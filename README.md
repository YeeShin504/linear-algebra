# linear-algebra

## About

This project builds on SymPy's Matrix class and is designed for students taking NUS MA1522 Linear Algebra for Computing. It has implementations of most of the algorithms taught during the course (as of Sem 1 AY24/25). Documentation of the code is still a work in progress. 

### Key Features

1. Import matrices from LaTeX directly. LaTeX expression can be obtained from Canvas by right-clicking on the equations.
2. Step-by-step workings for most algorithms (including LU Factorisation and SVD)

## Installation and Usage

### Installation

#### Prerequisites

This project is best supported in a Jupyter Notebook environment with Python 3.12+. You can download Python from [here](https://www.python.org/downloads/).

#### Clone the Repository

```bash
git clone https://github.com/YeeShin504/linear-algebra.git
cd linear-algebra
```

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

3. Install the required dependencies:
    ```bash
    pip install -r requirements.txt
    ```

### Usage

Create a Jupyter Notebook `test.ipynb`. Within the notebook, run the following code.
```python
from symbolic import *
A = Matrix([[1, 2, 3],
            [4, 5, 5],
            [7, 8, 9]])

b = Matrix([[1], 
            [2], 
            [3]])

augmented_matrix = A.aug_line().row_join(b)
display(augmented_matrix)
A.solve(rhs=b)
```

### FAQ

## Development

### Work in Progress

- [ ] Include better documentations for how to use the functions.
- [ ] Add conditions to generate specific matrices (either randomly or symbolically).
- [ ] Better `ref` algorithm to determine boundary conditions for a given matrix.
- [ ] Offline Optical Character Recognition (OCR) support to key in entries of matrices quickly and accurately during exam conditions.

### Credits

I would like to thank [@DenseLance](https://github.com/DenseLance) for his contributions.

### License

This project is licensed under the MIT License - see the [LICENSE](LICENSE.txt) file for details.