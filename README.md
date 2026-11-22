# linear-algebra (TI-84 Plus CE edition)

TI-BASIC programs for **NUS MA1522 Linear Algebra for Computing**, for the
**TI-84 Plus CE** (and TI-84 Plus CE Python / TI-83 Premium CE). This is a port of the
[`ma1522` Python library](https://github.com/YeeShin504/linear-algebra), which is built on
SymPy, into programs that run on the calculator itself.

Each program works on matrices you type into the calculator's matrix editor (`[A]`, `[B]`)
and stores its results in the other matrix variables, so you can view them or keep
calculating with them afterwards. Results are shown as fractions when possible.

**Contents**

* [Installing](#installing)
* [Basics](#basics): entering a matrix, running a program, moving through the screens
* [Quick reference](#quick-reference)
* Program guides
  * Systems and matrices: [RREF](#rref) · [SOLVE](#solve) · [LU](#lu) · [INVERSE](#inverse) · [CRAMER](#cramer)
  * Subspaces: [SPACES](#spaces) · [COLCON](#colcon) · [SUBSP](#subsp) · [COORD](#coord) · [LINTRAN](#lintran)
  * Orthogonality: [ORTHO](#ortho) · [GRAMSCH](#gramsch) · [LEASTSQ](#leastsq) · [POLYFIT](#polyfit)
  * Eigenvalues: [EIGEN](#eigen) · [ORTDIAG](#ortdiag) · [SVD](#svd) · [MARKOV](#markov)
  * Multiple choice: [TESTVEC](#testvec)
  * [The MA1522 menu](#the-ma1522-menu)
* [Things to know](#things-to-know)
* [Troubleshooting](#troubleshooting)
* [For developers](#for-developers)

## Installing

1. Install [TI Connect CE](https://education.ti.com/en/products/computer-software/ti-connect-ce-sw).
2. Connect the calculator, open [`bin/`](bin) in Finder, select **every** `.8xp` file (⌘A)
   and drag them onto TI Connect CE together (send to RAM). One drag sends them all.
   The `Z*` programs are helpers that the others call, so they must be on the calculator too.
3. On the calculator: `prgm` → `MA1522` → `enter` for the menu, or run any program directly.

## Basics

### Entering a matrix

1. Press `2nd` `x⁻¹` (MATRIX), move right to **EDIT**, and pick **1: [A]**.
2. Type the number of rows, `enter`, the number of columns, `enter`.
3. Type each entry followed by `enter`. Fractions such as `1/3` are fine. Use the `(-)` key
   for negative numbers, not the `−` (subtract) key.
4. Press `2nd` `mode` (QUIT) when you are done.

Do the same with **2: [B]** for programs that need a second matrix (a right-hand side, a
second set of vectors, an initial state).

In these guides a matrix is written row by row, the way the calculator writes it:
`[[1,2][3,4]]` is the 2×2 matrix with rows (1, 2) and (3, 4).

**Vectors go in as columns.** When a program works with a set of vectors (a basis, a
spanning set, inputs of a linear transformation), each vector is one **column** of the matrix.

### Running a program

Press `prgm`, choose the program from the list, and press `enter` twice (the first `enter`
pastes `prgmNAME` on the home screen, the second runs it). Or run `prgmMA1522` for a menu
of everything.

### Moving through the screens

* **`SHOW STEPS 1/0?`**: type `1` then `enter` to see every row operation, or `0` to skip
  straight to the answer. Other yes/no questions work the same way.

* **Most screens wait for you to press `enter`.** While a program is waiting, the small busy
  indicator in the top-right corner shows a dotted pattern. That is not a loading screen: the
  program is paused until you press `enter`.

* **Big matrices:** while a matrix is shown and paused, the arrow keys scroll it.
* **To stop a program** at any point, press `on`, then choose **1: Quit**.

### Looking at results afterwards

Every result is also saved, so you can look at it again or use it in later calculations.

* **Matrices:** `2nd` `x⁻¹` → **NAMES** → pick `[C]`, `[D]`, … → `enter`. To see fractions,
  add `►Frac` before `enter` (`math` → **1: ►Frac**).

* **Named lists** such as `ʟPIV` or `ʟEV`: `2nd` `stat` (LIST) → **NAMES** → pick it → `enter`.
* **Real variables** such as `R` (rank): `alpha` + letter → `enter`.
* **Reuse a result:** store it into `[A]` with the `sto→` key, e.g. `[C]` `sto→` `[A]` `enter`,
  then run another program on it.

## Quick reference

| Program   | What it does | Input | Output |
|-----------|--------------|-------|--------|
| `RREF`    | REF / RREF, every row operation shown if you ask for steps, rank, pivot columns | `[A]` | `[C]` = REF/RREF, `ʟPIV` |
| `SOLVE`   | General solution of `AX = B`: unique / infinitely many / inconsistent | `[A]`, `[B]` | `[C]` particular solution, `[D]` null space basis, `[E]` RREF of `[A\|B]` |
| `LU`      | LU factorisation by Gaussian elimination, `A = PLU` (with row swaps if needed) | `[A]` | `[C]` = L, `[D]` = U, `[E]` = P |
| `INVERSE` | Determinant, inverse from RREF of `[A\|I]` (with steps), adjoint and cofactor matrix | `[A]` | `[C]` = A⁻¹, `[D]` = adj A, `[E]` = cofactors |
| `CRAMER`  | Cramer's rule, shows each det(Aᵢ) | `[A]`, `[B]` | `[C]` = x |
| `SPACES`  | Rank, nullity, bases for column / row / null space, linear independence, spanning | `[A]` | `[C]` col space, `[D]` row space, `[E]` null space, `[F]` RREF |
| `COLCON`  | Conditions on **b** for `Ax = b` to be consistent | `[A]` | `[C]`: consistent ⇔ `[C]b = 0` |
| `SUBSP`   | span A ⊆ span B?, span A = span B?, intersection, extend to a basis of ℝᵐ, simplify a basis | `[A]`, `[B]` | `[C]` |
| `COORD`   | Coordinate vectors `[v]_S`, transition matrix from basis S to basis T | `[A]`, `[B]` | `[C]` |
| `LINTRAN` | Standard matrix of T from T(vᵢ) = wᵢ | `[A]` = vᵢ as columns, `[B]` = wᵢ as columns | `[C]` |
| `ORTHO`   | Orthogonal / orthonormal / orthogonal-matrix checks, orthogonal complement, normalise columns | `[A]` | `[C]` |
| `GRAMSCH` | Gram–Schmidt process, squared norms (so you can write uᵢ = vᵢ/√n exactly), QR factorisation with R also in exact form | `[A]` | `[C]` orthogonal, `[D]` orthonormal = Q, `[E]` = R, `[F]` = exact R numerators, `ʟNRM` |
| `LEASTSQ` | Least squares solution (`AᵀAx = Aᵀb`), projection of b onto col A, orthogonal part, distance | `[A]`, `[B]` | `[C]` x̂, `[D]` null basis, `[E]` projection, `[F]` b − projection |
| `POLYFIT` | Fit a polynomial with a Vandermonde matrix (exact interpolation or least squares) | x in `L₁`, y in `L₂`, degree | `[C]` coefficients a₀…a_d, `[E]` Vandermonde |
| `EIGEN`   | Characteristic polynomial, eigenvalues, algebraic and geometric multiplicities, eigenspace bases, diagonalisation, P⁻¹ | `[A]` | `ʟCP`, `ʟEV`, `ʟAM`, `ʟGM`, `[C]` = P, `[D]` = D, `[E]` = P⁻¹ |
| `ORTDIAG` | Orthogonal diagonalisation of a symmetric matrix, `A = PDPᵀ` | `[A]` | `[C]` = P, `[D]` = D, `[E]` = P before normalising |
| `SVD`     | Singular value decomposition `A = UΣVᵀ`, with the unnormalised vᵢ, Avᵢ and extra uᵢ shown in fractions for exact answers | `[A]` | `[C]` = U, `[D]` = Σ, `[E]` = V, `ʟSV` |
| `MARKOV`  | Stochastic matrix check, equilibrium vector(s), state after k steps | `[A]`, optional `[B]` = x₀ | `[C]` equilibrium, `[D]` = Aᵏx₀ |
| `TESTVEC` | Tests several candidate vectors at once: Av, in Null(A)?, eigenvector (and λ)?, in Col(A)? | `[A]`, `[B]` = the vectors as columns | screen |
| `MA1522`  | Menu that launches everything above | | |

Helpers (called by the programs above; you do not need to run them yourself):
`ZRREF` (row reduction with optional step display), `ZNULL` (null space from an RREF),
`ZGS` (Gram–Schmidt), `ZEIG` (characteristic polynomial and eigenvalues).

---

## Systems and matrices

### RREF

Row-reduces `[A]` to reduced row-echelon form (or stops at row-echelon form), and gives the
rank, nullity and pivot columns.

**Enter:** `[A]` = the matrix.

**Prompts:**

* `SHOW STEPS 1/0?`: `1` shows every row operation.
* `REF ONLY 1/0?`: `1` stops at row-echelon form (Gaussian elimination only), `0` goes on to RREF.

**Screens:** each step (if on), then the REF/RREF, then rank, nullity and pivot columns.

**Example:** `[A]` = `[[1,2,3][4,5,6][7,8,9]]`, steps on, REF only:

```text
R2 = R2 + C*R1         R3 = R3 + C*R1         R3 = R3 + C*R2
WHERE C = -4           WHERE C = -7           WHERE C = -2
[ 1  2  3]             [  1   2   3]          [ 1  2  3]
[ 0 -3 -6]             [  0  -3  -6]          [ 0 -3 -6]
[ 7  8  9]             [  0  -6 -12]          [ 0  0  0]

RANK = 2   NULLITY = 1   PIVOT COLUMNS {1,2}
```

With `REF ONLY` = `0` the answer is the RREF `[[1,0,-1][0,1,2][0,0,0]]`.

**Saved:** `[C]` = REF/RREF, `ʟPIV` = pivot columns, `R` = rank.

**Steps** are written as `Ri = Ri + C*Rj` (add C times row j to row i), `Ri = C*Ri`
(scale a row) and `Ri <-> Rj` (swap). Copy them into your working as row operations.

### SOLVE

Solves `AX = B` and says whether the system has a unique solution, infinitely many, or none.
`[B]` can have several columns (it solves for each column at once).

**Enter:** `[A]` = coefficient matrix, `[B]` = right-hand side (a column).

**Prompts:** `SHOW STEPS 1/0?` (row reduction of `[A|B]`).

**Screens:**

* `UNIQUE SOLUTION`: the solution `X = [C]`.
* `INFINITELY MANY`: `X = [C] + s1*D1 + ... + sF*DF`, with F (the number of parameters), the
  particular solution `[C]` and the null space basis `[D]` (its columns are D1…DF).

* `INCONSISTENT`: no solution. It shows the RREF of `[A|B]` so you can see the bad row.

**Example:** `[A]` = `[[1,2,3][2,4,7]]`, `[B]` = `[[1][2]]`, steps on:

```text
SHOW STEPS 1/0? 1
R2 = R2 + C*R1      WHERE C = -2
[1 2 3 1]
[0 0 1 0]
R1 = R1 + C*R2      WHERE C = -3
[1 2 0 1]
[0 0 1 0]
INFINITELY MANY
X=[C]+S1*D1+...+SF*DF       F = 1
PARTICULAR [C] = [1 0 0]ᵀ
NULL SPACE BASIS [D] = [-2 1 0]ᵀ
```

So x = (1, 0, 0) + s(−2, 1, 0), s ∈ ℝ.

**Saved:** `[C]` = particular solution (free variables set to 0), `[D]` = null space basis,
`[E]` = RREF of `[A|B]`, `F` = number of parameters.

### LU

LU factorisation by Gaussian elimination: `A = LU`, or `A = PLU` if rows have to be swapped.
`[A]` does not have to be square.

**Enter:** `[A]`.

**Prompts:** `SHOW STEPS 1/0?`

**Screens:** each elimination step (if on), then whether row swaps were needed (and P if
so), then L, then U.

**Example:** `[A]` = `[[2,4,1][4,9,3][-2,-1,7]]`, steps on:

```text
R2 = R2 - C*R1         R3 = R3 - C*R1         R3 = R3 - C*R2
C = L(2,1) = 2         C = L(3,1) = -1        C = L(3,2) = 3
[ 2  4  1]             [2 4 1]                [2 4 1]
[ 0  1  1]             [0 1 1]                [0 1 1]
[-2 -1  7]             [0 3 8]                [0 0 5]

NO ROW SWAPS   A = LU, P = I
L = [[1,0,0][2,1,0][-1,3,1]]      U = [[2,4,1][0,1,1][0,0,5]]
```

Each multiplier C is the entry of L it goes into, so you can read L off the steps.

**Saved:** `[C]` = L, `[D]` = U, `[E]` = P.

### INVERSE

Determinant, inverse (found by row-reducing `[A|I]`), adjoint and cofactor matrix.

**Enter:** `[A]` (square).

**Prompts:** `SHOW STEPS 1/0?` (row reduction of `[A|I]`).

**Screens:** det(A), then A⁻¹ (or `SINGULAR (RANK<N)`, `NO INVERSE` with the RREF of `[A|I]`), then
adj(A), then the cofactor matrix.

**Example:** `[A]` = `[[1,2][3,4]]`:

```text
DET(A) = -2
A^-1 = [C]          ADJ(A) = [D]       COFACTOR MATRIX [E]
[  -2    1]         [ 4 -2]            [ 4 -3]
[ 3/2 -1/2]         [-3  1]            [-2  1]
```

**Saved:** `O` = det(A), `[C]` = A⁻¹, `[D]` = adj(A), `[E]` = cofactor matrix
(adj(A) is its transpose, and A⁻¹ = adj(A)/det(A)).

### CRAMER

Cramer's rule: xᵢ = det(Aᵢ)/det(A), where Aᵢ is A with column i replaced by b.

**Enter:** `[A]` (square, det ≠ 0), `[B]` = b (one column).

**Screens:** det(A), then det(Aᵢ) and xᵢ for each i, then the solution.

**Example:** `[A]` = `[[2,1,-1][-3,-1,2][-2,1,2]]`, `[B]` = `[[8][-11][-3]]`:

```text
DET(A) = -1
DET(A1) = -2   X1 = 2
DET(A2) = -3   X2 = 3
DET(A3) = 1    X3 = -1
SOLUTION X = [C] = [2 3 -1]ᵀ
```

**Saved:** `[C]` = x.

---

## Subspaces

### SPACES

Everything about the four-subspaces picture of `[A]`: rank, nullity, bases for the column
space, row space and null space, and whether the columns are linearly independent / span ℝᵐ.

**Enter:** `[A]` (for a set of vectors, put them in as columns).

**Screens:** rank, nullity and pivot columns; independence / spanning; then the three bases.

**Example:** `[A]` = `[[1,2,0,1][2,4,1,4][3,6,2,7]]`:

```text
RANK = 2   NULLITY = 2   PIVOT COLUMNS {1,3}
COLUMNS ARE LINEARLY DEPENDENT
COLUMNS DO NOT SPAN R^M
COL SPACE BASIS [C]    ROW SPACE BASIS [D]    NULL SPACE BASIS [E]
[1 0]                  [1 2 0 1]              [-2 -1]
[2 1]                  [0 0 1 2]              [ 1  0]
[3 2]                                         [ 0 -2]
                                              [ 0  1]
```

The column space basis is the pivot columns of the **original** A (columns 1 and 3 here).
The row space basis is the non-zero rows of the RREF, written as rows.

**Saved:** `[C]` = column space basis, `[D]` = row space basis (rows), `[E]` = null space
basis (columns), `[F]` = RREF of A, `R` = rank, `F` = nullity.

### COLCON

Finds the conditions on b = (b₁, …, bₘ) for `Ax = b` to be consistent, i.e. the equations
that describe the column space of A.

**Enter:** `[A]`.

**Screens:** either `ALWAYS CONSISTENT` (when rank A = number of rows), or a matrix `[C]`:
the system is consistent exactly when `[C]·b = 0`, so each row of `[C]` is one condition.

**Example:** `[A]` = `[[1,2][2,4][3,6]]`:

```text
CONSISTENT IFF [C]*B = 0
[   1    0 -1/3]
[   0    1 -2/3]
```

Read each row as an equation: b₁ − b₃/3 = 0 and b₂ − 2b₃/3 = 0.

**Saved:** `[C]` = the conditions (one per row).

### SUBSP

Compares the subspaces V = span(columns of `[A]`) and W = span(columns of `[B]`).

**Enter:** `[A]` (vectors of V as columns), and `[B]` (vectors of W as columns) for options 1–3.

**Menu:**

| Option | Question | Result |
|--------|----------|--------|
| 1 `SPAN A IN SPAN B?` | Is V ⊆ W? | message |
| 2 `SPAN A = SPAN B?` | Is V = W? | message |
| 3 `INTERSECTION` | Basis for V ∩ W | `[C]` (columns), dimension |
| 4 `EXTEND TO R^M` | Extend the columns of A to a basis of ℝᵐ | `[C]`, which standard vectors eᵢ were added |
| 5 `SIMPLIFY BASIS` | A simpler basis for V (from the RREF of Aᵀ) | `[C]` (columns), dimension |

Options 1 and 2 also show the ranks of A, B, `[B|A]` and `[A|B]`, which is the reasoning you
write down: V ⊆ W exactly when rank `[B|A]` = rank B.

**Examples:**

```text
2: [A] = [[1,0][1,1][0,1]], [B] = [[1,1][2,1][1,0]]
   RANK A, B, [B|A], [A|B] = {2,2,2,2}   => SAME SUBSPACE

3: [A] = [[1,0][0,1][0,0]], [B] = [[0,1][1,0][0,1]]
   BASIS OF INTERSECTION = [0 1 0]ᵀ,  DIM = 1

4: [A] = [[1][1][0]]
   ADDED STD VECTORS E_I, I = 1, 3     [C] = [[1,1,0][1,0,0][0,0,1]]

5: [A] = [[1,2,3][1,2,4][1,2,5]]
   SIMPLER BASIS = columns (1,0,-1) and (0,1,2),  DIM = 2
```

For option 4, if the columns of A are not linearly independent you get a note, and only
their pivot columns are kept.

**Saved:** `[C]`.

### COORD

Coordinate vectors relative to a basis, and transition matrices between two bases.

**Enter:** `[A]` = basis S (vectors as columns), `[B]` = the vectors, or the second basis T.

**Menu:**

* **1 `COORDS OF B IN A`**: `[v]_S` for every column v of `[B]`, as the columns of `[C]`.
* **2 `TRANSITION A TO B`**: the transition matrix P from S to T, so `[v]_T = P [v]_S`.

**Examples:**

```text
1: [A] = [[1,1][1,-1]], [B] = [[3][1]]
   RESULT [C] = [2 1]ᵀ          (3,1) = 2(1,1) + 1(1,-1)

2: [A] = S = [[1,1][1,-1]], [B] = T = [[1,0][1,1]]
   RESULT [C] = [[1,1][0,-2]]   the transition matrix from S to T
```

For the transition matrix from T to S, swap `[A]` and `[B]` (or take the inverse of `[C]`).
You get an error message if the basis columns are not linearly independent, or if some
vector is not in their span.

**Saved:** `[C]`.

### LINTRAN

Finds the standard matrix of a linear transformation T from the values it takes on some
vectors: T(v₁) = w₁, …, T(vₖ) = wₖ.

**Enter:** `[A]` = v₁ … vₖ as columns, `[B]` = w₁ … wₖ as columns (same order).

**Screens:** the standard matrix `[C]`, with T(x) = `[C]`x.

**Example:** T(1,0) = (2,1) and T(1,1) = (3,3): `[A]` = `[[1,1][0,1]]`, `[B]` = `[[2,3][1,3]]`:

```text
STANDARD MATRIX [C]
[2 1]
[1 2]
```

If the vᵢ do not span the whole space, T is not unique: you get `NOT UNIQUE` and one
possible answer. If the vᵢ are dependent but the wᵢ do not match that dependence,
there is no such linear T, and the program says so.

**Saved:** `[C]` = standard matrix.

---

## Orthogonality

### ORTHO

Orthogonality checks and orthogonal complements for the columns of `[A]`.

**Enter:** `[A]` = the vectors as columns.

**Menu:**

* **1 `ORTHOGONAL CHECK`**: shows AᵀA (all the dot products uᵢ·uⱼ), then whether the
  columns are orthogonal, orthonormal, and (if A is square and orthonormal) whether A is an
  orthogonal matrix (A⁻¹ = Aᵀ).

* **2 `ORTH COMPLEMENT`**: a basis for the orthogonal complement of the column space,
  which is the null space of Aᵀ.

* **3 `NORMALISE COLS`**: each column divided by its length.

**Examples:**

```text
1: [A] = [[1,1][1,-1]]       A^T*A = [[2,0][0,2]]  COLS ARE ORTHOGONAL, NOT ORTHONORMAL
2: [A] = [[1,0][1,1][1,2]]   ORTH COMPLEMENT: DIM = 1, basis (1,-2,1)
3: [A] = [[3,0][4,2]]        NORMALISED COLS [C] = [[3/5,0][4/5,1]]
```

To find the orthogonal complement of a set of vectors, put them in as the columns of `[A]`
and use option 2.

**Saved:** `[C]`.

### GRAMSCH

The Gram–Schmidt process on the columns of `[A]`, and the QR factorisation A = QR.

**Enter:** `[A]` = the vectors u₁ … uₖ as columns, in order.

**Screens:** the orthogonal basis v₁ … vₖ, their squared lengths, the orthonormal basis
(= Q), then R, then R in exact form.

**Example:** `[A]` = `[[1,1][1,0][0,1]]`:

```text
ORTHOGONAL BASIS [C]       SQUARED NORMS |VI|^2
[   1  1/2]                {2, 3/2}
[   1 -1/2]
[   0    1]
ORTHONORMAL BASIS [D] (decimals)      R = [E] (decimals)
```

The orthonormal vectors involve square roots, so they show as decimals. Use the squared
norms to write them exactly: here u₁ = (1,1,0)/√2 and u₂ = (1/2,−1/2,1)/√(3/2).

The last screen gives R exactly. `[F]` = VᵀA holds the numerators vᵢ·aⱼ in fractions, and
each entry of R is rᵢⱼ = F(i,j)/√‖vᵢ‖², using the squared norms from before:

```text
EXACT FORM, NI=|VI|^2:
R(I,J)=F(I,J)/SQRT(NI)
F = [F]
[  2   1]                so R = [[2/√2, 1/√2][0, (3/2)/√(3/2)]] = [[√2, 1/√2][0, √(3/2)]]
[  0 3/2]
```

If a column depends on the earlier ones it is dropped (its vᵢ would be 0), and there is
then no QR factorisation.

**Saved:** `[C]` = orthogonal basis, `ʟNRM` = squared norms, `[D]` = orthonormal basis = Q,
`[E]` = R, `[F]` = VᵀA (exact numerators of R).

### LEASTSQ

Least squares solution of `Ax = b` (by solving the normal equation AᵀAx = Aᵀb), the
projection of b onto the column space of A, and the distance from b to it.

**Enter:** `[A]`, `[B]` = b (one column).

**Screens:** the least squares solution x̂ (unique, or a general solution with parameters),
the projection Ax̂, the part of b orthogonal to col A, and the distance ‖b − Ax̂‖ (and its
square).

**Example:** `[A]` = `[[1,0][1,1][1,2]]`, `[B]` = `[[6][0][0]]`:

```text
UNIQUE LSQ SOLUTION   X = [5 -3]ᵀ
PROJECTION OF B ONTO COL(A) = [5 2 -1]ᵀ
B - PROJ = [1 -2 1]ᵀ
DIST = 2.449489743   SQUARED = 6        (so the distance is √6)
```

**Saved:** `[C]` = x̂, `[D]` = null space basis of AᵀA (when x̂ is not unique),
`[E]` = projection of b, `[F]` = b − projection.

### POLYFIT

Fits a polynomial y = a₀ + a₁x + … + a_d xᵈ to data points, using the Vandermonde matrix
and least squares. With d = (number of points − 1) the curve goes through every point exactly.

**Enter:** the x values in `L₁` and the y values in `L₂`: `stat` → **1: Edit…**, then type
the numbers into the `L1` and `L2` columns.

**Prompts:** `DEGREE?` = d.

**Screens:** the coefficients a₀ … a_d, lowest power first.

**Example:** points (0,1), (1,3), (2,7), (3,13), degree 2:

```text
Y=A0+A1X+...+ADX^D
[A0,...,AD] = [C] = [1 1 1]ᵀ        so y = 1 + x + x²
```

**Saved:** `[C]` = coefficients a₀ … a_d, `[E]` = Vandermonde matrix.

---

## Eigenvalues

### EIGEN

Characteristic polynomial, eigenvalues, multiplicities, eigenspaces and diagonalisation.

**Enter:** `[A]` (square).

**Screens:**

1. The characteristic polynomial det(xI − A), as a list of coefficients, highest power first.
2. The eigenvalues and their algebraic multiplicities.
3. For each eigenvalue: its algebraic and geometric multiplicity and a basis for its eigenspace.
4. Either `DIAGONALIZABLE` with P and D (A = PDP⁻¹), or `NOT DIAGONALIZABLE`.
5. If diagonalisable, P⁻¹ in fractions (what you need for Aᵐ = PDᵐP⁻¹).

**Example:** `[A]` = `[[4,1][2,3]]`:

```text
CHAR POLY {1,-7,10}                  x² - 7x + 10
EIGENVALUES {5,2}   ALGEBRAIC MULT. {1,1}
EIGENVALUE 5   ALG MULT, GEO MULT {1,1}   EIGENSPACE BASIS [1 1]ᵀ
EIGENVALUE 2   ALG MULT, GEO MULT {1,1}   EIGENSPACE BASIS [-1/2 1]ᵀ
DIAGONALIZABLE   A = P*D*P^-1
P = [C] = [[1,-1/2][1,1]]      D = [D] = [[5,0][0,2]]
P^-1 = [E] = [[2/3,1/3][-2/3,2/3]]
```

For `[A]` = `[[1,1][0,1]]` you get eigenvalue 1 with ALG MULT, GEO MULT `{2,1}`, then
`NOT DIAGONALIZABLE`.

If the eigenvalues are complex, they are shown as complex numbers, followed by
`NOT DIAGONALIZABLE OVER THE REALS`. This leaves the calculator in `a+bi` mode; set it back with
`mode` → **REAL**.

**Saved:** `ʟCP` = characteristic polynomial, `ʟEV` = eigenvalues, `ʟAM` / `ʟGM` = algebraic /
geometric multiplicities, `[C]` = P, `[D]` = D, `[E]` = P⁻¹.

### ORTDIAG

Orthogonal diagonalisation of a symmetric matrix: A = PDPᵀ with P orthogonal. Each
eigenspace basis is made orthonormal with Gram–Schmidt.

**Enter:** `[A]` (symmetric; you get a message if it is not).

**Screens:** for each eigenvalue, an orthogonal basis of its eigenspace; then P, D, and P
before normalising.

**Example:** `[A]` = `[[2,1,1][1,2,1][1,1,2]]`:

```text
EIGENVALUE 4   ORTHOGONAL BASIS (1,1,1)
EIGENVALUE 1   ORTHOGONAL BASIS (-1,1,0), (-1/2,-1/2,1)
P = [C] (decimals)    D = [D] = diag(4,1,1)
P BEFORE NORMALISING = [E] = [[1,-1,-1/2][1,1,-1/2][1,0,1]]
```

P has square roots in it, so it shows as decimals. Write it exactly by dividing each
column of `[E]` by its length: here (1,1,1)/√3, (−1,1,0)/√2 and (−1/2,−1/2,1)/√(3/2).

**Saved:** `[C]` = P, `[D]` = D, `[E]` = P before normalising, `ʟEV` / `ʟAM` as in `EIGEN`.

### SVD

Singular value decomposition A = UΣVᵀ of any m×n matrix.

**Enter:** `[A]`.

**Screens:**

1. For each eigenvalue of AᵀA (largest first): its eigenvectors **before normalising**, and
   (if the eigenvalue is not 0) A times them, in fractions. Since Av = σu, each column of
   A·v points along a uᵢ: uᵢ = Avᵢ/‖Avᵢ‖.
2. If U needs extra columns (rank < m): the extra uᵢ before normalising, in fractions.
3. The singular values (and their squares, the eigenvalues of AᵀA), then U, Σ and V.

**Example:** `[A]` = `[[3,2,2][2,3,-2]]`:

```text
EIGENVALUE OF ATA 25   V VECTORS UNNORMALISED [1 1 0]ᵀ        A*V [5 5]ᵀ
EIGENVALUE OF ATA 9    V VECTORS UNNORMALISED [1/4 -1/4 1]ᵀ   A*V [9/4 -9/4]ᵀ
EIGENVALUE OF ATA 0    V VECTORS UNNORMALISED [-2 2 1]ᵀ
SINGULAR VALUES {5,3,0}    SQUARED {25,9,0}
U = [C] (2×2, decimals)    SIGMA = [D] = [[5,0,0][0,3,0]]    V = [E] (3×3, decimals)
```

V is made of orthonormal eigenvectors of AᵀA (largest eigenvalue first), and the columns
of U are uᵢ = Avᵢ/σᵢ, extended to an orthonormal basis of ℝᵐ when needed. Use the squared
singular values to write σᵢ exactly (here σ₁ = √25 = 5, σ₂ = √9 = 3), and the first screens
to write U and V exactly: here v₁ = (1,1,0)/√2, v₂ = (1,−1,4)/√18, v₃ = (−2,2,1)/3,
u₁ = (1,1)/√2, u₂ = (1,−1)/√2.

**Saved:** `[C]` = U, `[D]` = Σ, `[E]` = V, `ʟSV` = singular values.

### MARKOV

Markov chains: checks that the transition matrix is stochastic, finds the equilibrium
(steady-state) vector, and optionally the state after k steps.

**Enter:** `[A]` = transition matrix (each **column** sums to 1), and `[B]` = initial state x₀
(a column) if you want the k-step state.

**Prompts:** `K STEPS FOR A^K*B, 0=NO?`: type k to get Aᵏx₀, or `0` to skip (then `[B]`
is not needed).

**Example:** `[A]` = `[[0.7,0.2][0.3,0.8]]`, `[B]` = `[[1][0]]`, k = 2:

```text
A IS STOCHASTIC
EQUILIBRIUM VECTOR(S) (AX=X, SUM=1) = [C] = [2/5 3/5]ᵀ
K STEPS FOR A^K*B, 0=NO? 2
A^K*B = [D] = [11/20 9/20]ᵀ
```

If the matrix is not stochastic (a negative entry, or a column not summing to 1) you get a
warning, but the equilibrium is still computed.

**Saved:** `[C]` = equilibrium vector(s) (entries sum to 1), `[D]` = Aᵏx₀.

---

## Multiple choice

### TESTVEC

Checks several candidate vectors against one matrix in a single run, for questions like
"which of the following are eigenvectors of A / solutions of Ax = 0 / in the column space of A?".

**Enter:** `[A]` = the matrix, `[B]` = the candidate vectors as **columns** (one column per option).

**Screens:** for each vector in turn:

* If it has as many entries as `[A]` has columns: A·v, then `A*V=0: IN NULL(A)`, and for a
  square `[A]` either `EIGENVECTOR, LAMBDA=` with λ, or `NOT AN EIGENVECTOR`.
* If it has as many entries as `[A]` has rows: `IN COL(A)` or `NOT IN COL(A)`
  (rank of `[A|v]` compared with rank of `[A]`).

For a square `[A]` you get both.

**Example:** `[A]` = `[[1,1][1,1]]`, `[B]` = `[[1,1,1][1,-1,0]]` (options (1,1), (1,−1), (1,0)):

```text
VECTOR 1   A*V = [2 2]ᵀ   EIGENVECTOR, LAMBDA= 2      IN COL(A)
VECTOR 2   A*V = [0 0]ᵀ   A*V=0: IN NULL(A)
                          EIGENVECTOR, LAMBDA=0       NOT IN COL(A)
VECTOR 3   A*V = [1 1]ᵀ   NOT AN EIGENVECTOR          NOT IN COL(A)
```

The zero vector is never an eigenvector; it is flagged as `(ZERO VECTOR)`.

**Saved:** `B` = rank of `[A]` (for the column space test), `[E]` = A·v for the last vector.

---

## The MA1522 menu

`prgmMA1522` launches every program from one menu, so you do not have to scroll the
program list. The topics are spread over four pages:

| Page | Options |
|------|---------|
| 1/3 | REF / RREF · SOLVE AX=B · LU = PLU · DET/INVERSE/ADJ · CRAMER RULE · NEXT PAGE · QUIT |
| 2/3 | SUBSPACES OF A · COL CONSTRAINTS · SUBSPACE TOOLS · COORD/TRANSITION · LINEAR TRANSFORM · NEXT PAGE · BACK |
| 3/3 | ORTHOGONALITY · GRAM-SCHMIDT/QR · LEAST SQUARES · POLYNOMIAL FIT · EIGEN/DIAGONAL · MORE · BACK |
| MORE | ORTHOG DIAGONAL · SVD · MARKOV CHAIN · TEST VECTORS · BACK · QUIT |

Pick an option with the arrow keys and `enter`, or press its number. When the program
finishes you are back on the home screen; run `prgmMA1522` again for the next one
(`2nd` `enter` recalls the last command).

## Things to know

* **Results are decimals, shown as fractions where possible.** The calculator has no
  computer algebra system, so it works with 14-digit decimals and displays them with
  `►Frac`. Irrational values such as 1/√2 stay as decimals. For exact answers, `GRAMSCH`,
  `ORTDIAG` and `SVD` also show the vectors before normalising (and `GRAMSCH` the exact
  numerators of R), so you can write each unit vector as v/√‖v‖².

* **Eigenvalues** are found as roots of the characteristic polynomial. Repeated eigenvalues
  are detected and snapped to simple fractions (denominator ≤ 12). For large matrices
  (5×5, 6×6) `EIGEN`, `ORTDIAG` and `SVD` can take several seconds.

* **Variables are overwritten.** The programs use the real variables A–Z, matrices
  `[C]`–`[J]` and some named lists (`ʟPIV`, `ʟZD`, …). `[A]`, `[B]`, `L₁` and `L₂` are only read.
  `EIGEN` leaves the calculator in `a+bi` mode if it finds complex eigenvalues.

* Row operations in the step display are numbered `R1`–`R9`, so steps are meant for matrices
  with at most 9 rows.

* With steps on, the first non-zero entry is used as the pivot (like working by hand).
  With steps off, the largest entry is used, which is more accurate.

* **Check the answer.** Answers are numerical, so a quick check is worth it in an exam:
  multiply back (e.g. `[A]*[C]` for a solution, `[C]*[D]*[C]⁻¹` for a diagonalisation).

## Troubleshooting

| What you see | Cause and fix |
|--------------|---------------|
| Dotted busy indicator in the top-right, nothing happens | The program is waiting for you. Press `enter`. |
| `ERR:UNDEFINED` | `[A]` (or `[B]`) is empty: enter it in the matrix editor first. Otherwise a helper program is missing: send all the files in `bin/` again, including the `Z*` ones. |
| `ERR:DIM MISMATCH` or `ERR:INVALID DIM` | `[A]` and `[B]` do not fit together, e.g. `[B]` has a different number of rows from `[A]`. |
| `ERR:ARCHIVED` | The programs went to Archive. `2nd` `+` (MEM) → **2: Mem Management** → **Prgm**, select each one and press `enter` to move it to RAM. |
| An error (`ERR:DOMAIN` / `ERR:INVALID DIM`) only with steps on | The matrix has more than 9 rows, and steps number rows `R1`–`R9`. Run with steps off. |
| `A MUST BE SQUARE` | `INVERSE`, `CRAMER`, `EIGEN`, `ORTDIAG` and `MARKOV` need a square `[A]`. |
| Program takes long | `EIGEN`, `ORTDIAG` and `SVD` on 5×5 or bigger can take several seconds. The busy indicator moves while it works. Press `on` to stop. |
| Complex numbers everywhere afterwards | `EIGEN` found complex eigenvalues and switched to `a+bi` mode. `mode` → **REAL**. |

## For developers

```text
src/      TI-BASIC source, one program per .tib file (Unicode token names, // comments)
bin/      built .8xp files, ready to send to the calculator
tools/    build.py (source -> .8xp), sim.py (TI-BASIC interpreter for testing)
tests/    runs every program in the interpreter and checks the results
```

```bash
python3 tools/build.py                   # rebuild bin/*.8xp (downloads the TI-Toolkit token sheet once)
python3 -m unittest discover -s tests -v # run the programs in the simulator
```

`tools/sim.py` runs a subset of TI-BASIC with IEEE doubles. It is close to the calculator but
it is not an emulator, so test changes on a real calculator (or CEmu) as well.

## Credits

Based on the `ma1522` Python library by [@YeeShin504](https://github.com/YeeShin504), with
contributions from [@DenseLance](https://github.com/DenseLance). Token table from
[TI-Toolkit/tokens](https://github.com/TI-Toolkit/tokens).
