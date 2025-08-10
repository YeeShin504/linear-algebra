
# MA1522 Linear Algebra for Computing Tutorial 6

## Problem 1

**(a)** Let $\mathbf{u}_1 = \begin{pmatrix} 1 \\ 2 \\ -1 \end{pmatrix}$, $\mathbf{u}_2 = \begin{pmatrix} 0 \\ 2 \\ 1 \end{pmatrix}$, $\mathbf{u}_3 = \begin{pmatrix} 0 \\ -1 \\ 3 \end{pmatrix}$. Show that $S = \{\mathbf{u}_1, \mathbf{u}_2, \mathbf{u}_3\}$ forms a basis for $\mathbb{R}^3$.

**(b)** Suppose $\mathbf{w} = \begin{pmatrix} 1 \\ 1 \\ 1 \end{pmatrix}$. Find the coordinate vector of $\mathbf{w}$ relative to $S$.

**(c)** Let $T = \{\mathbf{v}_1, \mathbf{v}_2, \mathbf{v}_3\}$ be another basis for $\mathbb{R}^3$ where $\mathbf{v}_1 = \begin{pmatrix} 1 \\ 5 \\ 4 \end{pmatrix}$, $\mathbf{v}_2 = \begin{pmatrix} -1 \\ 3 \\ 7 \end{pmatrix}$, $\mathbf{v}_3 = \begin{pmatrix} 2 \\ 2 \\ 4 \end{pmatrix}$. Find the transition matrix from $T$ to $S$.

**(d)** Find the transition matrix from $S$ to $T$.

**(e)** Use the vector $\mathbf{w}$ in Part (b). Find the coordinate vector of $\mathbf{w}$ relative to $T$.

## Problem 2

Let $V$ be a subspace of $\mathbb{R}^n$ and $S = \{\mathbf{u}_1, \mathbf{u}_2, \mathbf{u}_3\}$ be a basis for a subspace $V$. Define $T = \{\mathbf{v}_1, \mathbf{v}_2, \mathbf{v}_3\}$ where $$\mathbf{v}_1 = \mathbf{u}_1 + \mathbf{u}_2 + \mathbf{u}_3, \quad \mathbf{v}_2 = \mathbf{u}_2 + \mathbf{u}_3 \quad \text{and} \quad \mathbf{v}_3 = \mathbf{u}_2 - \mathbf{u}_3$$

**(a)** Show that $T$ is a basis for $V$.

**(b)** Find the transition matrix from $S$ to $T$.

## Problem 3

**(a)** Let $\mathbf{A} = \begin{pmatrix} 1 & -1 & 1 \\ 1 & 1 & -1 \\ -1 & -1 & 1 \end{pmatrix}$ and $\mathbf{b} = \begin{pmatrix} 2 \\ 1 \\ 0 \end{pmatrix}$. Is $\mathbf{b}$ in the column space of $\mathbf{A}$? If it is, express it as a linear combination of the columns of $\mathbf{A}$.

**(b)** Let $\mathbf{A} = \begin{pmatrix} 1 & 9 & 1 \\ -1 & 3 & 1 \\ 1 & 1 & 1 \end{pmatrix}$ and $\mathbf{b} = \begin{pmatrix}5 & 1 & -1\end{pmatrix}$. Is $\mathbf{b}$ in the row space of $\mathbf{A}$? If it is, express it as a linear combination of the rows of $\mathbf{A}$.

**(c)** Let $\mathbf{A} = \begin{pmatrix} 1 & 2 & 0 & 1 \\ 0 & 1 & 2 & 1 \\ 1 & 2 & 1 & 3 \\ 0 & 1 & 2 & 2 \end{pmatrix}$. Is the row space and column space of $\mathbf{A}$ the whole $\mathbb{R}^4$?

## Problem 4

For each of the following matrices $\mathbf{A}$,

(i) Find a basis for the row space of $\mathbf{A}$. (ii) Find a basis for the column space of $\mathbf{A}$. (iii) Find a basis for the nullspace of $\mathbf{A}$. (iv) Hence determine $\text{rank}(\mathbf{A})$, $\text{nullity}(\mathbf{A})$ and verify the dimension theorem for matrices. (v) Is $\mathbf{A}$ full rank?

**(a)** $\mathbf{A} = \begin{pmatrix} 1 & 2 & 5 & 3 \\ 1 & -4 & -1 & -9 \\ -1 & 0 & -3 & 1 \\ 2 & 1 & 7 & 0 \\ 0 & 1 & 1 & 2 \end{pmatrix}$

**(b)** $\mathbf{A} = \begin{pmatrix} 1 & 3 & 7 \\ 2 & 1 & 8 \\ 3 & -5 & -1 \\ 2 & -2 & 2 \\ 1 & 1 & 5 \end{pmatrix}$

## Problem 5

Let $W$ be a subspace of $\mathbb{R}^5$ spanned by the following vectors $$\mathbf{u}_1 = \begin{pmatrix} 1 \\ -2 \\ 0 \\ 0 \\ 3 \end{pmatrix}, \quad \mathbf{u}_2 = \begin{pmatrix} 2 \\ -5 \\ -3 \\ -2 \\ 6 \end{pmatrix}, \quad \mathbf{u}_3 = \begin{pmatrix} 0 \\ 5 \\ 15 \\ 10 \\ 0 \end{pmatrix}, \quad \mathbf{u}_4 = \begin{pmatrix} 2 \\ 1 \\ 15 \\ 8 \\ 6 \end{pmatrix}$$

**(a)** Find a basis for $W$.

**(b)** What is $\dim(W)$?

**(c)** Extend the basis $W$ found in (a) to a basis for $\mathbb{R}^5$.

## Problem 6

Let $S = \left\{ \begin{pmatrix} 1 \\ 0 \\ 1 \\ 3 \end{pmatrix}, \begin{pmatrix} 2 \\ -1 \\ 0 \\ 1 \end{pmatrix}, \begin{pmatrix} -1 \\ 3 \\ 5 \\ 12 \end{pmatrix}, \begin{pmatrix} 0 \\ 1 \\ 2 \\ 5 \end{pmatrix}, \begin{pmatrix} 3 \\ -1 \\ 1 \\ 4 \end{pmatrix} \right\}$ and $V = \text{span}(S)$. Find a subset $S' \subseteq S$ such that $S'$ forms a basis for $V$.

# MA1522 Linear Algebra for Computing Tutorial 7

## Problem 1

**(a)** Let $a_1 x_1 + a_2 x_2 + \cdots + a_n x_n = b$ be a linear equation. Express this linear system as $\mathbf{a} \cdot \mathbf{x} = b$ for some (column) vectors $\mathbf{a}$ and $\mathbf{x}$.

**(b)** Find the solution set of the linear system $$\begin{cases}
x_1 &+& 3x_2 &-& 2x_3 &&  &=& 0 \\
2x_1 &+& 6x_2 &-& 5x_3 &-& 2x_4 &=& 0  \\
&& && 5x_3 &+& 10x_4 &=& 0  
\end{cases}$$


**(c)** Find a nonzero vector $\mathbf{v} \in \mathbb{R}^4$ such that $\mathbf{a}_1 \cdot \mathbf{v} = 0$, $\mathbf{a}_2 \cdot \mathbf{v} = 0$, and $\mathbf{a}_3 \cdot \mathbf{v} = 0$, where $$\mathbf{a}_1 = \begin{pmatrix} 1 \\ 3 \\ -2 \\ 0 \end{pmatrix}, \quad \mathbf{a}_2 = \begin{pmatrix} 2 \\ 6 \\ -5 \\ -2 \end{pmatrix}, \quad \mathbf{a}_3 = \begin{pmatrix} 0 \\ 0 \\ 5 \\ 10 \end{pmatrix}$$

This exercise demonstrates the fact that if $\mathbf{A}$ is a $m \times n$ matrix, then the solution set of the homogeneous linear system $\mathbf{A}\mathbf{x} = \mathbf{0}$ consists of all the vectors in $\mathbb{R}^n$ that are orthogonal to every row vector of $\mathbf{A}$.

## Problem 2

Let $\{\mathbf{v}_1, \mathbf{v}_2, \mathbf{v}_3\}$ be an orthonormal set. Suppose $$\mathbf{x} = \mathbf{v}_1 - 2\mathbf{v}_2 - 2\mathbf{v}_3 \quad \text{and} \quad \mathbf{y} = 2\mathbf{v}_1 - 3\mathbf{v}_2 + \mathbf{v}_3$$

Determine the value for each of the following:

**(a)** $\mathbf{x} \cdot \mathbf{y}$.

**(b)** $||\mathbf{x}||$ and $||\mathbf{y}||$.

**(c)** The angle $\theta$ between $\mathbf{x}$ and $\mathbf{y}$.

## Problem 3

Let $\mathbf{v}_1 = \begin{pmatrix} 1 \\ 2  \\ \ -1 \end{pmatrix}$, $\mathbf{v}_2 = \begin{pmatrix} 1  \\ 0 \\ 1 \end{pmatrix}$, and $\mathbf{V} = \begin{pmatrix} \mathbf{v}_1 & \mathbf{v}_2 \end{pmatrix}$.

**(a)** Compute $\mathbf{v}_1 \cdot \mathbf{v}_1$, $\mathbf{v}_1 \cdot \mathbf{v}_2$, $\mathbf{v}_2 \cdot \mathbf{v}_1$ and $\mathbf{v}_2 \cdot \mathbf{v}_2$.

**(b)** Compute $\mathbf{V}^T \mathbf{V}$. What do the entries of $\mathbf{V}^T \mathbf{V}$ represent?

## Problem 4

Let $W$ be a subspace of $\mathbb{R}^n$. The orthogonal complement of $W$, denoted as $W^{\perp}$, is defined to be $$W^{\perp} := \{ \mathbf{v} \in \mathbb{R}^n : \mathbf{v} \cdot \mathbf{w} = 0 \text{ for all } \mathbf{w} \in W \}$$

Let $\mathbf{w}_1 = \begin{pmatrix} 1 \\ 1 \\ 1 \\ 1 \\ 1 \end{pmatrix}$, $\mathbf{w}_2 = \begin{pmatrix} 1 \\ 2 \\ -1 \\ -2 \\ 0 \end{pmatrix}$, and $\mathbf{w}_3 = \begin{pmatrix} 1 \\ -1 \\ 1 \\ -1 \\ 0 \end{pmatrix}$, and $W = \text{span}\{\mathbf{w}_1, \mathbf{w}_2, \mathbf{w}_3\}$.

**(a)** Show that $S = \{\mathbf{w}_1, \mathbf{w}_2, \mathbf{w}_3\}$ is linearly independent.

**(b)** Show that $S$ is orthogonal.

**(c)** Show that $W^{\perp}$ is a subspace of $\mathbb{R}^5$ by showing that it is a span of a set. What is the dimension? (Hint: See Question 1.)

**(d)** Obtain an orthonormal set $T$ by normalizing $\mathbf{w}_1$, $\mathbf{w}_2$, $\mathbf{w}_3$.

**(e)** Let $\mathbf{v} = \begin{pmatrix} 2 \\ 0 \\ 1 \\ 1 \\ -1 \end{pmatrix}$. Find the projection of $\mathbf{v}$ onto $W$.

**(f)** Let $\mathbf{v}_W$ be the projection of $\mathbf{v}$ onto $W$. Show that $\mathbf{v} - \mathbf{v}_W$ is in $W^{\perp}$.

This exercise demonstrates the fact that every vector $\mathbf{v}$ in $\mathbb{R}^5$ can be written as $\mathbf{v} = \mathbf{v}_W + \mathbf{v}_{W^{\perp}}$, for some $\mathbf{v}_W$ in $W$ and $\mathbf{v}_{W^{\perp}}$ in $W^{\perp}$. In other words, $W + W^{\perp} = \mathbb{R}^5$.

## Problem 5

Let $S = \{\mathbf{u}_1, \mathbf{u}_2, \mathbf{u}_3, \mathbf{u}_4\}$ where $$\mathbf{u}_1 = \begin{pmatrix} 1 \\ 2 \\ 2 \\ -1 \end{pmatrix}, \quad \mathbf{u}_2 = \begin{pmatrix} 1 \\ 1 \\ -1 \\ 1 \end{pmatrix}, \quad \mathbf{u}_3 = \begin{pmatrix} -1 \\ 1 \\ -1 \\ -1 \end{pmatrix}, \quad \mathbf{u}_4 = \begin{pmatrix} -2 \\ 1 \\ 1 \\ 2 \end{pmatrix}$$

**(a)** Check that $S$ is an orthogonal basis for $\mathbb{R}^4$.

**(b)** Is it possible to find a nonzero vector $\mathbf{w}$ in $\mathbb{R}^4$ such that $S \cup {\mathbf{w}}$ is an orthogonal set?

**(c)** Obtain an orthonormal set $T$ by normalizing $\mathbf{u}_1$, $\mathbf{u}_2$, $\mathbf{u}_3$, $\mathbf{u}_4$.

**(d)** Let $\mathbf{v} = \begin{pmatrix} 0 \\ 1 \\ 2 \\ 3 \end{pmatrix}$. Find $[\mathbf{v}]_S$ and $[\mathbf{v}]_T$.

**(e)** Suppose $\mathbf{w}$ is a vector in $\mathbb{R}^4$ such that $[\mathbf{w}]_S = \begin{pmatrix} 1 \\ 2 \\ 1 \\ 1 \end{pmatrix}$. Find $[\mathbf{w}]_T$.

# MA1522 Linear Algebra for Computing Tutorial 8

## Problem 1

Apply Gram-Schmidt Process to convert

**(a)** $\left\{ \begin{pmatrix} 1 \\ 1 \\ 1 \\ 1 \end{pmatrix}, \begin{pmatrix} 1 \\ -1 \\ 1 \\ 0 \end{pmatrix}, \begin{pmatrix} 1 \\ 1 \\ -1 \\ -1 \end{pmatrix}, \begin{pmatrix} 1 \\ 2 \\ 0 \\ 1 \end{pmatrix} \right\}$ into an orthonormal basis for $\mathbb{R}^4$.

**(b)** $\left\{ \begin{pmatrix} 1 \\ 2 \\ 2 \\ 1 \end{pmatrix}, \begin{pmatrix} 1 \\ 2 \\ 1 \\ 0 \end{pmatrix}, \begin{pmatrix} 1 \\ 0 \\ 1 \\ 0 \end{pmatrix}, \begin{pmatrix} 1 \\ 0 \\ 2 \\ 1 \end{pmatrix} \right\}$ into an orthonormal set. Is the set obtained an orthonormal basis? Why?

## Problem 2

Let $\mathbf{A} = \begin{pmatrix} 0 & 1 & 1 & 0 \\ 1 & -1 & 1 & -1 \\ 1 & 0 & 1 & 0 \\ 1 & 1 & 1 & 1 \end{pmatrix}$ and $\mathbf{b} = \begin{pmatrix} 6 \\ 3 \\ -1 \\ 1 \end{pmatrix}$.

**(a)** Is the linear system $\mathbf{A}\mathbf{x} = \mathbf{b}$ inconsistent?

**(b)** Find a least squares solution to the system. Is the solution unique?

**(c)** Use your answer in (b), compute the projection of $\mathbf{b}$ onto the column space of $\mathbf{A}$. Is the solution unique?

## Problem 3 (Application)

A line $p(x) = a_1 x + a_0$ is said to be the least squares approximating line for a given set of data points $(x_1, y_1), (x_2, y_2), \ldots, (x_m, y_m)$ if the sum $$S = [y_1 - p(x_1)]^2 + [y_2 - p(x_2)]^2 + \cdots + [y_m - p(x_m)]^2$$ is minimized. Writing $$\mathbf{x} = \begin{pmatrix} x_1 \\ x_2 \\ \vdots \\ x_m \end{pmatrix}, \quad \mathbf{y} = \begin{pmatrix} y_1 \\ y_2 \\ \vdots \\ y_m \end{pmatrix}, \quad \text{and} \quad p(\mathbf{x}) = \begin{pmatrix} p(x_1) \\ p(x_2) \\ \vdots \\ p(x_m) \end{pmatrix} = \begin{pmatrix} a_1 x_1 + a_0 \\ a_1 x_2 + a_0 \\ \vdots \\ a_1 x_m + a_0 \end{pmatrix}$$

the problem is now rephrased as finding $a_0, a_1$ such that $$S = ||\mathbf{y} - p(\mathbf{x})||^2$$ is minimized. Observe that if we let $$\mathbf{N} = \begin{pmatrix} 1 & x_1 \\ 1 & x_2 \\ \vdots & \vdots \\ 1 & x_m \end{pmatrix} \quad \text{and} \quad \mathbf{a} = \begin{pmatrix} a_0 \\ a_1 \end{pmatrix}$$

then $\mathbf{N}\mathbf{a} = p(\mathbf{x})$. And so our aim is to find $\mathbf{a}$ that minimizes $||\mathbf{y} - \mathbf{N}\mathbf{a}||^2$.

It is known the equation representing the dependency of the resistance of a cylindrically shaped conductor (a wire) at $20°C$ is given by $$R = \rho \frac{L}{A}$$ where $R$ is the resistance measured in Ohms $\Omega$, $L$ is the length of the material in meters $m$, $A$ is the cross-sectional area of the material in meter squared $m^2$, and $\rho$ is the resistivity of the material in Ohm meters $\Omega m$.

A student wants to measure the resistivity of a certain material. Keeping the cross-sectional area constant at $0.002 m^2$, he connected the power sources along the material at various lengths and measured the resistance and obtained the following data.

|$L$|$0.01$|$0.012$|$0.015$|$0.02$|
|---|---|---|---|---|
|$R$|$2.75 \times 10^{-4}$|$3.31 \times 10^{-4}$|$3.92 \times 10^{-4}$|$4.95 \times 10^{-4}$|

It is known that the Ohm meter might not be calibrated. Taking that into account, the student wants to find a linear graph $R = \frac{\rho}{0.002} L + R_0$ from the data obtained to compute the resistivity of the material.

**(a)** Relabeling, we let $R = y$, $\frac{\rho}{0.002} = a_1$ and $R_0 = a_0$. Is it possible to find a graph $y = a_1 x + a_0$ satisfying the points?

**(b)** Find the least square approximating line for the data points and hence find the resistivity of the material. Would this material make a good wire?

## Problem 4 (Application)

Suppose the equation governing the relation between data pairs is not known. We may want to then find a polynomial $$p(x) = a_0 + a_1 x + a_2 x^2 + \cdots + a_n x^n$$ of degree $n$, $n \leq m - 1$, that best approximates the data pairs $(x_1, y_1), (x_2, y_2), \ldots, (x_m, y_m)$. A least square approximating polynomial of degree $n$ is such that $$||\mathbf{y} - p(\mathbf{x})||^2$$ is minimized. If we write $$\mathbf{x} = \begin{pmatrix} x_1 \\ x_2 \\ \vdots \\ x_m \end{pmatrix}, \quad \mathbf{y} = \begin{pmatrix} y_1 \\ y_2 \\ \vdots \\ y_m \end{pmatrix}, \quad \mathbf{N} = \begin{pmatrix} 1 & x_1 & x_1^2 & \cdots & x_1^n \\ 1 & x_2 & x_2^2 & \cdots & x_2^n \\ \vdots & \vdots & \vdots & \ddots & \vdots \\ 1 & x_m & x_m^2 & \cdots & x_m^n \end{pmatrix} \quad \text{and} \quad \mathbf{a} = \begin{pmatrix} a_0 \\ a_1 \\ \vdots \\ a_n \end{pmatrix}$$

then $p(\mathbf{x}) = \mathbf{N}\mathbf{a}$, and the task is to find $\mathbf{a}$ such that $||\mathbf{y} - \mathbf{N}\mathbf{a}||^2$ is minimized. Observe that $\mathbf{N}$ is a matrix minor of the Vandermonde matrix. If at least $n + 1$ of the $x$-values $x_1, x_2, \ldots, x_m$ are distinct, the columns of $\mathbf{N}$ are linearly independent, and thus $\mathbf{a}$ is uniquely determined by $$\mathbf{a} = (\mathbf{N}^T \mathbf{N})^{-1} \mathbf{N}^T \mathbf{y}$$

We shall now find a quartic polynomial $$p(x) = a_0 + a_1 x + a_2 x^2 + a_3 x^3 + a_4 x^4$$ that is a least square approximating polynomial for the following data points

|$x$|$4$|$4.5$|$5$|$5.5$|$6$|$6.5$|$7$|$8$|$8.5$|
|---|---|---|---|---|---|---|---|---|---|
|$y$|$0.8651$|$0.4828$|$2.590$|$-4.389$|$-7.858$|$3.103$|$7.456$|$0.0965$|$4.326$|

## Problem 5

Let $\mathbf{A} = \begin{pmatrix} 1 & 1 & 0 \\ 1 & 1 & 0 \\ 1 & 1 & 1 \\ 0 & 1 & 1 \end{pmatrix}$.

**(a)** Find a $QR$ factorization of $\mathbf{A}$.

**(b)** Use your answer in (a) to find the least square solution to $\mathbf{A}\mathbf{x} = \mathbf{b}$, where $\mathbf{b} = \begin{pmatrix} 1 \\ 1 \\ 0 \\ 0 \end{pmatrix}$.

# MA1522 Linear Algebra for Computing Tutorial 9

## Problem 1

A father wishes to distribute an amount of money among his three sons Jack, Jim, and John. He wishes to distribute such that the following conditions are all satisfied.

(i) The amount Jack receives plus twice the amount Jim receives is $300. (ii) The amount Jim receives plus the amount John receives is $300. (iii) Jack receives $300 more than twice of what John receives.

**(a)** Is it possible for the following conditions to all be satisfied?

**(b)** If it is not possible, find a least square solution. (Make sure that your least square solution is feasible. For example, one cannot give a negative amount of money to anybody.)

## Problem 2

**(a)** Suppose $\mathbf{A}$ is a $m \times n$ matrix where $m > n$. Let $\mathbf{A} = \mathbf{Q}\mathbf{R}$ be a $QR$ factorization of $\mathbf{A}$. Explain how you might use this to write $$\mathbf{A} = \mathbf{Q}'\mathbf{R}'$$ where $\mathbf{Q}'$ is an $m \times m$ orthogonal matrix, and $\mathbf{R}'$ a $m \times n$ matrix with $m - n$ zero rows at the bottom. This is known as the full $QR$ factorization of $\mathbf{A}$.

**(b)** In MATLAB, enter the following.

```matlab
>> A=sym([1 1 0;1 1 0;1 1 1;0 1 1])
>> [Q R]=qr(A)
```

What is $\mathbf{Q}$ and $\mathbf{R}$? Compare this with the answer in tutorial 8 question 5(a).

**(c)** Explain how you might use the command `qr` in MATLAB to find a $QR$ factorization of a $m \times n$ matrix $\mathbf{A}$?

## Problem 3 (Cayley-Hamilton theorem)

Consider $$p(X) = X^3 - 4X^2 - X + 4I$$

**(a)** Compute $p(X)$ for $X = \begin{pmatrix} 1 & 1 & 2 \\ 1 & 2 & 1 \\ 2 & 1 & 1 \end{pmatrix}$.

**(b)** Find the characteristic polynomial of $X$.

**(c)** Show that $X$ is invertible. Express the inverse of $X$ as a function of $X$.

This question demonstrates the Cayley-Hamilton theorem, which states that if $p(x)$ is the characteristic polynomial of $X$, then $p(X) = 0$. This also shows that if $0$ is not an eigenvalue of $X$, then the constant term of the characteristic polynomial $p(x)$ is nonzero, and we can use that to compute the inverse of $X$.

## Problem 4

For each of the following matrices $\mathbf{A}$, determine if $\mathbf{A}$ is diagonalizable. If $\mathbf{A}$ is diagonalizable, find an invertible $\mathbf{P}$ that diagonalizes $\mathbf{A}$ and determine $\mathbf{P}^{-1}\mathbf{A}\mathbf{P}$.

**(a)** $\mathbf{A} = \begin{pmatrix} 1 & -3 & 3 \\ 3 & -5 & 3 \\ 6 & -6 & 4 \end{pmatrix}$

**(b)** $\mathbf{A} = \begin{pmatrix} 9 & 8 & 6 & 3 \\ 0 & -1 & 3 & -4 \\ 0 & 0 & 2 & 0 \\ 0 & 0 & 0 & 3 \end{pmatrix}$

**(c)** $\mathbf{A} = \begin{pmatrix} 1 & 0 & 0 \\ 1 & 1 & 0 \\ 0 & 1 & 1 \end{pmatrix}$

**(d)** $\mathbf{A} = \begin{pmatrix} 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \\ 1 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 \end{pmatrix}$

**(e)** $\mathbf{A} = \begin{pmatrix} -1 & 1 & 1 \\ 1 & 1 & -1 \\ -4 & 2 & 3 \end{pmatrix}$

## Problem 5

**(a)** Show that $\lambda$ is an eigenvalue of $\mathbf{A}$ if and only if it is an eigenvalue of $\mathbf{A}^T$.

**(b)** Suppose $\mathbf{A}$ is diagonalizable. Is $\mathbf{A}^T$ diagonalizable? Justify your answer.

**(c)** Suppose $\mathbf{v}$ is an eigenvector of $\mathbf{A}$ associated to eigenvalue $\lambda$. Show that $\mathbf{v}$ is an eigenvector of $\mathbf{A}^k$ associated to eigenvalue $\lambda^k$ for any positive integer $k$.

**(d)** If $\mathbf{A}$ is invertible, show that $\mathbf{v}$ is an eigenvector of $\mathbf{A}^k$ associated to eigenvalue $\lambda^k$ for any negative integer $k$.

**(e)** A square matrix is said to be nilpotent if there is a positive integer $k$ such that $\mathbf{A}^k = \mathbf{0}$. Show that if $\mathbf{A}$ is nilpotent, then $0$ is the only eigenvalue.

**(f)** Let $\mathbf{A}$ be a $n \times n$ matrix with one eigenvalue $\lambda$ with algebraic multiplicity $n$. Show that $\mathbf{A}$ is diagonalizable if and only if $\mathbf{A}$ is a scalar matrix, $\mathbf{A} = \lambda \mathbf{I}$.

**(g)** Show that the only diagonalizable nilpotent matrix is the zero matrix.

# MA1522 Linear Algebra for Computing Tutorial 10

## Problem 1

A population of ants is put into a maze with 3 compartments labeled a, b, and c. If the ant is in compartment a, an hour later, there is a 20% chance it will go to compartment b, and a 40% chance it will go to compartment c. If it is in compartment b, an hour later, there is a 10% chance it will go to compartment a, and a 30% chance it will go to compartment c. If it is in compartment c, an hour later, there is a 50% chance it will go to compartment a, and a 20% chance it will go to compartment b. Suppose 100 ants have been placed in compartment a.

**(a)** Find the transition probability matrix $\mathbf{A}$. Show that it is a stochastic matrix.

**(b)** By diagonalizing $\mathbf{A}$, find the number of ants in each compartment after 3 hours.

**(c)** (MATLAB) We can use MATLAB to diagonalize the matrix $\mathbf{A}$. Type

```matlab
>> [P D]=eig(sym(A))
```

The matrix $\mathbf{P}$ will be an invertible matrix, and $\mathbf{D}$ will be a diagonal matrix. Compare the answer with what you have obtained in (b).

**(d)** In the long run (assuming no ants died), where will the majority of the ants be?

**(e)** Suppose initially the numbers of ants in compartments a, b and c are $\alpha$, $\beta$, and $\gamma$ respectively. What is the population distribution in the long run (assuming no ants died)?

## Problem 2

By diagonalizing $\mathbf{A} = \begin{pmatrix} 1 & 0 & 3 \\ 0 & 4 & 0 \\ 0 & 0 & 4 \end{pmatrix}$, find a matrix $\mathbf{B}$ such that $\mathbf{B}^2 = \mathbf{A}$.

## Problem 3

For each of the following symmetric matrices $\mathbf{A}$, find an orthogonal matrix $\mathbf{P}$ that orthogonally diagonalizes $\mathbf{A}$.

**(a)** $\mathbf{A} = \begin{pmatrix} 3 & 1 \\ 1 & 3 \end{pmatrix}$

**(b)** $\mathbf{A} = \begin{pmatrix} 2 & 2 & -2 \\ 2 & -1 & 4 \\ -2 & 4 & -1 \end{pmatrix}$

## Problem 4 (MATLAB)

Let $\mathbf{A} = \begin{pmatrix} 1 & -2 & 0 & 0 \\ -2 & 1 & 0 & 0 \\ 0 & 0 & 1 & -2 \\ 0 & 0 & -2 & 1 \end{pmatrix}$.

**(a)** Find an orthogonal matrix $\mathbf{P}$ that orthogonally diagonalizes $\mathbf{A}$, and compute $\mathbf{P}^T \mathbf{A} \mathbf{P}$.

**(b)** We will use MATLAB to orthogonally diagonalize $\mathbf{A}$. Type

```matlab
>> A=[1 -2 0 0;-2 1 0 0;0 0 1 -2;0 0 -2 1];
>> [P D]=eig(A);
>> sym(P), sym(D)
```

Compare the result with your answer in (a).

## Problem 5

Find the SVD of the following matrices $\mathbf{A}$.

**(a)** $\mathbf{A} = \begin{pmatrix} 3 & 2 \\ 2 & 3 \\ 2 & -2 \end{pmatrix}$

**(b)** $\mathbf{A} = \begin{pmatrix} 3 & 2 & 2 \\ 2 & 3 & -2 \end{pmatrix}$

**(c)** $\mathbf{A} = \begin{pmatrix} 1 & 0 & 1 \\ 0 & 1 & 1 \\ 1 & 1 & 2 \end{pmatrix}$

## Problem 6 (MATLAB)

Let $\mathbf{A} = \begin{pmatrix} -18 & 13 & -4 & 4 \\ 2 & 19 & -4 & 12 \\ -14 & 11 & -12 & 8 \\ -2 & 21 & 4 & 8 \end{pmatrix}$.

**(a)** Find a SVD of $\mathbf{A}$.

**(b)** In MATLAB, type

```matlab
>> [U S V]=svd(A)
```

Compare the result with your answer in (a).
# MA1522 Linear Algebra for Computing Tutorial 11

## Problem 1

For each of the following: (i) Determine whether the following are linear transformations. (ii) Write down the standard matrix for each of the linear transformations. (iii) Find a basis for the range for each of the linear transformations. (iv) Find a basis for the kernel for each of the linear transformations.

**(a)** $T_1 : \mathbb{R}^2 \to \mathbb{R}^2$ such that $T_1\begin{pmatrix} x \\ y \end{pmatrix} = \begin{pmatrix} x + y \\ y - x \end{pmatrix}$ for $\begin{pmatrix} x \\ y \end{pmatrix} \in \mathbb{R}^2$.

**(b)** $T_2 : \mathbb{R}^2 \to \mathbb{R}^2$ such that $T_2\begin{pmatrix} x \\ y \end{pmatrix} = \begin{pmatrix} 2x \\ 0 \end{pmatrix}$ for $\begin{pmatrix} x \\ y \end{pmatrix} \in \mathbb{R}^2$.

**(c)** $T_3 : \mathbb{R}^2 \to \mathbb{R}^3$ such that $T_3\begin{pmatrix} x \\ y \end{pmatrix} = \begin{pmatrix} x + y \\ 0 \\ 0 \end{pmatrix}$ for $\begin{pmatrix} x \\ y \end{pmatrix} \in \mathbb{R}^2$.

**(d)** $T_4 : \mathbb{R}^3 \to \mathbb{R}^3$ such that $T_4\begin{pmatrix} x \\ y \\ z \end{pmatrix} = \begin{pmatrix} 1 \\ y - x \\ y - z \end{pmatrix}$ for $\begin{pmatrix} x \ y \ z \end{pmatrix} \in \mathbb{R}^3$.

**(e)** $T_5 : \mathbb{R}^5 \to \mathbb{R}$ such that $T_5\begin{pmatrix} x_1 \\ \vdots \\ x_5 \end{pmatrix} = x_3 + 2x_4 - x_5$ for $\begin{pmatrix} x_1 \\ \vdots \\ x_5 \end{pmatrix} \in \mathbb{R}^5$.

**(f)** $T_6 : \mathbb{R}^n \to \mathbb{R}$ such that $T_6(\mathbf{x}) = \mathbf{x} \cdot \mathbf{x}$ for $\mathbf{x} \in \mathbb{R}^n$.

## Problem 2

Let $F : \mathbb{R}^3 \to \mathbb{R}^3$ and $G : \mathbb{R}^3 \to \mathbb{R}^3$ be linear transformations such that

$$F\left( \begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix} \right) = \begin{pmatrix} x_1 - 2x_2 \\ x_1 + x_2 - 3x_3 \\ 5x_2 - x_3 \end{pmatrix} \quad \text{and} \quad G\left( \begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix} \right) = \begin{pmatrix} x_3 - x_1 \\ x_2 + 5x_1 \\ x_1 + x_2 + x_3 \end{pmatrix}$$

and let $\mathbf{A}_F$ and $\mathbf{B}_G$ be the standard matrix of $F$ and $G$, respectively.

**(a)** Find $\mathbf{A}_F$ and $\mathbf{B}_G$.

**(b)** Define $(F + G)(\mathbf{x}) := F(\mathbf{x}) + G(\mathbf{x})$, for all $\mathbf{x} \in \mathbb{R}^3$. Is $(F + G)$ a linear transformation? If it is, find its standard matrix.

**(c)** Write down the formula for $F(G(\mathbf{x}))$ and find its standard matrix.

**(d)** Find a linear transformation $H : \mathbb{R}^3 \to \mathbb{R}^3$ such that $H(G(\mathbf{x})) = \mathbf{x}$, for all $\mathbf{x} \in \mathbb{R}^3$.

## Problem 3

For each of the following linear transformations, (i) determine whether there is enough information for us to find the formula of $T$; and (ii) find the formula and the standard matrix for $T$ if possible.

**(a)** $T : \mathbb{R}^3 \to \mathbb{R}^4$ such that $$T\left( \begin{pmatrix} 1 \\ 0 \\ 0 \end{pmatrix} \right) = \begin{pmatrix} 1 \\ 3 \\ 0 \\ 1 \end{pmatrix}, \quad T\left( \begin{pmatrix} 0 \\ 1 \\ 0 \end{pmatrix} \right) = \begin{pmatrix} 2 \\ 2 \\ -1 \\ 4 \end{pmatrix}, \quad T\left( \begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix} \right) = \begin{pmatrix} 0 \\ 4 \\ 1 \\ 6 \end{pmatrix}$$

**(b)** $T : \mathbb{R}^2 \to \mathbb{R}^2$ such that $$T\left( \begin{pmatrix} 1 \\ -1 \end{pmatrix} \right) = \begin{pmatrix} 2 \\ 0 \end{pmatrix}, \quad T\left( \begin{pmatrix} 1 \\ 1 \end{pmatrix} \right) = \begin{pmatrix} 0 \\ 2 \end{pmatrix}, \quad T\left( \begin{pmatrix} 2 \\ 0 \end{pmatrix} \right) = \begin{pmatrix} 2 \\ 2 \end{pmatrix}$$

**(c)** $T : \mathbb{R}^3 \to \mathbb{R}$ such that $$T\left( \begin{pmatrix} 1 \\ -1 \\ 0 \end{pmatrix} \right) = -1, \quad T\left( \begin{pmatrix} 0 \\ 1 \\ -1 \end{pmatrix} \right) = 1, \quad T\left( \begin{pmatrix} -1 \\ 0 \\ 1 \end{pmatrix} \right) = 0$$

## Problem 4

For each of the following linear transformations $T$, determine its rank and nullity, and whether it is one-to-one, and/or onto.

**(a)** $T : \mathbb{R}^4 \to \mathbb{R}^6$ such that the rank is 4.

**(b)** $T : \mathbb{R}^6 \to \mathbb{R}^4$ such that the nullity is 2.

**(c)** $T : \mathbb{R}^4 \to \mathbb{R}^6$ such that the reduced row-echelon form of its standard matrix has 3 nonzero rows.

**(d)** $T : \mathbb{R}^3 \to \mathbb{R}^3$ such that $T$ is one-to-one.

