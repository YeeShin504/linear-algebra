import sympy as sym
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from symbolic import Matrix


def _aug_line(self, pos: int = -1) -> "Matrix":
    """
    Inserts an augmented line at the specified position.

    This method adds an augmented line (i.e., a line for augmented matrix operations)
    to the matrix at the specified column position. If no position is provided (default: -1),
    the line is inserted at the last column.

    Parameters:
    - pos (int, optional): The position (column index) where the augmented line
        will be inserted. Default is -1, which means the augmented line is added
        at the end of the matrix.

    Returns:
    - Matrix: The current matrix with the augmented line added at the specified position.

    Example:
    >>> mat = Matrix([[1, 2], [3, 4]])
    >>> mat.aug_line(1)
    Matrix([[1, 2], [3, 4], [0, 0]])

    Notes:
    - The method updates the `__aug_pos` attribute to track the position of the inserted line.
    """

    new_pos = pos
    if new_pos < 0:
        new_pos += self.cols + 1

    if not 0 <= new_pos <= self.cols:
        raise IndexError(
            f"Position for augmented line ({pos}) out of range ({self.cols})."
        )

    self._aug_pos.add(new_pos)
    return self
