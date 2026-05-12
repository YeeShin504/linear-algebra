import multiprocessing as mp
from contextlib import suppress
from typing import TYPE_CHECKING
from warnings import warn

import numpy as np

if TYPE_CHECKING:
    from .symbolic import Matrix


def _orthogonal_completion_worker(U: "Matrix", verbosity: int, queue: mp.Queue) -> None:
    """Worker function to compute orthogonal complement and Gram-Schmidt."""
    try:
        complement = U.orthogonal_complement(verbosity=verbosity)
        orth = complement.gram_schmidt(factor=False, verbosity=verbosity)
        queue.put(("ok", orth))
    except Exception as error:
        queue.put(("err", repr(error)))


def _get_orth(U: "Matrix", timeout_seconds: int = 30, verbosity: int = 0) -> "Matrix":
    """Compute an orthonormal completion for U.

    Runs orthogonal-complement and Gram-Schmidt in a subprocess and falls back
    to a numerical QR-based completion on timeout or worker failure.

    Args:
        U: The input matrix to complete.
        timeout_seconds: The maximum time to wait for the worker process.
        verbosity: The level of verbosity for logging.

    Returns:
        An orthonormal completion of U.
    """

    def _numerical_fallback(reason: str, error: Exception | None = None) -> "Matrix":
        message = (
            f"Gram-Schmidt fallback ({reason}). "
            "Using numerical QR completion for orthogonal complement."
        )
        if error is not None:
            message = f"{message} Original error: {error}"
        warn(message, RuntimeWarning, stacklevel=2)

        # Numerical QR-based completion using U to produce orthonormal extra cols
        U_np = np.array(U.evalf().tolist(), dtype=np.complex128)
        q_np, _ = np.linalg.qr(np.column_stack([U_np, np.eye(U.rows)]))
        n_extra = U.rows - U.cols

        # Lazy import to avoid circular import at module import time
        from .symbolic import Matrix

        if n_extra <= 0:
            return Matrix()
        orth = Matrix(q_np[:, U.cols : U.cols + n_extra].tolist())
        with suppress(Exception):
            orth = orth.identify()
        return orth

    ctx = mp.get_context("spawn")
    queue: mp.Queue = ctx.Queue(maxsize=1)
    proc = ctx.Process(
        target=_orthogonal_completion_worker,
        args=(U, verbosity, queue),
        daemon=True,
    )
    proc.start()
    proc.join(timeout_seconds)

    if proc.is_alive():
        proc.terminate()
        proc.join(1)
        return _numerical_fallback("timeout")

    if queue.empty():
        return _numerical_fallback("worker produced no result")

    status, payload = queue.get()
    if status == "ok":
        if payload is None:
            return _numerical_fallback("worker returned None")
        return payload
    return _numerical_fallback("worker failure", Exception(payload))
