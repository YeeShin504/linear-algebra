from unittest.mock import MagicMock

import pytest

from ma1522 import Matrix
from ma1522._workers import _get_orth, _orthogonal_completion_worker


class _QueueStub:
    def __init__(self, items=None):
        self.items = list(items or [])
        self.put_calls = []

    def put(self, item):
        self.put_calls.append(item)

    def empty(self):
        return len(self.items) == 0

    def get(self):
        return self.items.pop(0)


class _ContextStub:
    def __init__(self, queue, process):
        self._queue = queue
        self._process = process

    def Queue(self, maxsize=1):
        return self._queue

    def Process(self, **kwargs):
        return self._process


class _ComplementStub:
    def __init__(self, result):
        self.result = result

    def gram_schmidt(self, factor, verbosity):
        return self.result


class _UWorkerSuccessStub:
    def __init__(self, result):
        self.result = result

    def orthogonal_complement(self, verbosity):
        return _ComplementStub(self.result)


class _UWorkerErrorStub:
    def orthogonal_complement(self, verbosity):
        raise ValueError("error")


def _build_context(queue_items, is_alive=False):
    queue = _QueueStub(items=queue_items)
    proc = MagicMock()
    proc.is_alive.return_value = is_alive
    ctx = _ContextStub(queue=queue, process=proc)
    return ctx, queue, proc


class TestOrthogonalCompletionWorker:
    def test_worker_puts_ok_payload(self):
        queue = _QueueStub()
        expected = Matrix([[1], [0]])
        _orthogonal_completion_worker(
            _UWorkerSuccessStub(expected), verbosity=0, queue=queue
        )

        assert queue.put_calls == [("ok", expected)]

    def test_worker_puts_error_payload(self):
        queue = _QueueStub()
        _orthogonal_completion_worker(_UWorkerErrorStub(), verbosity=0, queue=queue)

        assert queue.put_calls[0][0] == "err"
        assert "ValueError('error')" in queue.put_calls[0][1]


class TestGetOrth:
    def test_returns_worker_payload_on_success(self, monkeypatch):
        expected = Matrix([[1], [0]])
        ctx, _, proc = _build_context(queue_items=[("ok", expected)], is_alive=False)
        monkeypatch.setattr("ma1522._workers.mp.get_context", lambda _: ctx)

        result = _get_orth(Matrix([[1], [0]]), timeout_seconds=1, verbosity=0)

        assert result == expected
        proc.start.assert_called_once()
        proc.join.assert_called_once_with(1)
        proc.terminate.assert_not_called()

    def test_timeout_uses_numerical_fallback(self, monkeypatch):
        ctx, _, proc = _build_context(queue_items=[], is_alive=True)
        monkeypatch.setattr("ma1522._workers.mp.get_context", lambda _: ctx)

        u = Matrix([[1], [0], [0]])
        with pytest.warns(RuntimeWarning, match=r"Gram-Schmidt fallback \(timeout\)"):
            result = _get_orth(u, timeout_seconds=1, verbosity=0)

        assert result.shape == (3, 2)
        proc.terminate.assert_called_once()
        assert proc.join.call_args_list[-1][0] == (1,)

    def test_empty_queue_uses_numerical_fallback(self, monkeypatch):
        ctx, _, _ = _build_context(queue_items=[], is_alive=False)
        monkeypatch.setattr("ma1522._workers.mp.get_context", lambda _: ctx)

        with pytest.warns(
            RuntimeWarning,
            match=r"Gram-Schmidt fallback \(worker produced no result\)",
        ):
            result = _get_orth(Matrix([[1, 0], [0, 1]]), timeout_seconds=1, verbosity=0)

        assert result.shape == (0, 0)

    def test_ok_none_payload_uses_numerical_fallback(self, monkeypatch):
        ctx, _, _ = _build_context(queue_items=[("ok", None)], is_alive=False)
        monkeypatch.setattr("ma1522._workers.mp.get_context", lambda _: ctx)

        with pytest.warns(
            RuntimeWarning,
            match=r"Gram-Schmidt fallback \(worker returned None\)",
        ):
            result = _get_orth(Matrix([[1, 0], [0, 1]]), timeout_seconds=1, verbosity=0)

        assert result.shape == (0, 0)

    def test_worker_error_uses_numerical_fallback(self, monkeypatch):
        ctx, _, _ = _build_context(
            queue_items=[("err", "ValueError('x')")],
            is_alive=False,
        )
        monkeypatch.setattr("ma1522._workers.mp.get_context", lambda _: ctx)

        with pytest.warns(
            RuntimeWarning,
            match=r"Gram-Schmidt fallback \(worker failure\)",
        ):
            _get_orth(Matrix([[1, 0], [0, 1]]), timeout_seconds=1, verbosity=0)
