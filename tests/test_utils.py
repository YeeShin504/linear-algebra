import sympy as sym
from unittest.mock import patch

from ma1522.utils import _powerset, _is_zero


class TestUtils:
    def test_powerset(self):
        assert list(_powerset([1, 2])) == [(), (1,), (2,), (1, 2)]
        assert list(_powerset([])) == [()]

    def test_is_zero(self):
        a, b = sym.symbols("a b", real=True)
        assert _is_zero(-a * b / 2 - a + b) is True
        assert _is_zero(a - b) is True
        assert _is_zero(0) is True
        assert _is_zero(a * b) is True
        assert _is_zero(b) is True
        assert _is_zero(1) is False

        n, d = sym.fraction(1)
        assert _is_zero(d) is False

    def test_textify(self):
        from ma1522.utils import _textify
        assert _textify("hello_world") == r"\text{hello\_world}"
        assert _textify("Simple") == r"\text{Simple}"

    def test_latex_wrapping(self):
        from ma1522.utils import _wrap_latex, _unwrap_latex
        assert _wrap_latex("x") == "$x$"
        assert _unwrap_latex("$x$") == "x"
        assert _unwrap_latex("$$x$$") == "x"
        assert _unwrap_latex(None) == ""

    def test_ipython_detection(self):
        """Mock various environments to verify IPython detection logic."""
        from ma1522.utils import _is_IPython
        from unittest.mock import MagicMock

        # Standard Python (get_ipython doesn't exist)
        with patch("IPython.core.getipython.get_ipython", side_effect=NameError):
            assert _is_IPython() is False

        # Jupyter/ZMQ Shell (Success)
        with patch("IPython.core.getipython.get_ipython") as mock_get:
            mock_shell = MagicMock()
            mock_shell.__class__.__name__ = "ZMQInteractiveShell"
            mock_get.return_value = mock_shell
            assert _is_IPython() is True

        # Terminal IPython (Success)
        with patch("IPython.core.getipython.get_ipython") as mock_get:
            mock_shell = MagicMock()
            mock_shell.__class__.__name__ = "TerminalInteractiveShell"
            mock_get.return_value = mock_shell
            assert _is_IPython() is True

        # Unrecognized Shell (Failure)
        with patch("IPython.core.getipython.get_ipython") as mock_get:
            mock_shell = MagicMock()
            mock_shell.__class__.__name__ = "UnknownShell"
            mock_get.return_value = mock_shell
            assert _is_IPython() is False

    def test_ipython_detection_without_ipython_installed(self):
        from ma1522.utils import _is_IPython

        def fake_import(name, *args, **kwargs):
            if name == "IPython.core.getipython":
                raise ImportError
            return original_import(name, *args, **kwargs)

        original_import = __import__
        with patch("builtins.__import__", side_effect=fake_import):
            assert _is_IPython() is False

    def test_display_ipython_and_fallback(self):
        from ma1522.utils import display
        
        # IPython SUCCESS
        with patch("ma1522.utils._is_IPython", return_value=True):
            with patch("IPython.display.display") as mock_disp:
                display({"a": 1}, 5, opt="dict")
                assert mock_disp.called

                display("math expression", opt="math")
                assert mock_disp.called
                
                display("standard")
                assert mock_disp.called
        # IPython ABSENT (Fallback to print)
        with patch("ma1522.utils._is_IPython", return_value=False):
            with patch("builtins.print") as mock_print:
                display("hello world")
                assert mock_print.called
                # Check if it was called with the right argument
                mock_print.assert_any_call("hello world")

    def test_gen_latex_repr_dict_special(self):
        from ma1522.utils import _gen_latex_repr_dict
        
        class Special:
            def _repr_latex_(self):
                return "$special$"
            def __repr__(self):
                return "$special$"
        
        obj = {"key": Special()}
        res = _gen_latex_repr_dict(obj)
        assert "special" in res

    def test_standardise_symbol(self):
        from ma1522.utils import _standardise_symbol
        import sympy as sym
        symbols = {sym.Symbol("x1"), sym.Symbol("y2")}
        std = _standardise_symbol(symbols)
        std_names = [str(s) for s in std]
        assert "x_1" in std_names
        assert "y_2" in std_names
