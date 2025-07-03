import IPython.display
import sympy as sym

# from symbolic import Matrix


def is_IPython():
    # Adapted from https://stackoverflow.com/questions/15411967/how-can-i-check-if-code-is-executed-in-the-ipython-notebook
    try:
        shell = get_ipython().__class__.__name__
        if shell == "ZMQInteractiveShell" or shell == "TerminalInteractiveShell":
            return True  # Jupyter notebook, qtconsole or terminal running IPython
        else:
            return False  # Other type
    except NameError:
        return False  # Probably standard Python interpreter


def aug_print(matrix) -> None:
    # get latex representation of matrix
    raw = sym.latex(matrix, mat_str="array")

    # create formatting string s to insert augment line visually
    ls = [pos for pos in matrix._aug_pos if 0 < pos < matrix.cols]
    ls.sort()
    delta = [ls[0]]
    delta.extend([ls[i] - ls[i - 1] for i in range(1, len(ls))])
    remainder = matrix.cols - sum(delta)
    delta.append(remainder)
    s = "{" + "|".join(["c" * i for i in delta]) + "}"
    default_s = "{" + "c" * matrix.cols + "}"

    formatted = raw.replace(default_s, s)
    display(IPython.display.Math(formatted))


def display(input) -> None:
    if is_IPython():
        if hasattr(input, "_aug_pos") and len(input._aug_pos) > 0:
            aug_print(input)
        else:
            IPython.display.display(input)
    else:
        sym.pprint(input)
