#!/usr/bin/env python3
"""A small TI-BASIC interpreter for testing the programs in src/ on a computer.

It understands the subset of TI-84 Plus CE BASIC used by these programs (real/complex
numbers, lists, matrices, If/Then/Else, For/While/Repeat, Lbl/Goto, Menu, prgm calls,
matrix row operations, ...). It is a testing aid, not a full emulator: floating point is
IEEE double instead of the calculator's 14-digit BCD, and the screen is printed as text.

    from sim import Calc
    c = Calc("src")
    c.set_matrix("A", [[1, 2], [3, 4]])
    c.run("RREF", inputs=[0, 0])
    print(c.matrix("C"))
"""

from __future__ import annotations

import cmath
import math
import re
from fractions import Fraction
from pathlib import Path


class TIError(Exception):
    pass


class _Stop(Exception):
    pass


class _Return(Exception):
    pass


# --------------------------------------------------------------------------- lexer
FUNCS = [
    "identity(", "augment(", "rowSwap(", "*row+(", "*row(", "row+(", "Matr►list(",
    "List►matr(", "dim(", "max(", "min(", "abs(", "round(", "sum(", "sub(", "not(",
    "seq(", "real(", "imag(", "det(", "√(", "int(", "iPart(", "fPart(", "length(",
    "Fill(", "SortA(", "SortD(", "Menu(", "Output(", "For(",
]
WORDS = [
    "If ", "Then", "Else", "End", "While ", "Repeat ", "Return", "Stop", "Lbl ", "Goto ",
    "Disp ", "Pause ", "Pause", "Input ", "Prompt ", "ClrHome", "prgm", "DelVar ",
    "a+b𝑖", "Real", " and ", " or ", " xor ", "►Frac",
]
OPS = ["→", "⁻¹", "ᵀ", "²", "⁻", "≠", "≤", "≥", "=", "<", ">", "+", "-", "*", "/", "^",
       "(", ")", ",", "{", "}", "[", "]", "𝑖", "π", '"']
FIXED = sorted(FUNCS + WORDS + OPS, key=len, reverse=True)
NUM_RE = re.compile(r"(\d+\.?\d*|\.\d+)(ᴇ⁻?\d+)?")
MAT_RE = re.compile(r"\[([A-J])\]")
LIST_RE = re.compile(r"ʟ([A-Zθ][A-Z0-9θ]{0,4})|L([₁₂₃₄₅₆])")
STR_RE = re.compile(r"Str(\d)")


def lex(line: str) -> list[tuple[str, object]]:
    out = []
    i = 0
    while i < len(line):
        ch = line[i]
        if ch == '"':
            j = line.find('"', i + 1)
            if j < 0:
                j = len(line)
            out.append(("str", line[i + 1 : j]))
            i = j + 1
            continue
        m = MAT_RE.match(line, i)
        if m:
            out.append(("mat", m.group(1)))
            i = m.end()
            continue
        m = LIST_RE.match(line, i)
        if m:
            name = m.group(1) or ("L" + "₁₂₃₄₅₆".index(m.group(2)).__add__(1).__str__())
            out.append(("list", name))
            i = m.end()
            continue
        m = STR_RE.match(line, i)
        if m:
            out.append(("strvar", m.group(1)))
            i = m.end()
            continue
        for f in FIXED:
            if line.startswith(f, i):
                out.append(("op", f))
                i += len(f)
                break
        else:
            m = NUM_RE.match(line, i)
            if m:
                txt = m.group(0).replace("ᴇ", "e").replace("⁻", "-")
                out.append(("num", float(txt)))
                i = m.end()
            elif ch.isupper() or ch == "θ":
                out.append(("var", ch))
                i += 1
            elif ch == " ":
                i += 1
            else:
                raise TIError(f"lex error at {line[i:]!r}")
    return out


def split_statements(text: str) -> list[str]:
    stmts = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("//"):
            continue
        cur, in_str = "", False
        for ch in line:
            if ch == '"':
                in_str = not in_str
            if ch == ":" and not in_str:
                stmts.append(cur)
                cur = ""
            else:
                cur += ch
        stmts.append(cur)
    return [s for s in stmts if s != ""]


# --------------------------------------------------------------------------- values
def is_mat(x):
    return isinstance(x, list) and x and isinstance(x[0], list)


def is_list(x):
    return isinstance(x, list) and not is_mat(x)


def num(x):
    if isinstance(x, complex) and x.imag == 0:
        return x.real
    return x


def fmt_num(x, frac=False):
    x = num(x)
    if isinstance(x, complex):
        return f"{fmt_num(x.real, frac)}+{fmt_num(x.imag, frac)}i"
    if frac:
        f = Fraction(x).limit_denominator(10000)
        if abs(float(f) - x) < 1e-10 * max(1, abs(x)):
            return str(f)
    return f"{x:.10g}"


def fmt(v, frac=False):
    if isinstance(v, str):
        return v
    if is_mat(v):
        rows = [[fmt_num(e, frac) for e in r] for r in v]
        w = max(len(e) for r in rows for e in r)
        return "\n".join("[" + " ".join(e.rjust(w) for e in r) + "]" for r in rows)
    if isinstance(v, list):
        return "{" + ",".join(fmt_num(e, frac) for e in v) + "}"
    return fmt_num(v, frac)


# --------------------------------------------------------------------------- interpreter
def nonreal(x) -> bool:
    return isinstance(x, complex) and x.imag != 0


class Calc:
    def __init__(self, src_dir: str | Path, echo: bool = False):
        self.src_dir = Path(src_dir)
        self.reals: dict[str, object] = {c: 0.0 for c in "ABCDEFGHIJKLMNOPQRSTUVWXYZθ"}
        self.lists: dict[str, list] = {f"L{i}": [] for i in range(1, 7)}
        # On the calculator a list is either real or complex. Storing a nonreal number into
        # an element of a real list is ERR:DATA TYPE, so track which lists are complex.
        self.clists: set[str] = set()
        self.mats: dict[str, list] = {}
        self.strs: dict[str, str] = {}
        self.complex_mode = False
        self.inputs: list = []
        self.screen: list[str] = []
        self.echo = echo
        self._progs: dict[str, list[str]] = {}
        self.steps = 0

    # ---- public helpers
    def set_matrix(self, name, rows):
        self.mats[name] = [[float(e) for e in r] for r in rows]

    def matrix(self, name):
        return [[num(e) for e in r] for r in self.mats[name]]

    def list(self, name):
        return [num(e) for e in self.lists[name]]

    def run(self, prog: str, inputs=()):
        self.inputs = list(inputs)
        try:
            self._call(prog)
        except _Stop:
            pass
        return self

    def out(self, s):
        self.screen.append(s)
        if self.echo:
            print(s)

    # ---- program execution
    def _load(self, name):
        if name not in self._progs:
            path = self.src_dir / f"{name}.tib"
            self._progs[name] = [lex(s) for s in split_statements(path.read_text(encoding="utf-8"))]
        return self._progs[name]

    def _call(self, name):
        prog = self._load(name)
        try:
            self._exec(prog, name)
        except _Return:
            pass

    @staticmethod
    def _kind(toks):
        t = toks[0]
        return t[1] if t[0] == "op" else None

    def _skip_block(self, prog, pc, stop_at_else):
        """pc points just after a block opener. Return index of matching End (or Else)."""
        depth = 0
        while pc < len(prog):
            k = self._kind(prog[pc])
            if k in ("Then", "For(", "While ", "Repeat "):
                depth += 1
            elif k == "End":
                if depth == 0:
                    return pc
                depth -= 1
            elif k == "Else" and depth == 0 and stop_at_else:
                return pc
            pc += 1
        raise TIError("missing End")

    def _exec(self, prog, name):
        pc = 0
        stack: list = []
        while pc < len(prog):
            self.steps += 1
            if self.steps > 5_000_000:
                raise TIError("step limit")
            toks = prog[pc]
            k = self._kind(toks)
            try:
                if k == "If ":
                    cond = self._truth(self.eval(toks[1:]))
                    if pc + 1 < len(prog) and self._kind(prog[pc + 1]) == "Then":
                        if cond:
                            stack.append(("if",))
                            pc += 2
                        else:
                            j = self._skip_block(prog, pc + 2, True)
                            if self._kind(prog[j]) == "Else":
                                stack.append(("if",))
                            pc = j + 1
                    else:
                        pc += 1 if cond else 2
                    continue
                if k == "Else":
                    j = self._skip_block(prog, pc + 1, False)
                    stack.pop()
                    pc = j + 1
                    continue
                if k == "End":
                    fr = stack.pop()
                    if fr[0] == "for":
                        _, var, end, step, body = fr
                        v = self.reals[var] + step
                        self.reals[var] = v
                        if (step > 0 and v <= end) or (step < 0 and v >= end):
                            stack.append(fr)
                            pc = body
                            continue
                    elif fr[0] == "while":
                        pc = fr[1]
                        continue
                    elif fr[0] == "repeat":
                        if not self._truth(self.eval(fr[2])):
                            stack.append(fr)
                            pc = fr[1]
                            continue
                    pc += 1
                    continue
                if k == "For(":
                    args = self._args(toks[1:])
                    var = args[0][0][1]
                    a = self.eval(args[1])
                    b = self.eval(args[2])
                    s = self.eval(args[3]) if len(args) > 3 else 1.0
                    self.reals[var] = a
                    if (s > 0 and a > b) or (s < 0 and a < b):
                        pc = self._skip_block(prog, pc + 1, False) + 1
                    else:
                        stack.append(("for", var, b, s, pc + 1))
                        pc += 1
                    continue
                if k == "While ":
                    if self._truth(self.eval(toks[1:])):
                        stack.append(("while", pc))
                        pc += 1
                    else:
                        pc = self._skip_block(prog, pc + 1, False) + 1
                    continue
                if k == "Repeat ":
                    stack.append(("repeat", pc + 1, toks[1:]))
                    pc += 1
                    continue
                if k == "Goto ":
                    pc = self._find_label(prog, toks[1:])
                    continue
                if k == "Menu(":
                    args = self._args(toks[1:])
                    opts = [(self.eval(args[i]), args[i + 1]) for i in range(1, len(args), 2)]
                    self.out(f"MENU {self.eval(args[0])}: " + " | ".join(o[0] for o in opts))
                    choice = int(self.inputs.pop(0))
                    self.out(f"  -> {opts[choice - 1][0]}")
                    pc = self._find_label(prog, opts[choice - 1][1])
                    continue
                self._statement(toks)
            except (_Return, _Stop):
                raise
            except TIError as e:
                raise TIError(f"{name} stmt {pc + 1}: {e}") from None
            except Exception as e:
                raise TIError(f"{name} stmt {pc + 1}: {type(e).__name__}: {e}") from None
            pc += 1

    def _find_label(self, prog, lbl_toks):
        want = "".join(str(t[1]) for t in lbl_toks)
        for i, t in enumerate(prog):
            if self._kind(t) == "Lbl " and "".join(str(x[1]) for x in t[1:]) == want:
                return i + 1
        raise TIError(f"ERR:LABEL {want}")

    def _statement(self, toks):
        k = self._kind(toks)
        if k in ("Lbl ", "Then"):
            return
        if k == "Return":
            raise _Return
        if k == "Stop":
            raise _Stop
        if k == "ClrHome":
            self.out("-" * 26)
            return
        if k == "a+b𝑖":
            self.complex_mode = True
            return
        if k == "Real":
            self.complex_mode = False
            return
        if k == "prgm":
            self._call("".join(t[1] for t in toks[1:]))
            return
        if k in ("Disp ", "Pause ", "Pause"):
            if len(toks) == 1:
                return
            for a in self._args(toks[1:]):
                frac = a and a[-1] == ("op", "►Frac")
                self.out(fmt(self.eval(a[:-1] if frac else a), frac))
            return
        if k == "Output(":
            a = self._args(toks[1:])
            frac = a[2][-1] == ("op", "►Frac")
            self.out(fmt(self.eval(a[2][:-1] if frac else a[2]), frac))
            return
        if k in ("Input ", "Prompt "):
            a = self._args(toks[1:])
            target = a[-1]
            val = self.inputs.pop(0)
            if k == "Input " and len(a) > 1:
                self.out(f"{self.eval(a[0])}{val}")
            self._assign(target, float(val))
            return
        if k == "DelVar ":
            t = toks[1]
            if t[0] == "list":
                self.lists.pop(t[1], None)
                self.clists.discard(t[1])
            elif t[0] == "mat":
                self.mats.pop(t[1], None)
            return
        if k == "Fill(":
            a = self._args(toks[1:])
            v = self.eval(a[0])
            t = a[1][0]
            if t[0] == "list":
                self.lists[t[1]] = [v] * len(self.lists[t[1]])
            else:
                m = self.mats[t[1]]
                self.mats[t[1]] = [[v] * len(m[0]) for _ in m]
            return
        if k in ("SortA(", "SortD("):
            a = self._args(toks[1:])
            names = [x[0][1] for x in a]
            key = self.lists[names[0]]
            order = sorted(range(len(key)), key=lambda i: key[i], reverse=(k == "SortD("))
            for n in names:
                self.lists[n] = [self.lists[n][i] for i in order]
            return
        if k == "Matr►list(":
            a = self._args(toks[1:])
            m = self.eval(a[0])
            col = int(self.eval(a[1])) - 1
            self.lists[a[2][0][1]] = [r[col] for r in m]
            return
        if k == "List►matr(":
            a = self._args(toks[1:])
            lst = self.eval(a[0])
            self.mats[a[1][0][1]] = [[e] for e in lst]
            return
        # assignment or bare expression
        idx = [i for i, t in enumerate(toks) if t == ("op", "→")]
        if idx:
            i = idx[-1]
            self._assign(toks[i + 1 :], self.eval(toks[:i]))
        else:
            self.eval(toks)

    # ---- assignment
    def _assign(self, target, val):
        t = target[0]
        if t[0] == "var":
            if is_list(val) or is_mat(val) or isinstance(val, str):
                raise TIError("ERR:DATA TYPE")
            self.reals[t[1]] = self._chk(val)
            return
        if t[0] == "strvar":
            self.strs[t[1]] = val
            return
        if t == ("op", "dim("):
            obj = target[1]
            if obj[0] == "list":
                n = int(val)
                if obj[1] not in self.lists:
                    self.clists.discard(obj[1])  # dim( on a new list creates a real list
                old = self.lists.get(obj[1], [])
                self.lists[obj[1]] = (old + [0.0] * n)[:n]
            else:
                r, c = int(val[0]), int(val[1])
                old = self.mats.get(obj[1], [])
                self.mats[obj[1]] = [
                    [old[i][j] if i < len(old) and j < len(old[0]) else 0.0 for j in range(c)]
                    for i in range(r)
                ]
            return
        if t[0] == "list":
            if len(target) == 1:
                if not is_list(val):
                    raise TIError("ERR:DATA TYPE")
                self.lists[t[1]] = list(val)
                # Strict model: a stored list is complex only if some element is nonreal
                # (the calculator may turn e.g. 0*complex list into a real list).
                if any(nonreal(e) for e in val):
                    self.clists.add(t[1])
                else:
                    self.clists.discard(t[1])
                return
            i = int(self.eval(target[2:-1] if target[-1] == ("op", ")") else target[2:])) - 1
            lst = self.lists[t[1]]
            if i == len(lst):
                lst.append(0.0)
            if not 0 <= i < len(lst):
                raise TIError("ERR:INVALID DIM")
            if nonreal(val) and t[1] not in self.clists:
                raise TIError(f"ERR:DATA TYPE (nonreal value stored into real list ʟ{t[1]})")
            lst[i] = val
            return
        if t[0] == "mat":
            if len(target) == 1:
                if not is_mat(val):
                    raise TIError("ERR:DATA TYPE")
                self.mats[t[1]] = [list(r) for r in val]
                return
            inner = target[2:-1] if target[-1] == ("op", ")") else target[2:]
            a = self._args(inner)
            i, j = int(self.eval(a[0])) - 1, int(self.eval(a[1])) - 1
            m = self.mats[t[1]]
            if not (0 <= i < len(m) and 0 <= j < len(m[0])):
                raise TIError("ERR:INVALID DIM")
            if isinstance(val, complex) and val.imag != 0:
                raise TIError("ERR:DATA TYPE (complex in matrix)")
            m[i][j] = num(val)
            return
        raise TIError(f"bad assignment target {target}")

    def _chk(self, v):
        if isinstance(v, complex):
            if v.imag != 0 and not self.complex_mode:
                raise TIError("ERR:NONREAL ANS")
            if v.imag == 0:
                return v.real
        return v

    @staticmethod
    def _truth(v):
        if isinstance(v, list):
            raise TIError("ERR:DATA TYPE in condition")
        return v != 0

    # ---- expression parsing (Pratt-style recursive descent)
    def _args(self, toks):
        """Split tokens on top-level commas; drop a trailing ')' that closes the call."""
        if toks and toks[-1] == ("op", ")"):
            depth = 0
            for t in toks[:-1]:
                if t[0] == "op" and (t[1].endswith("(") or t[1] in "{["):
                    depth += 1
                elif t[0] == "op" and t[1] in (")", "}", "]"):
                    depth -= 1
            if depth == 0:
                toks = toks[:-1]
        parts, cur, depth = [], [], 0
        for t in toks:
            if t[0] == "op" and (t[1].endswith("(") or t[1] in ("{", "[")):
                depth += 1
            elif t[0] == "op" and t[1] in (")", "}", "]"):
                depth -= 1
            if t == ("op", ",") and depth == 0:
                parts.append(cur)
                cur = []
            else:
                cur.append(t)
        parts.append(cur)
        return parts

    def eval(self, toks):
        self._t = list(toks)
        self._p = 0
        v = self._or()
        if self._p != len(self._t):
            raise TIError(f"syntax: leftover {self._t[self._p:]}")
        return v

    def _peek(self):
        return self._t[self._p] if self._p < len(self._t) else (None, None)

    def _take(self):
        t = self._t[self._p]
        self._p += 1
        return t

    def _expect(self, s):
        t = self._peek()
        if t == ("op", s):
            self._p += 1
        elif t != (None, None):  # closing parens may be omitted at end of line
            raise TIError(f"expected {s}, got {t}")

    def _or(self):
        v = self._and()
        while self._peek() in (("op", " or "), ("op", " xor ")):
            op = self._take()[1]
            w = self._and()
            v = self._bin(v, w, (lambda a, b: float(bool(a) or bool(b))) if op == " or " else (lambda a, b: float(bool(a) != bool(b))))
        return v

    def _and(self):
        v = self._rel()
        while self._peek() == ("op", " and "):
            self._take()
            w = self._rel()
            v = self._bin(v, w, lambda a, b: float(bool(a) and bool(b)))
        return v

    RELS = {"=": lambda a, b: float(a == b), "≠": lambda a, b: float(a != b),
            "<": lambda a, b: float(a < b), ">": lambda a, b: float(a > b),
            "≤": lambda a, b: float(a <= b), "≥": lambda a, b: float(a >= b)}

    def _rel(self):
        v = self._add()
        while self._peek()[0] == "op" and self._peek()[1] in self.RELS:
            f = self.RELS[self._take()[1]]
            w = self._add()
            v = self._bin(v, w, f)
        return v

    def _add(self):
        v = self._mul()
        while self._peek() in (("op", "+"), ("op", "-")):
            op = self._take()[1]
            w = self._mul()
            if isinstance(v, str) and isinstance(w, str) and op == "+":
                v = v + w
            elif is_mat(v) or is_mat(w):
                v = self._matadd(v, w, 1 if op == "+" else -1)
            else:
                v = self._bin(v, w, (lambda a, b: a + b) if op == "+" else (lambda a, b: a - b))
        return v

    def _starts_primary(self):
        t = self._peek()
        if t[0] in ("num", "var", "list", "mat", "strvar"):
            return True
        return t[0] == "op" and (t[1] in ("(", "{", "𝑖", "π", "⁻") or t[1] in FUNCS) and t[1] != "⁻"

    def _mul(self):
        v = self._neg()
        while True:
            t = self._peek()
            if t in (("op", "*"), ("op", "/")):
                op = self._take()[1]
                w = self._neg()
            elif self._starts_primary():
                op, w = "*", self._neg()
            else:
                return v
            if op == "*" and (is_mat(v) or is_mat(w)):
                v = self._matmul(v, w)
            elif op == "/" and is_mat(v):
                v = [[self._div(e, w) for e in r] for r in v]
            else:
                v = self._bin(v, w, (lambda a, b: a * b) if op == "*" else self._div)
            v = self._fix(v)

    def _neg(self):
        if self._peek() == ("op", "⁻"):
            self._take()
            v = self._neg()
            if is_mat(v):
                return [[-e for e in r] for r in v]
            return self._map(v, lambda a: -a)
        return self._pow()

    def _pow(self):
        v = self._postfix()
        while self._peek() == ("op", "^"):
            self._take()
            if self._peek() == ("op", "⁻"):
                self._take()
                w = self._map(self._postfix(), lambda a: -a)
            else:
                w = self._postfix()
            if is_mat(v):
                n = int(w)
                r = self._identity(len(v))
                for _ in range(n):
                    r = self._matmul(r, v)
                v = r
            else:
                v = self._bin(v, w, self._power)
        return v

    def _power(self, a, b):
        if a == 0 and b == 0:
            raise TIError("ERR:DOMAIN 0^0")
        if isinstance(a, complex) or isinstance(b, complex):
            return self._chk(complex(a) ** b)
        if a < 0 and b != int(b):
            if not self.complex_mode:
                raise TIError("ERR:NONREAL ANS")
            return complex(a) ** b
        return float(a) ** b

    def _postfix(self):
        v = self._primary()
        while True:
            t = self._peek()
            if t == ("op", "²"):
                self._take()
                v = self._matmul(v, v) if is_mat(v) else self._map(v, lambda a: a * a)
            elif t == ("op", "ᵀ"):
                self._take()
                v = [list(c) for c in zip(*v)]
            elif t == ("op", "⁻¹"):
                self._take()
                v = self._inverse(v) if is_mat(v) else self._map(v, lambda a: self._div(1, a))
            else:
                return v

    def _fix(self, v):
        if isinstance(v, complex):
            return self._chk(v)
        return v

    def _primary(self):
        t = self._take()
        kind, val = t
        if kind == "num":
            return val
        if kind == "str":
            return val
        if kind == "var":
            return self.reals[val]
        if kind == "strvar":
            return self.strs[val]
        if kind == "list":
            if val not in self.lists:
                raise TIError(f"ERR:UNDEFINED ʟ{val}")
            lst = self.lists[val]
            if self._peek() == ("op", "("):
                self._take()
                i = int(self._or()) - 1
                self._expect(")")
                if not 0 <= i < len(lst):
                    raise TIError(f"ERR:INVALID DIM ʟ{val}({i + 1})")
                return lst[i]
            return list(lst)
        if kind == "mat":
            if val not in self.mats:
                raise TIError(f"ERR:UNDEFINED [{val}]")
            m = self.mats[val]
            if self._peek() == ("op", "("):
                self._take()
                i = int(self._or()) - 1
                self._expect(",")
                j = int(self._or()) - 1
                self._expect(")")
                if not (0 <= i < len(m) and 0 <= j < len(m[0])):
                    raise TIError(f"ERR:INVALID DIM [{val}]({i + 1},{j + 1})")
                return m[i][j]
            return [list(r) for r in m]
        if kind == "op":
            if val == "(":
                v = self._or()
                self._expect(")")
                return v
            if val == "𝑖":
                if not self.complex_mode:
                    pass  # literal i is allowed; results are checked on store
                return 1j
            if val == "π":
                return math.pi
            if val == "{":
                items = []
                while self._peek() != ("op", "}"):
                    items.append(self._or())
                    if self._peek() == ("op", ","):
                        self._take()
                    elif self._peek() == (None, None):
                        break
                self._expect("}")
                return items
            if val == "[":
                rows = []
                while self._peek() == ("op", "["):
                    self._take()
                    row = []
                    while self._peek() != ("op", "]"):
                        row.append(self._or())
                        if self._peek() == ("op", ","):
                            self._take()
                    self._take()
                    rows.append(row)
                self._expect("]")
                return rows
            if val in FUNCS:
                return self._func(val)
        raise TIError(f"syntax: unexpected {t}")

    def _func(self, name):
        if name == "seq(":
            # seq(expr,var,start,end[,step]) - parse lazily
            start = self._p
            depth = 0
            parts, cur = [], []
            while True:
                t = self._take()
                if t[0] == "op" and (t[1].endswith("(") or t[1] in "{["):
                    depth += 1
                if t[0] == "op" and t[1] in (")", "}", "]"):
                    if depth == 0:
                        break
                    depth -= 1
                if t == ("op", ",") and depth == 0:
                    parts.append(cur)
                    cur = []
                else:
                    cur.append(t)
            parts.append(cur)
            saved = (self._t, self._p)
            var = parts[1][0][1]
            a, b = self.eval(parts[2]), self.eval(parts[3])
            s = self.eval(parts[4]) if len(parts) > 4 else 1
            out = []
            x = a
            while (s > 0 and x <= b) or (s < 0 and x >= b):
                self.reals[var] = x
                out.append(self.eval(parts[0]))
                if nonreal(out[-1]) and not nonreal(out[0]):
                    # Strict model: the list type follows the first element.
                    raise TIError("ERR:DATA TYPE (seq( first element real, later one complex)")
                x += s
            self._t, self._p = saved
            return out
        args = []
        if self._peek() != ("op", ")"):
            while True:
                args.append(self._or())
                if self._peek() == ("op", ","):
                    self._take()
                    continue
                break
        self._expect(")")
        return self._call_func(name, args)

    def _call_func(self, name, a):
        f = name
        if f == "dim(":
            x = a[0]
            return [float(len(x)), float(len(x[0]))] if is_mat(x) else float(len(x))
        if f in ("max(", "min("):
            g = max if f == "max(" else min
            if len(a) == 1:
                if not a[0]:
                    raise TIError("ERR:INVALID DIM")
                return g(a[0])
            return self._bin(a[0], a[1], lambda x, y: g(x, y))
        if f == "abs(":
            if is_mat(a[0]):
                return [[abs(e) for e in r] for r in a[0]]
            return self._map(a[0], abs)
        if f == "round(":
            n = int(a[1]) if len(a) > 1 else 9
            rnd = lambda x: complex(round(x.real, n), round(x.imag, n)) if isinstance(x, complex) else float(round(x, n))
            if is_mat(a[0]):
                return [[rnd(e) for e in r] for r in a[0]]
            return self._map(a[0], rnd)
        if f == "sum(":
            if not a[0]:
                raise TIError("ERR:INVALID DIM")
            return sum(a[0])
        if f == "sub(":
            return a[0][int(a[1]) - 1 : int(a[1]) - 1 + int(a[2])]
        if f == "not(":
            return self._map(a[0], lambda x: float(x == 0))
        if f == "real(":
            return self._map(a[0], lambda x: float(x.real) if isinstance(x, complex) else x)
        if f == "imag(":
            return self._map(a[0], lambda x: float(x.imag) if isinstance(x, complex) else 0.0)
        if f == "√(":
            def sq(x):
                if isinstance(x, complex) or x < 0:
                    if not self.complex_mode and not isinstance(x, complex):
                        raise TIError("ERR:NONREAL ANS")
                    return cmath.sqrt(x)
                return math.sqrt(x)
            return self._map(a[0], sq)
        if f in ("int(", "iPart("):
            return self._map(a[0], lambda x: float(math.floor(x) if f == "int(" else int(x)))
        if f == "fPart(":
            return self._map(a[0], lambda x: x - int(x))
        if f == "length(":
            return float(len(a[0]))
        if f == "det(":
            return self._det(a[0])
        if f == "identity(":
            return self._identity(int(a[0]))
        if f == "augment(":
            x, y = a
            if is_mat(x):
                if len(x) != len(y):
                    raise TIError("ERR:DIM MISMATCH")
                return [list(r) + list(s) for r, s in zip(x, y)]
            return list(x) + list(y)
        if f == "rowSwap(":
            m = [list(r) for r in a[0]]
            i, j = int(a[1]) - 1, int(a[2]) - 1
            m[i], m[j] = m[j], m[i]
            return m
        if f == "*row(":
            m = [list(r) for r in a[1]]
            i = int(a[2]) - 1
            m[i] = [a[0] * e for e in m[i]]
            return m
        if f == "*row+(":
            m = [list(r) for r in a[1]]
            i, j = int(a[2]) - 1, int(a[3]) - 1
            m[j] = [e + a[0] * d for e, d in zip(m[j], m[i])]
            return m
        if f == "row+(":
            m = [list(r) for r in a[0]]
            i, j = int(a[1]) - 1, int(a[2]) - 1
            m[j] = [e + d for e, d in zip(m[j], m[i])]
            return m
        raise TIError(f"unsupported function {name}")

    # ---- arithmetic helpers
    def _map(self, v, f):
        if isinstance(v, list):
            return [self._fix(f(e)) for e in v]
        return self._fix(f(v))

    def _bin(self, a, b, f):
        if is_mat(a) or is_mat(b):
            raise TIError("ERR:DATA TYPE (matrix in scalar op)")
        if isinstance(a, list) and isinstance(b, list):
            if len(a) != len(b):
                raise TIError("ERR:DIM MISMATCH")
            return [self._fix(f(x, y)) for x, y in zip(a, b)]
        if isinstance(a, list):
            return [self._fix(f(x, b)) for x in a]
        if isinstance(b, list):
            return [self._fix(f(a, y)) for y in b]
        return self._fix(f(a, b))

    @staticmethod
    def _div(a, b):
        if b == 0:
            raise TIError("ERR:DIVIDE BY 0")
        return a / b

    def _matadd(self, a, b, s):
        if not (is_mat(a) and is_mat(b)) or len(a) != len(b) or len(a[0]) != len(b[0]):
            raise TIError("ERR:DIM MISMATCH")
        return [[x + s * y for x, y in zip(r, q)] for r, q in zip(a, b)]

    def _matmul(self, a, b):
        if is_mat(a) and is_mat(b):
            if len(a[0]) != len(b):
                raise TIError("ERR:DIM MISMATCH")
            return [[sum(a[i][k] * b[k][j] for k in range(len(b))) for j in range(len(b[0]))] for i in range(len(a))]
        if is_mat(a):
            a, b = b, a
        if isinstance(a, list):
            raise TIError("ERR:DATA TYPE")
        return [[a * e for e in r] for r in b]

    @staticmethod
    def _identity(n):
        return [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]

    def _det(self, m):
        n = len(m)
        if n != len(m[0]):
            raise TIError("ERR:INVALID DIM")
        a = [list(r) for r in m]
        d = 1.0
        for c in range(n):
            p = max(range(c, n), key=lambda r: abs(a[r][c]))
            if abs(a[p][c]) < 1e-300:
                return 0.0
            if p != c:
                a[c], a[p] = a[p], a[c]
                d = -d
            d *= a[c][c]
            for r in range(c + 1, n):
                f = a[r][c] / a[c][c]
                a[r] = [x - f * y for x, y in zip(a[r], a[c])]
        return d

    def _inverse(self, m):
        n = len(m)
        a = [list(r) + [1.0 if i == j else 0.0 for j in range(n)] for i, r in enumerate(m)]
        for c in range(n):
            p = max(range(c, n), key=lambda r: abs(a[r][c]))
            if abs(a[p][c]) < 1e-13:
                raise TIError("ERR:SINGULAR MAT")
            a[c], a[p] = a[p], a[c]
            piv = a[c][c]
            a[c] = [x / piv for x in a[c]]
            for r in range(n):
                if r != c:
                    f = a[r][c]
                    a[r] = [x - f * y for x, y in zip(a[r], a[c])]
        return [r[n:] for r in a]


if __name__ == "__main__":
    import sys

    c = Calc(Path(__file__).resolve().parent.parent / "src", echo=True)
    print("usage: import sim; see module docstring", file=sys.stderr)
