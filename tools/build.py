#!/usr/bin/env python3
"""Build TI-84 Plus CE program files (.8xp) from the plain-text TI-BASIC sources.

Usage:
    python3 tools/build.py            # builds every src/*.tib into bin/*.8xp

Source format (*.tib):
    * One TI-BASIC line per text line, written with the usual Unicode token names
      (→ ≠ ≤ ≥ ⁻ ⁻¹ ᵀ ² √( ʟ ►Frac 𝑖 θ [A] ...), exactly as shown on the calculator.
    * Lines starting with "//" are developer comments and are not sent to the calculator.
    * Leading indentation is stripped.
    * The program name is the file name (upper-cased, max 8 characters).

The token table comes from the TI-Toolkit token sheet (8X.xml,
https://github.com/TI-Toolkit/tokens). If it is not found next to this script it is
downloaded once and cached as tools/8X.xml.
"""

from __future__ import annotations

import sys
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TOKEN_URL = "https://raw.githubusercontent.com/TI-Toolkit/tokens/main/8X.xml"



def load_tokens(path: Path) -> dict[str, bytes]:
    if not path.exists():
        print(f"Downloading token sheet to {path} ...")
        urllib.request.urlretrieve(TOKEN_URL, path)
    root = ET.parse(path).getroot()
    table: dict[str, bytes] = {}

    def add(tok: ET.Element, prefix: bytes) -> None:
        value = bytes([int(tok.get("value")[1:], 16)])
        versions = tok.findall("version")
        if not versions:
            return
        lang = versions[-1].find("lang")
        if lang is None:
            return
        names = [lang.get("display")]
        names += [e.text for e in lang.findall("accessible")]
        names += [e.text for e in lang.findall("variant")]
        for n in names:
            if n and n not in table:
                table[n] = prefix + value

    for child in root:
        if child.tag == "token":
            add(child, b"")
        elif child.tag == "two-byte":
            prefix = bytes([int(child.get("value")[1:], 16)])
            for tok in child.findall("token"):
                add(tok, prefix)
                # Lowercase letters a-z: prefer the plain-letter tokens ($BB group) over
                # the statistics variables (a, b, c, d, e, r ...) that share the same glyph.
                lang = tok.findall("version")[-1].find("lang")
                disp = lang.get("display") if lang is not None else None
                if prefix == b"\xbb" and disp and len(disp) == 1 and disp.islower():
                    table[disp] = prefix + bytes([int(tok.get("value")[1:], 16)])
    return table


def tokenize(src: str, table: dict[str, bytes], name: str) -> bytes:
    maxlen = max(len(k) for k in table)
    out = bytearray()
    lines = []
    for raw in src.splitlines():
        line = raw.strip()
        if not line or line.startswith("//"):
            continue
        lines.append(line)
    for ln, line in enumerate(lines, 1):
        i = 0
        in_str = False
        while i < len(line):
            ch = line[i]
            if in_str:
                # Inside strings only use one-character tokens so text is shown verbatim.
                if ch == '"':
                    in_str = False
                if ch not in table:
                    raise SystemExit(f"{name}:{ln}: no token for {ch!r} in string")
                out += table[ch]
                i += 1
                continue
            for L in range(min(maxlen, len(line) - i), 0, -1):
                piece = line[i : i + L]
                if piece in table:
                    out += table[piece]
                    i += L
                    if piece == '"':
                        in_str = True
                    elif piece in ("prgm", "ʟ"):
                        # A program / list name is a run of single letter and digit tokens.
                        # Never let greedy matching turn e.g. "MS" in GRAMSCH into a stats token.
                        while i < len(line) and (line[i] in "θ" or (line[i].isascii() and line[i].isalnum())):
                            out += table[line[i]]
                            i += 1
                    break
            else:
                raise SystemExit(f"{name}:{ln}: cannot tokenize {line[i:]!r}")
        out.append(0x3F)  # newline
    if out and out[-1] == 0x3F:
        out = out[:-1]
    return bytes(out)


def make_entry(name: str, tokens: bytes) -> bytes:
    name_b = name.encode("ascii").ljust(8, b"\x00")
    var_data = len(tokens).to_bytes(2, "little") + tokens
    return (
        b"\x0d\x00"
        + len(var_data).to_bytes(2, "little")
        + b"\x05"  # type: program
        + name_b
        + b"\x00"  # version
        + b"\x00"  # flag (RAM)
        + len(var_data).to_bytes(2, "little")
        + var_data
    )


def make_file(entries: list[bytes]) -> bytes:
    """A computer-side TI file holding the given variable entries."""
    body = b"".join(entries)
    if len(body) > 0xFFFF:
        raise SystemExit("too much data for one file")
    comment = b"MA1522 linear algebra for TI-84 Plus CE".ljust(42, b"\x00")
    header = b"**TI83F*\x1a\x0a\x00" + comment + len(body).to_bytes(2, "little")
    checksum = (sum(body) & 0xFFFF).to_bytes(2, "little")
    return header + body + checksum


def make_8xp(name: str, tokens: bytes) -> bytes:
    return make_file([make_entry(name, tokens)])


def main() -> None:
    table = load_tokens(HERE / "8X.xml")
    src_dir = ROOT / "src"
    out_dir = ROOT / "bin"
    out_dir.mkdir(exist_ok=True)
    files = sorted(src_dir.glob("*.tib"))
    if not files:
        raise SystemExit("no sources found")
    for f in files:
        name = f.stem.upper()
        if len(name) > 8 or not name[0].isalpha() or not name.isalnum():
            raise SystemExit(f"bad program name {name!r}")
        toks = tokenize(f.read_text(encoding="utf-8"), table, f.name)
        (out_dir / f"{name}.8xp").write_bytes(make_8xp(name, toks))
        print(f"{name:8s} {len(toks):6d} bytes")


if __name__ == "__main__":
    sys.exit(main())
