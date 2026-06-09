#!/usr/bin/env python3
"""Reject non-UTF-8 text and common mojibake in tracked source files."""

from __future__ import annotations

import argparse
import ast
import re
import subprocess
import sys
from pathlib import Path


BINARY_EXTENSIONS = {
    ".a",
    ".dll",
    ".egg",
    ".gcda",
    ".gcno",
    ".h5",
    ".hdf5",
    ".jpg",
    ".jpeg",
    ".lib",
    ".mod",
    ".npy",
    ".npz",
    ".o",
    ".pdf",
    ".png",
    ".pyd",
    ".pyc",
    ".pyo",
    ".smod",
    ".so",
    ".tar",
    ".whl",
    ".zip",
}

TEXT_EXTENSIONS = {
    "",
    ".c",
    ".cfg",
    ".csv",
    ".f",
    ".f90",
    ".h",
    ".in",
    ".ini",
    ".json",
    ".md",
    ".patch",
    ".ps1",
    ".py",
    ".pyi",
    ".rst",
    ".sh",
    ".toml",
    ".txt",
    ".yaml",
    ".yml",
}

MOJIBAKE_MARKERS = ("".join(chr(code) for code in (0x951F, 0x65A4, 0x62F7)), chr(0xFFFD))
LATIN1_MOJIBAKE_CHARS = "".join(
    chr(code) for code in (0x00C3, 0x00C2, 0x00E2, 0x00E4, 0x00E5, 0x00E6, 0x00E7, 0x00E8, 0x00E9, 0x00EF, 0x00F0)
)
LATIN1_MOJIBAKE = re.compile(f"[{LATIN1_MOJIBAKE_CHARS}][\\u0080-\\u00ff]")
C1_CONTROL = re.compile(r"[\u0080-\u009f]")
PRIVATE_USE = re.compile(r"[\ue000-\uf8ff]")


class PythonTextIoVisitor(ast.NodeVisitor):
    def __init__(self) -> None:
        self.problems: list[tuple[int, str]] = []

    def visit_Call(self, node: ast.Call) -> None:
        call_name = self._call_name(node.func)
        if call_name in {"open", "Path.open"}:
            if not self._has_keyword(node, "encoding") and not self._is_binary_open(node):
                self.problems.append((node.lineno, f"{call_name} text IO without explicit encoding"))
        elif call_name in {"Path.read_text", "Path.write_text"}:
            if not self._has_keyword(node, "encoding"):
                self.problems.append((node.lineno, f"{call_name} without explicit encoding"))
        self.generic_visit(node)

    @staticmethod
    def _call_name(func: ast.AST) -> str | None:
        if isinstance(func, ast.Name) and func.id == "open":
            return "open"
        if isinstance(func, ast.Attribute):
            if func.attr in {"open", "read_text", "write_text"}:
                return f"Path.{func.attr}"
        return None

    @staticmethod
    def _has_keyword(node: ast.Call, name: str) -> bool:
        return any(keyword.arg == name for keyword in node.keywords)

    @staticmethod
    def _is_binary_open(node: ast.Call) -> bool:
        mode: ast.AST | None = None
        if len(node.args) >= 2:
            mode = node.args[1]
        for keyword in node.keywords:
            if keyword.arg == "mode":
                mode = keyword.value
        return isinstance(mode, ast.Constant) and isinstance(mode.value, str) and "b" in mode.value


def run_git(root: Path, args: list[str]) -> list[str]:
    result = subprocess.run(
        ["git", *args],
        cwd=root,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    names = result.stdout.decode("utf-8").split("\0")
    return [name for name in names if name]


def repo_root() -> Path:
    result = subprocess.run(
        ["git", "rev-parse", "--show-toplevel"],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    return Path(result.stdout.decode("utf-8").strip())


def candidate_files(root: Path, staged: bool) -> list[Path]:
    if staged:
        names = run_git(
            root,
            ["diff", "--cached", "--name-only", "-z", "--diff-filter=ACMR"],
        )
    else:
        names = run_git(root, ["ls-files", "-z"])
    return [root / name for name in names]


def is_probably_text(path: Path, data: bytes) -> bool:
    if path.suffix.lower() in BINARY_EXTENSIONS:
        return False
    if path.suffix.lower() in TEXT_EXTENSIONS:
        return True
    return b"\0" not in data[:8192]


def line_number(text: str, offset: int) -> int:
    return text.count("\n", 0, offset) + 1


def mojibake_hits(text: str) -> list[tuple[int, str]]:
    hits: list[tuple[int, str]] = []
    for marker in MOJIBAKE_MARKERS:
        offset = text.find(marker)
        if offset >= 0:
            hits.append((line_number(text, offset), f"mojibake marker {marker!r}"))
    for pattern, label in (
        (C1_CONTROL, "C1 control character"),
        (PRIVATE_USE, "private-use character"),
        (LATIN1_MOJIBAKE, "latin-1 decoded UTF-8 mojibake"),
    ):
        match = pattern.search(text)
        if match:
            hits.append((line_number(text, match.start()), label))
    return hits


def check_file(path: Path) -> list[str]:
    if not path.is_file():
        return []
    data = path.read_bytes()
    if not is_probably_text(path, data):
        return []
    label = path.as_posix()
    if data.startswith(b"\xef\xbb\xbf"):
        return [f"{label}:1: UTF-8 BOM is not allowed"]
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError as exc:
        return [f"{label}:{exc.start + 1}: invalid UTF-8 byte sequence"]
    problems = [f"{label}:{line}: {reason}" for line, reason in mojibake_hits(text)]
    if path.suffix == ".py":
        problems.extend(check_python_text_io(path, text))
    return problems


def check_python_text_io(path: Path, text: str) -> list[str]:
    try:
        tree = ast.parse(text, filename=str(path))
    except SyntaxError:
        return []
    visitor = PythonTextIoVisitor()
    visitor.visit(tree)
    return [f"{path.as_posix()}:{line}: {reason}" for line, reason in visitor.problems]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--staged", action="store_true", help="check staged paths only")
    args = parser.parse_args(argv)

    root = repo_root()
    problems: list[str] = []
    for path in candidate_files(root, args.staged):
        problems.extend(check_file(path))

    if problems:
        print("Text encoding check failed:", file=sys.stderr)
        for problem in problems:
            print(f"  {problem}", file=sys.stderr)
        return 1
    print("Text encoding check passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
