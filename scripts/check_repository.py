#!/usr/bin/env python3
"""Run lightweight publication checks on the repository's Python scripts."""

from __future__ import annotations

import compileall
import ast
import io
import os
import re
import subprocess
import sys
import tempfile
import tokenize
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
CYRILLIC = re.compile(r"[\u0400-\u04FF]")
LOCAL_PATH_MARKERS = ("/home/" + "asik/", "C:\\" + "Users\\", "/Users/")
CORRUPTION_MARKERS = (
    "Processing details",
    "data/Processing details",
    "Запуски тулов на",
)
MISSING_INPUT_SMOKE_TESTS = (
    "benchmarking_mhc_i.py",
    "benchmarking_mhc_ii.py",
    "benchmarking_mhc_i_per_allele.py",
    "benchmarking_mhc_ii_per_allele.py",
)


def main() -> int:
    failures: list[str] = []
    python_files = sorted(SCRIPTS.glob("*.py"))
    if not python_files:
        failures.append("No Python scripts were found.")

    for path in python_files:
        text = path.read_text(encoding="utf-8")
        try:
            tokens = tokenize.generate_tokens(io.StringIO(text).readline)
            for token in tokens:
                if token.type == tokenize.COMMENT and CYRILLIC.search(token.string):
                    failures.append(
                        f"Cyrillic comment remains: {path.relative_to(ROOT)}:{token.start[0]}"
                    )
        except (IndentationError, tokenize.TokenError) as exc:
            failures.append(f"Tokenization failed for {path.relative_to(ROOT)}: {exc}")

        try:
            module = ast.parse(text)
            if (
                module.body
                and isinstance(module.body[0], ast.Expr)
                and isinstance(module.body[0].value, ast.Constant)
                and isinstance(module.body[0].value.value, str)
                and CYRILLIC.search(module.body[0].value.value)
            ):
                failures.append(f"Cyrillic module docstring remains: {path.relative_to(ROOT)}")
        except SyntaxError as exc:
            failures.append(f"Parsing failed for {path.relative_to(ROOT)}: {exc}")

        if path.resolve() != Path(__file__).resolve():
            for marker in LOCAL_PATH_MARKERS:
                if marker in text:
                    failures.append(
                        f"Machine-specific path remains ({marker!r}): "
                        f"{path.relative_to(ROOT)}"
                    )
            for marker in CORRUPTION_MARKERS:
                if marker in text:
                    failures.append(
                        f"Known conversion-corruption marker remains ({marker!r}): "
                        f"{path.relative_to(ROOT)}"
                    )

    if not compileall.compile_dir(str(SCRIPTS), quiet=1):
        failures.append("One or more Python files failed to compile.")

    for script_name in MISSING_INPUT_SMOKE_TESTS:
        with tempfile.TemporaryDirectory() as empty_data_root:
            environment = os.environ.copy()
            environment["RM_ROOT"] = empty_data_root
            result = subprocess.run(
                [sys.executable, str(SCRIPTS / script_name), "--antigen", "p24"],
                cwd=ROOT,
                env=environment,
                capture_output=True,
                text=True,
                check=False,
            )
            if result.returncode == 0:
                failures.append(f"{script_name} accepted missing required inputs.")
            expected_path = str(Path(empty_data_root) / "p24" / "p24.xlsx")
            if expected_path not in result.stdout + result.stderr:
                failures.append(
                    f"{script_name} did not report the expected portable input path."
                )
            if any(Path(empty_data_root).rglob("*.xlsx")):
                failures.append(f"{script_name} wrote an Excel file from empty inputs.")

    if failures:
        print("Repository check failed:")
        for failure in failures:
            print(f"- {failure}")
        return 1

    print(f"Repository check passed for {len(python_files)} Python files.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
