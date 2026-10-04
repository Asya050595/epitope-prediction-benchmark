#!/usr/bin/env python3
"""Run lightweight publication checks on the repository's Python scripts."""

from __future__ import annotations

import compileall
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
CYRILLIC = re.compile(r"[\u0400-\u04FF]")
LOCAL_PATH_MARKERS = ("/home/" + "asik/", "C:\\" + "Users\\")


def main() -> int:
    failures: list[str] = []
    python_files = sorted(SCRIPTS.glob("*.py"))
    if not python_files:
        failures.append("No Python scripts were found.")

    for path in python_files:
        text = path.read_text(encoding="utf-8")
        if CYRILLIC.search(text):
            failures.append(f"Cyrillic text remains: {path.relative_to(ROOT)}")
        for marker in LOCAL_PATH_MARKERS:
            if marker in text:
                failures.append(
                    f"Machine-specific path remains ({marker!r}): "
                    f"{path.relative_to(ROOT)}"
                )

    if not compileall.compile_dir(str(SCRIPTS), quiet=1):
        failures.append("One or more Python files failed to compile.")

    if failures:
        print("Repository check failed:")
        for failure in failures:
            print(f"- {failure}")
        return 1

    print(f"Repository check passed for {len(python_files)} Python files.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
