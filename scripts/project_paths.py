#!/usr/bin/env python3
"""Shared, portable project paths."""

from __future__ import annotations

import os
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = Path(os.environ.get("RM_ROOT", PROJECT_ROOT / "data")).expanduser().resolve()
