#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pairwise exact McNemar analysis for MHC epitope-prediction tools.

Inputs for each antigen
-----------------------
1. One reference Excel dataset containing experimentally validated allele-peptide
   pairs. MHC I pairs are read from columns C-D; MHC II pairs from F-G.
2. RM Excel files for each prediction tool. If a workbook contains a Union_OR
   sheet, that sheet is used; otherwise all sheets with detectable allele and
   peptide columns are combined.

UP files are intentionally not read. The analysis compares sensitivity on the
same experimentally validated positive pairs:
  1 = exact allele-peptide pair found in a tool RM file
  0 = validated pair missed by the tool

For each antigen the script creates one human-readable Excel workbook with:
  - Pairwise_results
  - Pvalue_Holm_matrix
  - Sensitivity_diff_matrix
  - Tool_overview
  - Discordant_pairs
  - Input_QC
  - README

Statistical methods
-------------------
- Two-sided exact McNemar test (exact binomial test on discordant cells b and c)
- Holm correction (the primary and only multiple-testing correction)
- Wilson 95% confidence intervals for each tool sensitivity
- Paired bootstrap 95% confidence interval for sensitivity differences

Allele handling is strict: each cell must contain exactly one HLA allele.
Values containing /, +, ;, or , are rejected and reported in Input_QC.
"""

from __future__ import annotations

from project_paths import DATA_ROOT

import argparse
import itertools
import math
import os
import re
import sys
from pathlib import Path
from statistics import NormalDist
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from openpyxl import load_workbook
from openpyxl.formatting.rule import CellIsRule, ColorScaleRule, FormulaRule
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import get_column_letter
from openpyxl.worksheet.table import Table, TableStyleInfo

try:
    from scipy.stats import binomtest
except Exception:
    binomtest = None

try:
    from allele_sets import get_allele_sets
except ImportError as exc:
    sys.exit(
        "ERROR: place allele_sets.py in the same directory as this script.\n"
        f"Details: {exc}"
    )


# =============================================================================
# User settings
# =============================================================================

DEFAULT_ROOT = f"{DATA_ROOT}"
BASE = Path(os.environ.get("RM_ROOT", DEFAULT_ROOT))
ANTIGENS = ["p24", "pp65"]

MHC_LABEL = "MHC I"
MHC_CLASS = "I"
MHC_DIR_NAME = "Matches MHC I"
OUT_DIR_NAME = "McNemar test MHC I"
OUT_FILE_CLASS = "MHC_I"

ALPHA = 0.05
BOOTSTRAP_ITERATIONS = 10_000
RANDOM_SEED = 42

# Fixed positions in the antigen reference workbooks (0-based indices).
DATASET_COLUMNS = {
    "MHC I": {"allele_idx": 2, "peptide_idx": 3, "allele_excel": "C", "peptide_excel": "D"},
    "MHC II": {"allele_idx": 5, "peptide_idx": 6, "allele_excel": "F", "peptide_excel": "G"},
}

EMPTY_VALUES = {"", "NAN", "NONE", "NA", "N/A", "-", "—", "NULL"}
UNEXPECTED_MULTI_ALLELE_PATTERN = re.compile(r"[/+;,]")

# Excel palette.
DARK_BLUE = "17365D"
MID_BLUE = "5B9BD5"
LIGHT_BLUE = "D9EAF7"
PALE_BLUE = "EAF2F8"
GREEN = "E2F0D9"
DARK_GREEN = "548235"
YELLOW = "FFF2CC"
RED = "FCE4D6"
DARK_RED = "C00000"
GREY = "E7E6E6"
WHITE = "FFFFFF"
GRID = "D0D0D0"


# =============================================================================
# Basic normalization and validation
# =============================================================================

PairKey = Tuple[str, str]


def clean_cell(value) -> str:
    if value is None or pd.isna(value):
        return ""
    text = str(value).strip()
    if text.upper() in EMPTY_VALUES:
        return ""
    return text


def normalize_peptide(value) -> str:
    """Uppercase a peptide sequence and remove whitespace."""
    text = clean_cell(value)
    if not text:
        return ""
    return re.sub(r"\s+", "", text.upper())


def normalize_allele(value) -> Tuple[str, str]:
    """
    Normalize exactly one HLA allele.

    Returns
    -------
    normalized_allele, error_code

    Examples
    --------
    HLA-A*02:01       -> A0201
    A*02:01           -> A0201
    HLA-DRB1*01:01    -> DRB10101
    DRB1_0101         -> DRB10101

    A cell containing /, +, ;, or , is rejected because the reference datasets
    are expected to contain one allele per cell, not heterodimers or lists.
    """
    text = clean_cell(value)
    if not text:
        return "", "BLANK_ALLELE"

    if UNEXPECTED_MULTI_ALLELE_PATTERN.search(text):
        return "", "MULTIPLE_ALLELES_OR_UNEXPECTED_SEPARATOR"

    normalized = text.upper()
    normalized = re.sub(r"\s+", "", normalized)
    normalized = normalized.replace("HLA-", "")
    normalized = normalized.replace("HLA_", "")
    if normalized.startswith("HLA"):
        normalized = normalized[3:]
    normalized = re.sub(r"[^A-Z0-9]", "", normalized)

    if not normalized:
        return "", "INVALID_ALLELE"

    return normalized, ""


def make_pair_key(allele, peptide) -> Tuple[Optional[PairKey], str]:
    allele_norm, allele_error = normalize_allele(allele)
    peptide_norm = normalize_peptide(peptide)

    if allele_error:
        return None, allele_error
    if not peptide_norm:
        return None, "BLANK_OR_INVALID_PEPTIDE"

    return (allele_norm, peptide_norm), ""


def allele_norm_class(allele_norm: str) -> str:
    """Diagnostic HLA-class inference; fixed reference columns remain authoritative."""
    if re.match(r"^[ABCEFG]\d", allele_norm):
        return "MHC I"
    if re.match(r"^(DR|DQ|DP)", allele_norm):
        return "MHC II"
    return "Unknown"


# =============================================================================
# Reference dataset parser
# =============================================================================


def read_dataset_pairs(path: Path) -> Tuple[Dict[PairKey, dict], List[dict]]:
    """Read unique experimentally validated pairs from fixed class-specific columns."""
    pairs: Dict[PairKey, dict] = {}
    report: List[dict] = []
    spec = DATASET_COLUMNS[MHC_LABEL]

    try:
        df = pd.read_excel(path, sheet_name=0, dtype=str)
    except Exception as exc:
        return {}, [{
            "Antigen": path.stem,
            "MHC class": MHC_LABEL,
            "Source": "Reference dataset",
            "Tool": "",
            "File": str(path),
            "Sheet": "Sheet1",
            "Status": "ERROR",
            "Allele column": "",
            "Peptide column": "",
            "Rows read": 0,
            "Pairs detected": 0,
            "Unique pairs added": 0,
            "Duplicates removed": 0,
            "Blank/invalid rows": 0,
            "Multiple alleles rejected": 0,
            "Unexpected HLA class": 0,
            "Outside reference set": 0,
            "Message": str(exc),
        }]

    allele_idx = spec["allele_idx"]
    peptide_idx = spec["peptide_idx"]

    if df.shape[1] <= max(allele_idx, peptide_idx):
        return {}, [{
            "Antigen": path.stem,
            "MHC class": MHC_LABEL,
            "Source": "Reference dataset",
            "Tool": "",
            "File": str(path),
            "Sheet": "Sheet1",
            "Status": "ERROR",
            "Allele column": "",
            "Peptide column": "",
            "Rows read": len(df),
            "Pairs detected": 0,
            "Unique pairs added": 0,
            "Duplicates removed": 0,
            "Blank/invalid rows": 0,
            "Multiple alleles rejected": 0,
            "Unexpected HLA class": 0,
            "Outside reference set": 0,
            "Message": (
                f"Dataset has {df.shape[1]} columns, but {MHC_LABEL} requires "
                f"Excel columns {spec['allele_excel']}-{spec['peptide_excel']}."
            ),
        }]

    allele_col = df.columns[allele_idx]
    peptide_col = df.columns[peptide_idx]

    detected = 0
    empty_selected_rows = 0
    blank_or_invalid = 0
    multiple_alleles = 0
    duplicates = 0
    unexpected_class = 0

    for excel_row, row in df.iterrows():
        allele_raw = clean_cell(row.iloc[allele_idx])
        peptide_raw = clean_cell(row.iloc[peptide_idx])

        # The MHC I and MHC II lists can have different lengths in the same sheet.
        # Rows where both selected cells are empty are ordinary trailing empty rows,
        # not malformed biological records.
        if not allele_raw and not peptide_raw:
            empty_selected_rows += 1
            continue

        key, error_code = make_pair_key(allele_raw, peptide_raw)

        if key is None:
            if error_code == "MULTIPLE_ALLELES_OR_UNEXPECTED_SEPARATOR":
                multiple_alleles += 1
            else:
                blank_or_invalid += 1
            continue

        detected += 1
        if allele_norm_class(key[0]) != MHC_LABEL:
            unexpected_class += 1

        if key in pairs:
            duplicates += 1
            continue

        pairs[key] = {
            "Allele": allele_raw,
            "Peptide": peptide_raw,
            "Allele_norm": key[0],
            "Peptide_norm": key[1],
            "Dataset_excel_row": int(excel_row) + 2,
        }

    report.append({
        "Antigen": path.stem,
        "MHC class": MHC_LABEL,
        "Source": "Reference dataset",
        "Tool": "",
        "File": str(path),
        "Sheet": "Sheet1",
        "Status": "OK" if pairs else "ERROR",
        "Allele column": f"{spec['allele_excel']} / {allele_col}",
        "Peptide column": f"{spec['peptide_excel']} / {peptide_col}",
        "Rows read": len(df),
        "Pairs detected": detected,
        "Unique pairs added": len(pairs),
        "Duplicates removed": duplicates,
        "Blank/invalid rows": blank_or_invalid,
        "Multiple alleles rejected": multiple_alleles,
        "Unexpected HLA class": unexpected_class,
        "Outside reference set": 0,
        "Message": f"Parsed from fixed class-specific reference columns. Empty selected rows ignored: {empty_selected_rows}.",
    })
    return pairs, report


# =============================================================================
# RM parser
# =============================================================================


def norm_col_name(column) -> str:
    return re.sub(r"[^a-z0-9]+", " ", str(column).lower()).strip()


ALLELE_EXACT_PRIORITY = [
    "dataset allele",
    "dataset hla allele",
    "true allele",
    "matched allele",
    "match allele",
    "hla allele",
    "mhc allele",
    "allele",
    "hla",
]

PEPTIDE_EXACT_PRIORITY = [
    "dataset peptide",
    "dataset epitope",
    "true peptide",
    "matched peptide",
    "match peptide",
    "peptide",
    "epitope",
    "epitope sequence",
    "peptide sequence",
    "linear sequence",
    "ctl epitopes",
    "htl epitopes",
]


def looks_like_allele_value(value: str) -> bool:
    text = clean_cell(value).upper()
    if not text or UNEXPECTED_MULTI_ALLELE_PATTERN.search(text):
        return False
    text = re.sub(r"\s+", "", text)
    patterns = [
        r"HLA[-_]?[ABCEFG]\*?\d{2}:?\d{2}",
        r"^[ABCEFG]\*?\d{2}:?\d{2}",
        r"^[ABCEFG]_?\d{4}",
        r"HLA[-_]?DRB\d?\*?\d{2}:?\d{2}",
        r"^DRB\d?\*?\d{2}:?\d{2}",
        r"^DRB\d?_?\d{4}",
        r"HLA[-_]?[D][QP][AB]\d\*?\d{2}:?\d{2}",
        r"^[D][QP][AB]\d\*?\d{2}:?\d{2}",
        r"^[D][QP][AB]\d_?\d{4}",
    ]
    return any(re.search(pattern, text) for pattern in patterns)


def looks_like_peptide_value(value: str) -> bool:
    text = normalize_peptide(value)
    if not text:
        return False
    if not re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWYX]+", text):
        return False
    return 5 <= len(text) <= 45


def score_column_by_values(df: pd.DataFrame, column, kind: str, max_rows: int = 500) -> int:
    values = df[column].dropna().astype(str).head(max_rows)
    if kind == "allele":
        return sum(looks_like_allele_value(value) for value in values)
    if kind == "peptide":
        return sum(looks_like_peptide_value(value) for value in values)
    return 0


def find_col_by_priority(
    df: pd.DataFrame,
    priorities: Sequence[str],
    kind: str,
    avoid_col=None,
):
    columns = list(df.columns)
    normalized_names = {column: norm_col_name(column) for column in columns}

    for target in priorities:
        for column in columns:
            if column == avoid_col:
                continue
            if normalized_names[column] == target:
                return column

    for column in columns:
        if column == avoid_col:
            continue
        name = normalized_names[column]
        if kind == "allele" and ("allele" in name or name == "hla" or "hla " in name or "mhc" in name):
            return column
        if kind == "peptide" and (
            ("peptide" in name or "epitope" in name or "sequence" in name)
            and "core" not in name
            and "icore" not in name
        ):
            return column

    scored = []
    for column in columns:
        if column == avoid_col:
            continue
        score = score_column_by_values(df, column, kind)
        if score > 0:
            scored.append((score, column))

    if not scored:
        return None
    scored.sort(reverse=True, key=lambda item: item[0])
    return scored[0][1]


def find_pair_columns_for_tp(df: pd.DataFrame):
    allele_col = find_col_by_priority(df, ALLELE_EXACT_PRIORITY, "allele")
    peptide_col = find_col_by_priority(df, PEPTIDE_EXACT_PRIORITY, "peptide", avoid_col=allele_col)
    return allele_col, peptide_col


def choose_rm_sheets(sheet_names: Iterable[str]) -> List[str]:
    """Prefer a precomputed union sheet; otherwise inspect all available sheets."""
    names = list(sheet_names)
    union_names = [
        name for name in names
        if norm_col_name(name) in {"union or", "union", "rm", "recovered matches", "true positives"}
    ]
    if union_names:
        # Prefer Union_OR over other possible summary labels.
        union_names.sort(key=lambda name: (norm_col_name(name) != "union or", names.index(name)))
        return [union_names[0]]
    return names


def read_rm_pairs_from_excel(path: Path, antigen: str, tool: str) -> Tuple[Dict[PairKey, dict], List[dict]]:
    """Read unique RM pairs, preferring Union_OR when present."""
    pairs: Dict[PairKey, dict] = {}
    report: List[dict] = []

    try:
        excel = pd.ExcelFile(path)
        selected_sheets = choose_rm_sheets(excel.sheet_names)
    except Exception as exc:
        return {}, [{
            "Antigen": antigen,
            "MHC class": MHC_LABEL,
            "Source": "RM file",
            "Tool": tool,
            "File": str(path),
            "Sheet": "",
            "Status": "ERROR",
            "Allele column": "",
            "Peptide column": "",
            "Rows read": 0,
            "Pairs detected": 0,
            "Unique pairs added": 0,
            "Duplicates removed": 0,
            "Blank/invalid rows": 0,
            "Multiple alleles rejected": 0,
            "Unexpected HLA class": 0,
            "Outside reference set": 0,
            "Message": str(exc),
        }]

    for sheet_name in selected_sheets:
        try:
            df = pd.read_excel(path, sheet_name=sheet_name, dtype=str)
        except Exception as exc:
            report.append({
                "Antigen": antigen,
                "MHC class": MHC_LABEL,
                "Source": "RM file",
                "Tool": tool,
                "File": str(path),
                "Sheet": sheet_name,
                "Status": "ERROR",
                "Allele column": "",
                "Peptide column": "",
                "Rows read": 0,
                "Pairs detected": 0,
                "Unique pairs added": 0,
                "Duplicates removed": 0,
                "Blank/invalid rows": 0,
                "Multiple alleles rejected": 0,
                "Unexpected HLA class": 0,
                "Outside reference set": 0,
                "Message": str(exc),
            })
            continue

        if df.empty:
            report.append({
                "Antigen": antigen,
                "MHC class": MHC_LABEL,
                "Source": "RM file",
                "Tool": tool,
                "File": str(path),
                "Sheet": sheet_name,
                "Status": "SKIPPED",
                "Allele column": "",
                "Peptide column": "",
                "Rows read": 0,
                "Pairs detected": 0,
                "Unique pairs added": 0,
                "Duplicates removed": 0,
                "Blank/invalid rows": 0,
                "Multiple alleles rejected": 0,
                "Unexpected HLA class": 0,
                "Outside reference set": 0,
                "Message": "Empty sheet.",
            })
            continue

        allele_col, peptide_col = find_pair_columns_for_tp(df)
        if allele_col is None or peptide_col is None:
            report.append({
                "Antigen": antigen,
                "MHC class": MHC_LABEL,
                "Source": "RM file",
                "Tool": tool,
                "File": str(path),
                "Sheet": sheet_name,
                "Status": "SKIPPED",
                "Allele column": allele_col or "",
                "Peptide column": peptide_col or "",
                "Rows read": len(df),
                "Pairs detected": 0,
                "Unique pairs added": 0,
                "Duplicates removed": 0,
                "Blank/invalid rows": 0,
                "Multiple alleles rejected": 0,
                "Unexpected HLA class": 0,
                "Outside reference set": 0,
                "Message": "Could not detect allele and peptide columns.",
            })
            continue

        before = len(pairs)
        detected = 0
        duplicates = 0
        blank_or_invalid = 0
        multiple_alleles = 0
        unexpected_class = 0

        for _, row in df.iterrows():
            allele_raw = clean_cell(row.get(allele_col, ""))
            peptide_raw = clean_cell(row.get(peptide_col, ""))
            key, error_code = make_pair_key(allele_raw, peptide_raw)

            if key is None:
                if error_code == "MULTIPLE_ALLELES_OR_UNEXPECTED_SEPARATOR":
                    multiple_alleles += 1
                else:
                    blank_or_invalid += 1
                continue

            detected += 1
            if allele_norm_class(key[0]) != MHC_LABEL:
                unexpected_class += 1

            if key in pairs:
                duplicates += 1
                continue

            pairs[key] = {
                "Allele": allele_raw,
                "Peptide": peptide_raw,
                "Allele_norm": key[0],
                "Peptide_norm": key[1],
            }

        report.append({
            "Antigen": antigen,
            "MHC class": MHC_LABEL,
            "Source": "RM file",
            "Tool": tool,
            "File": str(path),
            "Sheet": sheet_name,
            "Status": "OK",
            "Allele column": str(allele_col),
            "Peptide column": str(peptide_col),
            "Rows read": len(df),
            "Pairs detected": detected,
            "Unique pairs added": len(pairs) - before,
            "Duplicates removed": duplicates,
            "Blank/invalid rows": blank_or_invalid,
            "Multiple alleles rejected": multiple_alleles,
            "Unexpected HLA class": unexpected_class,
            "Outside reference set": 0,
            "Message": "Preferred Union_OR/RM sheet used." if len(selected_sheets) == 1 else "All detectable sheets combined.",
        })

    return pairs, report


# =============================================================================
# File discovery and tool names
# =============================================================================


def strip_copy_suffix(stem: str) -> str:
    return re.sub(r"\s*\(\d+\)$", "", stem).strip()


def friendly_tool_name(raw_name: str) -> str:
    name = raw_name.replace("_", " ")
    name = re.sub(r"\s+", " ", name).strip()
    return name


def infer_tool_name(rm_file: Path, antigen: str) -> str:
    """Use the filename first so separate IEDB modes are not accidentally merged."""
    stem = strip_copy_suffix(rm_file.stem)
    stem = re.sub(rf"^RM_{re.escape(antigen)}_", "", stem, flags=re.IGNORECASE)
    if stem and stem != rm_file.stem:
        return friendly_tool_name(stem)

    parent = rm_file.parent.name
    if parent.lower().startswith("matches "):
        return friendly_tool_name(parent[len("Matches "):])
    return friendly_tool_name(stem.replace("RM_", ""))


def discover_rm_files(antigen: str, base: Path) -> Dict[str, List[Path]]:
    root = base / f"Processing details{antigen}" / MHC_DIR_NAME
    tool_to_files: Dict[str, List[Path]] = {}
    if not root.exists():
        return tool_to_files

    pattern = f"RM_{antigen}_*.xlsx"
    for path in sorted(root.rglob(pattern)):
        if path.name.startswith("~$"):
            continue
        tool = infer_tool_name(path, antigen)
        tool_key = re.sub(r"[^a-z0-9]+", "", tool.lower())
        if "mhcpred" in tool_key or "rankpep" in tool_key:
            continue
        tool_to_files.setdefault(tool, []).append(path)
    return tool_to_files


# =============================================================================
# Statistics
# =============================================================================


def exact_binomial_two_sided_p(k: int, n: int) -> float:
    """Two-sided exact binomial p-value for p=0.5; fallback when SciPy is absent."""
    if n == 0:
        return 1.0
    k = min(k, n - k)
    log_half = math.log(0.5)
    probabilities = []
    for i in range(k + 1):
        log_p = (
            math.lgamma(n + 1)
            - math.lgamma(i + 1)
            - math.lgamma(n - i + 1)
            + n * log_half
        )
        probabilities.append(math.exp(log_p))
    return min(1.0, 2.0 * sum(probabilities))


def mcnemar_exact_p_value(only_tool_1: int, only_tool_2: int) -> float:
    b = int(only_tool_1)
    c = int(only_tool_2)
    n = b + c
    if n == 0:
        return 1.0
    if binomtest is not None:
        return float(binomtest(b, n=n, p=0.5, alternative="two-sided").pvalue)
    return exact_binomial_two_sided_p(b, n)


def adjust_pvalues_holm(pvalues: Sequence[float]) -> List[float]:
    """Holm step-down FWER adjustment."""
    m = len(pvalues)
    if m == 0:
        return []
    order = sorted(range(m), key=lambda index: pvalues[index])
    adjusted_sorted = [0.0] * m
    running_max = 0.0
    for rank, original_index in enumerate(order):
        candidate = (m - rank) * float(pvalues[original_index])
        running_max = max(running_max, candidate)
        adjusted_sorted[rank] = min(1.0, running_max)

    adjusted = [0.0] * m
    for rank, original_index in enumerate(order):
        adjusted[original_index] = adjusted_sorted[rank]
    return adjusted


def wilson_ci(successes: int, total: int, confidence: float = 0.95) -> Tuple[float, float]:
    """Wilson score interval returned as percentages."""
    if total <= 0:
        return float("nan"), float("nan")
    z = NormalDist().inv_cdf(1 - (1 - confidence) / 2)
    p = successes / total
    denominator = 1 + z * z / total
    center = (p + z * z / (2 * total)) / denominator
    half_width = (
        z
        * math.sqrt(p * (1 - p) / total + z * z / (4 * total * total))
        / denominator
    )
    return max(0.0, center - half_width) * 100, min(1.0, center + half_width) * 100


def paired_bootstrap_difference_ci(
    both_found: int,
    only_tool_1: int,
    only_tool_2: int,
    neither: int,
    iterations: int,
    seed: int,
) -> Tuple[float, float]:
    """Percentile bootstrap CI for paired sensitivity difference in percentage points."""
    counts = np.array([both_found, only_tool_1, only_tool_2, neither], dtype=int)
    total = int(counts.sum())
    if total <= 0:
        return float("nan"), float("nan")
    probabilities = counts / total
    rng = np.random.default_rng(seed)
    sampled = rng.multinomial(total, probabilities, size=iterations)
    differences = (sampled[:, 1] - sampled[:, 2]) / total * 100.0
    low, high = np.quantile(differences, [0.025, 0.975])
    return float(low), float(high)


def format_ci(low: float, high: float, digits: int = 1, signed: bool = False) -> str:
    if math.isnan(low) or math.isnan(high):
        return "NA"
    if signed:
        return f"{low:+.{digits}f} to {high:+.{digits}f}"
    return f"{low:.{digits}f} to {high:.{digits}f}"


def interpret_result(row: dict) -> str:
    difference = row["Sensitivity difference (pp)"]
    if difference > 0:
        higher_tool = row["Tool 1"]
    elif difference < 0:
        higher_tool = row["Tool 2"]
    else:
        higher_tool = "Neither tool"

    if bool(row["Significant after Holm"]):
        return f"Statistically significant after Holm correction; {higher_tool} had higher sensitivity."
    return "Not statistically significant after Holm correction."


# =============================================================================
# Output tables
# =============================================================================


def build_discordant_rows(
    antigen: str,
    comparison_id: str,
    tool_1: str,
    tool_2: str,
    only_tool_1: Iterable[PairKey],
    only_tool_2: Iterable[PairKey],
    dataset_pairs: Dict[PairKey, dict],
) -> List[dict]:
    rows: List[dict] = []
    categories = [
        (only_tool_1, tool_1, tool_2, "Found only by Tool 1"),
        (only_tool_2, tool_2, tool_1, "Found only by Tool 2"),
    ]
    for pair_keys, found_by, missed_by, category in categories:
        for key in sorted(pair_keys, key=lambda item: (dataset_pairs[item]["Allele"], dataset_pairs[item]["Peptide"])):
            info = dataset_pairs[key]
            rows.append({
                "Antigen": antigen,
                "MHC class": MHC_LABEL,
                "Comparison ID": comparison_id,
                "Tool 1": tool_1,
                "Tool 2": tool_2,
                "Category": category,
                "Found by": found_by,
                "Missed by": missed_by,
                "Allele": info["Allele"],
                "Peptide": info["Peptide"],
                "Reference row": info.get("Dataset_excel_row", ""),
                "Allele normalized": info["Allele_norm"],
                "Peptide normalized": info["Peptide_norm"],
            })
    return rows


def make_matrix(tools: Sequence[str], rows: Sequence[dict], value_column: str, antisymmetric: bool = False) -> pd.DataFrame:
    matrix = pd.DataFrame(np.nan, index=tools, columns=tools, dtype=float)
    for tool in tools:
        if antisymmetric:
            matrix.loc[tool, tool] = 0.0
    for row in rows:
        tool_1 = row["Tool 1"]
        tool_2 = row["Tool 2"]
        value = float(row[value_column])
        matrix.loc[tool_1, tool_2] = value
        matrix.loc[tool_2, tool_1] = -value if antisymmetric else value
    matrix.index.name = "Tool"
    return matrix


# =============================================================================
# Excel formatting
# =============================================================================


def thin_border() -> Border:
    side = Side(style="thin", color=GRID)
    return Border(left=side, right=side, top=side, bottom=side)


def style_title(ws, title: str, subtitle: str, last_column: int):
    end_letter = get_column_letter(max(1, last_column))

    # Style the complete ranges before merging so Excel renderers show the fill
    # consistently across the full title width.
    for column in range(1, last_column + 1):
        title_cell = ws.cell(1, column)
        title_cell.fill = PatternFill("solid", fgColor=DARK_BLUE)
        title_cell.font = Font(bold=True, size=16, color=WHITE)
        title_cell.alignment = Alignment(horizontal="left", vertical="center")

        subtitle_cell = ws.cell(2, column)
        subtitle_cell.fill = PatternFill("solid", fgColor=PALE_BLUE)
        subtitle_cell.font = Font(italic=True, size=10, color="404040")
        subtitle_cell.alignment = Alignment(horizontal="left", vertical="center", wrap_text=True)

    ws.merge_cells(f"A1:{end_letter}1")
    ws["A1"] = title
    ws.row_dimensions[1].height = 26

    ws.merge_cells(f"A2:{end_letter}2")
    ws["A2"] = subtitle
    ws.row_dimensions[2].height = 32


def apply_header_style(ws, header_row: int, max_column: int):
    for cell in ws.iter_cols(min_col=1, max_col=max_column, min_row=header_row, max_row=header_row):
        item = cell[0]
        item.fill = PatternFill("solid", fgColor=MID_BLUE)
        item.font = Font(bold=True, color=WHITE)
        item.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
        item.border = thin_border()
    ws.row_dimensions[header_row].height = 42


def add_excel_table(ws, header_row: int, last_row: int, last_col: int, table_name: str):
    if last_row <= header_row:
        return
    reference = f"A{header_row}:{get_column_letter(last_col)}{last_row}"
    table = Table(displayName=table_name, ref=reference)
    table.tableStyleInfo = TableStyleInfo(
        name="TableStyleMedium2",
        showFirstColumn=False,
        showLastColumn=False,
        showRowStripes=True,
        showColumnStripes=False,
    )
    ws.add_table(table)


def set_column_widths(ws, widths: Dict[str, float]):
    for column, width in widths.items():
        ws.column_dimensions[column].width = width


def format_pvalue_cell(cell):
    cell.number_format = "0.0000E+00"
    cell.alignment = Alignment(horizontal="center")


def style_pairwise_sheet(ws, row_count: int, col_count: int):
    style_title(
        ws,
        f"Pairwise exact McNemar results — {MHC_LABEL}",
        "Holm-adjusted p-value is the primary confirmatory result. Raw exact p-values and diagonal counts are retained in hidden technical columns.",
        col_count,
    )

    # Group header row 4; DataFrame headers are in row 5.
    groups = [
        (1, 4, "Comparison"),
        (5, 8, "Tool 1"),
        (9, 12, "Tool 2"),
        (13, 14, "Effect size"),
        (15, 17, "Discordant results used by McNemar"),
        (18, 19, "Primary multiple-testing correction"),
        (20, 20, "Interpretation"),
        (21, 23, "Technical details"),
    ]
    for start_col, end_col, label in groups:
        start_letter = get_column_letter(start_col)
        end_letter = get_column_letter(end_col)
        if start_col != end_col:
            ws.merge_cells(f"{start_letter}4:{end_letter}4")
        cell = ws.cell(4, start_col, label)
        cell.fill = PatternFill("solid", fgColor=DARK_BLUE)
        cell.font = Font(bold=True, color=WHITE)
        cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
        cell.border = thin_border()
    ws.row_dimensions[4].height = 28

    apply_header_style(ws, 5, col_count)
    last_row = 5 + row_count
    add_excel_table(ws, 5, last_row, col_count, "PairwiseResultsTable")
    ws.freeze_panes = "D6"
    ws.auto_filter.ref = f"A5:{get_column_letter(col_count)}{last_row}"

    widths = {
        "A": 14, "B": 26, "C": 26, "D": 18,
        "E": 15, "F": 15, "G": 18, "H": 20,
        "I": 15, "J": 15, "K": 18, "L": 20,
        "M": 20, "N": 22,
        "O": 18, "P": 18, "Q": 18,
        "R": 18, "S": 18, "T": 58,
        "U": 18, "V": 16, "W": 16,
    }
    set_column_widths(ws, widths)

    for row in range(6, last_row + 1):
        for col in range(1, col_count + 1):
            ws.cell(row, col).alignment = Alignment(vertical="top", wrap_text=True)
        for col in [7, 11, 13]:
            ws.cell(row, col).number_format = "0.0"
        for col in [18, 21]:
            format_pvalue_cell(ws.cell(row, col))

    if last_row >= 6:
        ws.conditional_formatting.add(
            f"R6:R{last_row}",
            CellIsRule(operator="lessThan", formula=[str(ALPHA)], fill=PatternFill("solid", fgColor=GREEN)),
        )
        ws.conditional_formatting.add(
            f"S6:S{last_row}",
            FormulaRule(formula=['$S6="Yes"'], fill=PatternFill("solid", fgColor=GREEN), font=Font(color=DARK_GREEN, bold=True)),
        )
        ws.conditional_formatting.add(
            f"S6:S{last_row}",
            FormulaRule(formula=['$S6="No"'], fill=PatternFill("solid", fgColor=GREY)),
        )
        ws.conditional_formatting.add(
            f"M6:M{last_row}",
            ColorScaleRule(
                start_type="min", start_color="F8696B",
                mid_type="num", mid_value=0, mid_color="FFFFFF",
                end_type="max", end_color="63BE7B",
            ),
        )

    # Hide technical raw p-value and agreement cells, but keep them reproducibly available.
    for letter in ["U", "V", "W"]:
        ws.column_dimensions[letter].hidden = True


def style_matrix_sheet(ws, title: str, subtitle: str, size: int, kind: str):
    last_col = size + 1
    style_title(ws, title, subtitle, last_col)
    header_row = 4
    apply_header_style(ws, header_row, last_col)
    ws.freeze_panes = "B5"
    ws.column_dimensions["A"].width = 30
    for col in range(2, last_col + 1):
        ws.column_dimensions[get_column_letter(col)].width = 15

    last_row = header_row + size
    data_range = f"B5:{get_column_letter(last_col)}{last_row}"
    for row in ws.iter_rows(min_row=5, max_row=last_row, min_col=2, max_col=last_col):
        for cell in row:
            cell.alignment = Alignment(horizontal="center", vertical="center")
            cell.border = thin_border()
            if kind == "pvalue" and cell.value is not None:
                cell.number_format = "0.0000E+00"
            elif kind == "difference" and cell.value is not None:
                cell.number_format = "+0.0;-0.0;0.0"

    if size > 0:
        if kind == "pvalue":
            ws.conditional_formatting.add(
                data_range,
                ColorScaleRule(
                    start_type="num", start_value=0, start_color="63BE7B",
                    mid_type="num", mid_value=0.05, mid_color="FFEB84",
                    end_type="num", end_value=1, end_color="F8696B",
                ),
            )
            ws.conditional_formatting.add(
                data_range,
                CellIsRule(operator="lessThan", formula=[str(ALPHA)], font=Font(bold=True, color="006100")),
            )
        else:
            ws.conditional_formatting.add(
                data_range,
                ColorScaleRule(
                    start_type="min", start_color="F8696B",
                    mid_type="num", mid_value=0, mid_color="FFFFFF",
                    end_type="max", end_color="63BE7B",
                ),
            )


def style_overview_sheet(ws, row_count: int, col_count: int):
    style_title(
        ws,
        f"Tool sensitivity overview — {MHC_LABEL}",
        "Sensitivity is the proportion of experimentally validated allele-peptide pairs detected by each tool. Confidence intervals are Wilson 95% intervals.",
        col_count,
    )
    apply_header_style(ws, 4, col_count)
    last_row = 4 + row_count
    add_excel_table(ws, 4, last_row, col_count, "ToolOverviewTable")
    ws.freeze_panes = "A5"
    set_column_widths(ws, {
        "A": 30, "B": 22, "C": 15, "D": 15, "E": 18,
        "F": 18, "G": 18, "H": 12, "I": 42, "J": 18,
    })
    for row in range(5, last_row + 1):
        ws.cell(row, 5).number_format = "0.0"
        ws.cell(row, 6).number_format = "0.0"
        ws.cell(row, 7).number_format = "0.0"
        for col in range(1, col_count + 1):
            ws.cell(row, col).alignment = Alignment(vertical="top", wrap_text=True)
    if last_row >= 5:
        ws.conditional_formatting.add(
            f"E5:E{last_row}",
            ColorScaleRule(
                start_type="min", start_color="F8696B",
                mid_type="percentile", mid_value=50, mid_color="FFEB84",
                end_type="max", end_color="63BE7B",
            ),
        )
    # Technical provenance columns.
    for letter in ["I", "J"]:
        ws.column_dimensions[letter].hidden = True


def style_standard_table_sheet(
    ws,
    title: str,
    subtitle: str,
    header_row: int,
    row_count: int,
    col_count: int,
    table_name: str,
    widths: Optional[Dict[str, float]] = None,
):
    style_title(ws, title, subtitle, col_count)
    apply_header_style(ws, header_row, col_count)
    last_row = header_row + row_count
    add_excel_table(ws, header_row, last_row, col_count, table_name)
    ws.freeze_panes = f"A{header_row + 1}"
    if widths:
        set_column_widths(ws, widths)
    for row in ws.iter_rows(min_row=header_row + 1, max_row=last_row, min_col=1, max_col=col_count):
        for cell in row:
            cell.alignment = Alignment(vertical="top", wrap_text=True)


def style_qc_sheet(ws, row_count: int, col_count: int):
    style_standard_table_sheet(
        ws,
        f"Input quality control — {MHC_LABEL}",
        "Check ERROR, SKIPPED, multiple-allele, unexpected-class and outside-reference counts before interpreting the statistical results.",
        4,
        row_count,
        col_count,
        "InputQCTable",
        {
            "A": 14, "B": 12, "C": 20, "D": 30, "E": 55, "F": 20,
            "G": 12, "H": 22, "I": 22, "J": 12, "K": 15, "L": 18,
            "M": 18, "N": 18, "O": 24, "P": 20, "Q": 22, "R": 55,
        },
    )
    last_row = 4 + row_count
    if last_row >= 5:
        ws.conditional_formatting.add(
            f"G5:G{last_row}",
            FormulaRule(formula=['$G5="ERROR"'], fill=PatternFill("solid", fgColor=RED), font=Font(color=DARK_RED, bold=True)),
        )
        ws.conditional_formatting.add(
            f"G5:G{last_row}",
            FormulaRule(formula=['$G5="SKIPPED"'], fill=PatternFill("solid", fgColor=YELLOW)),
        )
        ws.conditional_formatting.add(
            f"G5:G{last_row}",
            FormulaRule(formula=['$G5="OK"'], fill=PatternFill("solid", fgColor=GREEN)),
        )


def style_readme_sheet(ws):
    ws.column_dimensions["A"].width = 34
    ws.column_dimensions["B"].width = 110
    ws.freeze_panes = "A2"
    for cell in ws[1]:
        cell.fill = PatternFill("solid", fgColor=DARK_BLUE)
        cell.font = Font(bold=True, color=WHITE)
        cell.alignment = Alignment(horizontal="center", vertical="center")
    for row in ws.iter_rows(min_row=2):
        for cell in row:
            cell.alignment = Alignment(vertical="top", wrap_text=True)
            cell.border = thin_border()


def style_workbook(path: Path, sheet_meta: dict):
    book = load_workbook(path)

    pairwise = book["Pairwise_results"]
    style_pairwise_sheet(pairwise, sheet_meta["pairwise_rows"], sheet_meta["pairwise_cols"])

    style_matrix_sheet(
        book["Pvalue_Holm_matrix"],
        f"Holm-adjusted exact McNemar p-values — {MHC_LABEL}",
        "Primary confirmatory correction. Values below 0.05 are highlighted.",
        sheet_meta["tool_count"],
        "pvalue",
    )
    style_matrix_sheet(
        book["Sensitivity_diff_matrix"],
        f"Pairwise sensitivity differences — {MHC_LABEL}",
        "Each cell is row-tool sensitivity minus column-tool sensitivity, in percentage points.",
        sheet_meta["tool_count"],
        "difference",
    )

    style_overview_sheet(book["Tool_overview"], sheet_meta["overview_rows"], sheet_meta["overview_cols"])

    discordant = book["Discordant_pairs"]
    style_standard_table_sheet(
        discordant,
        f"Discordant validated pairs — {MHC_LABEL}",
        "Only pairs detected by one tool and missed by the other are listed; these are the b and c cells used by the exact McNemar test.",
        4,
        sheet_meta["discordant_rows"],
        sheet_meta["discordant_cols"],
        "DiscordantPairsTable",
        {
            "A": 14, "B": 12, "C": 15, "D": 28, "E": 28, "F": 24,
            "G": 28, "H": 28, "I": 20, "J": 28, "K": 15, "L": 22, "M": 24,
        },
    )
    discordant.column_dimensions["L"].hidden = True
    discordant.column_dimensions["M"].hidden = True

    style_qc_sheet(book["Input_QC"], sheet_meta["qc_rows"], sheet_meta["qc_cols"])
    style_readme_sheet(book["README"])

    # Apply a light border to matrix row/column headers and keep all sheet zooms readable.
    for ws in book.worksheets:
        ws.sheet_view.zoomScale = 85

    book.save(path)


# =============================================================================
# Main analysis
# =============================================================================


def analyze_antigen(
    antigen: str,
    base: Path,
    bootstrap_iterations: int,
    random_seed: int,
) -> Optional[Path]:
    antigen_dir = base / f"Processing details{antigen}"
    dataset_path = antigen_dir / f"{antigen}.xlsx"
    out_dir = antigen_dir / OUT_DIR_NAME
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"McNemar_{antigen}_{OUT_FILE_CLASS}.xlsx"

    print("=" * 88)
    print(f"Antigen: {antigen} | {MHC_LABEL}")
    print("=" * 88)

    if not dataset_path.exists():
        print(f"[ERROR] Reference dataset not found: {dataset_path}")
        return None

    all_dataset_pairs, dataset_report = read_dataset_pairs(dataset_path)
    evaluable_alleles, _ = get_allele_sets(MHC_CLASS, antigen)
    evaluable_by_key = {
        normalize_allele(allele)[0]: allele for allele in evaluable_alleles
    }
    dataset_pairs = {
        key: info for key, info in all_dataset_pairs.items()
        if key[0] in evaluable_by_key
    }
    excluded_non_evaluable = len(all_dataset_pairs) - len(dataset_pairs)
    present_allele_keys = {allele for allele, _ in dataset_pairs}
    missing_evaluable = {
        canonical for key, canonical in evaluable_by_key.items()
        if key not in present_allele_keys
    }
    dataset_report.append({
        "Antigen": antigen,
        "MHC class": MHC_LABEL,
        "Source": "Evaluable-set filter",
        "Tool": "",
        "File": str(dataset_path),
        "Sheet": "Sheet1",
        "Status": "OK" if dataset_pairs else "ERROR",
        "Allele column": DATASET_COLUMNS[MHC_LABEL]["allele_excel"],
        "Peptide column": DATASET_COLUMNS[MHC_LABEL]["peptide_excel"],
        "Rows read": len(all_dataset_pairs),
        "Pairs detected": len(all_dataset_pairs),
        "Unique pairs added": len(dataset_pairs),
        "Duplicates removed": 0,
        "Blank/invalid rows": 0,
        "Multiple alleles rejected": 0,
        "Unexpected HLA class": 0,
        "Outside reference set": excluded_non_evaluable,
        "Message": (
            f"Only pairs whose alleles belong to the evaluable set were retained; "
            f"{excluded_non_evaluable} pair(s) were excluded. "
            f"Evaluable alleles absent from the experimental dataset: {sorted(missing_evaluable)}."
        ),
    })
    true_set = set(dataset_pairs)
    if not true_set:
        print(f"[ERROR] No validated {MHC_LABEL} pairs were detected in {dataset_path}")
        return None

    print(
        f"Evaluable set: {len(evaluable_alleles)} alleles; "
        f"validated pairs retained: {len(true_set)}; "
        f"non-evaluable pairs excluded: {excluded_non_evaluable}"
    )
    if missing_evaluable:
        print(
            f"[WARNING] {len(missing_evaluable)} evaluable allele(s) are absent "
            f"from the experimental dataset: {sorted(missing_evaluable)}"
        )

    tool_to_files = discover_rm_files(antigen, base)
    if not tool_to_files:
        print(f"[ERROR] No RM files found under {antigen_dir / MHC_DIR_NAME}")
        return None

    tool_pairs: Dict[str, set] = {}
    qc_rows: List[dict] = list(dataset_report)
    overview_rows: List[dict] = []

    for tool, files in sorted(tool_to_files.items()):
        merged_pairs: Dict[PairKey, dict] = {}
        tool_reports: List[dict] = []

        for file in files:
            pairs, report = read_rm_pairs_from_excel(file, antigen, tool)
            tool_reports.extend(report)
            qc_rows.extend(report)
            for key, info in pairs.items():
                merged_pairs.setdefault(key, info)

        all_rm_keys = set(merged_pairs)
        rm_in_reference = all_rm_keys & true_set
        outside_reference = all_rm_keys - true_set
        tool_pairs[tool] = rm_in_reference

        # Add a tool-level QC row because outside-reference pairs are computed only
        # after all files belonging to the same tool have been merged.
        qc_rows.append({
            "Antigen": antigen,
            "MHC class": MHC_LABEL,
            "Source": "Tool summary",
            "Tool": tool,
            "File": "\n".join(str(file) for file in files),
            "Sheet": "Merged tool RM set",
            "Status": "OK" if any(row.get("Status") == "OK" for row in tool_reports) else "ERROR",
            "Allele column": "",
            "Peptide column": "",
            "Rows read": sum(int(row.get("Rows read", 0) or 0) for row in tool_reports),
            "Pairs detected": len(all_rm_keys),
            "Unique pairs added": len(rm_in_reference),
            "Duplicates removed": 0,
            "Blank/invalid rows": sum(int(row.get("Blank/invalid rows", 0) or 0) for row in tool_reports),
            "Multiple alleles rejected": sum(int(row.get("Multiple alleles rejected", 0) or 0) for row in tool_reports),
            "Unexpected HLA class": sum(int(row.get("Unexpected HLA class", 0) or 0) for row in tool_reports),
            "Outside reference set": len(outside_reference),
            "Message": "Pairs outside the evaluable validated reference set were excluded; UP files were not read.",
        })

        found = len(rm_in_reference)
        missed = len(true_set) - found
        sensitivity = found / len(true_set) * 100
        ci_low, ci_high = wilson_ci(found, len(true_set))
        overview_rows.append({
            "Tool": tool,
            "Validated pairs (N)": len(true_set),
            "Recovered (RM)": found,
            "Unrecovered (UP)": missed,
            "Sensitivity (%)": sensitivity,
            "95% CI lower (%)": ci_low,
            "95% CI upper (%)": ci_high,
            "Sensitivity (95% CI)": f"{sensitivity:.1f}% ({ci_low:.1f} to {ci_high:.1f})",
            "RM input files": "\n".join(str(file) for file in files),
            "Outside reference set (excluded)": len(outside_reference),
        })
        print(f"  {tool}: RM={found}, UP={missed}, sensitivity={sensitivity:.1f}%")

    tools = sorted(tool_pairs)
    if len(tools) < 2:
        print("[ERROR] At least two tools are required for pairwise McNemar tests.")
        return None

    pairwise_rows: List[dict] = []
    discordant_rows: List[dict] = []

    overview_lookup = {row["Tool"]: row for row in overview_rows}

    for comparison_index, (tool_1, tool_2) in enumerate(itertools.combinations(tools, 2), start=1):
        comparison_id = f"C{comparison_index:03d}"
        set_1 = tool_pairs[tool_1]
        set_2 = tool_pairs[tool_2]

        both_found = set_1 & set_2
        only_tool_1 = set_1 - set_2
        only_tool_2 = set_2 - set_1
        neither = true_set - (set_1 | set_2)

        a = len(both_found)
        b = len(only_tool_1)
        c = len(only_tool_2)
        d = len(neither)

        sensitivity_1 = overview_lookup[tool_1]["Sensitivity (%)"]
        sensitivity_2 = overview_lookup[tool_2]["Sensitivity (%)"]
        difference = sensitivity_1 - sensitivity_2
        difference_low, difference_high = paired_bootstrap_difference_ci(
            a, b, c, d,
            iterations=bootstrap_iterations,
            seed=random_seed + comparison_index,
        )

        pairwise_rows.append({
            "Comparison ID": comparison_id,
            "Tool 1": tool_1,
            "Tool 2": tool_2,
            "Validated pairs (N)": len(true_set),
            "Tool 1 recovered (RM)": len(set_1),
            "Tool 1 unrecovered (UP)": len(true_set) - len(set_1),
            "Tool 1 sensitivity (%)": sensitivity_1,
            "Tool 1 95% CI": format_ci(
                overview_lookup[tool_1]["95% CI lower (%)"],
                overview_lookup[tool_1]["95% CI upper (%)"],
            ),
            "Tool 2 recovered (RM)": len(set_2),
            "Tool 2 unrecovered (UP)": len(true_set) - len(set_2),
            "Tool 2 sensitivity (%)": sensitivity_2,
            "Tool 2 95% CI": format_ci(
                overview_lookup[tool_2]["95% CI lower (%)"],
                overview_lookup[tool_2]["95% CI upper (%)"],
            ),
            "Sensitivity difference (pp)": difference,
            "Difference 95% CI": format_ci(difference_low, difference_high, signed=True),
            "Found only by Tool 1 (b)": b,
            "Found only by Tool 2 (c)": c,
            "Discordant pairs (b+c)": b + c,
            "Holm-adjusted p": np.nan,
            "Significant after Holm": False,
            "Interpretation": "",
            "Raw exact McNemar p": mcnemar_exact_p_value(b, c),
            "Both found (a)": a,
            "Neither found (d)": d,
        })

        discordant_rows.extend(build_discordant_rows(
            antigen,
            comparison_id,
            tool_1,
            tool_2,
            only_tool_1,
            only_tool_2,
            dataset_pairs,
        ))

    raw_pvalues = [float(row["Raw exact McNemar p"]) for row in pairwise_rows]
    holm_pvalues = adjust_pvalues_holm(raw_pvalues)

    for row, holm_p in zip(pairwise_rows, holm_pvalues):
        row["Holm-adjusted p"] = holm_p
        row["Significant after Holm"] = holm_p < ALPHA
        row["Interpretation"] = interpret_result(row)
        # User-facing Yes/No values are easier to scan in Excel.
        row["Significant after Holm"] = "Yes" if row["Significant after Holm"] else "No"

    # Restore bools temporarily for the interpretation helper was already applied.
    pairwise_df = pd.DataFrame(pairwise_rows)
    overview_df = pd.DataFrame(overview_rows).sort_values(
        ["Sensitivity (%)", "Tool"], ascending=[False, True]
    ).reset_index(drop=True)
    discordant_df = pd.DataFrame(discordant_rows)
    qc_df = pd.DataFrame(qc_rows)

    holm_matrix = make_matrix(tools, pairwise_rows, "Holm-adjusted p")
    difference_matrix = make_matrix(tools, pairwise_rows, "Sensitivity difference (pp)", antisymmetric=True)

    readme_rows = [
        ("Purpose", "Pairwise comparison of prediction-tool sensitivity on experimentally validated allele-peptide pairs restricted to the evaluable allele set."),
        ("Required inputs", "One class-specific reference positive dataset plus RM files for each tool. UP files are not read."),
        ("Reference columns", f"{MHC_LABEL}: Excel columns {DATASET_COLUMNS[MHC_LABEL]['allele_excel']}-{DATASET_COLUMNS[MHC_LABEL]['peptide_excel']} in Sheet1."),
        ("Reference set", "Only allele-peptide pairs whose alleles belong to the antigen- and HLA-class-specific evaluable set from allele_sets.py are analysed."),
        ("Allele rule", "Exactly one HLA allele per cell. Values containing /, +, ;, or , are rejected and reported in Input_QC; heterodimer normalization is not performed."),
        ("RM sheet rule", "Union_OR is preferred when present; otherwise an RM/Recovered Matches sheet or all sheets with detectable allele and peptide columns are combined."),
        ("Sensitivity", "Recovered matches (RM) / evaluable validated pairs N. Unrecovered pairs (UP) = N - RM."),
        ("Sensitivity 95% CI", "Wilson score interval."),
        ("Difference", "Tool 1 sensitivity minus Tool 2 sensitivity, expressed in percentage points."),
        ("Difference 95% CI", f"Paired percentile bootstrap with {bootstrap_iterations:,} resamples and deterministic seeds derived from {random_seed}."),
        ("Exact McNemar test", "Two-sided exact binomial test on discordant cells b (Tool 1 only) and c (Tool 2 only)."),
        ("Holm adjustment", "Primary and only multiple-testing correction within one antigen and MHC class; Holm-adjusted p-values determine statistical significance."),
        ("Raw p-value", "Retained in a hidden technical column for reproducibility; significance is interpreted using the Holm-adjusted p-value."),
        ("P-value matrix", "Symmetric matrix of Holm-adjusted p-values. Diagonal cells are blank."),
        ("Sensitivity difference matrix", "Cell = row-tool sensitivity minus column-tool sensitivity in percentage points."),
        ("Discordant_pairs", "Contains only validated pairs found by one tool and missed by the other; these pairs determine the McNemar test."),
        ("Input_QC", "Review errors, skipped sheets, rejected multi-allele cells, unexpected HLA classes, and pairs excluded because they were outside the reference positive set."),
    ]
    readme_df = pd.DataFrame(readme_rows, columns=["Item", "Explanation"])

    with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
        pairwise_df.to_excel(writer, sheet_name="Pairwise_results", index=False, startrow=4)
        holm_matrix.to_excel(writer, sheet_name="Pvalue_Holm_matrix", startrow=3)
        difference_matrix.to_excel(writer, sheet_name="Sensitivity_diff_matrix", startrow=3)
        overview_df.to_excel(writer, sheet_name="Tool_overview", index=False, startrow=3)
        discordant_df.to_excel(writer, sheet_name="Discordant_pairs", index=False, startrow=3)
        qc_df.to_excel(writer, sheet_name="Input_QC", index=False, startrow=3)
        readme_df.to_excel(writer, sheet_name="README", index=False)

    style_workbook(out_path, {
        "pairwise_rows": len(pairwise_df),
        "pairwise_cols": len(pairwise_df.columns),
        "tool_count": len(tools),
        "overview_rows": len(overview_df),
        "overview_cols": len(overview_df.columns),
        "discordant_rows": len(discordant_df),
        "discordant_cols": len(discordant_df.columns),
        "qc_rows": len(qc_df),
        "qc_cols": len(qc_df.columns),
    })

    print(f"Saved: {out_path}")
    print()
    return out_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=f"Exact McNemar sensitivity analysis for {MHC_LABEL} epitope prediction tools."
    )
    parser.add_argument(
        "--base",
        type=Path,
        default=BASE,
        help=f"Base Control directory. Default: {BASE}",
    )
    parser.add_argument(
        "--antigens",
        nargs="+",
        default=ANTIGENS,
        choices=ANTIGENS,
        help="Antigens to analyse. Default: p24 pp65",
    )
    parser.add_argument(
        "--bootstrap-iterations",
        type=int,
        default=BOOTSTRAP_ITERATIONS,
        help=f"Paired bootstrap resamples. Default: {BOOTSTRAP_ITERATIONS}",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=RANDOM_SEED,
        help=f"Base random seed. Default: {RANDOM_SEED}",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    generated = []
    for antigen in args.antigens:
        result = analyze_antigen(
            antigen=antigen,
            base=args.base,
            bootstrap_iterations=args.bootstrap_iterations,
            random_seed=args.seed,
        )
        if result is not None:
            generated.append(result)

    expected = len(args.antigens)
    print(f"Generated {len(generated)} of {expected} requested {MHC_LABEL} workbooks.")
    return 0 if len(generated) == expected else 1


if __name__ == "__main__":
    sys.exit(main())
