#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""English documentation for this module or helper is provided in the repository README and in the surrounding code."""

from __future__ import annotations

from project_paths import DATA_ROOT

import math
import os
import re
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import get_column_letter

try:
    from allele_sets import get_allele_sets
except ImportError as exc:
    sys.exit(
        "Processing details"
        f"Processing details{exc}"
    )


# ──────────────────────────────────────────────────────────────────────────────
# Configuration.
# ──────────────────────────────────────────────────────────────────────────────

MHC_CLASS = "I"
MHC_LABEL = "MHC I"

DEFAULT_ROOT = f"{DATA_ROOT}"
CONTROL_DIR = Path(os.environ.get("RM_ROOT", DEFAULT_ROOT))
ANTIGENS = ["p24", "pp65"]

TOOL_SUFFIXES = [
    "IEDB_I_Consensus",
    "IEDB_I_NetMHCpan_4.1_BA",
    "IEDB_I_NetMHCpan_4.1_EL",
    "NetMHCpan_4.1",
    "NetMHC_4.0",
    "NetCTL_1.2",
]

DATASET_PATHS = {
    "p24": CONTROL_DIR / "p24" / "p24.xlsx",
    "pp65": CONTROL_DIR / "pp65" / "pp65.xlsx",
}

# Implementation detail; see the repository documentation.
# Implementation detail; see the repository documentation.
MHC_I_PREFIXES = ("A*", "B*", "C*")
MHC_II_PREFIXES = ("DR", "DQ", "DP")

Pair = Tuple[str, str]  # (canonical_allele, peptide)


# ──────────────────────────────────────────────────────────────────────────────
# Normalization.
# ──────────────────────────────────────────────────────────────────────────────

def normalize_column_name(name: object) -> str:
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    text = str(name).strip().lower()
    text = re.sub(r"[\s\-./()]+", "_", text)
    text = re.sub(r"[^a-z0-9_]+", "", text)
    text = re.sub(r"_+", "_", text).strip("_")
    return text


def clean_peptide(value: object) -> Optional[str]:
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    if pd.isna(value):
        return None
    peptide = str(value).strip().upper()
    peptide = re.sub(r"\s+", "", peptide)
    if peptide in {"", "NAN", "NONE", "NULL", "-"}:
        return None
    return peptide


def canonicalize_allele(value: object) -> Optional[str]:
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    if pd.isna(value):
        return None

    allele = str(value).strip().upper()
    if allele in {"", "NAN", "NONE", "NULL", "-"}:
        return None

    allele = allele.replace(" ", "")
    allele = allele.replace("HLA-", "")
    allele = allele.replace("HLA_", "")
    allele = allele.replace("HLA", "")
    allele = allele.replace("/", "_")

    # DRB1_0101 -> DRB1*01:01; A_0201 -> A*02:01
    allele = allele.replace("_", "*") if "*" not in allele and "_" in allele else allele

    # Normalization.
    m = re.match(r"^(A|B|C|E|F|G)\*(\d{2}):(\d{2,3})(?::\d{2,3})?$", allele)
    if m:
        return f"{m.group(1)}*{m.group(2)}:{m.group(3)}"

    m = re.match(r"^(DRB\d|DQA\d|DQB\d|DPA\d|DPB\d)\*(\d{2}):(\d{2,3})(?::\d{2,3})?$", allele)
    if m:
        return f"{m.group(1)}*{m.group(2)}:{m.group(3)}"

    # Implementation detail; see the repository documentation.
    m = re.match(r"^(A|B|C|E|F|G)\*(\d{2})(\d{2,3})$", allele)
    if m:
        return f"{m.group(1)}*{m.group(2)}:{m.group(3)}"

    m = re.match(r"^(DRB\d|DQA\d|DQB\d|DPA\d|DPB\d)\*(\d{2})(\d{2,3})$", allele)
    if m:
        return f"{m.group(1)}*{m.group(2)}:{m.group(3)}"

    # A02:01 -> A*02:01
    m = re.match(r"^(A|B|C|E|F|G)(\d{2}):(\d{2,3})(?::\d{2,3})?$", allele)
    if m:
        return f"{m.group(1)}*{m.group(2)}:{m.group(3)}"

    # A0201 -> A*02:01
    m = re.match(r"^(A|B|C|E|F|G)(\d{2})(\d{2,3})$", allele)
    if m:
        return f"{m.group(1)}*{m.group(2)}:{m.group(3)}"

    # DRB10101 -> DRB1*01:01; DQB10302 -> DQB1*03:02
    m = re.match(r"^(DRB\d|DQA\d|DQB\d|DPA\d|DPB\d)(\d{2})(\d{2,3})$", allele)
    if m:
        return f"{m.group(1)}*{m.group(2)}:{m.group(3)}"

    return allele




def extract_canonical_alleles(value: object) -> List[str]:
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    if pd.isna(value):
        return []

    text = str(value).strip().upper()
    if text in {"", "NAN", "NONE", "NULL", "-"}:
        return []

    found: List[str] = []

    # MHC I: HLA-A*02:01, HLA-A02:01, HLA-A0201, A_0201
    for m in re.finditer(
        r"(?:HLA[-_\s]?)?(A|B|C|E|F|G)[*_\-\s]?(\d{2}):?(\d{2,3})",
        text,
    ):
        found.append(f"{m.group(1)}*{m.group(2)}:{m.group(3)}")

    # Implementation detail; see the repository documentation.
    for m in re.finditer(
        r"(?:HLA[-_\s]?)?(DRB\d|DQA\d|DQB\d|DPA\d|DPB\d)[*_\-\s]?(\d{2}):?(\d{2,3})",
        text,
    ):
        found.append(f"{m.group(1)}*{m.group(2)}:{m.group(3)}")

    if not found:
        allele = canonicalize_allele(value)
        if allele is not None:
            found.append(allele)

    # Output generation.
    unique: List[str] = []
    seen = set()
    for allele in found:
        if allele not in seen:
            seen.add(allele)
            unique.append(allele)
    return unique


def belongs_to_mhc_class(allele: str, mhc_class: str) -> bool:
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    if mhc_class == "I":
        return allele.startswith(MHC_I_PREFIXES)
    if mhc_class == "II":
        return allele.startswith(MHC_II_PREFIXES)
    raise ValueError(f"Unknown MHC class: {mhc_class}")


def find_allele_column(df: pd.DataFrame) -> Optional[str]:
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    exact_priority = [
        "allele",
        "hla",
        "hla_allele",
        "mhc_allele",
        "mapped_allele",
        "allele_mapped",
        "matched_allele",
        "input_allele",
        "experimental_allele",
        "tool_allele",
        "mhc_restriction",
        "mhc_restriction_name",
        "mhc_restriction_allele",
        "hla_restriction",
        "hla_restriction_name",
        "restricting_hla_allele",
        "restricting_hla_alleles",
        "restricting_allele",
        "restriction_allele",
        "restriction_name",
    ]
    norm_to_original = {normalize_column_name(c): c for c in df.columns}
    for key in exact_priority:
        if key in norm_to_original:
            return norm_to_original[key]

    for col in df.columns:
        norm = normalize_column_name(col)
        if any(bad in norm for bad in ["count", "number", "num", "total", "class", "length"]):
            continue
        if (
            "allele" in norm
            or norm in {"hla", "mhc"}
            or norm.startswith("hla_")
            or ("restriction" in norm and ("mhc" in norm or "hla" in norm or "allele" in norm))
        ):
            return col
    return None


def find_peptide_column(df: pd.DataFrame) -> Optional[str]:
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    exact_priority = [
        "peptide",
        "epitope",
        "epitope_sequence",
        "peptide_sequence",
        "sequence",
        "matched_peptide",
        "match_peptide",
        "input_peptide",
        "experimental_peptide",
        "target_peptide",
    ]
    norm_to_original = {normalize_column_name(c): c for c in df.columns}
    for key in exact_priority:
        if key in norm_to_original:
            return norm_to_original[key]

    for col in df.columns:
        norm = normalize_column_name(col)
        if any(bad in norm for bad in ["core", "icore", "length", "count", "number", "num", "total"]):
            continue
        if "peptide" in norm or "epitope" in norm or norm == "sequence":
            return col
    return None



def get_column_by_normalized_name(df: pd.DataFrame, keys: Sequence[str]) -> Optional[str]:
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    norm_to_original = {normalize_column_name(c): c for c in df.columns}
    for key in keys:
        if key in norm_to_original:
            return norm_to_original[key]
    return None


def find_dataset_style_pair_columns(df: pd.DataFrame, mhc_class: str) -> Optional[Tuple[str, str]]:
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    if mhc_class == "I":
        allele_col = get_column_by_normalized_name(
            df,
            [
                "hla_allele",
                "mhc_i_hla_allele",
                "mhc_i_allele",
                "ctl_hla_allele",
                "ctl_allele",
            ],
        )
        peptide_col = get_column_by_normalized_name(
            df,
            [
                "ctl_epitopes",
                "ctl_epitope",
                "ctl_peptides",
                "ctl_peptide",
                "mhc_i_epitopes",
                "mhc_i_epitope",
                "mhc_i_peptides",
                "mhc_i_peptide",
            ],
        )
    elif mhc_class == "II":
        allele_col = get_column_by_normalized_name(
            df,
            [
                "hla_allele_1",      # Implementation detail; see the repository documentation.
                "mhc_ii_hla_allele",
                "mhc_ii_allele",
                "htl_hla_allele",
                "htl_allele",
            ],
        )
        peptide_col = get_column_by_normalized_name(
            df,
            [
                "htl_epitopes",
                "htl_epitope",
                "htl_peptides",
                "htl_peptide",
                "mhc_ii_epitopes",
                "mhc_ii_epitope",
                "mhc_ii_peptides",
                "mhc_ii_peptide",
            ],
        )
    else:
        raise ValueError(f"Unknown MHC class: {mhc_class}")

    if allele_col is not None and peptide_col is not None:
        return allele_col, peptide_col
    return None


# ──────────────────────────────────────────────────────────────────────────────
# Input loading.
# ──────────────────────────────────────────────────────────────────────────────

def read_pairs_from_excel(path: Path, mhc_class: str) -> Set[Pair]:
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    pairs: Set[Pair] = set()

    if not path.exists():
        print(f"Required input or value was not found{path}")
        return pairs

    try:
        xls = pd.ExcelFile(path)
    except Exception as exc:
        print(f"Processing details{path}\n      {exc}")
        return pairs

    for sheet_name in xls.sheet_names:
        try:
            df = pd.read_excel(xls, sheet_name=sheet_name, dtype=str)
        except Exception as exc:
            print(f"Processing details{sheet_name}Processing details{path.name}: {exc}")
            continue

        if df.empty:
            continue

        pair_cols = find_dataset_style_pair_columns(df, mhc_class)
        if pair_cols is not None:
            allele_col, peptide_col = pair_cols
        else:
            allele_col = find_allele_column(df)
            peptide_col = find_peptide_column(df)
            if allele_col is None or peptide_col is None:
                continue

        for allele_raw, peptide_raw in zip(df[allele_col], df[peptide_col]):
            alleles = extract_canonical_alleles(allele_raw)
            peptide = clean_peptide(peptide_raw)
            if not alleles or peptide is None:
                continue
            for allele in alleles:
                if not belongs_to_mhc_class(allele, mhc_class):
                    continue
                pairs.add((allele, peptide))

    return pairs


def discover_tool_files(antigen: str, mhc_class: str) -> Dict[str, Dict[str, Optional[Path]]]:
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    matches_dir = CONTROL_DIR / f"Processing details{antigen}" / f"Matches MHC {mhc_class}"
    result: Dict[str, Dict[str, Optional[Path]]] = {}

    if not matches_dir.exists():
        print(f"Required input or value was not found{matches_dir}")
        return result

    for tool_label in TOOL_SUFFIXES:
        folder = matches_dir / f"Matches {tool_label}"
        rm_files = sorted(
            p for p in folder.glob(f"RM_{antigen}_{tool_label}*.xls*")
            if not p.name.startswith("~$")
        )
        up_files = sorted(
            p for p in folder.glob(f"UP_{antigen}_{tool_label}*.xls*")
            if not p.name.startswith("~$")
        )

        rm_file = choose_best_file(rm_files, "RM", antigen, tool_label)
        up_file = choose_best_file(up_files, "UP", antigen, tool_label)

        if rm_file is None and up_file is None:
            print(f"    ⚠ {tool_label}Required input or value was not found{folder}")
            continue

        result[tool_label] = {"RM": rm_file, "UP": up_file}

        if rm_file is None:
            print(f"    ⚠ {tool_label}Required input or value was not found")
        if up_file is None:
            print(f"    ⚠ {tool_label}Required input or value was not found")

    return result


def choose_best_file(files: Sequence[Path], kind: str, antigen: str, tool_label: str) -> Optional[Path]:
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    if not files:
        return None

    expected_names = [
        f"{kind}_{antigen}_{tool_label}.xlsx",
        f"{kind}_{antigen}_{tool_label}.xls",
    ]
    for expected in expected_names:
        for path in files:
            if path.name == expected:
                return path

    # Implementation detail; see the repository documentation.
    if len(files) == 1:
        return files[0]

    # Path configuration.
    chosen = sorted(files, key=lambda p: (len(p.name), p.name.lower()))[0]
    print(f"    ⚠ {tool_label}Processing details{kind}Processing details{chosen.name}")
    return chosen


# ──────────────────────────────────────────────────────────────────────────────
# Implementation detail; see the repository documentation.
# ──────────────────────────────────────────────────────────────────────────────

def jaccard_index(a: Set[Pair], b: Set[Pair]) -> Optional[float]:
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    union = a | b
    if len(union) == 0:
        return None
    return len(a & b) / len(union)


def build_jaccard_matrix(pair_sets: Dict[str, Set[Pair]]) -> pd.DataFrame:
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    tools = sorted(pair_sets.keys(), key=str.lower)
    matrix = pd.DataFrame(index=tools, columns=tools, dtype=object)

    for tool_a in tools:
        for tool_b in tools:
            if tool_a == tool_b:
                matrix.loc[tool_a, tool_b] = 1.0 if len(pair_sets[tool_a]) > 0 else None
            else:
                matrix.loc[tool_a, tool_b] = jaccard_index(pair_sets[tool_a], pair_sets[tool_b])

    return matrix


def write_excel(output_path: Path, rm_matrix: pd.DataFrame, all_matrix: pd.DataFrame) -> None:
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        rm_matrix.to_excel(writer, sheet_name="Jaccard_RM")
        all_matrix.to_excel(writer, sheet_name="Jaccard_RM_plus_UP")

    format_workbook(output_path)


def format_workbook(path: Path) -> None:
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    wb = load_workbook(path)

    header_fill = PatternFill("solid", fgColor="D9EAF7")
    index_fill = PatternFill("solid", fgColor="F2F2F2")
    thin = Side(style="thin", color="BFBFBF")
    border = Border(left=thin, right=thin, top=thin, bottom=thin)

    for ws in wb.worksheets:
        ws.freeze_panes = "B2"

        max_row = ws.max_row
        max_col = ws.max_column

        for cell in ws[1]:
            cell.font = Font(bold=True)
            cell.fill = header_fill
            cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
            cell.border = border

        for row in range(2, max_row + 1):
            index_cell = ws.cell(row=row, column=1)
            index_cell.font = Font(bold=True)
            index_cell.fill = index_fill
            index_cell.alignment = Alignment(horizontal="left", vertical="center", wrap_text=True)
            index_cell.border = border

        for row in ws.iter_rows(min_row=2, min_col=2, max_row=max_row, max_col=max_col):
            for cell in row:
                cell.alignment = Alignment(horizontal="center", vertical="center")
                cell.border = border
                if isinstance(cell.value, (int, float)) and not math.isnan(float(cell.value)):
                    cell.number_format = "0.000"

        # Implementation detail; see the repository documentation.
        ws.column_dimensions["A"].width = 36
        for col_idx in range(2, max_col + 1):
            ws.column_dimensions[get_column_letter(col_idx)].width = 16

    wb.save(path)


# ──────────────────────────────────────────────────────────────────────────────
# Main entry point.
# ──────────────────────────────────────────────────────────────────────────────

def process_antigen(antigen: str) -> None:
    print("=" * 80)
    print(f"Processing details{antigen} | {MHC_LABEL}")
    print("=" * 80)

    dataset_path = DATASET_PATHS[antigen]
    evaluable_alleles, _ = get_allele_sets(MHC_CLASS, antigen)
    dataset_pairs_all = read_pairs_from_excel(dataset_path, MHC_CLASS)
    dataset_pairs = {pair for pair in dataset_pairs_all if pair[0] in evaluable_alleles}
    print(
        f"  Evaluable set: {len(evaluable_alleles)}Allele status"
        f"Processing details{len(dataset_pairs)}Processing details{len(dataset_pairs_all)}"
    )

    tool_files = discover_tool_files(antigen, MHC_CLASS)
    if not tool_files:
        print(f"Processing details{antigen} {MHC_LABEL}Processing details")
        return

    rm_sets: Dict[str, Set[Pair]] = {}
    all_sets: Dict[str, Set[Pair]] = {}

    print("Processing details")
    for tool_label, files in tool_files.items():
        rm_file = files.get("RM")
        up_file = files.get("UP")

        rm_all = read_pairs_from_excel(rm_file, MHC_CLASS) if rm_file is not None else set()
        up_all = read_pairs_from_excel(up_file, MHC_CLASS) if up_file is not None else set()
        rm_pairs = {pair for pair in rm_all if pair[0] in evaluable_alleles}
        up_pairs = {pair for pair in up_all if pair[0] in evaluable_alleles}

        # Implementation detail; see the repository documentation.
        if dataset_pairs:
            outside_dataset = rm_pairs - dataset_pairs
            if outside_dataset:
                print(
                    f"    ⚠ {tool_label}: {len(outside_dataset)}Required input or value was not found"
                    f"Processing details"
                )

        rm_sets[tool_label] = rm_pairs
        all_sets[tool_label] = rm_pairs | up_pairs

        print(
            f"    {tool_label}: RM={len(rm_pairs)} ({len(rm_all) - len(rm_pairs)}Processing details"
            f"UP={len(up_pairs)} ({len(up_all) - len(up_pairs)}Processing details"
            f"RM+UP unique={len(all_sets[tool_label])}"
        )

    rm_matrix = build_jaccard_matrix(rm_sets)
    all_matrix = build_jaccard_matrix(all_sets)

    output_path = (
        CONTROL_DIR
        / f"Processing details{antigen}"
        / f"Jaccard Index MHC {MHC_CLASS}"
        / f"Jaccard{antigen}_MHC_{MHC_CLASS}.xlsx"
    )
    write_excel(output_path, rm_matrix, all_matrix)
    print(f"Saved output{output_path}\n")


def main() -> None:
    for antigen in ANTIGENS:
        process_antigen(antigen)


if __name__ == "__main__":
    main()
