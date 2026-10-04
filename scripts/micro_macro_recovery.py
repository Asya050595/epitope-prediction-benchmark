#!/usr/bin/env python3
"""English documentation for this module or helper is provided in the repository README and in the surrounding code."""

from __future__ import annotations

from project_paths import DATA_ROOT

import importlib.util
import os
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterable

import pandas as pd


# ---------------------------------------------------------------------------
# Configuration.
# ---------------------------------------------------------------------------

DEFAULT_ROOT = Path(f"{DATA_ROOT}")
ROOT = Path(os.environ.get("RM_ROOT", DEFAULT_ROOT)).expanduser()

DEFAULT_OUT_DIR = (
    ROOT / "p24" / "Micro-Macro"
)
OUT_DIR = Path(
    os.environ.get("MICRO_MACRO_OUT_DIR", DEFAULT_OUT_DIR)
).expanduser()
OUT_FILE = OUT_DIR / "micro_macro_recovery.xlsx"

ANTIGENS = ("p24", "pp65", "PtxS1")

REFERENCE_FILES = {
    "p24": ROOT / "p24" / "p24.xlsx",
    "pp65": ROOT / "pp65" / "pp65.xlsx",
    "PtxS1": ROOT / "PtxS1" / "PtxS1.xlsx",
}


# ---------------------------------------------------------------------------
# Implementation detail; see the repository documentation.
# ---------------------------------------------------------------------------

def build_tool_paths(mhc_class: str, antigen: str) -> dict[str, list[Path]]:
    base = ROOT / f"Processing details{antigen}"
    if mhc_class == "I":
        matches = base / "Matches MHC I"
        return {
            "IEDB": [
                matches / "Matches IEDB_I_Consensus"
                / f"RM_{antigen}_IEDB_I_Consensus.xlsx",
                matches / "Matches IEDB_I_NetMHCpan_4.1_BA"
                / f"RM_{antigen}_IEDB_I_NetMHCpan_4.1_BA.xlsx",
                matches / "Matches IEDB_I_NetMHCpan_4.1_EL"
                / f"RM_{antigen}_IEDB_I_NetMHCpan_4.1_EL.xlsx",
            ],
            "NetCTL 1.2": [
                matches / "Matches NetCTL_1.2"
                / f"RM_{antigen}_NetCTL_1.2.xlsx"
            ],
            "NetMHC 4.0": [
                matches / "Matches NetMHC_4.0"
                / f"RM_{antigen}_NetMHC_4.0.xlsx"
            ],
            "NetMHCpan 4.1": [
                matches / "Matches NetMHCpan_4.1"
                / f"RM_{antigen}_NetMHCpan_4.1.xlsx"
            ],
        }

    if mhc_class == "II":
        matches = base / "Matches MHC II"
        return {
            "IEDB II": [
                matches / "Matches IEDB_II_Consensus"
                / f"RM_{antigen}_IEDB_II_Consensus.xlsx",
                matches / "Matches IEDB_II_NetMHCIIpan_4.1_BA"
                / f"RM_{antigen}_IEDB_II_NetMHCIIpan_4.1_BA.xlsx",
                matches / "Matches IEDB_II_NetMHCIIpan_4.1_EL"
                / f"RM_{antigen}_IEDB_II_NetMHCIIpan_4.1_EL.xlsx",
            ],
            "NetMHCII 2.3": [
                matches / "Matches NetMHC_II_2.3"
                / f"RM_{antigen}_NetMHC_II_2.3.xlsx"
            ],
            "NetMHCIIpan 4.1": [
                matches / "Matches NetMHCIIpan_4.1"
                / f"RM_{antigen}_NetMHCIIpan_4.1.xlsx"
            ],
        }

    raise ValueError(f"Processing details{mhc_class!r}")


# ---------------------------------------------------------------------------
# Input loading.
# ---------------------------------------------------------------------------

def find_allele_sets_file() -> Path:
    explicit = os.environ.get("ALLELE_SETS_PATH")
    if explicit:
        path = Path(explicit).expanduser()
        if not path.is_file():
            raise FileNotFoundError(f"Required input or value was not found{path}")
        return path

    candidates = [
        Path(__file__).resolve().with_name("allele_sets.py"),
        Path(__file__).resolve().with_name("allele_sets(2).py"),
        Path(__file__).resolve().with_name("allele_sets(1).py"),
        ROOT / "allele_sets.py",
        ROOT / "allele_sets(2).py",
        ROOT / "allele_sets(1).py",
    ]
    found = [path for path in candidates if path.is_file()]
    if not found:
        raise FileNotFoundError(
            "Required input or value was not found"
            "Processing details"
        )
    return found[0]


def load_evaluable_sets() -> dict[str, dict[str, set[str]]]:
    path = find_allele_sets_file()
    spec = importlib.util.spec_from_file_location("micro_macro_allele_sets", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Processing details{path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    required = ("EVALUABLE_SET_MHC_I", "EVALUABLE_SET_MHC_II")
    missing = [name for name in required if not hasattr(module, name)]
    if missing:
        raise AttributeError(f"Processing details{path}Required input or value was not found{missing}")

    return {
        "I": {
            antigen: {normalize_allele(a) for a in alleles}
            for antigen, alleles in module.EVALUABLE_SET_MHC_I.items()
        },
        "II": {
            antigen: {normalize_allele(a) for a in alleles}
            for antigen, alleles in module.EVALUABLE_SET_MHC_II.items()
        },
    }


# ---------------------------------------------------------------------------
# Normalization.
# ---------------------------------------------------------------------------

ALLELE_ALIASES = {
    "allele", "hlaallele", "mhcallele", "allelename", "hla",
    "restriction", "restrictionallele",
}
PEPTIDE_ALIASES = {
    "peptide", "peptidesequence", "epitope", "epitopesequence", "sequence",
}
def compact_name(value: object) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value).strip().lower())


def find_column(columns: Iterable[object], aliases: set[str]) -> object | None:
    for column in columns:
        if compact_name(column) in aliases:
            return column
    return None


def normalize_peptide(value: object) -> str:
    if pd.isna(value):
        return ""
    return re.sub(r"\s+", "", str(value)).upper()


def normalize_allele(value: object) -> str:
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    if pd.isna(value):
        return ""
    raw = str(value).strip().upper()
    raw = re.sub(r"^HLA[-_ ]?", "", raw)
    raw = raw.replace("_", "*").replace(" ", "")

    match = re.search(
        r"(DRB[1-5]|DQA1|DQB1|DPA1|DPB1|[ABCE])\*?"
        r"(\d{2})[:]?([0-9]{2})",
        raw,
    )
    if not match:
        return raw
    gene, field1, field2 = match.groups()
    return f"{gene}*{field1}:{field2}"


# ---------------------------------------------------------------------------
# Input loading.
# ---------------------------------------------------------------------------

def load_reference_pairs(
    path: Path,
    mhc_class: str,
) -> set[tuple[str, str]]:
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    if not path.is_file():
        raise FileNotFoundError(f"Required input or value was not found{path}")

    df = pd.read_excel(path, sheet_name=0)
    expected_columns = {
        "I": ("HLA allele", "CTL epitopes"),
        "II": ("HLA allele.1", "HTL epitopes"),
    }
    try:
        allele_col, peptide_col = expected_columns[mhc_class]
    except KeyError as exc:
        raise ValueError(f"Processing details{mhc_class!r}") from exc

    missing_columns = [
        column for column in (allele_col, peptide_col) if column not in df.columns
    ]
    if missing_columns:
        raise ValueError(
            f"Processing details{path}Required input or value was not found{missing_columns}. "
            f"Processing details{list(df.columns)}"
        )

    pairs = {
        (normalize_allele(allele), normalize_peptide(peptide))
        for allele, peptide in zip(df[allele_col], df[peptide_col])
        if normalize_allele(allele) and normalize_peptide(peptide)
    }
    if not pairs:
        raise ValueError(
            f"Processing details{path}Required input or value was not found{mhc_class}"
        )
    return pairs


def read_pairs_from_excel(path: Path) -> set[tuple[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(f"Required input or value was not found{path}")
    excel = pd.ExcelFile(path)
    sheet = "Union_OR" if "Union_OR" in excel.sheet_names else excel.sheet_names[0]
    df = pd.read_excel(path, sheet_name=sheet)
    allele_col = find_column(df.columns, ALLELE_ALIASES)
    peptide_col = find_column(df.columns, PEPTIDE_ALIASES)
    if allele_col is None or peptide_col is None:
        raise ValueError(
            f"Processing details{path}Processing details{sheet}Required input or value was not found"
            f"Processing details{list(df.columns)}"
        )
    return {
        (normalize_allele(a), normalize_peptide(p))
        for a, p in zip(df[allele_col], df[peptide_col])
        if normalize_allele(a) and normalize_peptide(p)
    }


def load_tool_rm_pairs(paths: Iterable[Path]) -> set[tuple[str, str]]:
    pairs: set[tuple[str, str]] = set()
    for path in paths:
        pairs |= read_pairs_from_excel(path)
    return pairs


# ---------------------------------------------------------------------------
# Metric calculation.
# ---------------------------------------------------------------------------

def calculate_recovery(
    reference_pairs: set[tuple[str, str]],
    recovered_pairs: set[tuple[str, str]],
    evaluable_alleles: set[str],
) -> tuple[dict[str, object], list[dict[str, object]]]:
    evaluable_reference = {
        pair for pair in reference_pairs if pair[0] in evaluable_alleles
    }
    if not evaluable_reference:
        raise ValueError("Processing details")

    # Implementation detail; see the repository documentation.
    # Implementation detail; see the repository documentation.
    recovered_reference = recovered_pairs & evaluable_reference

    ref_by_allele: dict[str, set[tuple[str, str]]] = defaultdict(set)
    rec_by_allele: dict[str, set[tuple[str, str]]] = defaultdict(set)
    for pair in evaluable_reference:
        ref_by_allele[pair[0]].add(pair)
    for pair in recovered_reference:
        rec_by_allele[pair[0]].add(pair)

    represented_alleles = sorted(ref_by_allele)
    absent_alleles = sorted(evaluable_alleles - set(represented_alleles))
    if absent_alleles:
        raise ValueError(
            "Allele status"
            + ", ".join(absent_alleles)
        )
    allele_rows: list[dict[str, object]] = []
    allele_recoveries: list[float] = []
    for allele in represented_alleles:
        denominator = len(ref_by_allele[allele])
        numerator = len(rec_by_allele[allele])
        recovery = numerator / denominator
        allele_recoveries.append(recovery)
        allele_rows.append({
            "Allele": allele,
            "Recovered pairs": numerator,
            "Reference pairs": denominator,
            "Allele-specific recovery": recovery,
        })

    summary = {
        "Recovered pairs": len(recovered_reference),
        "Evaluable reference pairs": len(evaluable_reference),
        "Micro recovery": len(recovered_reference) / len(evaluable_reference),
        "HLA alleles in macro": len(represented_alleles),
        "Macro recovery": sum(allele_recoveries) / len(allele_recoveries),
        "Evaluable alleles absent from reference": ", ".join(absent_alleles),
        "RM pairs outside evaluable reference": len(
            recovered_pairs - evaluable_reference
        ),
    }
    return summary, allele_rows


def run_analysis() -> tuple[pd.DataFrame, pd.DataFrame]:
    evaluable_sets = load_evaluable_sets()

    summary_rows: list[dict[str, object]] = []
    allele_rows_all: list[dict[str, object]] = []

    for mhc_class in ("I", "II"):
        for antigen in ANTIGENS:
            evaluable = evaluable_sets[mhc_class][antigen]
            reference_file = REFERENCE_FILES[antigen]
            reference_pairs = load_reference_pairs(reference_file, mhc_class)
            for tool, paths in build_tool_paths(mhc_class, antigen).items():
                recovered_pairs = load_tool_rm_pairs(paths)
                summary, allele_rows = calculate_recovery(
                    reference_pairs, recovered_pairs, evaluable
                )
                summary_rows.append({
                    "HLA class": mhc_class,
                    "Antigen": antigen,
                    "Prediction tool": tool,
                    **summary,
                })
                for row in allele_rows:
                    allele_rows_all.append({
                        "HLA class": mhc_class,
                        "Antigen": antigen,
                        "Prediction tool": tool,
                        **row,
                    })
                print(
                    f"HLA {mhc_class} | {antigen:<5} | {tool:<11} "
                    f"micro={summary['Micro recovery']:.1%}, "
                    f"macro={summary['Macro recovery']:.1%}"
                )

    return pd.DataFrame(summary_rows), pd.DataFrame(allele_rows_all)


# ---------------------------------------------------------------------------
# Implementation detail; see the repository documentation.
# ---------------------------------------------------------------------------

TOOL_ORDER = {
    "I": ["IEDB", "NetMHCpan 4.1", "NetMHC 4.0", "NetCTL 1.2"],
    "II": ["IEDB II", "NetMHCIIpan 4.1", "NetMHCII 2.3"],
}


def write_excel(summary: pd.DataFrame, allele_detail: pd.DataFrame) -> None:
    try:
        from openpyxl import Workbook
        from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
        from openpyxl.utils import get_column_letter
    except ImportError as exc:
        raise RuntimeError("Processing details") from exc

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    workbook = Workbook()
    sheet = workbook.active
    sheet.title = "Micro-Macro"

    dark = "1F4E78"
    light = "D9EAF7"
    section = "B4C6E7"
    white = "FFFFFF"
    thin_gray = Side(style="thin", color="B7B7B7")

    sheet.merge_cells("A1:G1")
    sheet["A1"] = (
        "Micro- and macro-averaged recovery of epitope prediction workflows"
    )
    sheet["A1"].font = Font(bold=True, size=12, color=white)
    sheet["A1"].fill = PatternFill("solid", fgColor=dark)
    sheet["A1"].alignment = Alignment(horizontal="center", vertical="center")
    sheet.row_dimensions[1].height = 24

    row = 3
    for mhc_class in ("I", "II"):
        sheet.merge_cells(start_row=row, start_column=1, end_row=row, end_column=7)
        cell = sheet.cell(row=row, column=1, value=f"HLA class {mhc_class}")
        cell.font = Font(bold=True)
        cell.fill = PatternFill("solid", fgColor=section)
        row += 1

        sheet.merge_cells(start_row=row, start_column=1, end_row=row + 1, end_column=1)
        sheet.cell(row=row, column=1, value="Prediction tool")
        for index, antigen in enumerate(ANTIGENS):
            start_col = 2 + index * 2
            sheet.merge_cells(
                start_row=row, start_column=start_col,
                end_row=row, end_column=start_col + 1,
            )
            sheet.cell(row=row, column=start_col, value=antigen)
            sheet.cell(row=row + 1, column=start_col, value="Micro")
            sheet.cell(row=row + 1, column=start_col + 1, value="Macro")

        for cells in sheet.iter_rows(min_row=row, max_row=row + 1, min_col=1, max_col=7):
            for header in cells:
                header.font = Font(bold=True, color=white)
                header.fill = PatternFill("solid", fgColor=dark)
                header.alignment = Alignment(horizontal="center", vertical="center")
                header.border = Border(
                    left=thin_gray, right=thin_gray,
                    top=thin_gray, bottom=thin_gray,
                )
        row += 2

        for tool in TOOL_ORDER[mhc_class]:
            sheet.cell(row=row, column=1, value=tool)
            for index, antigen in enumerate(ANTIGENS):
                selected = summary[
                    (summary["HLA class"] == mhc_class)
                    & (summary["Antigen"] == antigen)
                    & (summary["Prediction tool"] == tool)
                ]
                if len(selected) != 1:
                    raise RuntimeError(
                        f"Validation status{mhc_class}, {antigen}, {tool}; "
                        f"Processing details{len(selected)}"
                    )
                record = selected.iloc[0]
                start_col = 2 + index * 2
                sheet.cell(row=row, column=start_col, value=record["Micro recovery"])
                sheet.cell(row=row, column=start_col + 1, value=record["Macro recovery"])

            for col in range(1, 8):
                cell = sheet.cell(row=row, column=col)
                cell.border = Border(
                    left=thin_gray, right=thin_gray,
                    top=thin_gray, bottom=thin_gray,
                )
                if col > 1:
                    cell.number_format = "0.0%"
                    cell.alignment = Alignment(horizontal="center")
            if (row % 2) == 0:
                for col in range(1, 8):
                    sheet.cell(row=row, column=col).fill = PatternFill(
                        "solid", fgColor=light
                    )
            row += 1
        row += 2

    note_row = row
    sheet.merge_cells(start_row=note_row, start_column=1, end_row=note_row, end_column=7)
    sheet.cell(
        row=note_row,
        column=1,
        value=(
            "Values are percentages. Micro-averaged recovery was calculated "
            "across all evaluable allele–peptide pairs; macro-averaged recovery "
            "is the unweighted mean of allele-specific recovery estimates."
        ),
    )
    sheet.cell(row=note_row, column=1).font = Font(italic=True, size=9)
    sheet.cell(row=note_row, column=1).alignment = Alignment(wrap_text=True)
    sheet.row_dimensions[note_row].height = 34
    sheet.freeze_panes = "B5"
    sheet.column_dimensions["A"].width = 22
    for col in range(2, 8):
        sheet.column_dimensions[get_column_letter(col)].width = 12

    # Validation.
    details = workbook.create_sheet("Calculation details")
    detail_columns = list(summary.columns)
    details.append(detail_columns)
    for record in summary.itertuples(index=False, name=None):
        details.append(record)
    for cell in details[1]:
        cell.font = Font(bold=True, color=white)
        cell.fill = PatternFill("solid", fgColor=dark)
        cell.alignment = Alignment(horizontal="center", wrap_text=True)
    for name in ("Micro recovery", "Macro recovery"):
        col = detail_columns.index(name) + 1
        for data_row in range(2, details.max_row + 1):
            details.cell(data_row, col).number_format = "0.0%"
    details.freeze_panes = "A2"
    details.auto_filter.ref = details.dimensions
    for col, width in enumerate((11, 11, 18, 16, 24, 15, 20, 15, 31, 28), start=1):
        details.column_dimensions[get_column_letter(col)].width = width

    # Allele-level values used for macro averaging.
    allele_sheet = workbook.create_sheet("Allele-level recovery")
    allele_columns = list(allele_detail.columns)
    allele_sheet.append(allele_columns)
    for record in allele_detail.itertuples(index=False, name=None):
        allele_sheet.append(record)
    for cell in allele_sheet[1]:
        cell.font = Font(bold=True, color=white)
        cell.fill = PatternFill("solid", fgColor=dark)
        cell.alignment = Alignment(horizontal="center", wrap_text=True)
    recovery_col = allele_columns.index("Allele-specific recovery") + 1
    for data_row in range(2, allele_sheet.max_row + 1):
        allele_sheet.cell(data_row, recovery_col).number_format = "0.0%"
    allele_sheet.freeze_panes = "A2"
    allele_sheet.auto_filter.ref = allele_sheet.dimensions
    for col, width in enumerate((11, 11, 18, 15, 17, 17, 24), start=1):
        allele_sheet.column_dimensions[get_column_letter(col)].width = width

    workbook.save(OUT_FILE)
    print(f"Completed successfully{OUT_FILE}")


def main() -> int:
    try:
        summary, allele_detail = run_analysis()
        write_excel(summary, allele_detail)
    except Exception as exc:
        print(f"Error{exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
