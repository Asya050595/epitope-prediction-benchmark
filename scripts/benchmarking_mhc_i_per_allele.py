#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""English documentation for this module or helper is provided in the repository README and in the surrounding code."""

from __future__ import annotations

from project_paths import DATA_ROOT

import argparse
import os
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

try:
    import pandas as pd
except ImportError as e:
    sys.exit(
        "Processing details"
        "  pip install pandas openpyxl\n"
        f"Processing details{e}"
    )

try:
    from allele_sets import get_allele_sets
except ImportError as e:
    sys.exit(
        "Processing details"
        f"Processing details{e}"
    )

try:
    from openpyxl import load_workbook
    from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
    from openpyxl.worksheet.table import Table, TableStyleInfo
except ImportError as e:
    sys.exit(
        "Processing details"
        "  pip install openpyxl\n"
        f"Processing details{e}"
    )


# ══════════════════════════════════════════════════════════════════════════════
# Configuration.
# ══════════════════════════════════════════════════════════════════════════════

DEFAULT_ROOT = f"{DATA_ROOT}"
BASE = Path(os.environ.get("RM_ROOT", DEFAULT_ROOT))
ANTIGENS = ["p24", "pp65", "PtxS1"]
MHC_CLASS = "I"
SCRIPT_STEM = "benchmarking_mhc_i_per_allele"

TOOLS = [
    {
        "label": "IEDB I",
        "sources": [
            ("Matches IEDB_I_Consensus", "IEDB_I_Consensus"),
            ("Matches IEDB_I_NetMHCpan_4.1_BA", "IEDB_I_NetMHCpan_4.1_BA"),
            ("Matches IEDB_I_NetMHCpan_4.1_EL", "IEDB_I_NetMHCpan_4.1_EL"),
        ],
    },
    {"label": "NetCTL", "sources": [("Matches NetCTL_1.2", "NetCTL_1.2")]},
    {"label": "NetMHC", "sources": [("Matches NetMHC_4.0", "NetMHC_4.0")]},
    {"label": "NetMHCpan", "sources": [("Matches NetMHCpan_4.1", "NetMHCpan_4.1")]},
]

DATASET_ALLELE_CANDIDATES = ['HLA allele', 'Allele', 'MHC allele']
DATASET_PEPTIDE_CANDIDATES = ['CTL epitopes', 'CTL epitope', 'CTL Epitope', 'Epitope', 'Peptide']

SHEET_PRIORITY = {
    "RM": ["Union_OR", "RM", "Recovered Matches", "Recovered_Matches", "Matches"],
    "UP": ["Union_OR", "UP", "Unrecovered Pairs", "Unrecovered_Pairs", "Matches"],
}

AA_RE = re.compile(r"^[ACDEFGHIKLMNPQRSTVWYXBZUO]+$", re.IGNORECASE)
KNOWN_HLA_GENES = [
    "DPA1", "DPB1", "DQA1", "DQB1", "DRB1", "DRB3", "DRB4", "DRB5",
    "A", "B", "C", "E", "G",
]
GENE_ORDER = {
    "A": 1, "B": 2, "C": 3, "E": 4, "G": 5,
    "DRB1": 11, "DRB3": 12, "DRB4": 13, "DRB5": 14,
    "DQA1": 21, "DQB1": 22, "DPA1": 23, "DPB1": 24,
}


# ══════════════════════════════════════════════════════════════════════════════
# Normalization.
# ══════════════════════════════════════════════════════════════════════════════

def clean_colname(c: object) -> str:
    return str(c).strip().lower().replace("_", " ").replace("-", " ")


def guess_allele_col(cols: Sequence[object]) -> Optional[object]:
    for c in cols:
        cc = clean_colname(c)
        if cc in {"allele", "hla allele"}:
            return c
    for c in cols:
        cc = clean_colname(c)
        if "allele" in cc and "source" not in cc:
            return c
    return None


def guess_peptide_col(cols: Sequence[object]) -> Optional[object]:
    for c in cols:
        cc = clean_colname(c)
        if cc == "peptide":
            return c
    for c in cols:
        cc = clean_colname(c)
        if "peptide" in cc and "source" not in cc:
            return c
    return None


def find_named_col(cols: Sequence[object], candidates: Sequence[str]) -> Optional[object]:
    by_exact = {str(c).strip(): c for c in cols}
    for name in candidates:
        if name in by_exact:
            return by_exact[name]

    normalized = {clean_colname(c): c for c in cols}
    for name in candidates:
        key = clean_colname(name)
        if key in normalized:
            return normalized[key]
    return None


def normalize_allele(value: object) -> str:
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    s = str(value).strip().upper()
    if not s or s == "NAN":
        return ""
    s = s.replace("HLA-", "").replace("HLA_", "").replace("HLA", "")
    s = s.replace(" ", "")

    if "*" in s and ":" in s:
        gene, rest = s.split("*", 1)
        rest = rest.replace("_", "")
        return f"{gene}*{rest}"

    raw = s.replace("*", "").replace("_", "").replace(":", "")
    for gene in sorted(KNOWN_HLA_GENES, key=len, reverse=True):
        if raw.startswith(gene):
            rest = raw[len(gene):]
            if len(rest) >= 4 and rest[:4].isdigit():
                return f"{gene}*{rest[:2]}:{rest[2:4]}"

    return s


def normalize_peptide(value: object) -> str:
    p = str(value).strip().upper()
    if not p or p == "NAN":
        return ""
    p = re.sub(r"\s+", "", p)
    p = re.sub(r"[^A-Z]", "", p)
    return p


def allele_sort_key(allele: str) -> Tuple[int, str, int, int, str]:
    m = re.match(r"^([A-Z0-9]+)\*(\d+):(\d+)", allele)
    if not m:
        return (999, allele, 999, 999, allele)
    gene, field1, field2 = m.group(1), m.group(2), m.group(3)
    return (GENE_ORDER.get(gene, 998), gene, int(field1), int(field2), allele)


# ══════════════════════════════════════════════════════════════════════════════
# Input loading.
# ══════════════════════════════════════════════════════════════════════════════

def dataset_path_for_antigen(base: Path, antigen: str) -> Path:
    return base / f"Processing details{antigen}" / f"{antigen}.xlsx"


def find_dataset_sheet_and_cols(path: Path) -> Tuple[str, object, object]:
    xls = pd.ExcelFile(path)
    for sheet in xls.sheet_names:
        header = pd.read_excel(xls, sheet_name=sheet, nrows=0)
        a_col = find_named_col(header.columns, DATASET_ALLELE_CANDIDATES)
        p_col = find_named_col(header.columns, DATASET_PEPTIDE_CANDIDATES)
        if a_col is not None and p_col is not None:
            return sheet, a_col, p_col
    raise ValueError(
        f"Required input or value was not found{MHC_CLASS}. "
        f"Processing details{DATASET_ALLELE_CANDIDATES}, peptide={DATASET_PEPTIDE_CANDIDATES}; "
        f"Processing details{path}"
    )


def load_experimental_dataset(
    base: Path,
    antigen: str,
    *,
    strict: bool = False,
) -> Tuple[pd.DataFrame, List[str], Dict[str, object]]:
    path = dataset_path_for_antigen(base, antigen)
    log: Dict[str, object] = {
        "Antigen": antigen,
        "Kind": "DATASET",
        "Path": str(path),
        "Sheet": "",
        "Allele_col": "",
        "Peptide_col": "",
        "Status": "OK",
        "Pairs_after_file_dedup": 0,
        "Unique_alleles": 0,
        "Message": "",
    }

    if not path.exists():
        msg = f"missing dataset: {path}"
        log.update({"Status": "MISSING", "Message": msg})
        if strict:
            raise FileNotFoundError(msg)
        print(f"  [!] DATASET | {antigen:5s} | missing")
        return pd.DataFrame(columns=["antigen", "allele", "peptide"]), [], log

    try:
        sheet, a_col, p_col = find_dataset_sheet_and_cols(path)
        df = pd.read_excel(path, sheet_name=sheet, keep_default_na=False)
        log.update({"Sheet": sheet, "Allele_col": str(a_col), "Peptide_col": str(p_col)})
    except Exception as e:
        msg = f"cannot read dataset: {e}"
        log.update({"Status": "ERROR", "Message": msg})
        if strict:
            raise RuntimeError(msg) from e
        print(f"  [!] DATASET | {antigen:5s} | {msg}")
        return pd.DataFrame(columns=["antigen", "allele", "peptide"]), [], log

    rows = []
    skipped_bad_peptides = 0
    for _, row in df.iterrows():
        allele = normalize_allele(row[a_col])
        peptide = normalize_peptide(row[p_col])
        if not allele or not peptide:
            continue
        if not AA_RE.match(peptide):
            skipped_bad_peptides += 1
            continue
        rows.append((antigen, allele, peptide))

    out = pd.DataFrame(rows, columns=["antigen", "allele", "peptide"])
    if not out.empty:
        out = out.drop_duplicates(subset=["antigen", "allele", "peptide"])

    alleles = sorted(out["allele"].unique().tolist(), key=allele_sort_key) if not out.empty else []
    log["Pairs_after_file_dedup"] = int(len(out))
    log["Unique_alleles"] = int(len(alleles))
    if skipped_bad_peptides:
        log["Message"] = f"skipped suspicious dataset peptides: {skipped_bad_peptides}"

    print(
        f"  DATASET | {antigen:5s} | sheet={log['Sheet']} | "
        f"alleles={len(alleles)} | pairs={len(out)}"
    )
    return out, alleles, log


# ══════════════════════════════════════════════════════════════════════════════
# Input loading.
# ══════════════════════════════════════════════════════════════════════════════

def sheet_has_pair_cols(xls: pd.ExcelFile, sheet_name: str) -> bool:
    try:
        header = pd.read_excel(xls, sheet_name=sheet_name, nrows=0)
    except Exception:
        return False
    return guess_allele_col(header.columns) is not None and guess_peptide_col(header.columns) is not None


def find_pair_sheet(xls: pd.ExcelFile, kind: str) -> Optional[str]:
    for sh in SHEET_PRIORITY[kind]:
        if sh in xls.sheet_names and sheet_has_pair_cols(xls, sh):
            return sh
    for sh in xls.sheet_names:
        if sheet_has_pair_cols(xls, sh):
            return sh
    return None


def find_input_file(base: Path, antigen: str, kind: str, source_dir: str, suffix: str) -> Path:
    folder = base / f"Processing details{antigen}" / f"Matches MHC {MHC_CLASS}" / source_dir
    exact = folder / f"{kind}_{antigen}_{suffix}.xlsx"
    if exact.exists():
        return exact

    candidates = sorted(
        p for p in folder.glob(f"{kind}_{antigen}_{suffix}*.xlsx")
        if not p.name.startswith("~$")
    )
    return candidates[0] if candidates else exact


def load_pair_file(
    path: Path,
    *,
    kind: str,
    antigen: str,
    tool_label: str,
    source_suffix: str,
    strict: bool = False,
) -> Tuple[pd.DataFrame, Dict[str, object]]:
    columns = ["antigen", "tool", "source", "kind", "allele", "peptide", "input_path", "sheet"]
    log: Dict[str, object] = {
        "Antigen": antigen,
        "Tool": tool_label,
        "Source": source_suffix,
        "Kind": kind,
        "Path": str(path),
        "Sheet": "",
        "Status": "OK",
        "Pairs_after_file_dedup": 0,
        "Message": "",
    }

    if not path.exists():
        msg = f"missing: {path}"
        log.update({"Status": "MISSING", "Message": msg})
        if strict:
            raise FileNotFoundError(msg)
        print(f"  [!] {kind:2s} | {antigen:5s} | {tool_label:12s} | missing")
        return pd.DataFrame(columns=columns), log

    try:
        xls = pd.ExcelFile(path)
        sheet = find_pair_sheet(xls, kind)
        if sheet is None:
            raise ValueError(f"Required input or value was not found{xls.sheet_names}")
        df = pd.read_excel(xls, sheet_name=sheet, keep_default_na=False)
        log["Sheet"] = sheet
    except Exception as e:
        msg = f"cannot read: {e}"
        log.update({"Status": "ERROR", "Message": msg})
        if strict:
            raise RuntimeError(f"Processing details{path}: {e}") from e
        print(f"  [!] {kind:2s} | {antigen:5s} | {tool_label:12s} | {path.name} | {msg}")
        return pd.DataFrame(columns=columns), log

    a_col = guess_allele_col(df.columns)
    p_col = guess_peptide_col(df.columns)
    if a_col is None or p_col is None:
        msg = f"no allele/peptide columns; columns={list(df.columns)}"
        log.update({"Status": "ERROR", "Message": msg})
        if strict:
            raise ValueError(f"{path}: {msg}")
        print(f"  [!] {kind:2s} | {antigen:5s} | {tool_label:12s} | {path.name} | {msg}")
        return pd.DataFrame(columns=columns), log

    rows = []
    skipped_bad_peptides = 0
    for _, row in df.iterrows():
        allele = normalize_allele(row[a_col])
        peptide = normalize_peptide(row[p_col])
        if not allele or not peptide:
            continue
        if not AA_RE.match(peptide):
            skipped_bad_peptides += 1
            continue
        rows.append((antigen, tool_label, source_suffix, kind, allele, peptide, str(path), log["Sheet"]))

    out = pd.DataFrame(rows, columns=columns)
    if not out.empty:
        out = out.drop_duplicates(subset=["antigen", "tool", "source", "kind", "allele", "peptide"])

    if skipped_bad_peptides:
        log["Message"] = f"skipped suspicious peptides: {skipped_bad_peptides}"
    log["Pairs_after_file_dedup"] = int(len(out))

    print(
        f"  {kind:2s} | {antigen:5s} | {tool_label:12s} | "
        f"{path.name} | sheet={log['Sheet']} | pairs={len(out)}"
    )
    return out, log


def load_pairs_for_antigen(
    base: Path,
    antigen: str,
    strict: bool = False,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, List[str], List[Dict[str, object]]]:
    logs: List[Dict[str, object]] = []

    print(f"\n=== {SCRIPT_STEM} | antigen={antigen} ===")
    expected, expected_alleles, dataset_log = load_experimental_dataset(base, antigen, strict=strict)
    logs.append(dataset_log)

    rm_logical_frames = []
    up_logical_frames = []

    for tool in TOOLS:
        label = tool["label"]
        source_rm_frames = []
        source_up_frames = []

        for source_dir, suffix in tool["sources"]:
            rm_path = find_input_file(base, antigen, "RM", source_dir, suffix)
            up_path = find_input_file(base, antigen, "UP", source_dir, suffix)

            rm_df, rm_log = load_pair_file(
                rm_path, kind="RM", antigen=antigen, tool_label=label,
                source_suffix=suffix, strict=strict,
            )
            up_df, up_log = load_pair_file(
                up_path, kind="UP", antigen=antigen, tool_label=label,
                source_suffix=suffix, strict=strict,
            )

            source_rm_frames.append(rm_df)
            source_up_frames.append(up_df)
            logs.extend([rm_log, up_log])

        rm_tool = pd.concat(source_rm_frames, ignore_index=True) if source_rm_frames else pd.DataFrame()
        up_tool = pd.concat(source_up_frames, ignore_index=True) if source_up_frames else pd.DataFrame()

        # Implementation detail; see the repository documentation.
        if not rm_tool.empty:
            rm_tool = rm_tool.drop_duplicates(subset=["antigen", "tool", "allele", "peptide"])
        if not up_tool.empty:
            up_tool = up_tool.drop_duplicates(subset=["antigen", "tool", "allele", "peptide"])

        # Implementation detail; see the repository documentation.
        if not rm_tool.empty and not up_tool.empty:
            rm_keys = set(zip(rm_tool["antigen"], rm_tool["tool"], rm_tool["allele"], rm_tool["peptide"]))
            up_keys = list(zip(up_tool["antigen"], up_tool["tool"], up_tool["allele"], up_tool["peptide"]))
            overlap_mask = [key in rm_keys for key in up_keys]
            n_overlap = int(sum(overlap_mask))
            if n_overlap:
                print(f"  [!] {label}: removed {n_overlap} UP pairs that also appear as RM")
                up_tool = up_tool.loc[[not x for x in overlap_mask]].copy()

        rm_logical_frames.append(rm_tool)
        up_logical_frames.append(up_tool)

    rm = pd.concat(rm_logical_frames, ignore_index=True) if rm_logical_frames else pd.DataFrame()
    up = pd.concat(up_logical_frames, ignore_index=True) if up_logical_frames else pd.DataFrame()

    return expected, rm, up, expected_alleles, logs


def reconcile_reference_set(
    expected: pd.DataFrame,
    rm: pd.DataFrame,
    up: pd.DataFrame,
    allowed_alleles: set[str],
    set_label: str,
    strict: bool = False,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, List[str]]:
    reference = expected[expected["allele"].isin(allowed_alleles)].copy()
    reference = reference.drop_duplicates(subset=["antigen", "allele", "peptide"])
    reference_keys = set(zip(reference["allele"], reference["peptide"]))
    reference_alleles = sorted(reference["allele"].unique().tolist(), key=allele_sort_key)

    missing_alleles = allowed_alleles - set(reference_alleles)
    if missing_alleles:
        message = f"{set_label}: {len(missing_alleles)} alleles are absent from the experimental dataset"
        if strict:
            raise ValueError(f"{message}: {sorted(missing_alleles)}")
        print(f"  [!] {message}")

    rm_rows = []
    up_rows = []
    for tool in TOOLS:
        label = tool["label"]
        tool_rm = rm[rm["tool"] == label] if not rm.empty else pd.DataFrame()
        tool_up = up[up["tool"] == label] if not up.empty else pd.DataFrame()
        rm_keys = set(zip(tool_rm["allele"], tool_rm["peptide"])) & reference_keys if not tool_rm.empty else set()
        input_up_keys = set(zip(tool_up["allele"], tool_up["peptide"])) & reference_keys if not tool_up.empty else set()
        added_up_keys = reference_keys - rm_keys - input_up_keys
        final_up_keys = (input_up_keys | added_up_keys) - rm_keys

        rm_rows.extend((label, allele, peptide) for allele, peptide in rm_keys)
        up_rows.extend((label, allele, peptide) for allele, peptide in final_up_keys)
        print(
            f"  {set_label:18s} | {label:12s} | reference={len(reference_keys):3d} "
            f"RM={len(rm_keys):3d} UP={len(final_up_keys):3d} added_UP={len(added_up_keys):3d}"
        )

    columns = ["tool", "allele", "peptide"]
    return (
        reference,
        pd.DataFrame(rm_rows, columns=columns),
        pd.DataFrame(up_rows, columns=columns),
        reference_alleles,
    )


# ══════════════════════════════════════════════════════════════════════════════
# BENCHMARKING PER ALLELE
# ══════════════════════════════════════════════════════════════════════════════

def count_by_tool_allele(df: pd.DataFrame, tool_label: str, allele: str) -> int:
    if df.empty:
        return 0
    return int(((df["tool"] == tool_label) & (df["allele"] == allele)).sum())


def count_expected_by_allele(expected: pd.DataFrame, allele: str) -> int:
    if expected.empty:
        return 0
    return int((expected["allele"] == allele).sum())


def build_per_allele_table(
    expected: pd.DataFrame,
    rm: pd.DataFrame,
    up: pd.DataFrame,
    expected_alleles: List[str],
) -> pd.DataFrame:
    rows = []
    for tool in TOOLS:
        label = tool["label"]
        for allele in expected_alleles:
            exp_n = count_expected_by_allele(expected, allele)
            rm_n = count_by_tool_allele(rm, label, allele)
            up_n = count_by_tool_allele(up, label, allele)
            denom = rm_n + up_n
            recovery_rate = (rm_n / denom * 100) if denom else None
            rows.append({
                "Tool": label,
                "Allele": allele,
                "Experimental_pairs": exp_n,
                "RM per allele": rm_n,
                "UP per allele": up_n,
                "Recovery rate per allele (%)": recovery_rate,
            })
    return pd.DataFrame(rows)


def build_alleles_expected_table(
    expected: pd.DataFrame,
    expected_alleles: List[str],
    set_label: str,
) -> pd.DataFrame:
    rows = []
    for allele in expected_alleles:
        rows.append({
            "Set": set_label,
            "Allele": allele,
            "Experimental_pairs": count_expected_by_allele(expected, allele),
        })
    return pd.DataFrame(rows)


def build_validation_table(
    benchmark: pd.DataFrame,
    expected_alleles: List[str],
    antigen: str,
    set_label: str,
) -> pd.DataFrame:
    rows = []
    expected_count = len(expected_alleles)
    for tool in TOOLS:
        label = tool["label"]
        sub = benchmark[benchmark["Tool"] == label]
        rows.append({
            "Set": set_label,
            "Antigen": antigen,
            "Tool": label,
            "Expected_unique_alleles_in_dataset": expected_count,
            "Rows_for_tool_in_output": int(len(sub)),
            "Allele_count_OK": bool(len(sub) == expected_count),
        })
    return pd.DataFrame(rows)


def output_path_for_antigen(base: Path, antigen: str) -> Path:
    return (
        base
        / f"Processing details{antigen}"
        / f"Benchmarking MHC {MHC_CLASS}"
        / f"Benchmarking_{antigen}_MHC_{MHC_CLASS}_per_allele.xlsx"
    )


def safe_table_name(prefix: str, antigen: str) -> str:
    raw = f"{prefix}_{MHC_CLASS}_{antigen}"
    safe = re.sub(r"[^A-Za-z0-9_]", "_", raw)
    if safe and safe[0].isdigit():
        safe = "T_" + safe
    return safe[:240]


def style_sheet_as_table(ws, table_name: str, title: str, *, percent_cols: Sequence[str] = ()) -> None:
    header_fill = PatternFill("solid", fgColor="1F4E78")
    header_font = Font(color="FFFFFF", bold=True)
    title_font = Font(bold=True, size=14, color="1F4E78")
    thin = Side(style="thin", color="D9E2F3")
    border = Border(left=thin, right=thin, top=thin, bottom=thin)

    ws.insert_rows(1)
    max_col = ws.max_column
    last_col_letter = ws.cell(row=2, column=max_col).column_letter
    ws["A1"] = title
    ws.merge_cells(start_row=1, start_column=1, end_row=1, end_column=max_col)
    ws["A1"].font = title_font
    ws["A1"].alignment = Alignment(horizontal="center")

    for cell in ws[2]:
        cell.fill = header_fill
        cell.font = header_font
        cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
        cell.border = border

    for row in ws.iter_rows(min_row=3, max_row=ws.max_row, min_col=1, max_col=max_col):
        for cell in row:
            cell.border = border
            cell.alignment = Alignment(
                horizontal="center" if cell.column > 2 else "left",
                vertical="top",
                wrap_text=True,
            )

    ws.freeze_panes = "A3"
    table_ref = f"A2:{last_col_letter}{ws.max_row}"
    table = Table(displayName=table_name, ref=table_ref)
    table.tableStyleInfo = TableStyleInfo(
        name="TableStyleMedium2",
        showFirstColumn=False,
        showLastColumn=False,
        showRowStripes=True,
        showColumnStripes=False,
    )
    ws.add_table(table)

    # Implementation detail; see the repository documentation.
    for i in range(1, max_col + 1):
        ws.column_dimensions[ws.cell(row=2, column=i).column_letter].width = 18

    headers = {ws.cell(row=2, column=i).value: i for i in range(1, max_col + 1)}
    width_map = {
        "Tool": 16,
        "Allele": 15,
        "Experimental_pairs": 18,
        "RM per allele": 15,
        "UP per allele": 15,
        "Recovery rate per allele (%)": 25,
        "Path": 80,
        "Message": 45,
    }
    for header, width in width_map.items():
        if header in headers:
            ws.column_dimensions[ws.cell(row=2, column=headers[header]).column_letter].width = width

    for header in percent_cols:
        col_idx = headers.get(header)
        if col_idx:
            col_letter = ws.cell(row=2, column=col_idx).column_letter
            for cell in ws[col_letter][2:]:
                cell.number_format = "0.00"


def save_per_allele_excel(
    evaluable_benchmark: pd.DataFrame,
    full_reference_benchmark: pd.DataFrame,
    alleles_expected: pd.DataFrame,
    validation: pd.DataFrame,
    logs: List[Dict[str, object]],
    out_path: Path,
    antigen: str,
) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
        evaluable_benchmark.to_excel(writer, sheet_name="Evaluable_set_per_allele", index=False)
        full_reference_benchmark.to_excel(writer, sheet_name="Full_reference_per_allele", index=False)
        alleles_expected.to_excel(writer, sheet_name="Alleles_expected", index=False)
        validation.to_excel(writer, sheet_name="Validation", index=False)
        pd.DataFrame(logs).to_excel(writer, sheet_name="Input_files", index=False)

    wb = load_workbook(out_path)

    for sheet_name, table_prefix, set_title in (
        ("Evaluable_set_per_allele", "Evaluable_per_allele", "Evaluable set"),
        ("Full_reference_per_allele", "Full_reference_per_allele", "Full reference set"),
    ):
        style_sheet_as_table(
            wb[sheet_name],
            safe_table_name(table_prefix, antigen),
            f"{set_title} per allele — {antigen} — MHC {MHC_CLASS}",
            percent_cols=["Recovery rate per allele (%)"],
        )
        ws = wb[sheet_name]
        headers = {ws.cell(row=2, column=i).value: i for i in range(1, ws.max_column + 1)}
        rm_col = headers["RM per allele"]
        up_col = headers["UP per allele"]
        recovery_rate_col = headers["Recovery rate per allele (%)"]
        for row in range(3, ws.max_row + 1):
            rm_cell = ws.cell(row=row, column=rm_col).coordinate
            up_cell = ws.cell(row=row, column=up_col).coordinate
            ws.cell(row=row, column=recovery_rate_col).value = f'=IF(({rm_cell}+{up_cell})=0,"",{rm_cell}/({rm_cell}+{up_cell})*100)'
    style_sheet_as_table(
        wb["Alleles_expected"],
        safe_table_name("Alleles_expected", antigen),
        f"Expected alleles from experimental dataset — {antigen} — MHC {MHC_CLASS}",
    )
    style_sheet_as_table(
        wb["Validation"],
        safe_table_name("Validation", antigen),
        f"Validation — {antigen} — MHC {MHC_CLASS}",
    )
    style_sheet_as_table(
        wb["Input_files"],
        safe_table_name("Input_files", antigen),
        f"Input files — {antigen} — MHC {MHC_CLASS}",
    )

    wb.save(out_path)
    print(f"[saved] {out_path}")


def print_benchmark_to_console(
    benchmark: pd.DataFrame,
    expected_alleles: List[str],
    set_label: str,
) -> None:
    printable = benchmark.copy()
    for col in ["Recovery rate per allele (%)"]:
        printable[col] = printable[col].map(lambda x: "" if pd.isna(x) else f"{x:.2f}")
    print(f"\n{set_label} — unique alleles: {len(expected_alleles)}")
    print(printable.to_string(index=False))


# ══════════════════════════════════════════════════════════════════════════════
# CLI
# ══════════════════════════════════════════════════════════════════════════════

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=f"Create MHC {MHC_CLASS} per-allele benchmarking Excel tables")
    parser.add_argument("--base", default=str(BASE), help="Processing details")
    parser.add_argument("--antigen", default="all", choices=ANTIGENS + ["all"], help="Antigen to process; default: all")
    parser.add_argument("--strict", action="store_true", help="Fail if an expected input file/sheet is missing or unreadable")
    parser.add_argument("--dry-run", action="store_true", help="Read inputs and print tables, but do not save Excel files")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    base = Path(args.base)
    antigens = ANTIGENS if args.antigen == "all" else [args.antigen]

    for antigen in antigens:
        expected, rm, up, _, logs = load_pairs_for_antigen(base, antigen, strict=args.strict)
        evaluable_alleles, full_reference_alleles = get_allele_sets(MHC_CLASS, antigen)
        eval_expected, eval_rm, eval_up, eval_alleles = reconcile_reference_set(
            expected, rm, up, evaluable_alleles, "Evaluable set", strict=args.strict,
        )
        full_expected, full_rm, full_up, full_alleles = reconcile_reference_set(
            expected, rm, up, full_reference_alleles, "Full reference set", strict=args.strict,
        )
        evaluable_benchmark = build_per_allele_table(eval_expected, eval_rm, eval_up, eval_alleles)
        full_reference_benchmark = build_per_allele_table(full_expected, full_rm, full_up, full_alleles)
        alleles_expected = pd.concat([
            build_alleles_expected_table(eval_expected, eval_alleles, "Evaluable set"),
            build_alleles_expected_table(full_expected, full_alleles, "Full reference set"),
        ], ignore_index=True)
        validation = pd.concat([
            build_validation_table(evaluable_benchmark, eval_alleles, antigen, "Evaluable set"),
            build_validation_table(full_reference_benchmark, full_alleles, antigen, "Full reference set"),
        ], ignore_index=True)

        print_benchmark_to_console(evaluable_benchmark, eval_alleles, "Evaluable set")
        print_benchmark_to_console(full_reference_benchmark, full_alleles, "Full reference set")

        out_path = output_path_for_antigen(base, antigen)
        if args.dry_run:
            print(f"[dry-run] would save: {out_path}")
        else:
            save_per_allele_excel(
                evaluable_benchmark, full_reference_benchmark,
                alleles_expected, validation, logs, out_path, antigen,
            )


if __name__ == "__main__":
    main()
