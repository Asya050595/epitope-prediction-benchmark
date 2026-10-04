#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Portable implementation of benchmarking_mhc_ii. See the repository README and data/README.md for inputs, outputs, and execution instructions."""

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
        "ERROR: установите pandas и openpyxl:\n"
        "  pip install pandas openpyxl\n"
        f"Детали: {e}"
    )

try:
    from allele_sets import get_allele_sets
except ImportError as e:
    sys.exit(
        "ERROR: положите allele_sets.py рядом с benchmarking-скриптом.\n"
        f"Детали: {e}"
    )

try:
    from openpyxl import load_workbook
    from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
    from openpyxl.worksheet.table import Table, TableStyleInfo
except ImportError as e:
    sys.exit(
        "ERROR: установите openpyxl:\n"
        "  pip install openpyxl\n"
        f"Детали: {e}"
    )


# ══════════════════════════════════════════════════════════════════════════════
# Configuration.
# ══════════════════════════════════════════════════════════════════════════════

DEFAULT_ROOT = f"{DATA_ROOT}"
BASE = Path(os.environ.get("RM_ROOT", DEFAULT_ROOT))
ANTIGENS = ["p24", "pp65", "PtxS1"]
MHC_CLASS = "II"
SCRIPT_STEM = "benchmarking_mhc_ii"

# Implementation detail; see the repository documentation.
TOOLS = [
    {
        "label": "IEDB II",
        "sources": [
            ("Matches IEDB_II_Consensus", "IEDB_II_Consensus"),
            ("Matches IEDB_II_NetMHCIIpan_4.1_BA", "IEDB_II_NetMHCIIpan_4.1_BA"),
            ("Matches IEDB_II_NetMHCIIpan_4.1_EL", "IEDB_II_NetMHCIIpan_4.1_EL"),
        ],
    },
    {"label": "NetMHC II", "sources": [("Matches NetMHC_II_2.3", "NetMHC_II_2.3")]},
    {"label": "NetMHCIIpan", "sources": [("Matches NetMHCIIpan_4.1", "NetMHCIIpan_4.1")]},
]

DATASET_ALLELE_CANDIDATES = ["HLA allele.1", "HTL HLA allele", "HLA class II allele", "HLA allele", "Allele", "MHC allele"]
DATASET_PEPTIDE_CANDIDATES = ["HTL epitopes", "HTL epitope", "Helper epitopes", "Epitope", "Peptide"]

SHEET_PRIORITY = {
    "RM": ["Union_OR", "RM", "Recovered Matches", "Recovered_Matches", "Matches"],
    "UP": ["Union_OR", "UP", "Unrecovered Pairs", "Unrecovered_Pairs", "Matches"],
}

AA_RE = re.compile(r"^[ACDEFGHIKLMNPQRSTVWYXBZUO]+$", re.IGNORECASE)
KNOWN_HLA_GENES = [
    "DPA1", "DPB1", "DQA1", "DQB1", "DRB1", "DRB3", "DRB4", "DRB5",
    "A", "B", "C", "E", "G",
]


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
    """Portable implementation of benchmarking_mhc_ii. See the repository README and data/README.md for inputs, outputs, and execution instructions."""
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


def dataset_path_for_antigen(base: Path, antigen: str) -> Path:
    return base / f"{antigen}" / f"{antigen}.xlsx"


def find_dataset_sheet_and_cols(path: Path) -> Tuple[str, object, object]:
    xls = pd.ExcelFile(path)
    for sheet in xls.sheet_names:
        header = pd.read_excel(xls, sheet_name=sheet, nrows=0)
        allele_col = find_named_col(header.columns, DATASET_ALLELE_CANDIDATES)
        peptide_col = find_named_col(header.columns, DATASET_PEPTIDE_CANDIDATES)
        if allele_col is not None and peptide_col is not None:
            return sheet, allele_col, peptide_col
    raise ValueError(
        f"Не найдены колонки экспериментального датасета для MHC {MHC_CLASS}. "
        f"Искала allele={DATASET_ALLELE_CANDIDATES}, peptide={DATASET_PEPTIDE_CANDIDATES}; "
        f"файл={path}"
    )


def load_experimental_dataset(
    base: Path,
    antigen: str,
    *,
    strict: bool = False,
) -> Tuple[pd.DataFrame, Dict[str, object]]:
    path = dataset_path_for_antigen(base, antigen)
    log: Dict[str, object] = {
        "Antigen": antigen,
        "Tool": "REFERENCE_DATASET",
        "Source": "",
        "Kind": "DATASET",
        "Path": str(path),
        "Sheet": "",
        "Status": "OK",
        "Pairs_after_file_dedup": 0,
        "Message": "",
    }
    empty = pd.DataFrame(columns=["antigen", "allele", "peptide"])

    if not path.exists():
        msg = f"missing dataset: {path}"
        log.update({"Status": "MISSING", "Message": msg})
        if strict:
            raise FileNotFoundError(msg)
        print(f"  [!] DATASET | {antigen:5s} | missing")
        return empty, log

    try:
        sheet, allele_col, peptide_col = find_dataset_sheet_and_cols(path)
        df = pd.read_excel(path, sheet_name=sheet, keep_default_na=False)
        log["Sheet"] = sheet
    except Exception as e:
        msg = f"cannot read dataset: {e}"
        log.update({"Status": "ERROR", "Message": msg})
        if strict:
            raise RuntimeError(msg) from e
        print(f"  [!] DATASET | {antigen:5s} | {msg}")
        return empty, log

    rows = []
    for _, row in df.iterrows():
        allele = normalize_allele(row[allele_col])
        peptide = normalize_peptide(row[peptide_col])
        if allele and peptide and AA_RE.match(peptide):
            rows.append((antigen, allele, peptide))

    expected = pd.DataFrame(rows, columns=["antigen", "allele", "peptide"])
    if not expected.empty:
        expected = expected.drop_duplicates(subset=["antigen", "allele", "peptide"])
    log["Pairs_after_file_dedup"] = int(len(expected))
    print(f"  DATASET | {antigen:5s} | sheet={log['Sheet']} | pairs={len(expected)}")
    return expected, log


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
    # Threshold handling.
    for sh in SHEET_PRIORITY[kind]:
        if sh in xls.sheet_names and sheet_has_pair_cols(xls, sh):
            return sh

    # Implementation detail; see the repository documentation.
    for sh in xls.sheet_names:
        if sheet_has_pair_cols(xls, sh):
            return sh
    return None


def find_input_file(base: Path, antigen: str, kind: str, source_dir: str, suffix: str) -> Path:
    """Portable implementation of benchmarking_mhc_ii. See the repository README and data/README.md for inputs, outputs, and execution instructions."""
    folder = base / f"{antigen}" / f"Matches MHC {MHC_CLASS}" / source_dir
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
    """Portable implementation of benchmarking_mhc_ii. See the repository README and data/README.md for inputs, outputs, and execution instructions."""
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
            raise ValueError(f"не найден лист с колонками allele/peptide; sheets={xls.sheet_names}")
        df = pd.read_excel(xls, sheet_name=sheet, keep_default_na=False)
        log["Sheet"] = sheet
    except Exception as e:
        msg = f"cannot read: {e}"
        log.update({"Status": "ERROR", "Message": msg})
        if strict:
            raise RuntimeError(f"Не удалось прочитать {path}: {e}") from e
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
        out = out.drop_duplicates(subset=["antigen", "tool", "kind", "allele", "peptide"])

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
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, List[Dict[str, object]]]:
    rm_frames = []
    up_frames = []
    logs: List[Dict[str, object]] = []

    print(f"\n=== {SCRIPT_STEM} | antigen={antigen} ===")
    expected, dataset_log = load_experimental_dataset(base, antigen, strict=strict)
    logs.append(dataset_log)
    for tool in TOOLS:
        label = tool["label"]
        for source_dir, suffix in tool["sources"]:
            rm_path = find_input_file(base, antigen, "RM", source_dir, suffix)
            up_path = find_input_file(base, antigen, "UP", source_dir, suffix)

            rm_df, rm_log = load_pair_file(
                rm_path,
                kind="RM",
                antigen=antigen,
                tool_label=label,
                source_suffix=suffix,
                strict=strict,
            )
            up_df, up_log = load_pair_file(
                up_path,
                kind="UP",
                antigen=antigen,
                tool_label=label,
                source_suffix=suffix,
                strict=strict,
            )

            rm_frames.append(rm_df)
            up_frames.append(up_df)
            logs.extend([rm_log, up_log])

    rm = pd.concat(rm_frames, ignore_index=True) if rm_frames else pd.DataFrame()
    up = pd.concat(up_frames, ignore_index=True) if up_frames else pd.DataFrame()

    # Implementation detail; see the repository documentation.
    if not rm.empty:
        rm = rm.drop_duplicates(subset=["antigen", "tool", "allele", "peptide"])
    if not up.empty:
        up = up.drop_duplicates(subset=["antigen", "tool", "allele", "peptide"])

    # Deduplication.
    if not rm.empty and not up.empty:
        rm_keys = set(zip(rm["antigen"], rm["tool"], rm["allele"], rm["peptide"]))
        up_keys = list(zip(up["antigen"], up["tool"], up["allele"], up["peptide"]))
        overlap_mask = [key in rm_keys for key in up_keys]
        n_overlap = int(sum(overlap_mask))
        if n_overlap:
            print(f"  [!] Removed {n_overlap} UP pairs that also appear as RM within the same logical tool")
            up = up.loc[[not x for x in overlap_mask]].copy()

    return expected, rm, up, logs


# ══════════════════════════════════════════════════════════════════════════════
# BENCHMARKING
# ══════════════════════════════════════════════════════════════════════════════

def reconcile_reference_set(
    expected: pd.DataFrame,
    rm: pd.DataFrame,
    up: pd.DataFrame,
    allowed_alleles: set[str],
    set_label: str,
    strict: bool = False,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    reference = expected[expected["allele"].isin(allowed_alleles)].copy()
    reference = reference.drop_duplicates(subset=["allele", "peptide"])
    reference_keys = set(zip(reference["allele"], reference["peptide"]))

    missing_alleles = allowed_alleles - set(reference["allele"])
    if missing_alleles:
        message = f"{set_label}: {len(missing_alleles)} alleles are absent from the experimental dataset"
        print(f"  [!] {message}: {sorted(missing_alleles)}")

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
    return pd.DataFrame(rm_rows, columns=columns), pd.DataFrame(up_rows, columns=columns)


def build_benchmark_table(rm: pd.DataFrame, up: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for tool in TOOLS:
        label = tool["label"]
        rm_n = int((rm["tool"] == label).sum()) if not rm.empty else 0
        up_n = int((up["tool"] == label).sum()) if not up.empty else 0
        denom = rm_n + up_n
        recovery_rate = (rm_n / denom * 100) if denom else None
        rows.append({"Tool": label, "Reference pairs": denom, "RM": rm_n, "UP": up_n, "Recovery rate (%)": recovery_rate})
    return pd.DataFrame(rows, columns=["Tool", "Reference pairs", "RM", "UP", "Recovery rate (%)"])


def output_path_for_antigen(base: Path, antigen: str) -> Path:
    return (
        base
        / f"{antigen}"
        / f"Benchmarking MHC {MHC_CLASS}"
        / f"Benchmarking_{antigen}_MHC_{MHC_CLASS}.xlsx"
    )


def save_benchmark_excel(
    evaluable_benchmark: pd.DataFrame,
    full_reference_benchmark: pd.DataFrame,
    logs: List[Dict[str, object]],
    out_path: Path,
    antigen: str,
) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
        for sheet_name, benchmark in (
            ("Evaluable_set", evaluable_benchmark),
            ("Full_reference_set", full_reference_benchmark),
        ):
            main = benchmark[["Tool", "Reference pairs", "RM", "UP"]].copy()
            main.to_excel(writer, sheet_name=sheet_name, index=False)
        pd.DataFrame(logs).to_excel(writer, sheet_name="Input_files", index=False)

    wb = load_workbook(out_path)
    details_ws = wb["Input_files"]
    header_fill = PatternFill("solid", fgColor="1F4E78")
    header_font = Font(color="FFFFFF", bold=True)
    thin = Side(style="thin", color="D9E2F3")
    border = Border(left=thin, right=thin, top=thin, bottom=thin)

    for sheet_name, set_title, table_prefix in (
        ("Evaluable_set", "Evaluable set", "Evaluable"),
        ("Full_reference_set", "Full reference set", "FullReference"),
    ):
        ws = wb[sheet_name]
        ws["E1"] = "Recovery rate (%)"
        ws.insert_rows(1)
        ws["A1"] = f"{set_title} — {antigen} — MHC {MHC_CLASS}"
        ws.merge_cells("A1:E1")
        ws["A1"].font = Font(bold=True, size=14, color="1F4E78")
        ws["A1"].alignment = Alignment(horizontal="center")

        for row in range(3, ws.max_row + 1):
            ws[f"E{row}"] = f'=IF((C{row}+D{row})=0,"",C{row}/(C{row}+D{row})*100)'
        for cell in ws[2]:
            cell.fill = header_fill
            cell.font = header_font
            cell.alignment = Alignment(horizontal="center")
            cell.border = border
        for row in ws.iter_rows(min_row=3, max_row=ws.max_row, min_col=1, max_col=5):
            for cell in row:
                cell.border = border
                cell.alignment = Alignment(horizontal="center" if cell.column > 1 else "left")

        for col_letter, width in {"A": 18, "B": 18, "C": 12, "D": 12, "E": 14}.items():
            ws.column_dimensions[col_letter].width = width
        for cell in ws["E"][2:]:
            cell.number_format = "0.00"
        ws.freeze_panes = "A3"

        table_ref = f"A2:E{ws.max_row}"
        table_name = f"{table_prefix}_{MHC_CLASS}_{antigen}".replace("-", "_")
        table = Table(displayName=table_name, ref=table_ref)
        table.tableStyleInfo = TableStyleInfo(
            name="TableStyleMedium2",
            showFirstColumn=False,
            showLastColumn=False,
            showRowStripes=True,
            showColumnStripes=False,
        )
        ws.add_table(table)

    # Formatting.
    details_ws.freeze_panes = "A2"
    for cell in details_ws[1]:
        cell.fill = header_fill
        cell.font = header_font
        cell.alignment = Alignment(horizontal="center")
        cell.border = border
    for row in details_ws.iter_rows(min_row=2):
        for cell in row:
            cell.border = border
            cell.alignment = Alignment(vertical="top", wrap_text=True)
    for col_letter, width in {
        "A": 10, "B": 16, "C": 34, "D": 8, "E": 80,
        "F": 18, "G": 12, "H": 20, "I": 40,
    }.items():
        details_ws.column_dimensions[col_letter].width = width

    wb.save(out_path)
    print(f"[saved] {out_path}")


def print_benchmark_to_console(benchmark: pd.DataFrame, set_label: str) -> None:
    printable = benchmark.copy()
    for col in ["Recovery rate (%)"]:
        printable[col] = printable[col].map(lambda x: "" if pd.isna(x) else f"{x:.2f}")
    print(f"\n{set_label}\n" + printable.to_string(index=False))


# ══════════════════════════════════════════════════════════════════════════════
# CLI
# ══════════════════════════════════════════════════════════════════════════════

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create MHC II benchmarking Excel tables from RM/UP match files")
    parser.add_argument("--base", default=str(BASE), help="Base directory with '...'")
    parser.add_argument("--antigen", default="all", choices=ANTIGENS + ["all"], help="Antigen to process; default: all")
    missing_group = parser.add_mutually_exclusive_group()
    missing_group.add_argument(
        "--strict", dest="strict", action="store_true", default=True,
        help="Fail if an expected input file or sheet is missing (default)",
    )
    missing_group.add_argument(
        "--allow-missing", dest="strict", action="store_false",
        help="Continue with missing inputs; intended only for diagnostics",
    )
    parser.add_argument("--dry-run", action="store_true", help="Read inputs and print tables, but do not save Excel files")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    base = Path(args.base)
    antigens = ANTIGENS if args.antigen == "all" else [args.antigen]

    for antigen in antigens:
        expected, rm, up, logs = load_pairs_for_antigen(base, antigen, strict=args.strict)
        evaluable_alleles, full_reference_alleles = get_allele_sets(MHC_CLASS, antigen)
        evaluable_rm, evaluable_up = reconcile_reference_set(
            expected, rm, up, evaluable_alleles, "Evaluable set", strict=args.strict,
        )
        full_rm, full_up = reconcile_reference_set(
            expected, rm, up, full_reference_alleles, "Full reference set", strict=args.strict,
        )
        evaluable_benchmark = build_benchmark_table(evaluable_rm, evaluable_up)
        full_reference_benchmark = build_benchmark_table(full_rm, full_up)
        print_benchmark_to_console(evaluable_benchmark, "Evaluable set")
        print_benchmark_to_console(full_reference_benchmark, "Full reference set")

        out_path = output_path_for_antigen(base, antigen)
        if args.dry_run:
            print(f"[dry-run] would save: {out_path}")
        else:
            save_benchmark_excel(
                evaluable_benchmark, full_reference_benchmark, logs, out_path, antigen,
            )


if __name__ == "__main__":
    main()
