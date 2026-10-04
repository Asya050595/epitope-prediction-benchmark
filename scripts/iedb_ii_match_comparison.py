"""English documentation for this module or helper is provided in the repository README and in the surrounding code."""

from project_paths import DATA_ROOT

import os
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter

# Implementation detail; see the repository documentation.

BASE = f"{DATA_ROOT}"
ANTIGENS = ["p24", "pp65", "PtxS1"]

# Implementation detail; see the repository documentation.
MODES = {
    "EL":        ("Matches IEDB_II_NetMHCIIpan_4.1_EL",  "IEDB_II_NetMHCIIpan_4.1_EL"),
    "BA":        ("Matches IEDB_II_NetMHCIIpan_4.1_BA",  "IEDB_II_NetMHCIIpan_4.1_BA"),
    "Consensus": ("Matches IEDB_II_Consensus",             "IEDB_II_Consensus"),
}

# Formatting.

HEADER_FILL  = PatternFill("solid", start_color="4472C4")
SUBHEAD_FILL = PatternFill("solid", start_color="D9E1F2")
SINGLE_FILL  = PatternFill("solid", start_color="E2EFDA")
UNION_FILL   = PatternFill("solid", start_color="FFF2CC")
VENN_FILL    = PatternFill("solid", start_color="FCE4D6")

HEADER_FONT  = Font(name="Arial", bold=True, color="FFFFFF", size=11)
SUBHEAD_FONT = Font(name="Arial", bold=True, size=10)
BODY_FONT    = Font(name="Arial", size=10)
BOLD_FONT    = Font(name="Arial", bold=True, size=10)

CENTER = Alignment(horizontal="center", vertical="center", wrap_text=True)
LEFT   = Alignment(horizontal="left",   vertical="center")

THIN_SIDE   = Side(style="thin",   color="BFBFBF")
THIN_BORDER = Border(left=THIN_SIDE, right=THIN_SIDE,
                     top=THIN_SIDE,  bottom=THIN_SIDE)
MED_SIDE    = Side(style="medium", color="4472C4")
MED_BORDER  = Border(left=MED_SIDE, right=MED_SIDE,
                     top=MED_SIDE,  bottom=MED_SIDE)

PCT_FMT = "0.0%"
NUM_FMT = "#,##0"


def style(cell, font=None, fill=None, align=None, border=None, num_format=None):
    if font:       cell.font          = font
    if fill:       cell.fill          = fill
    if align:      cell.alignment     = align
    if border:     cell.border        = border
    if num_format: cell.number_format = num_format


def write_header_row(ws, row, values, fills, fonts, col_start=1):
    for i, (val, fill, font) in enumerate(zip(values, fills, fonts)):
        c = ws.cell(row=row, column=col_start + i, value=val)
        style(c, font=font, fill=fill, align=CENTER, border=THIN_BORDER)


# Input loading.

MATCH_SHEET = "Union_OR"

def read_matches(path: str) -> set[tuple]:
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    try:
        xls = pd.ExcelFile(path)
    except Exception as exc:
        raise RuntimeError(f"Processing details{path}: {exc}")

    if MATCH_SHEET not in xls.sheet_names:
        raise ValueError(
            f"Processing details{path}Processing details{MATCH_SHEET}'. "
            f"Processing details{xls.sheet_names}"
        )

    df = pd.read_excel(xls, sheet_name=MATCH_SHEET)
    df.columns = df.columns.str.strip().str.lower()

    missing = {"allele", "peptide"} - set(df.columns)
    if missing:
        raise ValueError(
            f"Processing details{MATCH_SHEET}Processing details{path}Processing details{sorted(missing)}"
        )

    # Implementation detail; see the repository documentation.
    if df.empty:
        return set()

    # Implementation detail; see the repository documentation.
    df = df.dropna(subset=["allele", "peptide"])

    return set(zip(
        df["allele"].astype(str).str.strip(),
        df["peptide"].astype(str).str.strip(),
    ))


# Metric calculation.

def compute_stats(sets: dict[str, set]) -> dict:
    el  = sets["EL"]
    ba  = sets["BA"]
    con = sets["Consensus"]

    union   = el | ba | con
    n_union = len(union)

    single = {
        "EL":        len(el),
        "BA":        len(ba),
        "Consensus": len(con),
    }
    single_pct = {k: (v / n_union if n_union else 0) for k, v in single.items()}

    only_el  = el  - ba  - con
    only_ba  = ba  - el  - con
    only_con = con - el  - ba
    el_ba    = (el  & ba)  - con
    el_con   = (el  & con) - ba
    ba_con   = (ba  & con) - el
    all_3    = el  & ba  & con

    venn = {
        "Processing details":           len(only_el),
        "Processing details":           len(only_ba),
        "Processing details":    len(only_con),
        "EL + BA":             len(el_ba),
        "EL + Consensus":      len(el_con),
        "BA + Consensus":      len(ba_con),
        "EL + BA + Consensus": len(all_3),
    }
    venn_pct = {k: (v / n_union if n_union else 0) for k, v in venn.items()}

    return {
        "union":      n_union,
        "single":     single,
        "single_pct": single_pct,
        "venn":       venn,
        "venn_pct":   venn_pct,
    }


# Output generation.

def write_antigen_sheet(ws, antigen: str, stats: dict):
    ws.column_dimensions["A"].width = 30
    ws.column_dimensions["B"].width = 14
    ws.column_dimensions["C"].width = 14

    n_union = stats["union"]
    row = 1

    # Implementation detail; see the repository documentation.
    ws.merge_cells(start_row=row, start_column=1, end_row=row, end_column=3)
    c = ws.cell(row=row, column=1, value=f"Processing details{antigen}")
    style(c, font=HEADER_FONT, fill=HEADER_FILL, align=CENTER, border=MED_BORDER)
    row += 2

    # Implementation detail; see the repository documentation.
    ws.merge_cells(start_row=row, start_column=1, end_row=row, end_column=3)
    c = ws.cell(row=row, column=1,
                value="Processing details")
    style(c, font=SUBHEAD_FONT, fill=SUBHEAD_FILL, align=LEFT, border=THIN_BORDER)
    row += 1

    write_header_row(ws, row,
                     ["Mode", "Processing details", "Processing details"],
                     [SINGLE_FILL] * 3, [BOLD_FONT] * 3)
    row += 1

    mode_labels = {
        "EL":        "NetMHCIIpan 4.1 EL",
        "BA":        "NetMHCIIpan 4.1 BA",
        "Consensus": "Consensus",
    }
    for mode in ("EL", "BA", "Consensus"):
        n   = stats["single"][mode]
        pct = stats["single_pct"][mode]
        for col, (val, fmt) in enumerate(
                [(mode_labels[mode], None), (n, NUM_FMT), (pct, PCT_FMT)], start=1):
            c = ws.cell(row=row, column=col, value=val)
            style(c, font=BODY_FONT,
                  align=CENTER if col > 1 else LEFT,
                  border=THIN_BORDER,
                  num_format=fmt if fmt else "General")
        row += 1
    row += 1

    # Implementation detail; see the repository documentation.
    ws.merge_cells(start_row=row, start_column=1, end_row=row, end_column=3)
    c = ws.cell(row=row, column=1,
                value="Processing details")
    style(c, font=SUBHEAD_FONT, fill=SUBHEAD_FILL, align=LEFT, border=THIN_BORDER)
    row += 1

    c1 = ws.cell(row=row, column=1, value="Processing details")
    c2 = ws.cell(row=row, column=2, value=n_union)
    ws.merge_cells(start_row=row, start_column=2, end_row=row, end_column=3)
    style(c1, font=BOLD_FONT, fill=UNION_FILL, align=LEFT,   border=THIN_BORDER)
    style(c2, font=BOLD_FONT, fill=UNION_FILL, align=CENTER, border=THIN_BORDER,
          num_format=NUM_FMT)
    row += 2

    # Implementation detail; see the repository documentation.
    ws.merge_cells(start_row=row, start_column=1, end_row=row, end_column=3)
    c = ws.cell(row=row, column=1,
                value="Processing details")
    style(c, font=SUBHEAD_FONT, fill=SUBHEAD_FILL, align=LEFT, border=THIN_BORDER)
    row += 1

    write_header_row(ws, row,
                     ["Processing details", "Processing details", "Processing details"],
                     [VENN_FILL] * 3, [BOLD_FONT] * 3)
    row += 1

    for label, n in stats["venn"].items():
        pct = stats["venn_pct"][label]
        for col, (val, fmt) in enumerate(
                [(label, None), (n, NUM_FMT), (pct, PCT_FMT)], start=1):
            c = ws.cell(row=row, column=col, value=val)
            style(c, font=BODY_FONT,
                  align=CENTER if col > 1 else LEFT,
                  border=THIN_BORDER,
                  num_format=fmt if fmt else "General")
        row += 1

    c1 = ws.cell(row=row, column=1, value="Processing details")
    c2 = ws.cell(row=row, column=2, value=n_union)
    c3 = ws.cell(row=row, column=3, value=1.0 if n_union else 0)
    for c, fmt in [(c1, None), (c2, NUM_FMT), (c3, PCT_FMT)]:
        style(c, font=BOLD_FONT, fill=VENN_FILL,
              align=CENTER if c != c1 else LEFT,
              border=THIN_BORDER,
              num_format=fmt if fmt else "General")


# Implementation detail; see the repository documentation.

def write_summary_sheet(ws, all_stats: dict):
    ws.column_dimensions["A"].width = 10

    modes_display = ["NetMHCIIpan 4.1 EL", "NetMHCIIpan 4.1 BA", "Consensus"]
    modes_keys    = ["EL", "BA", "Consensus"]

    venn_labels = [
        "Processing details", "Processing details", "Processing details",
        "EL + BA", "EL + Consensus", "BA + Consensus",
        "EL + BA + Consensus",
    ]

    row = 1
    total_cols = 1 + len(modes_keys) * 2 + 1 + len(venn_labels) * 2

    ws.merge_cells(start_row=row, start_column=1,
                   end_row=row, end_column=total_cols)
    c = ws.cell(row=row, column=1,
                value="Processing details")
    style(c, font=HEADER_FONT, fill=HEADER_FILL, align=CENTER, border=MED_BORDER)
    row += 1

    col = 1
    ws.merge_cells(start_row=row, start_column=col, end_row=row+1, end_column=col)
    c = ws.cell(row=row, column=col, value="Antigen")
    style(c, font=BOLD_FONT, fill=SUBHEAD_FILL, align=CENTER, border=THIN_BORDER)
    col += 1

    ws.merge_cells(start_row=row, start_column=col,
                   end_row=row, end_column=col + len(modes_keys) * 2 - 1)
    c = ws.cell(row=row, column=col,
                value="Processing details")
    style(c, font=BOLD_FONT, fill=SINGLE_FILL, align=CENTER, border=THIN_BORDER)
    col += len(modes_keys) * 2

    ws.merge_cells(start_row=row, start_column=col, end_row=row+1, end_column=col)
    c = ws.cell(row=row, column=col, value="Union")
    style(c, font=BOLD_FONT, fill=UNION_FILL, align=CENTER, border=THIN_BORDER)
    col += 1

    ws.merge_cells(start_row=row, start_column=col,
                   end_row=row, end_column=col + len(venn_labels) * 2 - 1)
    c = ws.cell(row=row, column=col,
                value="Processing details")
    style(c, font=BOLD_FONT, fill=VENN_FILL, align=CENTER, border=THIN_BORDER)
    row += 1

    col = 2
    for m in modes_display:
        for suffix in (f"{m}\nn", f"{m}\n%"):
            c = ws.cell(row=row, column=col, value=suffix)
            style(c, font=BOLD_FONT, fill=SINGLE_FILL,
                  align=CENTER, border=THIN_BORDER)
            ws.column_dimensions[get_column_letter(col)].width = 16
            col += 1
    col += 1
    ws.column_dimensions[get_column_letter(col - 1)].width = 8
    for label in venn_labels:
        for suffix in (f"{label}\nn", f"{label}\n%"):
            c = ws.cell(row=row, column=col, value=suffix)
            style(c, font=BOLD_FONT, fill=VENN_FILL,
                  align=CENTER, border=THIN_BORDER)
            ws.column_dimensions[get_column_letter(col)].width = 16
            col += 1
    ws.row_dimensions[row].height = 30
    row += 1

    for antigen in ANTIGENS:
        if antigen not in all_stats:
            continue
        stats = all_stats[antigen]
        col = 1

        c = ws.cell(row=row, column=col, value=antigen)
        style(c, font=BOLD_FONT, align=CENTER, border=THIN_BORDER)
        col += 1

        for mk in modes_keys:
            n   = stats["single"][mk]
            pct = stats["single_pct"][mk]
            for val, fmt in [(n, NUM_FMT), (pct, PCT_FMT)]:
                c = ws.cell(row=row, column=col, value=val)
                style(c, font=BODY_FONT, fill=SINGLE_FILL,
                      align=CENTER, border=THIN_BORDER, num_format=fmt)
                col += 1

        c = ws.cell(row=row, column=col, value=stats["union"])
        style(c, font=BOLD_FONT, fill=UNION_FILL,
              align=CENTER, border=THIN_BORDER, num_format=NUM_FMT)
        col += 1

        for label in venn_labels:
            n   = stats["venn"][label]
            pct = stats["venn_pct"][label]
            for val, fmt in [(n, NUM_FMT), (pct, PCT_FMT)]:
                c = ws.cell(row=row, column=col, value=val)
                style(c, font=BODY_FONT, fill=VENN_FILL,
                      align=CENTER, border=THIN_BORDER, num_format=fmt)
                col += 1

        row += 1


# Main entry point.

def main():
    print("=" * 60)
    print("IEDB II match comparison (EL / BA / Consensus)")
    print("=" * 60)

    all_stats: dict[str, dict] = {}

    for antigen in ANTIGENS:
        print(f"Processing details{antigen}")

        antigen_base = os.path.join(
            BASE, f"Processing details{antigen}", "Matches MHC II"
        )
        out_dir = os.path.join(
            BASE, f"Processing details{antigen}", "Matches MHC II comparison"
        )
        os.makedirs(out_dir, exist_ok=True)

        sets: dict[str, set] = {}
        any_missing = False

        for mode_key, (subdir, file_suffix) in MODES.items():
            path = os.path.join(
                antigen_base, subdir,
                f"RM_{antigen}_{file_suffix}.xlsx"
            )
            if not os.path.isfile(path):
                print(f"Required input or value was not found{path}")
                any_missing = True
                continue
            sets[mode_key] = read_matches(path)
            print(f"  {mode_key}: {len(sets[mode_key])}Processing details")

        if any_missing or len(sets) < 3:
            print(f"Processing details{antigen}Required input or value was not found")
            continue

        stats = compute_stats(sets)
        all_stats[antigen] = stats

        print(f"  Union: {stats['union']}Processing details")
        print(f"  EL+BA+Consensus: {stats['venn']['EL + BA + Consensus']}")

        out_path = os.path.join(out_dir, f"IEDB_II_match_comparison_{antigen}.xlsx")
        wb = Workbook()
        ws = wb.active
        ws.title = antigen
        ws.sheet_view.showGridLines = False

        write_antigen_sheet(ws, antigen, stats)

        wb.save(out_path)
        print(f"Saved output{out_path}")

    if all_stats:
        first_antigen = next(iter(all_stats))
        summary_path = os.path.join(
            BASE,
            f"Processing details{first_antigen}",
            "Matches MHC II comparison",
            "IEDB_II_match_comparison_SUMMARY.xlsx"
        )
        wb_sum = Workbook()
        ws_sum = wb_sum.active
        ws_sum.title = "Summary"
        ws_sum.sheet_view.showGridLines = False
        write_summary_sheet(ws_sum, all_stats)
        wb_sum.save(summary_path)
        print(f"Processing details{summary_path}")

    print("Completed successfully")


if __name__ == "__main__":
    main()
