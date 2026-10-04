#!/usr/bin/env python3
"""English documentation for this module or helper is provided in the repository README and in the surrounding code."""

from project_paths import DATA_ROOT

import os
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter

# Threshold handling.
THRESH_EL = 2.0   # Threshold handling.
THRESH_BA = 2.0   # Threshold handling.

ANTIGENS = ["p24", "pp65", "PtxS1"]

SHEET_EL    = "Rank_EL_lte2"
SHEET_BA    = "Rank_BA_lte2"
SHEET_UNION = "Union_OR"

BASE_INPUT  = f"{DATA_ROOT}/{{ag}}/Matches MHC I/Matches NetMHCpan_4.1"
BASE_OUTPUT = f"{DATA_ROOT}/{{ag}}/Matches MHC I comparison"
INPUT_NAME  = "RM_{ag}_NetMHCpan_4.1.xlsx"
OUTPUT_NAME = "NetMHCpan_4.1_match_comparison.xlsx"

# Formatting.
HEADER_FILL  = PatternFill("solid", fgColor="1F4E79")
SUBHDR_FILL  = PatternFill("solid", fgColor="2E75B6")
SECTION_FILL = PatternFill("solid", fgColor="D6E4F0")
TOTAL_FILL   = PatternFill("solid", fgColor="E2EFDA")
WHITE_FILL   = PatternFill("solid", fgColor="FFFFFF")
ALT_FILL     = PatternFill("solid", fgColor="F2F7FC")

HEADER_FONT  = Font(bold=True, color="FFFFFF", name="Calibri", size=11)
SUBHDR_FONT  = Font(bold=True, color="FFFFFF", name="Calibri", size=10)
BODY_FONT    = Font(name="Calibri", size=10)
BOLD_FONT    = Font(bold=True, name="Calibri", size=10)
SECTION_FONT = Font(bold=True, color="1F4E79", name="Calibri", size=10)

CENTER = Alignment(horizontal="center", vertical="center", wrap_text=True)
LEFT   = Alignment(horizontal="left",   vertical="center", wrap_text=True)


def thin_border():
    s = Side(style="thin", color="BDD7EE")
    return Border(left=s, right=s, top=s, bottom=s)


def style_cell(cell, fill=None, font=None, alignment=None):
    if fill:      cell.fill      = fill
    if font:      cell.font      = font
    if alignment: cell.alignment = alignment
    cell.border = thin_border()


# Input loading.
def load_rm(ag):
    path = os.path.join(
        BASE_INPUT.replace("{ag}", ag),
        INPUT_NAME.replace("{ag}", ag)
    )
    print(f"Processing details{path}")
    if not os.path.exists(path):
        print(f"Required input or value was not found{path}")
        return None

    sheets = {}
    for sheet_name in (SHEET_EL, SHEET_BA, SHEET_UNION):
        df = pd.read_excel(path, sheet_name=sheet_name)
        sheets[sheet_name] = df
        print(f"Processing details{sheet_name}»: {len(df)}Processing details")
    return sheets


# Threshold handling.
def analyse(sheets):
    df_el  = sheets[SHEET_EL]
    df_ba  = sheets[SHEET_BA]
    df_uni = sheets[SHEET_UNION]

    pairs = lambda df: set(zip(df["Allele"], df["Peptide"]))

    pass_el  = pairs(df_el)
    pass_ba  = pairs(df_ba)
    pass_any = pairs(df_uni)

    # Implementation detail; see the repository documentation.
    # Deduplication.
    # Implementation detail; see the repository documentation.
    n_total = len(pass_any)

    if n_total == 0:
        raise ValueError(
            "Processing details"
        )

    # Implementation detail; see the repository documentation.
    # Implementation detail; see the repository documentation.
    # Implementation detail; see the repository documentation.
    # Implementation detail; see the repository documentation.
    # Implementation detail; see the repository documentation.
    computed_union = pass_el | pass_ba
    if computed_union != pass_any:
        only_in_computed = computed_union - pass_any
        only_in_union = pass_any - computed_union
        raise ValueError(
            f"Validation status"
            f"Required input or value was not found{len(only_in_computed)}; "
            f"Processing details{len(only_in_union)}."
        )

    pct = lambda k: round(k / n_total * 100, 1) if n_total else 0.0

    n_both = len(pass_el & pass_ba)
    n_el_only = len(pass_el - pass_ba)
    n_ba_only = len(pass_ba - pass_el)
    n_single = n_el_only + n_ba_only
    n_two_plus = n_both
    n_any = len(pass_any)

    return {
        "n_total":    n_total,              "pct_total":    100.0,
        "n_single":   n_single,            "pct_single":   pct(n_single),
        "n_two_plus": n_two_plus,          "pct_two_plus": pct(n_two_plus),
        "n_any":      n_any,               "pct_any":      pct(n_any),
        "n_pass_el":  len(pass_el),        "pct_pass_el":  pct(len(pass_el)),
        "n_el_only":  n_el_only,           "pct_el_only":  pct(n_el_only),
        "n_pass_ba":  len(pass_ba),        "pct_pass_ba":  pct(len(pass_ba)),
        "n_ba_only":  n_ba_only,           "pct_ba_only":  pct(n_ba_only),
        "n_both":     n_both,              "pct_both":     pct(n_both),
    }


# Output generation.
def write_excel(results, out_path):
    wb = Workbook()
    ws = wb.active
    ws.title = "Threshold comparison"

    # Output generation.
    ags = [ag for ag in ANTIGENS if ag in results]
    n_ags = len(ags)

    # Implementation detail; see the repository documentation.
    ws.column_dimensions["A"].width = 42
    for i in range(n_ags * 2):
        ws.column_dimensions[get_column_letter(2 + i)].width = 10

    # Implementation detail; see the repository documentation.
    last_col = get_column_letter(1 + n_ags * 2)

    row = 1

    # Implementation detail; see the repository documentation.
    ws.merge_cells(f"A{row}:{last_col}{row}")
    c = ws.cell(row, 1, "NetMHCpan 4.1 match comparison across thresholds")
    style_cell(c, fill=HEADER_FILL,
               font=Font(bold=True, color="FFFFFF", name="Calibri", size=13),
               alignment=CENTER)
    ws.row_dimensions[row].height = 24
    row += 1

    # Threshold handling.
    ws.merge_cells(f"A{row}:{last_col}{row}")
    c = ws.cell(row, 1,
        f"Processing details{THRESH_EL}  |  BA %Rank <= {THRESH_BA}")
    style_cell(c, fill=SUBHDR_FILL,
               font=Font(italic=True, color="FFFFFF", name="Calibri", size=10),
               alignment=CENTER)
    row += 1

    # Implementation detail; see the repository documentation.
    hdr1 = row
    hdr2 = row + 1

    # Implementation detail; see the repository documentation.
    ws.merge_cells(f"A{hdr1}:A{hdr2}")
    c = ws.cell(hdr1, 1, "Processing details")
    style_cell(c, fill=SUBHDR_FILL, font=SUBHDR_FONT, alignment=CENTER)

    # Implementation detail; see the repository documentation.
    for col_i, ag in enumerate(ags):
        ag_col = 2 + col_i * 2
        ws.merge_cells(
            f"{get_column_letter(ag_col)}{hdr1}:{get_column_letter(ag_col+1)}{hdr1}"
        )
        c = ws.cell(hdr1, ag_col, ag)
        style_cell(c, fill=SUBHDR_FILL, font=SUBHDR_FONT, alignment=CENTER)

    # Implementation detail; see the repository documentation.
    for col_i in range(n_ags):
        ag_col = 2 + col_i * 2
        for lbl, ci in [("n", ag_col), ("%", ag_col + 1)]:
            c = ws.cell(hdr2, ci, lbl)
            style_cell(c, fill=SUBHDR_FILL, font=SUBHDR_FONT, alignment=CENTER)

    row = hdr2 + 1

    # Helper functions.
    def data_row(label, key_n, key_pct, section=False, total=False, alt=False):
        nonlocal row
        fill = SECTION_FILL if section else (TOTAL_FILL if total else
                                             (ALT_FILL if alt else WHITE_FILL))
        font = SECTION_FONT if section else (BOLD_FONT if total else BODY_FONT)
        c = ws.cell(row, 1, label)
        style_cell(c, fill=fill, font=font, alignment=LEFT)
        for col_i, ag in enumerate(ags):
            ag_col = 2 + col_i * 2
            r = results.get(ag, {})
            vn = r.get(key_n,  "—")
            vp = r.get(key_pct, "—")
            cn = ws.cell(row, ag_col,     vn)
            cp = ws.cell(row, ag_col + 1, f"{vp}%" if isinstance(vp, float) else vp)
            style_cell(cn, fill=fill, font=font, alignment=CENTER)
            style_cell(cp, fill=fill, font=font, alignment=CENTER)
        row += 1

    def section_header(label):
        nonlocal row
        ws.merge_cells(f"A{row}:{last_col}{row}")
        c = ws.cell(row, 1, label)
        style_cell(c, fill=SECTION_FILL, font=SECTION_FONT, alignment=LEFT)
        row += 1

    # Implementation detail; see the repository documentation.
    data_row("Processing details", "n_total", "pct_total", total=True)

    # Implementation detail; see the repository documentation.
    section_header("Processing details")
    data_row("Processing details", "n_single", "pct_single", alt=True)

    # Threshold handling.
    section_header("Processing details")
    data_row("Processing details", "n_two_plus", "pct_two_plus", alt=True)

    # Implementation detail; see the repository documentation.
    section_header("Processing details")
    data_row("Processing details", "n_any", "pct_any", total=True)

    # Implementation detail; see the repository documentation.
    section_header("Processing details")
    data_row("   NetMHCpan 4.1 EL  (%Rank_EL <= 2.0)",
             "n_pass_el", "pct_pass_el")
    data_row("Processing details",
             "n_el_only", "pct_el_only", alt=True)
    data_row("   NetMHCpan 4.1 BA  (%Rank_BA <= 2.0)",
             "n_pass_ba", "pct_pass_ba")
    data_row("Processing details",
             "n_ba_only", "pct_ba_only", alt=True)
    data_row("Processing details",
             "n_both", "pct_both", total=True)

    # Implementation detail; see the repository documentation.
    row += 1
    ws.merge_cells(f"A{row}:{last_col}{row}")
    c = ws.cell(row, 1, "Processing details")
    style_cell(c, fill=HEADER_FILL,
               font=Font(bold=True, color="FFFFFF", name="Calibri", size=10),
               alignment=LEFT)
    row += 1
    notes = [
        "Processing details",
        "Processing details",
        "Processing details",
        "Processing details",
        "Required input or value was not found",
        "Required input or value was not found",
        "Processing details",
    ]
    for note in notes:
        ws.merge_cells(f"A{row}:{last_col}{row}")
        c = ws.cell(row, 1, f"* {note}")
        style_cell(c, fill=ALT_FILL,
                   font=Font(italic=True, name="Calibri", size=9, color="1F4E79"),
                   alignment=LEFT)
        row += 1

    # Implementation detail; see the repository documentation.
    for r in range(4, row):
        ws.row_dimensions[r].height = 18

    wb.save(out_path)
    print(f"Saved output{out_path}")


# ──────────────────────────── main ────────────────────────────────────────────
def main():
    print("=" * 60)
    print("Processing details")
    print("=" * 60)

    all_results = {}

    for ag in ANTIGENS:
        print(f"\n[{ag}Loaded input")
        sheets = load_rm(ag)
        if sheets is None:
            print(f"Processing details{ag}Required input or value was not found")
            continue
        try:
            result = analyse(sheets)
        except ValueError as e:
            print(f"Error{e}")
            print(f"Processing details{ag}Validation status")
            continue

        all_results[ag] = result

        print(f"Processing details{result['n_total']}")
        print(f"Processing details{result['n_single']} ({result['pct_single']}%)")
        print(f"Processing details{result['n_two_plus']} ({result['pct_two_plus']}%)")
        print(f"Processing details{result['n_any']} ({result['pct_any']}%)")
        print(f"Processing details{result['n_pass_el']} ({result['pct_pass_el']}%)")
        print(f"Processing details{result['n_pass_ba']} ({result['pct_pass_ba']}%)")
        print(f"Processing details{result['n_both']} ({result['pct_both']}%)")

        out_dir = BASE_OUTPUT.replace("{ag}", ag)
        os.makedirs(out_dir, exist_ok=True)
        out_path = os.path.join(out_dir, OUTPUT_NAME)
        write_excel({ag: result}, out_path)

    # Implementation detail; see the repository documentation.
    if all_results:
        print(f"Saved output")
        # Output generation.
        first_ag = next(ag for ag in ANTIGENS if ag in all_results)
        summary_dir = BASE_OUTPUT.replace("{ag}", first_ag)
        os.makedirs(summary_dir, exist_ok=True)
        summary_path = os.path.join(summary_dir, "NetMHCpan_4.1_match_comparison_ALL.xlsx")
        write_excel(all_results, summary_path)

    print(f"\n{'='*60}")
    print("Completed successfully")


if __name__ == "__main__":
    main()
