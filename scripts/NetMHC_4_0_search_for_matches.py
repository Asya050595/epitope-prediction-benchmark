#!/usr/bin/env python3
"""
NetMHC 4.0 — exact allele-peptide RM/UP.

Thresholds: Affinity < 500 nM, %Rank <= 2, and Union_OR.
RM = passed prediction present in the positive reference dataset.
UP = passed prediction absent from the positive reference dataset.
No FN calculation or export is performed.

Union_OR reports the best affinity and best rank for each pair independently,
so the displayed metrics are consistent with Passed_thresholds.
"""

from project_paths import DATA_ROOT

import os
import glob
import re
import sys

try:
    import openpyxl
    from openpyxl.styles import Font, PatternFill, Alignment
except ImportError:
    sys.exit("Processing details")


# ──────────────────────────────────────────────────────────────────────────────
# Configuration.
# ──────────────────────────────────────────────────────────────────────────────

BASE = f"{DATA_ROOT}"
TOOL = "NetMHC_4.0"
ALLELE_CHECK_FATAL = True   # Allele handling.

ANTIGENS = {
    'p24': {
        "excel": os.path.join(BASE, "p24", "p24.xlsx"),
        "xls_pattern": os.path.join(
            BASE, "p24", "MHC I", "NetMHC",
            "NetMHC_p24*8-14.xls"
        ),
        "out": os.path.join(
            BASE, "p24", "Matches MHC I", f"Matches {TOOL}"
        ),
        "expected_allele_set": {'HLA-A0101', 'HLA-A0201', 'HLA-A0207', 'HLA-A1101', 'HLA-A2402', 'HLA-A2501', 'HLA-A2601', 'HLA-A2602', 'HLA-A2603', 'HLA-A3002', 'HLA-A6801', 'HLA-A6802', 'HLA-B0702', 'HLA-B0801', 'HLA-B1402', 'HLA-B1501', 'HLA-B1502', 'HLA-B1503', 'HLA-B1517', 'HLA-B2705', 'HLA-B3501', 'HLA-B3503', 'HLA-B3901', 'HLA-B4001', 'HLA-B4002', 'HLA-B4201', 'HLA-B4402', 'HLA-B4403', 'HLA-B4501', 'HLA-B4801', 'HLA-B5101', 'HLA-B5301', 'HLA-B5701', 'HLA-B5801', 'HLA-B5802', 'HLA-C0303', 'HLA-C0401', 'HLA-C0602', 'HLA-C0701', 'HLA-C0802', 'HLA-E0101'},
    },
    'pp65': {
        "excel": os.path.join(BASE, "pp65", "pp65.xlsx"),
        "xls_pattern": os.path.join(
            BASE, "pp65", "MHC I", "NetMHC",
            "NetMHC_pp65*8-14.xls"
        ),
        "out": os.path.join(
            BASE, "pp65", "Matches MHC I", f"Matches {TOOL}"
        ),
        "expected_allele_set": {'HLA-A0101', 'HLA-A0201', 'HLA-A0207', 'HLA-A1101', 'HLA-A2402', 'HLA-A2601', 'HLA-A3001', 'HLA-A3201', 'HLA-A6801', 'HLA-B0702', 'HLA-B1501', 'HLA-B3501', 'HLA-B3503', 'HLA-B4001', 'HLA-B4002', 'HLA-B4201', 'HLA-B4402', 'HLA-B4403', 'HLA-B5101', 'HLA-B5301', 'HLA-B5701', 'HLA-B5801', 'HLA-C0401', 'HLA-C1502'},
    },
    'PtxS1': {
        "excel": os.path.join(BASE, "PtxS1", "PtxS1.xlsx"),
        "xls_pattern": os.path.join(
            BASE, "PtxS1", "MHC I", "NetMHC",
            "NetMHC_PtxS1*8-14.xls"
        ),
        "out": os.path.join(
            BASE, "PtxS1", "Matches MHC I", f"Matches {TOOL}"
        ),
        "expected_allele_set": {'HLA-A0201', 'HLA-B0702', 'HLA-B3501', 'HLA-B4403', 'HLA-B5301', 'HLA-B5701'},
    },
}

# Deduplication.
THRESHOLDS = [
    ("Affinity_lt500", "Affinity(nM) < 500", "affinity_nM", "lt",  500.0, "affinity_nM"),
    ("Rank_lte2",      "%Rank <= 2",          "rank_pct",    "lte", 2.0,   "rank_pct"),
]


# ──────────────────────────────────────────────────────────────────────────────
# Allele handling.
# ──────────────────────────────────────────────────────────────────────────────

def xls_allele_to_excel(allele: str) -> str:
    """HLA-B0702 → B*07:02"""
    core = allele.replace("HLA-", "")
    m = re.match(r"^([A-Za-z]+)(\d{4,})$", core)
    if not m:
        return allele
    gene, digits = m.group(1), m.group(2)
    return f"{gene}*{digits[:2]}:{digits[2:]}"


def excel_allele_to_xls(allele: str) -> str:
    """B*07:02 → HLA-B0702"""
    return "HLA-" + allele.strip().replace("*", "").replace(":", "")


# ──────────────────────────────────────────────────────────────────────────────
# Input loading.
# ──────────────────────────────────────────────────────────────────────────────


def load_experimental_pairs(excel_path: str) -> set:
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    wb = openpyxl.load_workbook(excel_path, read_only=True, data_only=True)
    ws = wb.active
    rows = ws.iter_rows(values_only=True)
    header = next(rows, None)
    if header is None:
        wb.close(); raise ValueError(f"Processing details{excel_path}")
    headers = [str(x).strip() if x is not None else "" for x in header]
    try:
        allele_i = headers.index("HLA allele")
        peptide_i = headers.index("CTL epitopes")
    except ValueError:
        wb.close(); raise ValueError(f"Processing details{excel_path}Processing details{headers}")
    pairs=set()
    for row in rows:
        allele=row[allele_i] if allele_i < len(row) else None
        peptide=row[peptide_i] if peptide_i < len(row) else None
        if allele and peptide:
            pairs.add((str(allele).strip().replace("HLA-", ""), str(peptide).strip().upper()))
    wb.close()
    return pairs


# ──────────────────────────────────────────────────────────────────────────────
# Implementation detail; see the repository documentation.
# ──────────────────────────────────────────────────────────────────────────────

_PEPTIDE_RE = re.compile(r"^[A-Z]{5,}$")


def _validate_allele_columns(row0: list, row1: list, src: str) -> list:
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    allele_col_map = []
    for i, v in enumerate(row0):
        if not (isinstance(v, str) and v.startswith("HLA-")):
            continue
        nm_col, rank_col = i, i + 1
        h_nm   = row1[nm_col]   if nm_col   < len(row1) else ""
        h_rank = row1[rank_col] if rank_col < len(row1) else ""
        if h_nm != "nM" or h_rank != "Rank":
            print(
                f"  [!] {src}Allele status{v} (col {i}) — "
                f"Validation status"
                f"Processing details{h_nm}'/'{h_rank}Allele status"
            )
            continue
        allele_col_map.append((v.strip(), nm_col, rank_col))
    return allele_col_map


def parse_xls_predictions(xls_path: str) -> list:
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    records = []
    src = os.path.basename(xls_path)

    with open(xls_path, "r", encoding="utf-8", errors="replace") as fh:
        lines = fh.readlines()

    if len(lines) < 3:
        print(f"  [!] {src}Table status")
        return records

    row0 = lines[0].rstrip("\n").split("\t")
    row1 = lines[1].rstrip("\n").split("\t")

    allele_col_map = _validate_allele_columns(row0, row1, src)
    if not allele_col_map:
        print(f"  [!] {src}Allele status")
        return records

    for line_no, line in enumerate(lines[2:], start=3):
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 3:
            continue
        peptide_raw = cols[1].strip()
        if not peptide_raw or not _PEPTIDE_RE.match(peptide_raw.upper()):
            continue
        peptide = peptide_raw.upper()

        for allele_str, nm_col, rank_col in allele_col_map:
            if nm_col >= len(cols) or rank_col >= len(cols):
                print(
                    f"  [!] {src}Table status{line_no}Allele status{allele_str} — "
                    f"Validation status{nm_col},{rank_col}, "
                    f"Table status{len(cols)}Table status"
                )
                continue
            try:
                nm_val   = float(cols[nm_col])
                rank_val = float(cols[rank_col])
            except ValueError:
                print(
                    f"  [!] {src}Table status{line_no}Allele status{allele_str}, "
                    f"Peptide status{peptide}Processing details"
                    f"nM='{cols[nm_col]}Processing details{cols[rank_col]}Processing details"
                )
                continue

            records.append({
                "allele_xls":   allele_str,
                "allele_excel": xls_allele_to_excel(allele_str),
                "peptide":      peptide,
                "affinity_nM":  nm_val,
                "rank_pct":     rank_val,
                "source_file":  src,
            })

    return records


# ──────────────────────────────────────────────────────────────────────────────
# Threshold handling.
# ──────────────────────────────────────────────────────────────────────────────

def build_best(all_records: list, dedup_field: str) -> dict:
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    best: dict = {}
    for r in all_records:
        key = (r["allele_xls"], r["peptide"])
        if key not in best or r[dedup_field] < best[key][dedup_field]:
            best[key] = r
    return best


# ──────────────────────────────────────────────────────────────────────────────
# Validation.
# ──────────────────────────────────────────────────────────────────────────────

def check_alleles(found: set, expected: set, antigen: str) -> bool:
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    sep = "-" * 60
    missing = sorted(expected - found)
    extra   = sorted(found - expected)
    ok      = not missing and not extra

    print(f"\n  {sep}")
    print(f"Validation status{antigen}]")
    print(f"  {sep}")
    print(f"Processing details{len(found)}")
    print(f"Validation status{len(expected)}  "
          f"{'✓ OK' if len(found) == len(expected) else 'Processing details'}")

    print(f"Validation status{len(found & expected)}):")
    for a in sorted(found & expected):
        print(f"    ✓ {a}")

    if missing:
        print(f"Required input or value was not found{len(missing)}):")
        for a in missing:
            print(f"    ✗ {a}")

    if extra:
        print(f"Validation status{len(extra)}):")
        for a in extra:
            print(f"    ? {a}")

    if ok:
        print(f"Validation status")
    print(f"  {sep}")
    return ok


# ──────────────────────────────────────────────────────────────────────────────
# Threshold handling.
# ──────────────────────────────────────────────────────────────────────────────


def classify_one_threshold(best: dict, thr: tuple, exp_xls_keys: set) -> tuple[list, list]:
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    _, _, field, op, threshold_val, _ = thr

    def passes(r):
        v = r[field]
        return v < threshold_val if op == "lt" else v <= threshold_val

    rm_rows, up_rows = [], []
    for key, r in best.items():
        if not passes(r):
            continue
        row = {
            "HLA allele": r["allele_excel"],
            "Peptide": r["peptide"],
            "Affinity(nM)": r["affinity_nM"],
            "%Rank": r["rank_pct"],
            "Source File": r["source_file"],
        }
        (rm_rows if key in exp_xls_keys else up_rows).append(row)
    return rm_rows, up_rows


# ──────────────────────────────────────────────────────────────────────────────
# Output generation.
# ──────────────────────────────────────────────────────────────────────────────

HEADER_FONT  = Font(bold=True, color="FFFFFF")
HEADER_FILL  = PatternFill(fill_type="solid", fgColor="2E5F8A")
HEADER_ALIGN = Alignment(horizontal="center")

SUMMARY_HEADER_FILL = PatternFill(fill_type="solid", fgColor="4A4A4A")


def _write_data_sheet(ws, columns: list, rows: list) -> None:
    ws.append(columns)
    for cell in ws[1]:
        cell.font      = HEADER_FONT
        cell.fill      = HEADER_FILL
        cell.alignment = HEADER_ALIGN
    for row in rows:
        ws.append([row.get(c, "") for c in columns])
    for col in ws.columns:
        max_len = max(
            (len(str(cell.value)) if cell.value is not None else 0)
            for cell in col
        )
        ws.column_dimensions[col[0].column_letter].width = min(max_len + 4, 45)



def _write_summary_sheet(ws, summary_rows: list) -> None:
    """summary_rows: list of (threshold_label, rm_count, up_count)."""
    columns = ["Threshold", "RM", "UP"]
    ws.append(columns)
    for cell in ws[1]:
        cell.font = Font(bold=True, color="FFFFFF")
        cell.fill = SUMMARY_HEADER_FILL
        cell.alignment = HEADER_ALIGN
    for label, rm_count, up_count in summary_rows:
        ws.append([label, rm_count, up_count])
    for col in ws.columns:
        ws.column_dimensions[col[0].column_letter].width = 18


def save_multisheet_excel(
    out_path: str,
    sheet_data: list,         # list of (sheet_name, columns, rows)
    summary_rows: list,       # for Summary sheet
) -> None:
    wb = openpyxl.Workbook()

    # Implementation detail; see the repository documentation.
    ws_summary = wb.active
    ws_summary.title = "Summary"
    _write_summary_sheet(ws_summary, summary_rows)

    # Implementation detail; see the repository documentation.
    for sheet_name, columns, rows in sheet_data:
        ws = wb.create_sheet(title=sheet_name)
        _write_data_sheet(ws, columns, rows)

    wb.save(out_path)
    total_data = sum(len(rows) for _, _, rows in sheet_data)
    print(f"Saved output{os.path.basename(out_path)}  "
          f"({len(sheet_data)}Processing details{total_data}Table status")


# ──────────────────────────────────────────────────────────────────────────────
# Main entry point.
# ──────────────────────────────────────────────────────────────────────────────

COLS_RM_UP = ["HLA allele", "Peptide", "Affinity(nM)", "%Rank", "Source File"]
COLS_RM_UP_UNION = [
    "HLA allele", "Peptide", "Affinity(nM)", "%Rank",
    "Passed_thresholds", "Source File",
]




def process_antigen(antigen_name: str, cfg: dict) -> str:
    sep = "=" * 62
    print(f"\n{sep}Processing details{antigen_name}\n{sep}")
    if not os.path.isfile(cfg["excel"]):
        return "Required input or value was not found"

    exp_pairs = load_experimental_pairs(cfg["excel"])
    exp_xls_keys = {(excel_allele_to_xls(a), p) for a, p in exp_pairs}
    print(f"  Reference allele-peptide pairs: {len(exp_pairs)}")
    xls_files = sorted(glob.glob(cfg["xls_pattern"]))
    if not xls_files:
        return "Required input or value was not found"

    all_records, alleles_seen = [], set()
    for xls_file in xls_files:
        recs = parse_xls_predictions(xls_file)
        alleles_seen |= {r["allele_xls"] for r in recs}
        all_records.extend(recs)
    if not all_records:
        return "Processing details"

    expected_set = set(cfg.get("expected_allele_set", set()))
    alleles_ok = check_alleles(alleles_seen, expected_set, antigen_name)
    if not alleles_ok and ALLELE_CHECK_FATAL:
        return "Allele status"
    if not alleles_ok:
        all_records = [r for r in all_records if r["allele_xls"] in expected_set]

    best_nm = build_best(all_records, "affinity_nM")
    best_rank = build_best(all_records, "rank_pct")
    results = {
        "Affinity_lt500": classify_one_threshold(best_nm, THRESHOLDS[0], exp_xls_keys),
        "Rank_lte2": classify_one_threshold(best_rank, THRESHOLDS[1], exp_xls_keys),
    }

    pool_nm = {k for k, r in best_nm.items() if r["affinity_nM"] < 500.0}
    pool_rank = {k for k, r in best_rank.items() if r["rank_pct"] <= 2.0}
    rm_union, up_union = [], []
    for key in sorted(pool_nm | pool_rank):
        r_nm, r_rank = best_nm.get(key), best_rank.get(key)
        base = r_rank or r_nm
        sources = []
        for r in (r_nm, r_rank):
            if r and r["source_file"] not in sources:
                sources.append(r["source_file"])
        in_nm, in_rank = key in pool_nm, key in pool_rank
        row = {
            "HLA allele": base["allele_excel"], "Peptide": base["peptide"],
            "Affinity(nM)": r_nm["affinity_nM"] if r_nm else "",
            "%Rank": r_rank["rank_pct"] if r_rank else "",
            "Passed_thresholds": "both" if in_nm and in_rank else ("nm" if in_nm else "rank"),
            "Source File": "; ".join(sources),
        }
        (rm_union if key in exp_xls_keys else up_union).append(row)
    results["Union_OR"] = (rm_union, up_union)

    print(f"\n  {'Threshold':<18} {'RM':>8} {'UP':>8}")
    for label in ("Affinity_lt500", "Rank_lte2", "Union_OR"):
        rm, up = results[label]
        print(f"  {label:<18} {len(rm):>8} {len(up):>8}")

    os.makedirs(cfg["out"], exist_ok=True)
    summary_rows = [(label, len(rm), len(up)) for label, (rm, up) in results.items()]
    for file_type, idx in (("RM", 0), ("UP", 1)):
        sheet_data = []
        for label in ("Affinity_lt500", "Rank_lte2", "Union_OR"):
            rows = results[label][idx]
            cols = COLS_RM_UP_UNION if label == "Union_OR" else COLS_RM_UP
            sheet_data.append((label, cols, rows))
        save_multisheet_excel(os.path.join(cfg["out"], f"{file_type}_{antigen_name}_{TOOL}.xlsx"), sheet_data, summary_rows)
    return "success"


def main():
    print("Validation status")
    print("Processing details")
    print(f"Allele status{ALLELE_CHECK_FATAL}\n")

    statuses: dict = {}
    for antigen_name, cfg in ANTIGENS.items():
        statuses[antigen_name] = process_antigen(antigen_name, cfg)

    # Implementation detail; see the repository documentation.
    all_ok = all(s == "success" for s in statuses.values())
    sep = "=" * 62
    print(f"\n\n{sep}")
    if all_ok:
        print("Processing status")
    else:
        print("Error")
    print(sep)
    for antigen_name, status in statuses.items():
        icon = "✓" if status == "success" else "✗"
        print(f"  {icon}  {antigen_name}: {status}")
    print(sep)


if __name__ == "__main__":
    main()
