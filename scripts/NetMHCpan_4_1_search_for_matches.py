#!/usr/bin/env python3
"""
NetMHCpan 4.1 — exact allele-peptide RM/UP.

Thresholds: %Rank_EL <= 2, %Rank_BA <= 2, and Union_OR.
RM = passed prediction present in the positive reference dataset.
UP = passed prediction absent from the positive reference dataset.
No FN calculation or export is performed.
"""

from project_paths import DATA_ROOT

import os
import sys
from pathlib import Path

import openpyxl
from openpyxl.styles import Font, PatternFill, Alignment


# ─────────────────────────────────────────────────────────────────────────────
# Configuration.
# ─────────────────────────────────────────────────────────────────────────────

BASE = f"{DATA_ROOT}"
TOOL = "NetMHCpan_4.1"

ANTIGENS = {
    "p24": {
        "xls_parts": [
            f"{BASE}/p24/MHC I/NetMHCpan/NetMHCpan_p24_part_1_8-14.xls",
            f"{BASE}/p24/MHC I/NetMHCpan/NetMHCpan_p24_part_2_8-14.xls",
            f"{BASE}/p24/MHC I/NetMHCpan/NetMHCpan_p24_part_3_8-14.xls",
            f"{BASE}/p24/MHC I/NetMHCpan/NetMHCpan_p24_part_4_8-14.xls",
        ],
        "xlsx": f"{BASE}/p24/p24.xlsx",
        "out":  f"{BASE}/p24/Matches MHC I/Matches NetMHCpan_4.1",
    },
    "pp65": {
        "xls_parts": [
            f"{BASE}/pp65/MHC I/NetMHCpan/NetMHCpan_pp65_part_1_8-14.xls",
            f"{BASE}/pp65/MHC I/NetMHCpan/NetMHCpan_pp65_part_2_8-14.xls",
        ],
        "xlsx": f"{BASE}/pp65/pp65.xlsx",
        "out":  f"{BASE}/pp65/Matches MHC I/Matches NetMHCpan_4.1",
    },
    "PtxS1": {
        "xls_parts": [
            f"{BASE}/PtxS1/MHC I/NetMHCpan/NetMHCpan_PtxS1_8-14.xls",
        ],
        "xlsx": f"{BASE}/PtxS1/PtxS1.xlsx",
        "out":  f"{BASE}/PtxS1/Matches MHC I/Matches NetMHCpan_4.1",
    },
}

# Normalization.
EXPECTED_ALLELES: dict[str, list[str]] = {
    'p24': [
        'A01:01', 'A02:01', 'A02:07', 'A11:01', 'A11:03', 'A24:02',
        'A25:01', 'A26:01', 'A26:02', 'A26:03', 'A26:38', 'A30:02',
        'A33:03', 'A68:01', 'A68:02', 'A74:01', 'B07:02', 'B08:01',
        'B13:02', 'B14:01', 'B14:02', 'B14:03', 'B15:01', 'B15:02',
        'B15:03', 'B15:10', 'B15:16', 'B15:17', 'B15:24', 'B15:40',
        'B27:05', 'B35:01', 'B35:02', 'B35:03', 'B35:05', 'B35:08',
        'B39:01', 'B39:10', 'B40:01', 'B40:02', 'B40:06', 'B42:01',
        'B42:02', 'B44:02', 'B44:03', 'B44:15', 'B45:01', 'B48:01',
        'B50:01', 'B51:01', 'B53:01', 'B57:01', 'B57:02', 'B57:03',
        'B58:01', 'B58:02', 'B67:01', 'C01:02', 'C03:03', 'C03:04',
        'C04:01', 'C06:02', 'C07:01', 'C08:02', 'C18:01', 'E01:01',
        'E01:03',
    ],
    'pp65': [
        'A01:01', 'A02:01', 'A02:07', 'A11:01', 'A24:02', 'A24:07',
        'A26:01', 'A30:01', 'A32:01', 'A68:01', 'B07:02', 'B15:01',
        'B35:01', 'B35:02', 'B35:03', 'B35:08', 'B35:11', 'B40:01',
        'B40:02', 'B40:06', 'B42:01', 'B42:02', 'B44:02', 'B44:03',
        'B51:01', 'B52:01', 'B53:01', 'B57:01', 'B58:01', 'C01:02',
        'C04:01', 'C08:01', 'C12:02', 'C15:02',
    ],
    'PtxS1': [
        'A02:01', 'B07:02', 'B35:01', 'B44:03', 'B53:01', 'B57:01',
    ],
}

_COLS_PER_ALLELE = 6   # core icore EL-score EL_Rank BA-score BA_Rank
_FIXED_COLS      = 3   # Pos Peptide ID
_TAIL_COLS       = 2   # Ave NB

# Threshold handling.
THRESHOLD_SETS = [
    ("Rank_EL_lte2", [("rank_el", "<=", 2.0)]),
    ("Rank_BA_lte2", [("rank_ba", "<=", 2.0)]),
    ("Union_OR",     [("rank_el", "<=", 2.0), ("rank_ba", "<=", 2.0)]),
]


# ─────────────────────────────────────────────────────────────────────────────
# Helper functions.
# ─────────────────────────────────────────────────────────────────────────────

def normalize_allele(raw: str) -> str:
    """HLA-A*01:01 / A*01:01 / HLA-A01:01  →  A01:01"""
    return raw.strip().replace("HLA-", "").replace("*", "")

def format_allele_for_output(allele: str) -> str:
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    allele = normalize_allele(allele)

    if len(allele) >= 2:
        return f"{allele[0]}*{allele[1:]}"

    return allele

# ─────────────────────────────────────────────────────────────────────────────
# Implementation detail; see the repository documentation.
# ─────────────────────────────────────────────────────────────────────────────


def parse_experimental_epitopes(xlsx_path: str) -> set[tuple[str, str]]:
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    wb=openpyxl.load_workbook(xlsx_path,read_only=True,data_only=True); ws=wb.active
    rows=ws.iter_rows(values_only=True); header=next(rows,None)
    if header is None: wb.close(); raise ValueError(f"Processing details{xlsx_path}")
    headers=[str(x).strip() if x is not None else "" for x in header]
    try: allele_i=headers.index("HLA allele"); peptide_i=headers.index("CTL epitopes")
    except ValueError:
        wb.close(); raise ValueError(f"Processing details{xlsx_path}Processing details{headers}")
    pairs=set()
    for row in rows:
        a=row[allele_i] if allele_i<len(row) else None; pep=row[peptide_i] if peptide_i<len(row) else None
        if a and pep: pairs.add((normalize_allele(str(a)),str(pep).strip().upper()))
    wb.close(); return pairs


# ─────────────────────────────────────────────────────────────────────────────
# Implementation detail; see the repository documentation.
# ─────────────────────────────────────────────────────────────────────────────

def _validate_xls_structure(
    lines: list[str], path: str
) -> tuple[list[str], list[int]] | None:
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    fname = os.path.basename(path)

    if len(lines) < 3:
        print(f"Error{fname}Table status")
        return None

    # Allele handling.
    allele_row  = lines[0].rstrip("\n").split("\t")
    raw_alleles = [tok.strip() for tok in allele_row if tok.strip()]
    col_starts  = [i for i, tok in enumerate(allele_row) if tok.strip()]

    if not raw_alleles:
        print(f"Error{fname}Allele status")
        return None

    # Deduplication.
    alleles_norm_list = [normalize_allele(a) for a in raw_alleles]
    if len(alleles_norm_list) != len(set(alleles_norm_list)):
        from collections import Counter
        counts = Counter(alleles_norm_list)
        dups = [a for a, c in counts.items() if c > 1]
        print(f"Error{fname}Allele status{dups}")
        return None

    n = len(raw_alleles)

    # Implementation detail; see the repository documentation.
    col_row = lines[1].rstrip("\n").split("\t")
    expected_total = _FIXED_COLS + _COLS_PER_ALLELE * n + _TAIL_COLS
    if len(col_row) != expected_total:
        print(
            f"Error{fname}Validation status{expected_total}Table status"
            f"({_FIXED_COLS}+{_COLS_PER_ALLELE}×{n}+{_TAIL_COLS}), "
            f"Processing details{len(col_row)}."
        )
        return None

    per_allele_expected = ["core", "icore", "EL-score", "EL_Rank", "BA-score", "BA_Rank"]
    for allele_raw, cs in zip(raw_alleles, col_starts):
        for offset, exp_name in enumerate(per_allele_expected):
            idx = cs + offset
            if idx >= len(col_row):
                print(
                    f"Error{fname}Allele status{allele_raw!r} — "
                    f"Table status{idx}Processing details"
                )
                return None
            actual = col_row[idx].strip()
            if actual != exp_name:
                print(
                    f"Error{fname}Allele status{allele_raw!r} (+{offset}) — "
                    f"Validation status{exp_name}Processing details{actual}'."
                )
                return None

    return alleles_norm_list, col_starts


def parse_netmhcpan_xls(
    xls_paths: list[str],
) -> tuple[list[dict], dict[str, set[str]], list[str]]:
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    records:          list[dict]           = []
    alleles_per_file: dict[str, set[str]]  = {}
    parse_errors:     list[str]            = []

    # Validation.
    # Deduplication.
    reference_signature: list[tuple[str, str, str]] | None = None
    reference_fname:     str                               = ""

    for xls_path in xls_paths:
        fname = os.path.basename(xls_path)
        with open(xls_path, "r", encoding="utf-8", errors="replace") as fh:
            lines = fh.readlines()

        result = _validate_xls_structure(lines, xls_path)
        if result is None:
            parse_errors.append(f"{fname}Error")
            continue

        alleles_norm, col_starts = result
        alleles_per_file[fname]  = set(alleles_norm)

        # Validation.
        current_signature: list[tuple[str, str, str]] = []
        signature_seen:    set[tuple[str, str, str]]  = set()
        for lineno, line in enumerate(lines[2:], start=3):
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 3:
                parse_errors.append(
                    f"{fname}:{lineno}Table status"
                )
                continue
            if not cols[1].strip():
                continue
            sig_key = (cols[0].strip(), cols[1].strip(), cols[2].strip())
            if sig_key in signature_seen:
                parse_errors.append(
                    f"{fname}:{lineno}Table status{sig_key}"
                )
                continue
            signature_seen.add(sig_key)
            current_signature.append(sig_key)

        if reference_signature is None:
            reference_signature = current_signature
            reference_fname     = fname
        else:
            if current_signature != reference_signature:
                parse_errors.append(
                    f"{fname}Validation status"
                    f"Processing details{reference_fname!r}Processing details"
                )

        # Implementation detail; see the repository documentation.
        for lineno, line in enumerate(lines[2:], start=3):
            line = line.rstrip("\n")
            if not line.strip():
                continue
            cols = line.split("\t")

            # Peptide handling.
            peptide_raw = cols[1].strip() if len(cols) > 1 else ""
            if not peptide_raw:
                parse_errors.append(f"{fname}:{lineno}Required input or value was not found")
                continue
            if not peptide_raw.replace("-", "").isalpha():
                parse_errors.append(
                    f"{fname}:{lineno}Processing details{peptide_raw!r}"
                )
                continue
            peptide = peptide_raw.upper()

            # Main entry point.
            if not (8 <= len(peptide) <= 14):
                parse_errors.append(
                    f"{fname}:{lineno}Peptide status{len(peptide)}Processing details{peptide!r}"
                )
                continue

            # Allele handling.
            min_len = col_starts[-1] + _COLS_PER_ALLELE if col_starts else 0
            if len(cols) < min_len:
                parse_errors.append(
                    f"{fname}:{lineno}Validation status"
                    f"({len(cols)} < {min_len}Table status"
                )
                continue

            for allele_norm, cs in zip(alleles_norm, col_starts):
                el_rank_idx = cs + 3
                ba_rank_idx = cs + 5

                if el_rank_idx >= len(cols) or ba_rank_idx >= len(cols):
                    parse_errors.append(
                        f"{fname}:{lineno}Allele status{allele_norm}: "
                        f"Table status"
                    )
                    continue

                el_str = cols[el_rank_idx].strip()
                ba_str = cols[ba_rank_idx].strip()
                try:
                    rank_el = float(el_str)
                    rank_ba = float(ba_str)
                except ValueError:
                    parse_errors.append(
                        f"{fname}:{lineno}Allele status{allele_norm}: "
                        f"Processing details{el_str!r}, BA={ba_str!r})"
                    )
                    continue

                if not (0.0 <= rank_el <= 100.0) or not (0.0 <= rank_ba <= 100.0):
                    parse_errors.append(
                        f"{fname}:{lineno}Allele status{allele_norm}: "
                        f"Processing details{rank_el}, BA={rank_ba}"
                    )
                    continue

                records.append({
                    "allele":      allele_norm,
                    "peptide":     peptide,
                    "rank_el":     rank_el,
                    "rank_ba":     rank_ba,
                    "source_file": fname,
                })

    return records, alleles_per_file, parse_errors


# ─────────────────────────────────────────────────────────────────────────────
# Validation.
# ─────────────────────────────────────────────────────────────────────────────

def validate_alleles_xls(
    alleles_per_file: dict[str, set[str]],
    expected: list[str],
    antigen: str,
) -> bool:
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    print(f"\n  {'─'*58}")
    print(f"Validation status{antigen}]:")
    print(f"Validation status{len(expected)}Allele status")

    ok = True

    # Deduplication.
    from collections import Counter
    all_seen: list[str] = []
    for alleles in alleles_per_file.values():
        all_seen.extend(alleles)
    counts = Counter(all_seen)
    duplicates = {a: [] for a in counts if counts[a] > 1}
    for fname, alleles in alleles_per_file.items():
        for a in alleles:
            if a in duplicates:
                duplicates[a].append(fname)
    if duplicates:
        print(f"Error")
        for a, fnames in sorted(duplicates.items()):
            print(f"    {a}: {', '.join(fnames)}")
        ok = False

    found_all = set().union(*alleles_per_file.values()) if alleles_per_file else set()
    expected_set = set(expected)
    missing = sorted(expected_set - found_all)
    extra   = sorted(found_all - expected_set)
    matched = sorted(expected_set & found_all)

    print(f"Processing details{len(found_all)}Allele status")
    if not missing and not extra:
        print(f"Processing details{len(matched)}Validation status")
    else:
        print(f"Validation status{len(matched)}")
        if missing:
            print(f"Required input or value was not found{len(missing)}):")
            for a in missing:
                print(f"      - {a}")
        if extra:
            print(f"Processing details{len(extra)}):")
            for a in extra:
                print(f"      + {a}")
        ok = False

    print(f"  {'─'*58}")
    return ok


def validate_alleles_experimental(
    pairs: set[tuple[str, str]],
    expected: list[str],
    antigen: str,
) -> bool:
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    print(f"\n  {'─'*58}")
    print(f"Validation status{antigen}]:")

    found = {allele for allele, _ in pairs}
    expected_set = set(expected)
    missing = sorted(expected_set - found)
    extra   = sorted(found - expected_set)
    matched = sorted(expected_set & found)

    print(f"Validation status{len(expected)}Allele status")
    print(f"Processing details{len(found)}Allele status")

    if not missing and not extra:
        print(f"Processing details{len(matched)}Validation status")
        ok = True
    else:
        print(f"Validation status{len(matched)}")
        if missing:
            print(f"Required input or value was not found{len(missing)}):")
            for a in missing:
                print(f"      - {a}")
        if extra:
            print(f"Processing details{len(extra)}):")
            for a in extra:
                print(f"      + {a}")
        ok = False

    print(f"  {'─'*58}")
    return ok


# ─────────────────────────────────────────────────────────────────────────────
# Output generation.
# ─────────────────────────────────────────────────────────────────────────────

_HEADER_FILL = PatternFill("solid", fgColor="D9D9D9")
_BOLD        = Font(bold=True)
_CENTER      = Alignment(horizontal="center")


RESULT_COLS = ["Allele", "Peptide", "%Rank_EL", "%Rank_BA", "Source File"]

_HEADER_FILL = PatternFill("solid", fgColor="1F4E78")
_HEADER_FONT = Font(color="FFFFFF", bold=True)
_CENTER = Alignment(horizontal="center")


def _write_header(ws, columns: list[str]) -> None:
    ws.append(columns)
    for cell in ws[1]:
        cell.font = _HEADER_FONT
        cell.fill = _HEADER_FILL
        cell.alignment = _CENTER


def _autofit(ws) -> None:
    for col in ws.columns:
        width = max((len(str(cell.value)) for cell in col if cell.value is not None), default=8)
        ws.column_dimensions[col[0].column_letter].width = min(width + 4, 45)


def _add_result_sheet(wb: openpyxl.Workbook, sheet_name: str,
                      records: list[dict], first: bool) -> None:
    ws = wb.active if first else wb.create_sheet(sheet_name)
    if first:
        ws.title = sheet_name
    _write_header(ws, RESULT_COLS)
    for r in sorted(records, key=lambda x: (x["allele"], x["peptide"])):
        ws.append([
            format_allele_for_output(r["allele"]), r["peptide"],
            r["rank_el"], r["rank_ba"], r.get("source_file", ""),
        ])
    _autofit(ws)


def save_result_workbook(path: str, results_by_threshold: list, result_index: int) -> None:
    """result_index: 1 for RM, 2 for UP in (sheet_name, RM, UP)."""
    wb = openpyxl.Workbook()
    for i, item in enumerate(results_by_threshold):
        sheet_name = item[0]
        records = item[result_index]
        _add_result_sheet(wb, sheet_name, records, i == 0)
    wb.save(path)
    last_ws = wb.worksheets[-1]
    print(f"Saved output{os.path.basename(path)} ({last_ws.max_row - 1}Table status")


# ─────────────────────────────────────────────────────────────────────────────
# Implementation detail; see the repository documentation.
# ─────────────────────────────────────────────────────────────────────────────



def process_antigen(name: str, cfg: dict) -> None:
    xls_parts = cfg["xls_parts"]
    xlsx_path = cfg["xlsx"]
    out_dir = Path(cfg["out"])
    print(f"\n{'='*62}Processing details{name}\n{'='*62}")

    missing = [p for p in xls_parts if not Path(p).exists()]
    if not Path(xlsx_path).exists():
        missing.append(xlsx_path)
    if missing:
        raise FileNotFoundError("Required input or value was not found" + "\n  ".join(missing))

    expected = EXPECTED_ALLELES[name]
    experimental = parse_experimental_epitopes(xlsx_path)
    print(f"  Reference allele-peptide pairs: {len(experimental)}")
    all_records, alleles_per_file, parse_errors = parse_netmhcpan_xls(xls_parts)
    if parse_errors:
        raise ValueError("Error" + "\n  ".join(parse_errors))
    if not all_records:
        raise ValueError(f"Processing details{name}")

    if not validate_alleles_xls(alleles_per_file, expected, name):
        raise ValueError(f"Validation status{name}")

    results_by_threshold = []
    ops_map = {"<=": lambda v, x: v <= x, "<": lambda v, x: v < x}
    print(f"\n  {'Mode':<18}{'RM':>8}{'UP':>10}")
    for set_name, thresholds in THRESHOLD_SETS:
        def record_passes(r: dict) -> bool:
            return any(ops_map[op](r[field], val) for field, op, val in thresholds)
        passing = [r for r in all_records if record_passes(r)]
        best_by_pair = {}
        for r in passing:
            key = (r["allele"], r["peptide"])
            if key not in best_by_pair:
                best_by_pair[key] = r; continue
            cur = best_by_pair[key]
            if set_name == "Rank_EL_lte2":
                better = r["rank_el"] < cur["rank_el"]
            elif set_name == "Rank_BA_lte2":
                better = r["rank_ba"] < cur["rank_ba"]
            else:
                better = min(r["rank_el"], r["rank_ba"]) < min(cur["rank_el"], cur["rank_ba"])
            if better: best_by_pair[key] = r
        binders = list(best_by_pair.values())
        rm = [r for r in binders if (r["allele"], r["peptide"]) in experimental]
        up = [r for r in binders if (r["allele"], r["peptide"]) not in experimental]
        print(f"  {set_name:<18}{len(rm):>8}{len(up):>10}")
        results_by_threshold.append((set_name, rm, up))

    out_dir.mkdir(parents=True, exist_ok=True)
    save_result_workbook(str(out_dir / f"RM_{name}_{TOOL}.xlsx"), results_by_threshold, 1)
    save_result_workbook(str(out_dir / f"UP_{name}_{TOOL}.xlsx"), results_by_threshold, 2)


# ─────────────────────────────────────────────────────────────────────────────
# Implementation detail; see the repository documentation.
# ─────────────────────────────────────────────────────────────────────────────

def main() -> None:
    print("NetMHCpan 4.1 — Search for Matches")
    print("Processing details", ", ".join(ANTIGENS))

    failed_antigens: list[str] = []

    for antigen_name, config in ANTIGENS.items():
        try:
            process_antigen(antigen_name, config)
        except (FileNotFoundError, ValueError) as error:
            print(f"Error{antigen_name}: {error}")
            failed_antigens.append(antigen_name)

    print("\n" + "=" * 62)
    if failed_antigens:
        print(f"Error{', '.join(failed_antigens)}")
        sys.exit(1)
    else:
        print("Completed successfully")
    print("=" * 62)


if __name__ == "__main__":
    main()
