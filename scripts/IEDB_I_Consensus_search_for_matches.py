#!/usr/bin/env python3
"""English documentation for this module or helper is provided in the repository README and in the surrounding code."""

from project_paths import DATA_ROOT

import os
import sys
import pandas as pd
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter
import openpyxl

# ---------------------------------------------------------------------------
# Configuration.
# ---------------------------------------------------------------------------
BASE = f"{DATA_ROOT}"
TOOL = "IEDB_I_Consensus"

ANTIGENS = {
    "p24": {
        "tsv":  f"{BASE}/p24/MHC I/IEDB I/IEDB_I_Consensus_p24_8-14.tsv",
        "xlsx": f"{BASE}/p24/p24.xlsx",
        "out":  f"{BASE}/p24/Matches MHC I/Matches {TOOL}",
    },
    "pp65": {
        "tsv":  f"{BASE}/pp65/MHC I/IEDB I/IEDB_I_Consensus_pp65_8-14.tsv",
        "xlsx": f"{BASE}/pp65/pp65.xlsx",
        "out":  f"{BASE}/pp65/Matches MHC I/Matches {TOOL}",
    },
    "PtxS1": {
        "tsv":  f"{BASE}/PtxS1/MHC I/IEDB I/IEDB_I_Consensus_PtxS1_8-14.tsv",
        "xlsx": f"{BASE}/PtxS1/PtxS1.xlsx",
        "out":  f"{BASE}/PtxS1/Matches MHC I/Matches {TOOL}",
    },
}

# Threshold handling.
THRESHOLDS = [
    ("median binding percentile", "<=", 1.0,   "Median_lte1"),
    ("consensus percentile",      "<=", 1.0,   "Consensus_lte1"),
    ("smm ic50",                  "<",  500.0, "SMM_IC50_lt500"),
    ("ann ic50",                  "<",  500.0, "ANN_IC50_lt500"),
]

# Implementation detail; see the repository documentation.
REQUIRED_TSV_COLS = {
    "allele", "peptide", "peptide length",
    "median binding percentile", "consensus percentile",
    "smm ic50", "ann ic50",
}

# Error handling.
ALLOWED_NAN_COLS = {"smm ic50"}

EXPECTED_ALLELES = {
    'p24': {
        'A*01:01', 'A*02:01', 'A*11:01', 'A*24:02', 'A*26:01', 'A*30:02',
        'A*68:01', 'A*68:02', 'B*07:02', 'B*08:01', 'B*15:01', 'B*35:01',
        'B*40:01', 'B*44:02', 'B*44:03', 'B*51:01', 'B*53:01', 'B*57:01',
        'B*58:01',
    },
    'pp65': {
        'A*01:01', 'A*02:01', 'A*11:01', 'A*24:02', 'A*26:01', 'A*30:01',
        'A*32:01', 'A*68:01', 'B*07:02', 'B*15:01', 'B*35:01', 'B*40:01',
        'B*44:02', 'B*44:03', 'B*51:01', 'B*53:01', 'B*57:01', 'B*58:01',
    },
    'PtxS1': {
        'A*02:01', 'B*07:02', 'B*35:01', 'B*44:03', 'B*53:01', 'B*57:01',
    },
}

PEPTIDE_LEN_MIN = 8
PEPTIDE_LEN_MAX = 14

# ---------------------------------------------------------------------------
# Helper functions.
# ---------------------------------------------------------------------------

def normalize_allele(allele: str) -> str:
    allele = str(allele).strip()
    return allele[4:] if allele.upper().startswith("HLA-") else allele


def load_validated(xlsx_path: str) -> set:
    df = pd.read_excel(xlsx_path)
    required = {"HLA allele", "CTL epitopes"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Processing details{xlsx_path}Required input or value was not found{missing}")
    pairs = set()
    for _, row in df[["HLA allele", "CTL epitopes"]].dropna().iterrows():
        a = normalize_allele(str(row["HLA allele"]))
        p = str(row["CTL epitopes"]).strip().upper()
        if a and p:
            pairs.add((a, p))
    return pairs


def load_predictions_raw(tsv_path: str) -> pd.DataFrame:
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    df = pd.read_csv(tsv_path, sep="\t")
    df.columns = [c.strip() for c in df.columns]
    return df


def check_tsv_structure(df: pd.DataFrame) -> tuple[bool, list]:
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    critical = False
    warnings = []
    print(f"Validation status")

    # Implementation detail; see the repository documentation.
    missing_cols = REQUIRED_TSV_COLS - set(df.columns)
    if missing_cols:
        msg = f"Required input or value was not found{', '.join(sorted(missing_cols))}"
        warnings.append(("Processing details", msg))
        print(f"Processing details{msg}")
        critical = True
    else:
        print(f"Processing details")

    if critical:
        print(f"  ─────────────────────────────────────────────────────")
        return critical, warnings

    # Implementation detail; see the repository documentation.
    score_cols = [col for col, _, _, _ in THRESHOLDS]
    for col in score_cols:
        converted = pd.to_numeric(df[col], errors="coerce")
        n_missing = converted.isna().sum()  # Implementation detail; see the repository documentation.
        if n_missing > 0:
            if col in ALLOWED_NAN_COLS:
                msg = f"Processing details{col}': {n_missing}Processing details"
                warnings.append(("[!]", msg))
                print(f"  [!] {msg}")
            else:
                msg = f"Processing details{col}': {n_missing}Error"
                warnings.append(("Processing details", msg))
                print(f"Processing details{msg}")
                critical = True
        else:
            print(f"Processing details{col}Processing details")

    # Validation.
    lengths_col = pd.to_numeric(df["peptide length"], errors="coerce")
    n_nan_len = lengths_col.isna().sum()
    if n_nan_len > 0:
        msg = f"Processing details{n_nan_len}Processing details"
        warnings.append(("Processing details", msg))
        print(f"Processing details{msg}")
        critical = True
    else:
        # Implementation detail; see the repository documentation.
        actual_lens = df["peptide"].str.len()
        mismatch = (actual_lens != lengths_col).sum()
        if mismatch > 0:
            msg = f"Validation status{mismatch}Table status"
            warnings.append(("Processing details", msg))
            print(f"Processing details{msg}")
            critical = True
        else:
            print(f"Peptide status")

        # Implementation detail; see the repository documentation.
        out_of_range = df[(lengths_col < PEPTIDE_LEN_MIN) | (lengths_col > PEPTIDE_LEN_MAX)]
        if not out_of_range.empty:
            bad_lens = sorted(lengths_col[
                (lengths_col < PEPTIDE_LEN_MIN) | (lengths_col > PEPTIDE_LEN_MAX)
            ].dropna().unique().tolist())
            msg = (f"Peptide status{PEPTIDE_LEN_MIN}-{PEPTIDE_LEN_MAX}: "
                   f"{len(out_of_range)}Table status{bad_lens}")
            warnings.append(("Processing details", msg))
            print(f"Processing details{msg}")
            critical = True
        else:
            print(f"Peptide status{PEPTIDE_LEN_MIN}-{PEPTIDE_LEN_MAX}")

    # Allele handling.
    peptides_by_allele = {
        allele: set(group["peptide"])
        for allele, group in df.groupby("allele")
    }
    if len(peptides_by_allele) > 1:
        allele_list = list(peptides_by_allele.keys())
        ref_allele  = allele_list[0]
        ref_set     = peptides_by_allele[ref_allele]
        differing   = [a for a in allele_list[1:] if peptides_by_allele[a] != ref_set]
        if differing:
            short_list = ", ".join(sorted(differing)[:5]) + ("..." if len(differing) > 5 else "")
            msg = (f"Allele status{ref_allele}): {short_list}")
            warnings.append(("Processing details", msg))
            print(f"Processing details{msg}")
            critical = True
        else:
            print(f"Allele status{len(ref_set)}Processing details")

    # Allele handling.
    n_dupes = df.duplicated(subset=["allele", "peptide"]).sum()
    if n_dupes > 0:
        msg = f"Processing details{n_dupes}Allele status"
        warnings.append(("[!]", msg))
        print(f"  [!] {msg}")
    else:
        print(f"Allele status")

    print(f"  ─────────────────────────────────────────────────────")
    return critical, warnings


def apply_types(df: pd.DataFrame) -> pd.DataFrame:
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    df = df.copy()
    df["allele"]  = df["allele"].apply(normalize_allele)
    df["peptide"] = df["peptide"].astype(str).str.strip()
    for col, _, _, _ in THRESHOLDS:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def check_alleles(antigen: str, found_alleles: list) -> bool:
    expected = EXPECTED_ALLELES.get(antigen)
    if expected is None:
        print(f"Validation status{antigen}Processing details")
        return True
    found_set = set(found_alleles)
    missing   = expected - found_set
    extra     = found_set - expected
    set_ok    = not missing and not extra
    print(f"Validation status")
    print(f"Validation status{len(expected)}Processing details{len(found_set)}")
    print(f"Validation status{'OK' if set_ok else 'Processing details'}")
    if missing:
        print(f"Required input or value was not found{len(missing)}): {', '.join(sorted(missing))}")
    if extra:
        print(f"Processing details{len(extra)}): {', '.join(sorted(extra))}")
    if set_ok:
        print(f"Validation status")
    print(f"  ─────────────────────────────────────────────────────")
    return set_ok


# ---------------------------------------------------------------------------
# Implementation detail; see the repository documentation.
# ---------------------------------------------------------------------------

HEADER_FILL   = PatternFill("solid", start_color="2F5496", end_color="2F5496")
SUMMARY_FILL  = PatternFill("solid", start_color="1F3864", end_color="1F3864")
HEADER_FONT   = Font(name="Arial", bold=True, color="FFFFFF", size=11)
BODY_FONT     = Font(name="Arial", size=10)
ALIGN_CENTER  = Alignment(horizontal="center", vertical="center", wrap_text=True)
ALIGN_LEFT    = Alignment(horizontal="left",   vertical="center", wrap_text=True)


def _autofit(ws):
    for col in ws.columns:
        max_len   = 0
        col_letter = get_column_letter(col[0].column)
        for cell in col:
            try:
                val = str(cell.value) if cell.value is not None else ""
                max_len = max(max_len, len(val))
            except Exception:
                pass
        ws.column_dimensions[col_letter].width = min(max_len + 4, 60)


def _format_sheet(ws, header_fill=None):
    fill = header_fill or HEADER_FILL
    for cell in ws[1]:
        cell.fill      = fill
        cell.font      = HEADER_FONT
        cell.alignment = ALIGN_CENTER
    for row in ws.iter_rows(min_row=2):
        for cell in row:
            cell.font      = BODY_FONT
            cell.alignment = ALIGN_LEFT
    ws.freeze_panes = "A2"
    _autofit(ws)


def save_result_excel(sheets: dict, path: str):
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        for sheet_name, df in sheets.items():
            df.to_excel(writer, index=False, sheet_name=sheet_name)

    wb = openpyxl.load_workbook(path)
    sheet_names = list(sheets.keys())
    for i, sheet_name in enumerate(sheet_names):
        ws   = wb[sheet_name]
        fill = SUMMARY_FILL if i == 0 else HEADER_FILL
        _format_sheet(ws, fill)

    wb.save(path)
    total = sum(len(df) for name, df in sheets.items() if name != "Summary")
    print(f"Saved output{os.path.basename(path)}  "
          f"({len(sheets)-1}Processing details")


# ---------------------------------------------------------------------------
# Main entry point.
# ---------------------------------------------------------------------------


def compute_per_threshold(pred_df: pd.DataFrame, validated_pairs: set, src: str) -> dict:
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    score_cols = [col for col, _, _, _ in THRESHOLDS]
    results = {}
    union_mask = pd.Series(False, index=pred_df.index)

    def subset_for_pairs(passed_df: pd.DataFrame, pairs: set) -> pd.DataFrame:
        cols = ["allele", "peptide"] + score_cols + ["source_file"]
        if not pairs:
            return pd.DataFrame(columns=cols)
        idx = pd.MultiIndex.from_frame(passed_df[["allele", "peptide"]])
        wanted = pd.MultiIndex.from_tuples(sorted(pairs), names=["allele", "peptide"])
        sub = passed_df[idx.isin(wanted)].copy()
        out = sub[["allele", "peptide"] + score_cols].copy()
        out["source_file"] = src
        return out.drop_duplicates(subset=["allele", "peptide"]).reset_index(drop=True)

    for col, op, val, label in THRESHOLDS:
        mask = pred_df[col] <= val if op == "<=" else pred_df[col] < val
        union_mask |= mask
        passed = pred_df[mask].copy()
        passed_pairs = set(zip(passed["allele"], passed["peptide"]))
        rm_pairs = validated_pairs & passed_pairs
        up_pairs = passed_pairs - validated_pairs
        results[label] = {
            "rm": subset_for_pairs(passed, rm_pairs),
            "up": subset_for_pairs(passed, up_pairs),
            "counts": {"RM": len(rm_pairs), "UP": len(up_pairs)},
        }
        print(f"  [{label}] RM={len(rm_pairs)}  UP={len(up_pairs)}")

    passed = pred_df[union_mask].copy()
    passed_pairs = set(zip(passed["allele"], passed["peptide"]))
    rm_pairs = validated_pairs & passed_pairs
    up_pairs = passed_pairs - validated_pairs
    results["Union_OR"] = {
        "rm": subset_for_pairs(passed, rm_pairs),
        "up": subset_for_pairs(passed, up_pairs),
        "counts": {"RM": len(rm_pairs), "UP": len(up_pairs)},
    }
    print(f"  [Union_OR] RM={len(rm_pairs)}  UP={len(up_pairs)}")
    return results


def build_summary(results: dict) -> pd.DataFrame:
    return pd.DataFrame([
        {"Threshold": label, "RM": data["counts"]["RM"], "UP": data["counts"]["UP"]}
        for label, data in results.items()
    ])


def process_antigen(name: str, tsv_path: str, xlsx_path: str, out_dir: str) -> bool:
    print(f"\n{'='*60}")
    print(f"Processing details{name}")
    print(f"  TSV : {tsv_path}")
    print(f"  XLSX: {xlsx_path}")

    # Validation.
    for label, path in [("TSV", tsv_path), ("XLSX", xlsx_path)]:
        if not os.path.isfile(path):
            print(f"Error{label}Required input or value was not found{path}")
            return False

    # Implementation detail; see the repository documentation.
    df_raw = load_predictions_raw(tsv_path)

    # Validation.
    critical, _ = check_tsv_structure(df_raw)
    if critical:
        print(f"Error")
        print(f"Processing details{name}Processing details")
        return False

    # Normalization.
    pred_df = apply_types(df_raw)

    # Deduplication.
    score_cols = [col for col, _, _, _ in THRESHOLDS]
    agg_dict   = {col: "min" for col in score_cols if col in pred_df.columns}
    pred_df    = pred_df.groupby(["allele", "peptide"], as_index=False).agg(agg_dict)

    # Validation.
    tool_alleles = sorted(pred_df["allele"].dropna().unique().tolist())
    print(f"Allele status{len(tool_alleles)}):")
    for a in tool_alleles:
        print(f"    {a}")

    if not check_alleles(name, tool_alleles):
        print(f"Error")
        print(f"Processing details{name}Processing details")
        return False

    validated_pairs = load_validated(xlsx_path)
    print(f"Processing details{len(validated_pairs)}")
    print(f"Table status{len(pred_df)}")

    os.makedirs(out_dir, exist_ok=True)
    src              = os.path.basename(tsv_path)

    print(f"Processing details")
    results = compute_per_threshold(pred_df, validated_pairs, src)

    summary_df = build_summary(results)
    sheet_order = ["Summary"] + [col for _, _, _, col in THRESHOLDS] + ["Union_OR"]

    for file_type, key in [("RM", "rm"), ("UP", "up")]:
        sheets = {"Summary": summary_df}
        for label in sheet_order[1:]:
            sheets[label] = results[label][key]
        save_result_excel(
            sheets,
            os.path.join(out_dir, f"{file_type}_{name}_{TOOL}.xlsx")
        )

    return True


def main():
    failed = []
    for antigen, paths in ANTIGENS.items():
        try:
            ok = process_antigen(
                name      = antigen,
                tsv_path  = paths["tsv"],
                xlsx_path = paths["xlsx"],
                out_dir   = paths["out"],
            )
            if not ok:
                failed.append(antigen)
        except Exception as exc:
            import traceback
            print(f"Processing details{antigen}: {exc}")
            traceback.print_exc()
            failed.append(antigen)

    print(f"\n{'='*60}")
    if failed:
        print(f"Error{', '.join(failed)}")
        sys.exit(1)
    else:
        print("Processing status")


if __name__ == "__main__":
    main()
