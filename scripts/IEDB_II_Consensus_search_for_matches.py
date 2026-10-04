"""English documentation for this module or helper is provided in the repository README and in the surrounding code."""

from project_paths import DATA_ROOT

import os
import re
import sys
import pandas as pd

# Configuration.

BASE = f"{DATA_ROOT}"
TOOL_NAME = "IEDB_II_Consensus"

ANTIGENS = {
    "p24": {
        "tsv":   f"{BASE}/p24/MHC II/IEDB II/IEDB_II_Consensus_p24_11-18_pep1.tsv",
        "excel": f"{BASE}/p24/p24.xlsx",
    },
    "pp65": {
        "tsv":   f"{BASE}/pp65/MHC II/IEDB II/IEDB_II_Consensus_pp65_11-18_pep1.tsv",
        "excel": f"{BASE}/pp65/pp65.xlsx",
    },
    "PtxS1": {
        "tsv":   f"{BASE}/PtxS1/MHC II/IEDB II/IEDB_II_Consensus_PtxS1_11-18_pep1.tsv",
        "excel": f"{BASE}/PtxS1/PtxS1.xlsx",
    },
}

# Threshold handling.
THRESHOLDS = [
    ("Median_percentile_lte10",    "median binding percentile", "le", 10.0),
    ("Consensus_percentile_lte10", "consensus percentile",      "le", 10.0),
    ("IC50_SMM_lt5000",            "smm_align ic50",            "lt", 5000.0),
    ("IC50_NN_lt5000",             "nn_align ic50",             "lt", 5000.0),
]
SCORE_COLS = [col for _, col, _, _ in THRESHOLDS]

EXCEL_ALLELE_COL  = "HLA allele.1"
EXCEL_PEPTIDE_COL = "HTL epitopes"

# Validation.
# Normalization.
EXPECTED_ALLELES = {
    'p24': {
        'DRB1*01:01', 'DRB1*03:01', 'DRB1*04:01', 'DRB1*04:05', 'DRB1*07:01', 'DRB1*09:01',
        'DRB1*11:01', 'DRB1*13:02', 'DRB1*15:01', 'DRB3*01:01', 'DRB4*01:01', 'DRB5*01:01',
    },
    'pp65': {
        'DRB1*03:01', 'DRB1*04:01', 'DRB1*07:01', 'DRB1*11:01', 'DRB1*15:01', 'DRB3*01:01',
        'DRB3*02:02',
    },
    'PtxS1': {
        'DRB1*01:01', 'DRB1*11:01',
    },
}


class AntigenSkipped(Exception):
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""


# Normalization.

# Implementation detail; see the repository documentation.
_DASH_LOOKALIKES = {
    "\u2010": "-", "\u2011": "-", "\u2012": "-", "\u2013": "-",
    "\u2014": "-", "\u2212": "-", "\uFE63": "-", "\uFF0D": "-",
    "\u00AD": "",  # Implementation detail; see the repository documentation.
}


def _normalize_dashes(s: str, raw_for_report) -> str:
    fixed = s
    for bad_char, repl in _DASH_LOOKALIKES.items():
        if bad_char in fixed:
            print(f"Allele status{raw_for_report!r}Processing details"
                  f"(U+{ord(bad_char):04X}Processing details")
            fixed = fixed.replace(bad_char, repl)
    return fixed


def normalize_allele(raw) -> str:
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    if raw is None:
        raise ValueError("Required input or value was not found")
    s = str(raw).strip()
    if not s:
        raise ValueError("Allele status")

    s = _normalize_dashes(s, raw)
    s = re.sub(r'^HLA-', '', s, flags=re.IGNORECASE)

    m = re.match(r'^([A-Z]+\d)[\*_]?(\d{2})[:_]?(\d{2})$', s, flags=re.IGNORECASE)
    if not m:
        raise ValueError(f"Allele status{raw!r}")

    gene, group, protein = m.groups()
    return f"{gene.upper()}*{group}:{protein}"


def normalize_alleles_column(series: pd.Series, label: str) -> tuple[pd.Series, list]:
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    normed, errors = [], []
    for i, raw in series.items():
        try:
            normed.append(normalize_allele(raw))
        except ValueError as exc:
            errors.append(f"{label}Table status{i}: {exc}")
            normed.append(None)
    return pd.Series(normed, index=series.index), errors


def normalize_peptide(p) -> str:
    return str(p).strip().upper()


# Input loading.

def load_experimental(excel_path: str) -> pd.DataFrame:
    df = pd.read_excel(excel_path, header=0)
    if EXCEL_ALLELE_COL not in df.columns or EXCEL_PEPTIDE_COL not in df.columns:
        raise AntigenSkipped(
            f"Processing details{excel_path}Required input or value was not found"
            f"'{EXCEL_ALLELE_COL}Processing details{EXCEL_PEPTIDE_COL}'. "
            f"Table status{list(df.columns)}"
        )
    htl = (
        df[[EXCEL_ALLELE_COL, EXCEL_PEPTIDE_COL]]
        .dropna(subset=[EXCEL_ALLELE_COL, EXCEL_PEPTIDE_COL])
        .copy()
    )
    htl.columns = ["allele_raw", "peptide_raw"]

    allele_norm, errors = normalize_alleles_column(htl["allele_raw"], "Excel")
    if errors:
        raise AntigenSkipped(
            "Allele status" + "\n      ".join(errors)
        )

    htl["allele"]  = allele_norm
    htl["peptide"] = htl["peptide_raw"].apply(normalize_peptide)
    htl = htl[["allele", "peptide"]].drop_duplicates().reset_index(drop=True)
    print(f"Processing details{len(htl)}")
    return htl


def load_predictions(tsv_path: str) -> pd.DataFrame:
    df = pd.read_csv(tsv_path, sep="\t", header=0, low_memory=False)

    required = ["allele", "peptide"] + SCORE_COLS
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise AntigenSkipped(
            f"Required input or value was not found{missing}. "
            f"Table status{list(df.columns)}"
        )

    allele_norm, errors = normalize_alleles_column(df["allele"], "TSV")
    if errors:
        raise AntigenSkipped(
            "Allele status" + "\n      ".join(errors)
        )

    df["allele_norm"] = allele_norm
    df["peptide"] = df["peptide"].apply(normalize_peptide)
    for col in SCORE_COLS:
        before_na = df[col].isna().sum()
        df[col] = pd.to_numeric(df[col], errors="coerce")
        after_na = df[col].isna().sum()
        if after_na:
            print(f"Table status{col}Processing details{after_na}")

    # Metric calculation.
    # Implementation detail; see the repository documentation.
    # Implementation detail; see the repository documentation.
    # Implementation detail; see the repository documentation.
    required_numeric = {"median binding percentile"}
    if TOOL_NAME == "IEDB_II_Consensus":
        required_numeric.add("consensus percentile")
    elif TOOL_NAME.endswith("_BA"):
        required_numeric.add("netmhciipan_ba ic50")
    bad_required = [c for c in required_numeric if df[c].isna().any()]
    if bad_required:
        raise AntigenSkipped(
            "Processing details" + ", ".join(sorted(bad_required))
        )

    print(f"Table status{len(df)}")
    return df


# Validation.

def validate_alleles(antigen: str, tool_alleles: set) -> None:
    print(f"Allele status{len(tool_alleles)}):")
    for a in sorted(tool_alleles):
        print(f"    {a}")

    expected = EXPECTED_ALLELES[antigen]
    extra   = tool_alleles - expected
    missing = expected - tool_alleles

    if len(tool_alleles) == len(expected) and not extra and not missing:
        print(f"Validation status{len(tool_alleles)}")
        print(f"Validation status")
        return

    print(f"Validation status{len(expected)}Processing details{len(tool_alleles)}")
    if extra:
        print(f"Processing details{len(extra)}): {', '.join(sorted(extra))}")
    if missing:
        print(f"Required input or value was not found{len(missing)}): {', '.join(sorted(missing))}")
    raise AntigenSkipped("Validation status")


# ═══════════════════════════════ RM / UP ═══════════════════════════════

def make_mask(df_pred: pd.DataFrame, col: str, op: str, value: float) -> pd.Series:
    if op == "le":
        return df_pred[col] <= value
    if op == "lt":
        return df_pred[col] < value
    raise ValueError(f"Processing details{op}")



def build_rm_up(df_passed: pd.DataFrame, exp_set: set, src: str):
    """English documentation for this module or helper is provided in the repository README and in the surrounding code."""
    if df_passed.empty:
        cols = ["allele", "peptide"] + SCORE_COLS + ["source_file"]
        empty = pd.DataFrame(columns=cols)
        return empty.copy(), empty.copy()

    # Output generation.
    agg = {c: "min" for c in SCORE_COLS}
    passed = (
        df_passed[["allele_norm", "peptide"] + SCORE_COLS]
        .groupby(["allele_norm", "peptide"], as_index=False)
        .agg(agg)
    )
    passed_set = set(zip(passed["allele_norm"], passed["peptide"]))
    rm_pairs = exp_set & passed_set
    up_pairs = passed_set - exp_set

    def slim(pairs: set) -> pd.DataFrame:
        cols = ["allele", "peptide"] + SCORE_COLS + ["source_file"]
        if not pairs:
            return pd.DataFrame(columns=cols)
        idx = pd.MultiIndex.from_frame(passed[["allele_norm", "peptide"]])
        wanted = pd.MultiIndex.from_tuples(sorted(pairs), names=["allele_norm", "peptide"])
        sub = passed[idx.isin(wanted)].copy()
        sub = sub.rename(columns={"allele_norm": "allele"})
        sub["source_file"] = src
        return sub[cols].reset_index(drop=True)

    return slim(rm_pairs), slim(up_pairs)


# Output generation.

def get_output_dir(antigen: str) -> str:
    return os.path.join(
        BASE,
        f"Processing details{antigen}",
        "Matches MHC II",
        f"Matches {TOOL_NAME}",
    )


def save_multisheet(sheets: dict, path: str) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        for sheet_name, df in sheets.items():
            df.to_excel(writer, sheet_name=sheet_name[:31], index=False)


# Implementation detail; see the repository documentation.

def process_antigen(antigen: str, paths: dict) -> bool:
    print(f"\n{'='*60}")
    print(f"Processing details{antigen}")
    print(f"{'='*60}")

    try:
        for key in ("tsv", "excel"):
            if not os.path.isfile(paths[key]):
                raise AntigenSkipped(f"Required input or value was not found{paths[key]}")

        df_exp  = load_experimental(paths["excel"])
        df_pred = load_predictions(paths["tsv"])

        tool_alleles = set(df_pred["allele_norm"].dropna().unique())
        validate_alleles(antigen, tool_alleles)

        exp_set = set(zip(df_exp["allele"], df_exp["peptide"]))
        src = os.path.basename(paths["tsv"])

        sheets_rm, sheets_up = {}, {}
        union_mask = pd.Series(False, index=df_pred.index)

        print()
        for sheet_name, col, op, value in THRESHOLDS:
            mask = make_mask(df_pred, col, op, value)
            union_mask = union_mask | mask
            df_passed = df_pred[mask]
            rm_df, up_df = build_rm_up(df_passed, exp_set, src)
            sheets_rm[sheet_name] = rm_df
            sheets_up[sheet_name] = up_df
            print(f"  [{sheet_name}Processing details{len(df_passed)}  |  RM: {len(rm_df)}, UP: {len(up_df)}")

        df_passed_union = df_pred[union_mask]
        rm_df, up_df = build_rm_up(df_passed_union, exp_set, src)
        sheets_rm["Union_OR"] = rm_df
        sheets_up["Union_OR"] = up_df
        print(f"Processing details{len(df_passed_union)}  |  RM: {len(rm_df)}, UP: {len(up_df)}")

        out_dir = get_output_dir(antigen)
        save_multisheet(sheets_rm, os.path.join(out_dir, f"RM_{antigen}_{TOOL_NAME}.xlsx"))
        save_multisheet(sheets_up, os.path.join(out_dir, f"UP_{antigen}_{TOOL_NAME}.xlsx"))
        print(f"Saved output{out_dir}")
        return True

    except AntigenSkipped as exc:
        print(f"\n  [ERROR] {exc}")
        print(f"Processing details{antigen}Processing details")
        return False


# ═══════════════════════════════ MAIN ═══════════════════════════════

def main():
    print("Validation status")
    print("Processing details"
          "IC50_SMM_lt5000, IC50_NN_lt5000, Union_OR")
    print("Processing details", ", ".join(ANTIGENS.keys()))

    results = {}
    for antigen, paths in ANTIGENS.items():
        try:
            results[antigen] = process_antigen(antigen, paths)
        except Exception as exc:  # Error handling.
            print(f"Error{antigen}: {exc}")
            results[antigen] = False

    ok     = [a for a, s in results.items() if s]
    failed = [a for a, s in results.items() if not s]

    print(f"\n{'='*60}")
    print(f"Processing status{len(ok)}/{len(ANTIGENS)} ({', '.join(ok) if ok else '—'})")
    if failed:
        print(f"Error{len(failed)}/{len(ANTIGENS)} ({', '.join(failed)})")
    print(f"{'='*60}\n")

    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
