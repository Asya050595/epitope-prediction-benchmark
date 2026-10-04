#!/usr/bin/env python3
"""
NetMHCII 2.3 — exact allele-peptide RM/UP.

Thresholds: %Rank <= 10, Affinity < 5000 nM, and Union_OR.
Reference columns are read by name ('HLA allele.1', 'HTL epitopes').
RM = passed prediction present in the positive reference dataset.
UP = passed prediction absent from the positive reference dataset.
No FN calculation or export is performed.
"""

from project_paths import DATA_ROOT

import os
import re
import glob
import openpyxl
import pandas as pd

BASE      = f"{DATA_ROOT}"
TOOL_NAME = "NetMHC_II_2.3"

ANTIGENS = {
    "p24": {
        "excel":     os.path.join(BASE, "p24",   "p24.xlsx"),
        "html_glob": os.path.join(BASE, "p24",   "MHC II", "NetMHC II", "NetMHCII_p24_part_*_*.html"),
    },
    "pp65": {
        "excel":     os.path.join(BASE, "pp65",  "pp65.xlsx"),
        "html_glob": os.path.join(BASE, "pp65",  "MHC II", "NetMHC II", "NetMHCII_pp65_*.html"),
    },
    "PtxS1": {
        "excel":     os.path.join(BASE, "PtxS1", "PtxS1.xlsx"),
        "html_glob": os.path.join(BASE, "PtxS1", "MHC II", "NetMHC II", "NetMHCII_PtxS1_*.html"),
    },
}

RANK_THRESHOLD     = 10.0
AFFINITY_THRESHOLD = 5000.0

# Output generation.
# Implementation detail; see the repository documentation.
# Allele handling.
EXPECTED_ALLELES = {
    'p24': [
        'DRB1_0101', 'DRB1_0301', 'DRB1_0401', 'DRB1_0404', 'DRB1_0405', 'DRB1_0701',
        'DRB1_0801', 'DRB1_0901', 'DRB1_1001', 'DRB1_1101', 'DRB1_1301', 'DRB1_1302',
        'DRB1_1501', 'DRB3_0101', 'DRB3_0301', 'DRB4_0101', 'DRB5_0101',
    ],
    'pp65': [
        'DRB1_0301', 'DRB1_0401', 'DRB1_0402', 'DRB1_0404', 'DRB1_0701', 'DRB1_1101',
        'DRB1_1501', 'DRB3_0101', 'DRB3_0202',
    ],
    'PtxS1': [
        'DRB1_0101', 'DRB1_1101',
    ],
}

# Threshold handling.
MODES = ["rank", "affinity", "union_or"]

MODE_SHEET_NAMES = {
    "rank":     "Rank_lte10",
    "affinity": "Affinity_lt5000",
    "union_or": "Union_OR",
}

MODE_LABELS = {
    "rank":     f"%Rank <= {RANK_THRESHOLD:g}",
    "affinity": f"Affinity < {AFFINITY_THRESHOLD:g}",
    "union_or": "Union OR",
}


# ──────────────────────────────────────────────────────────────────────────────

def get_output_dir(antigen: str) -> str:
    return os.path.join(
        BASE, f"Processing details{antigen}",
        "Matches MHC II", f"Matches {TOOL_NAME}",
    )


def html_allele_to_std(allele_html: str) -> str:
    m = re.match(r'^([A-Z0-9]+)_(\d{2})(\d{2})$', allele_html.strip())
    if m:
        return f"{m.group(1)}*{m.group(2)}:{m.group(3)}"
    return allele_html.strip()


def parse_html_predictions(html_path: str) -> list[dict]:
    records = []
    pattern = re.compile(
        r'^\s+'
        r'((?:DRB[0-9]|DP[AB][0-9]|DQ[AB][0-9])_\d+)'
        r'\s+(\d+)'
        r'\s+([A-Z]+)'
        r'\s+([A-Z]+)'
        r'\s+(-?\d+)'
        r'\s+(\S+)'
        r'\s+(\S+)'
        r'\s+(\S+)'
        r'\s+(\S+)'
        r'\s+(\S+)'
        r'(?:\s+(SB|WB))?'
        r'\s*$'
    )
    try:
        with open(html_path, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                m = pattern.match(line)
                if not m:
                    continue
                try:
                    affinity = float(m.group(7))
                    rank     = float(m.group(8))
                except ValueError:
                    continue
                records.append({
                    "allele_html": m.group(1),
                    "allele_std":  html_allele_to_std(m.group(1)),
                    "peptide":     m.group(3).upper(),
                    "affinity":    affinity,
                    "rank":        rank,
                    "bind_level":  m.group(11) or "",
                    "source_file": os.path.basename(html_path),
                })
    except FileNotFoundError:
        print(f"Required input or value was not found{html_path}")
    return records




def normalize_std_allele(raw: str) -> str:
    s = str(raw).strip().upper().replace("HLA-", "")
    m = re.match(r'^([A-Z]+\d)[*_]?(\d{2})[:_]?(\d{2})$', s)
    if not m:
        raise ValueError(f"Processing details{raw!r}")
    return f"{m.group(1)}*{m.group(2)}:{m.group(3)}"


def load_validated_htl(excel_path: str) -> set[tuple[str, str]]:
    df = pd.read_excel(excel_path, header=0)
    allele_col, peptide_col = "HLA allele.1", "HTL epitopes"
    missing = [c for c in (allele_col, peptide_col) if c not in df.columns]
    if missing:
        raise ValueError(f"Processing details{excel_path}Required input or value was not found{missing}")
    pairs = set()
    for _, row in df[[allele_col, peptide_col]].dropna().iterrows():
        pairs.add((normalize_std_allele(row[allele_col]), str(row[peptide_col]).strip().upper()))
    if not pairs:
        raise ValueError(f"Processing details{excel_path}Processing details")
    return pairs


def validate_alleles(antigen_name: str, found_std: set[str]) -> None:
    expected_std = {html_allele_to_std(a) for a in EXPECTED_ALLELES[antigen_name]}
    missing, extra = expected_std - found_std, found_std - expected_std
    print(f"\n  Prediction alleles: expected={len(expected_std)}, found={len(found_std)}")
    if missing: print(f"Required input or value was not found{', '.join(sorted(missing))}")
    if extra: print(f"Processing details{', '.join(sorted(extra))}")
    if missing or extra:
        raise ValueError("Validation status")


def aggregate_predictions(all_records: list[dict]) -> pd.DataFrame:
    df = pd.DataFrame(all_records).rename(columns={"allele_std":"Allele","peptide":"Peptide","rank":"%Rank","affinity":"Affinity(nM)","source_file":"Source File"})
    return df.groupby(["Allele","Peptide"], as_index=False).agg({"%Rank":"min","Affinity(nM)":"min","Source File":lambda s:"; ".join(dict.fromkeys(map(str,s)))})


def build_mode_frames(mode: str, pred_df: pd.DataFrame, validated_set: set) -> dict:
    if mode == "rank": mask = pred_df["%Rank"] <= RANK_THRESHOLD
    elif mode == "affinity": mask = pred_df["Affinity(nM)"] < AFFINITY_THRESHOLD
    else: mask = (pred_df["%Rank"] <= RANK_THRESHOLD) | (pred_df["Affinity(nM)"] < AFFINITY_THRESHOLD)
    passed = pred_df[mask].copy()
    pairs = set(zip(passed["Allele"], passed["Peptide"]))
    rm_pairs, up_pairs = pairs & validated_set, pairs - validated_set
    if mode == "rank": cols=["Allele","Peptide","%Rank","Source File"]
    elif mode == "affinity": cols=["Allele","Peptide","Affinity(nM)","Source File"]
    else:
        passed["Pass_Rank"] = passed["%Rank"] <= RANK_THRESHOLD
        passed["Pass_Affinity"] = passed["Affinity(nM)"] < AFFINITY_THRESHOLD
        passed["Passed_Thresholds_Count"] = passed[["Pass_Rank","Pass_Affinity"]].sum(axis=1)
        cols=["Allele","Peptide","%Rank","Affinity(nM)","Pass_Rank","Pass_Affinity","Passed_Thresholds_Count","Source File"]
    idx=pd.MultiIndex.from_frame(passed[["Allele","Peptide"]])
    def subset(ps):
        if not ps: return pd.DataFrame(columns=cols)
        wanted=pd.MultiIndex.from_tuples(sorted(ps),names=["Allele","Peptide"])
        return passed[idx.isin(wanted)][cols].reset_index(drop=True)
    return {"rm_df":subset(rm_pairs),"up_df":subset(up_pairs),"rm_pairs":rm_pairs,"up_pairs":up_pairs}


def save_excel_multi(sheets: dict[str, pd.DataFrame], path: str) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        for sheet_name, df in sheets.items():
            df.to_excel(writer, sheet_name=sheet_name, index=False)
    row_counts = ", ".join(f"{name}={len(df)}" for name, df in sheets.items())
    print(f"Saved output{path}\n    ({row_counts})")


# ──────────────────────────────────────────────────────────────────────────────



def process_antigen(antigen_name: str, cfg: dict) -> str:
    print(f"\n{'='*60}Processing details{antigen_name}\n{'='*60}")
    if not os.path.isfile(cfg["excel"]): return "Required input or value was not found"
    html_files=sorted(glob.glob(cfg["html_glob"]))
    if not html_files: return "Required input or value was not found"
    validated_set=load_validated_htl(cfg["excel"])
    all_records=[]
    for html_file in html_files:
        recs=parse_html_predictions(html_file)
        if not recs: raise ValueError(f"Processing details{html_file}")
        all_records.extend(recs)
    validate_alleles(antigen_name,{r["allele_std"] for r in all_records})
    pred_df=aggregate_predictions(all_records)
    results={m:build_mode_frames(m,pred_df,validated_set) for m in MODES}
    print(f"\n  {'Mode':<18}{'RM':>8}{'UP':>10}")
    for m in MODES: print(f"  {MODE_LABELS[m]:<18}{len(results[m]['rm_pairs']):>8}{len(results[m]['up_pairs']):>10}")
    out_dir=get_output_dir(antigen_name)
    for kind,key in (("RM","rm_df"),("UP","up_df")):
        save_excel_multi({MODE_SHEET_NAMES[m]:results[m][key] for m in MODES},os.path.join(out_dir,f"{kind}_{antigen_name}_{TOOL_NAME}.xlsx"))
    return "success"


def main():
    statuses={}
    for antigen,cfg in ANTIGENS.items():
        try: statuses[antigen]=process_antigen(antigen,cfg)
        except Exception as exc: print(f"Error{antigen}: {exc}"); statuses[antigen]=f"Error{exc}"
    for a,s in statuses.items(): print(f"{a}: {s}")
    if any(s!="success" for s in statuses.values()): raise SystemExit(1)


if __name__ == "__main__":
    print("Validation status")
    print(f"Processing details{RANK_THRESHOLD:g} / Affinity < {AFFINITY_THRESHOLD:g} nM) "
          f"Processing details")

    for antigen, cfg in ANTIGENS.items():
        process_antigen(antigen, cfg)

    print("Completed successfully")
